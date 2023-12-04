/*
Author: Mateusz Frejlich
Index: ******
Subject: Optyczne techniki skanowania i analizy danych 3D (OPTS)
Semester: 2022L
Task: 2
Object to find: 'glosnik'
*/

#include <ogx/Plugins/EasyPlugin.h>
#include <ogx/Data/Clouds/CloudHelpers.h>
#include <random>

using namespace ogx;
using namespace Data;

struct LoudspeakerSearching : public Plugin::EasyMethod 
{
	/*HARDCODED PARAMETERS AND DATA*/
	const double HAUSDORF_THRESHOLD = 50000;	// number of points which change version of Hausdorf algorithm

	/*Structure store data of a 3D point with vector of pointers to its neighbours, which make it an element of a graph structure.*/
	struct Point
	{
		Clouds::Point3D xyz;
		Clouds::Color color;
		std::vector<Point*> neighbours;
		int group_number = 0;
		bool is_calculated = 0;
		bool is_deleted = 0;
		Point(Clouds::Point3D p_coordinates, Clouds::Color p_color) :
			xyz(p_coordinates),
			color(p_color)
		{
			// point parameters have been set beyond
		}
	};

	/*Vector of graph elements with points parameters and pointers to their neighbours inside a spherical search kernel. It is used to save the time for searching neighbours in each operation. It might be made one time for many layers with diffrent parameters such as:
	* max color gradients in HSV
	* min points in a group
	* min neighbous of a point
	* dimensions of the loudspeaker
	if only search kernel radius is not	changed and operations are done for the same cloud.*/
	std::vector<Point> points_graph_structure;

	// parameters
	ResourceID m_node_id;
	float max_distance = 0.020;
	float points_to_remove_ratio = 0.9;
	int min_neighbours = 10;
	int min_points_group = 100;
	int neighbours_to_smooth = 10;
	int median_filter_size = 10;
	float max_brightness_gradient = 7.5;
	float max_saturation_gradient = 2.5;
	float max_hue_gradient = 180;
	float coverage_rate = 0.4;					// additional tolerance of hue and saturation differences for mean HSV=(*,1,1), 10 times smaller for HSV=(*,10,10) etc.
	float join_coefficient = 0.8;				// maximum ratio of ungrouped point among the neighbours to join the ungrouped point to any of the clusters
	float merge_coefficient = 0.02;				// the higher it is the more clusters will be merged
	float speaker_dim_x = 0.18;
	float speaker_dim_y = 0.16;
	float speaker_dim_z = 0.12;

	// constructor
	LoudspeakerSearching() : EasyMethod(L"Mateusz Frejlich", L"The program segments the cloud and recognizes the loudspeaker inside it.\nThe default values of all parameters are best fitted to the tested cloud.")
	{
	}

	// add input/output parameters
	virtual void DefineParameters(ParameterBank& bank) 
	{
		bank.Add(L"Node ID", m_node_id).AsNode();
		bank.Add(L"Points to remove ratio", points_to_remove_ratio).Min(0.0).Max(1.0);
		bank.Add(L"Max. distance to neighbour in a cluster [u]", max_distance).Describe(L"In distance unit used in the cloud.");
		bank.Add(L"Min. number of neighbours", min_neighbours).Describe(L"Points below this value will be deleted.");
		bank.Add(L"Min. number of points in a cluster", min_points_group).Describe(L"Clusters below this value will be canceled.");
		bank.Add(L"Number of neigbours to smooth", neighbours_to_smooth).Describe(L"Number of neighbours used to smoothing each point by fitting a surface. If there are less neighbours, smoothing will not be done.");
		bank.Add(L"Median filter size", median_filter_size).Describe(L"Number of neighbours used to smoothing the color of each point by median filter. If there are less neighbours, smoothing will not be done.");
		bank.Add(L"Max. hue gradient in a cluster [\260/u]", max_hue_gradient).Describe(L"One degree per distance unit.");
		bank.Add(L"Max. saturation gradient in a cluster [1/u]", max_saturation_gradient).Describe(L"Max. scope value per distance unit.");
		bank.Add(L"Max. brightness gradient in a cluster [1/u]", max_brightness_gradient).Describe(L"Max. scope value per distance unit.");
		bank.Add(L"Coverage rate", coverage_rate).Describe(L"The higher it is, the higher tolerance based on measurement inaccuracy of hue and saturation will be considered.");
		bank.Add(L"Join coefficient", join_coefficient).Describe(L"The higher it is, the more ungrouped points will be joined to any cluster.");
		bank.Add(L"Merge coefficient", merge_coefficient).Describe(L"The higher it is, the more clusters will be merged.");
		bank.Add(L"Loudspeaker dimension X [u]", speaker_dim_x).Describe(L"The largest dimension");
		bank.Add(L"Loudspeaker dimension Y [u]", speaker_dim_y).Describe(L"Medium dimention");
		bank.Add(L"Loudspeaker dimension Z [u]", speaker_dim_z).Describe(L"The smallest dimension");
	}

	virtual void Run(Context& context)
	{
		// get cloud data
		Clouds::ICloud* original_cloud = GetValidCloud(context, m_node_id);
		OGX_LINE.Msg(Info, L"All data is correct. Processing has started.");

		// reduce the number of points
		ResourceID new_node_id;
		Clouds::ICloud* cloud = CreateRandomlyReducedCloud(context, original_cloud, new_node_id);
		OGX_LINE.Msg(Info, L"Reduced cloud has been done.");

		// create points structure from created cloud
		BuildPointsGraphStructure(context, new_node_id);
		OGX_LINE.Msg(Info, L"Points structure has been done.");

		// create a new layer
		const String LAYER_NAME = L"Segmentation";
		auto layers = cloud->FindLayers(LAYER_NAME);
		auto segments_layer = layers.empty() ? cloud->CreateLayer(LAYER_NAME, 0.0) : layers[0];

		// reserve memory for results vector
		Clouds::PointsRange range;
		cloud->GetAccess().GetAllPoints(range);
		std::vector<StoredReal> segment_numbers;
		segment_numbers.reserve(points_graph_structure.size());

		// remove lonely points and smooth the cloud
		RemoveNoise(context);

		// reduce noise in colors using median filter
		MedianFilterColors(context);
		OGX_LINE.Msg(Info, L"Noise filtration has been done.");

		// visualize the cloud smoothing
		CopyGraphStructuctureToCloud(range);

		// make segmentation using Hausdorf algorithm
		Hausdorf(context);
		OGX_LINE.Msg(Info, L"Segmentation by Hausdorf algorithm has been done.");

		// fill holes with ungrouped points in clusters
		FillHolesInClusters(context);

		// remove too small groups
		RemoveTooSmallGroups(context);
		OGX_LINE.Msg(Info, L"Removing too small groups has been done.");

		// fill holes made by removing small groups
		FillHolesInClusters(context);
		OGX_LINE.Msg(Info, L"Filling holes in clusters has been done.");

		// merge clusters
		MergeClusters(context);
		OGX_LINE.Msg(Info, L"Merging clusters has been done.");

		// change group numbers to subsequent numbers with '0' assigned to ungrouped points
		TidyUpGroupNumbers(context, segment_numbers);
		OGX_LINE.Msg(Info, L"The ordering of group numbers has been done.");

		// assign all cluster numbers to created layer
		range.SetLayerVals(segment_numbers, *segments_layer);
		OGX_LINE.Msg(Info, L"Whole segmentation process has been finished.");

		// create next layer
		const String NEW_LAYER_NAME = L"Loudspeaker";
		auto new_layers = cloud->FindLayers(NEW_LAYER_NAME);
		auto speaker_layer = new_layers.empty() ? cloud->CreateLayer(NEW_LAYER_NAME, 0.0) : new_layers[0];
		std::vector<StoredReal> is_speaker;
		is_speaker.reserve(points_graph_structure.size());

		// find the best fitted cluster to the speaker
		FindTheSpeaker(context, is_speaker);
		OGX_LINE.Msg(Info, L"Loudspeaker identification has been done.");

		// mark found cluster in created layer
		range.SetLayerVals(is_speaker, *speaker_layer);
		OGX_LINE.Msg(Info, L"Processing has been finished. No errors.");

		//// check time
		//auto start = std::chrono::high_resolution_clock::now();
		//auto stop = std::chrono::high_resolution_clock::now();
		//auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		//auto time = duration.count();

		// release memory ???
	}

	/*Function reads data of all points from given node, check their correctness and returns a cloud of these points.*/
	Clouds::ICloud* GetValidCloud(Context& context, ResourceID node_id)
	{
		auto node = context.m_project->TransTreeFindNode(node_id);
		if (!node)
		{
			ReportError(L"Invalid node!");
		}
		auto element = node->GetElement();
		if (!element)
		{
			ReportError(L"Invalid element!");
		}
		auto cloud = element->GetData<Clouds::ICloud>();
		if (!cloud)
		{
			ReportError(L"Invalid cloud!");
		}
		if (points_to_remove_ratio < 0.0 || points_to_remove_ratio >= 1.0)
		{
			ReportError(L"Invalid points to remove ratio!");
		}
		if (max_distance <= 0.0)
		{
			ReportError(L"Invalid maximum distance to neighbour!");
		}
		if (min_neighbours < 0)
		{
			ReportError(L"Invalid minimum neigbours number!");
		}
		if (min_points_group < 0)
		{
			ReportError(L"Invalid minimum points cluster!");
		}
		if (neighbours_to_smooth < 5)
		{
			ReportError(L"Invalid number of neighbours to smooth!");
		}
		if (median_filter_size < 2)
		{
			ReportError(L"Invalid median filter size!");
		}
		if (max_brightness_gradient < 0.0)
		{
			ReportError(L"Invalid maximum brightness gradient!");
		}
		if (max_saturation_gradient < 0.0)
		{
			ReportError(L"Invalid maximum saturation gradient!");
		}
		if (max_hue_gradient < 0.0)
		{
			ReportError(L"Invalid maximum hue gradient!");
		}
		if (coverage_rate < 0.0)
		{
			ReportError(L"Invalid coverage rate!");
		}
		if (join_coefficient < 0.0 || join_coefficient > 1.0)
		{
			ReportError(L"Invalid join coefficient!");
		}
		if (merge_coefficient < 0.0)
		{
			ReportError(L"Invalid merge coefficient!");
		}
		if (speaker_dim_x <= 0.0 || speaker_dim_x < speaker_dim_y || speaker_dim_x < speaker_dim_z)
		{
			ReportError(L"Invalid speaker dimention X!");
		}
		if (speaker_dim_y <= 0.0 || speaker_dim_y < speaker_dim_z)
		{
			ReportError(L"Invalid speaker dimention Y!");
		}
		if (speaker_dim_z <= 0.0)
		{
			ReportError(L"Invalid speaker dimention Z!");
		}
		return cloud;
	}

	/*Function makes a new node with given ID, assign to them the cloud with random points from a given cloud with ratio given by a user.*/
	Clouds::ICloud* CreateRandomlyReducedCloud(Context& context, Clouds::ICloud*& cloud, ResourceID& node_id)
	{
		// check if cloud must be reduced
		if (points_to_remove_ratio == 0.0)
			return cloud;
		// create reduced cloud
		auto reduced_node = context.m_project->TransTreeFindNode(m_node_id)->CreateChild();
		reduced_node->Rename(L"Reduced cloud");
		auto reduced_element = context.Project().ElementCreate<Clouds::ICloud>();
		reduced_node->SetElement(reduced_element);
		auto reduced_cloud = reduced_element->GetData<Clouds::ICloud>();
		// get points parameters
		Clouds::PointsRange range;
		cloud->GetAccess().GetAllPoints(range);
		std::vector<Clouds::Point3D> coordinates;
		range.GetXYZ(coordinates);
		std::vector<Clouds::Color> colors;
		range.GetColors(colors);
		std::vector<Clouds::Point3D> reduced_coordinates;
		std::vector<Clouds::Color> reduced_colors;
		// get random points from cloud
		std::mt19937 random_numbers_engine;
		std::uniform_real_distribution<float> distribution(0.0, 1.0);
		float progress = 0.0;
		auto itr_rgb = colors.begin();
		for (auto itr_xyz = coordinates.begin(), end = coordinates.end(); itr_xyz != end; itr_xyz++, itr_rgb++)
		{
			if (distribution(random_numbers_engine) > points_to_remove_ratio)
			{
				reduced_coordinates.push_back(*itr_xyz);
				reduced_colors.push_back(*itr_rgb);
			}
			context.Feedback().Update(++progress / coordinates.size());
		}
		// add random points to the new cloud
		Clouds::PointsRange reduced_range;
		reduced_cloud->GetAccess().AllocPoints(reduced_coordinates.size(), &reduced_range);
		reduced_range.SetXYZ(reduced_coordinates);
		reduced_range.SetColors(reduced_colors);
		node_id = reduced_node->GetID();
		return reduced_cloud;
	}

	/*Function finds neighbours for all the points from the cloud in the given node and makes the lists of graph elements assigned to the vector in this structure.*/
	void BuildPointsGraphStructure(Context& context, ResourceID child_node_id)
	{
		// get data from a cloud
		auto cloud = context.m_project->TransTreeFindNode(child_node_id)->GetElement()->GetData<Clouds::ICloud>();
		Clouds::PointsRange range;
		cloud->GetAccess().GetAllPoints(range);
		std::vector<Clouds::Point3D> coordinates;
		std::vector<Clouds::Color> colors;
		range.GetXYZ(coordinates);
		range.GetColors(colors);
		// assign points parameters to the structure
		points_graph_structure.reserve(range.size());
		for (int i = 0; i < range.size(); i++)
			points_graph_structure.push_back(Point(coordinates[i], colors[i]));
		// sort all the points by 'x' coordinate to reduce the time of finding neighbours
		std::vector<Point*> sorted_by_x_ptrs;
		sorted_by_x_ptrs.reserve(range.size());
		for (auto itr = points_graph_structure.begin(), end = points_graph_structure.end(); itr != end; ++itr)
			sorted_by_x_ptrs.push_back(&*itr);
		std::sort(sorted_by_x_ptrs.begin(), sorted_by_x_ptrs.end(), [](Point* a, Point* b) {
			return a->xyz.x() < b->xyz.x();
			});
		// make neighbours lists
		float progress = 0.0;
		for (auto itr = sorted_by_x_ptrs.begin(), end = sorted_by_x_ptrs.end(); itr != end; ++itr)
		{
			auto first_point_xyz = (*itr)->xyz.cast<double>();
			float x_max = max_distance + (*itr)->xyz.x();
			// iterate through the rest of neighbourhood
			for (auto n_itr = itr + 1; n_itr != end; ++n_itr)
			{
				if ((*n_itr)->xyz.x() > x_max)
					break;
				if (Math::CalcPointToPointDistance3D(first_point_xyz, (*n_itr)->xyz.cast<double>()) <= max_distance)
				{
					// add neighbourhood between these points
					(*itr)->neighbours.push_back(*n_itr);
					(*n_itr)->neighbours.push_back(*itr);
				}
			}
			// sort neighbours by distance
			std::sort((*itr)->neighbours.begin(), (*itr)->neighbours.end(), [itr](Point* a, Point* b)
				{
					return Math::CalcPointToPointDistance3D(a->xyz.cast<double>(), (*itr)->xyz.cast<double>())
						< Math::CalcPointToPointDistance3D(b->xyz.cast<double>(), (*itr)->xyz.cast<double>());
				});
			context.Feedback().Update(++progress / points_graph_structure.size());
		}
	}
	
	/*Function copies all points' data from the graph structure to the given range of the cloud.*/
	void CopyGraphStructuctureToCloud(Clouds::PointsRange& range)
	{
		std::vector<Clouds::Point3D> coordinates;
		std::vector<Clouds::Color> colors;
		coordinates.reserve(points_graph_structure.size());
		colors.reserve(points_graph_structure.size());
		for (auto itr = points_graph_structure.begin(), end = points_graph_structure.end(); itr != end; ++itr)
		{
			coordinates.push_back(itr->xyz);
			colors.push_back(itr->color);
		}
		range.SetXYZ(coordinates);
		range.SetColors(colors);
	}

	/*Function marks lonely points as 'deleted' and smooths the cloud by projecting all points in the structure to a local fitted sphere or plane.*/
	void RemoveNoise(Context& context)
	{
		float progress = 0.0;
		for (auto itr = points_graph_structure.begin(), end = points_graph_structure.end(); itr != end; ++itr)
		{
			// remove lonely points
			if (itr->neighbours.size() < min_neighbours)
				itr->is_deleted = 1;
		}
		// calculate new coordinates to smooth the cloud
		std::vector<Clouds::Point3D> new_coordinates;
		new_coordinates.reserve(points_graph_structure.size());
		for (auto itr = points_graph_structure.begin(), end = points_graph_structure.end(); itr != end; ++itr)
		{
			// check if there are enough neighbours to precise fitting
			if (itr->neighbours.size() < neighbours_to_smooth || itr->is_deleted)
			{
				new_coordinates.push_back(itr->xyz);
				continue;
			}
			// calculate best fitted sphere and plane
			std::vector<Clouds::Point3D> n_xyz = { itr->xyz };
			for (auto n_itr = itr->neighbours.begin(), n_end = n_itr + neighbours_to_smooth; n_itr != n_end; ++n_itr)
			{
				n_xyz.push_back((*n_itr)->xyz);
			}
			auto sphere = Math::CalcBestSphere3D(n_xyz.begin(), n_xyz.end());
			auto plane = Math::CalcBestPlane3D(n_xyz.begin(), n_xyz.end());
			// calculate RMSE of fitting
			double sq_sum_sphere = 0.0;
			double sq_sum_plane = 0.0;
			for (auto n_itr = n_xyz.begin(), n_end = n_xyz.end(); n_itr != n_end; ++n_itr)
			{
				auto point_onto_sphere = sphere.Project(n_itr->cast<double>());
				auto point_onto_plane = Math::ProjectPointOntoPlane(plane, n_itr->cast<double>());
				sq_sum_sphere += pow(Math::CalcPointToPointDistance3D(point_onto_sphere.cast<double>(), n_itr->cast<double>()), 2);
				sq_sum_plane += pow(Math::CalcPointToPointDistance3D(point_onto_plane.cast<double>(), n_itr->cast<double>()), 2);
			}
			// choose better projection
			Clouds::Point3D new_xyz;
			if (sq_sum_sphere < sq_sum_plane)
			{
				new_xyz = Clouds::Point3D(
					sphere.Project(itr->xyz.cast<double>()).x(),
					sphere.Project(itr->xyz.cast<double>()).y(),
					sphere.Project(itr->xyz.cast<double>()).z()
				);
			}
			else
			{
				new_xyz = Clouds::Point3D(
					Math::ProjectPointOntoPlane(plane, itr->xyz.cast<double>()).x(),
					Math::ProjectPointOntoPlane(plane, itr->xyz.cast<double>()).y(),
					Math::ProjectPointOntoPlane(plane, itr->xyz.cast<double>()).z()
				);
			}
			new_coordinates.push_back(new_xyz);
			context.Feedback().Update(++progress / points_graph_structure.size());
		}
		// assign new coordinates to the structure
		for (int i = 0; i < points_graph_structure.size(); i++)
		{
			points_graph_structure[i].xyz = new_coordinates[i];
		}
		// sort neighbours by distance again
		for (auto itr = points_graph_structure.begin(), end = points_graph_structure.end(); itr != end; ++itr)
		{
			std::sort(itr->neighbours.begin(), itr->neighbours.end(), [itr](Point* a, Point* b)
				{
					return Math::CalcPointToPointDistance3D(a->xyz.cast<double>(), itr->xyz.cast<double>())
						< Math::CalcPointToPointDistance3D(b->xyz.cast<double>(), itr->xyz.cast<double>());
				});
		}
	}
	
	/*Function filters all points using the median colors of neighbours.*/
	void MedianFilterColors(Context& context)
	{
		float progress = 0.0;
		std::vector<Clouds::Color> colors;
		colors.reserve(points_graph_structure.size());
		for (auto itr = points_graph_structure.begin(), end = points_graph_structure.end(); itr != end; ++itr)
		{
			// set the previous color value when the point is deleted or it has not enough neighbours to precise filtration
			if (itr->is_deleted || itr->neighbours.size() < median_filter_size)
			{
				colors.push_back(itr->color);
				context.Feedback().Update(++progress / points_graph_structure.size());
				continue;
			}
			// get color values and the weight of each neighbour
			std::vector<std::pair<double, double>> red, green, blue;
			double current_weight = 0.0;
			for (auto n_itr = itr->neighbours.begin(), n_end = n_itr + median_filter_size; n_itr != n_end; ++n_itr)
			{
				double weight = pow(Math::CalcPointToPointDistance3D(itr->xyz.cast<double>(), (*n_itr)->xyz.cast<double>()), -2.0);
				red.push_back(std::pair<double, double>((*n_itr)->color.x(), weight));
				green.push_back(std::pair<double, double>((*n_itr)->color.y(), weight));
				blue.push_back(std::pair<double, double>((*n_itr)->color.z(), weight));
				current_weight += weight;
			}
			const double HALF_WEIGHT = current_weight / 2;
			// calculate median value of color for each channel
			std::vector<double> new_color;
			std::vector<std::vector<std::pair<double, double>>*> rgb = { &red,&green,&blue };
			std::for_each(rgb.begin(), rgb.end(), [&new_color, HALF_WEIGHT](std::vector<std::pair<double, double>>* channel) {
				std::sort(channel->begin(), channel->end(), [](std::pair<double, double> left, std::pair<double, double> right) {
					return left.first > right.first;
					});
				double current_weight = 0.0;
				double channel_value = 0.0;
				for (auto itr = channel->begin(), end = channel->end(); itr != end; ++itr)
				{
					current_weight += itr->second;
					if (current_weight > HALF_WEIGHT)
					{
						channel_value = itr->first;
						break;
					}
				}
				new_color.push_back(channel_value);
				});
			colors.push_back(Clouds::Color(new_color[0], new_color[1], new_color[2], 0));
			context.Feedback().Update(++progress / points_graph_structure.size());
		}
		// write the new color values to the structure
		for (int i = 0; i < points_graph_structure.size(); i++)
		{
			points_graph_structure[i].color = colors[i];
		}
	}

	/*Function groups the points into clusters choosing one of the implemented recursive Hausdorf algorithm depending on the number of points.*/
	void Hausdorf(Context& context)
	{
		float progress = 0.0;
		// set initial values
		for (auto itr = points_graph_structure.begin(), end = points_graph_structure.end(); itr != end; ++itr)
		{
			itr->is_calculated = 0;
			itr->group_number = 0;
		}
		int group_counter = 0;
		// realize Hausdorf algorithm for each point
		for (auto itr = points_graph_structure.begin(), end = points_graph_structure.end(); itr != end; ++itr)
		{
			if (itr->group_number != 0 || itr->is_deleted)
				continue;
			group_counter++;
			// check if recursion with this number of points results stack overflow
			if (points_graph_structure.size() < HAUSDORF_THRESHOLD)
				RecursiveHausdorf(context, *itr, group_counter, progress);
			else
				HalfRecursiveHausdorf(context, *itr, group_counter, progress);
		}
	}
	
	/*Function groups points into a cluster with given number, starting from a given point using recursive Hausdorf's algorithm adapted to the graph structure of data. It should not be used for large clouds due to stack overflow.*/
	void RecursiveHausdorf(Context& context, Point& start_point, const int& group_counter, float& progress)
	{
		// relize the algorithm for the point
		start_point.group_number = group_counter;
		for (auto n_itr = start_point.neighbours.begin(), n_end = start_point.neighbours.end(); n_itr != n_end; ++n_itr)
		{
			if ((*n_itr)->is_deleted || (*n_itr)->is_calculated)
				continue;

			double distance = Math::CalcPointToPointDistance3D(start_point.xyz.cast<double>(), (*n_itr)->xyz.cast<double>());
			if (CheckColorGradient(start_point.color, (*n_itr)->color, distance))
				(*n_itr)->group_number = group_counter;
		}
		start_point.is_calculated = 1;
		context.Feedback().Update(++progress / points_graph_structure.size());
		// repeate for the neighbours recursively
		for (auto n_itr = start_point.neighbours.begin(), n_end = start_point.neighbours.end(); n_itr != n_end; ++n_itr)
		{
			if ((*n_itr)->is_calculated == 0 && (*n_itr)->is_deleted == 0 && (*n_itr)->group_number == start_point.group_number)
				RecursiveHausdorf(context, **n_itr, group_counter, progress);
		}
	}
	
	/*Function groups points into a cluster with given number, starting from a given point using recursive Hausdorf's algorithm with 2 steps in one iteration of the function. This modification is less clear but prevents stack overflow.*/
	void HalfRecursiveHausdorf(Context& context, Point& start_point, const int& group_counter, float& progress)
	{
		// relize the algorithm for the point
		start_point.group_number = group_counter;
		for (auto n_itr = start_point.neighbours.begin(), n_end = start_point.neighbours.end(); n_itr != n_end; ++n_itr)
		{
			if ((*n_itr)->is_deleted || (*n_itr)->is_calculated)
				continue;

			double distance = Math::CalcPointToPointDistance3D(start_point.xyz.cast<double>(), (*n_itr)->xyz.cast<double>());
			if (CheckColorGradient(start_point.color, (*n_itr)->color, distance))
				(*n_itr)->group_number = group_counter;
		}
		start_point.is_calculated = 1;
		context.Feedback().Update(++progress / points_graph_structure.size());
		// repeate for the neighbours in a loop
		std::vector<Point*> points;
		for (auto itr = start_point.neighbours.begin(), end = start_point.neighbours.end(); itr != end; ++itr)
		{
			if ((*itr)->is_calculated || (*itr)->is_deleted || (*itr)->group_number != start_point.group_number)
				continue;

			for (auto n_itr = (*itr)->neighbours.begin(), n_end = (*itr)->neighbours.end(); n_itr != n_end; ++n_itr)
			{
				if ((*n_itr)->is_deleted || (*n_itr)->is_calculated)
					continue;

				double distance = Math::CalcPointToPointDistance3D(start_point.xyz.cast<double>(), (*n_itr)->xyz.cast<double>());
				if (CheckColorGradient((*itr)->color, (*n_itr)->color, distance))
					(*n_itr)->group_number = group_counter;
			}
			(*itr)->is_calculated = 1;
			points.push_back(*itr);
			context.Feedback().Update(++progress / points_graph_structure.size());
		}
		// repeate for the neighbours of neighbours recursively
		for (auto itr = points.begin(), end = points.end(); itr != end; ++itr)
		{
			for (auto n_itr = (*itr)->neighbours.begin(), n_end = (*itr)->neighbours.end(); n_itr != n_end; ++n_itr)
			{
				if ((*n_itr)->is_calculated == 0 && (*n_itr)->is_deleted == 0 && (*n_itr)->group_number == (*itr)->group_number)
					HalfRecursiveHausdorf(context, **n_itr, group_counter, progress);
			}
		}
	}
	
	/*Function returns 'true' if the color gradient of given data converted to HSV are within tolerance and 'false' otherwise.*/	
	bool CheckColorGradient(Clouds::Color color_1, Clouds::Color color_2, double distance)
	{
		// get RGB colors
		const double R1 = color_1.x();
		const double G1 = color_1.y();
		const double B1 = color_1.z();
		const double R2 = color_2.x();
		const double G2 = color_2.y();
		const double B2 = color_2.z();

		// convert to HSV
		const double V1 = std::max({ R1, G1, B1 }) / 255.0;
		const double S1 = V1 == 0.0 ? 0.0 : 1.0 - std::min({ R1, G1, B1 }) / 255.0 / V1;
		const double H1 = G1 >= B1 ? 180.0 / PI * acos((R1 - (G1 + B1) / 2.0) / sqrt(R1 * R1 + G1 * G1 + B1 * B1 - R1 * G1 - R1 * B1 - G1 * B1)) :
			360.0 - 180.0 / PI * acos((R1 - G1 / 2.0 - B1 / 2.0) / sqrt(R1 * R1 + G1 * G1 + B1 * B1 - R1 * G1 - R1 * B1 - G1 * B1));
		const double V2 = std::max({ R2, G2, B2 }) / 255;
		const double S2 = V2 == 0.0 ? 0.0 : 1.0 - std::min({ R2, G2, B2 }) / 255.0 / V2;
		const double H2 = G2 >= B2 ? 180.0 / PI * acos((R2 - (G2 + B2) / 2.0) / sqrt(R2 * R2 + G2 * G2 + B2 * B2 - R2 * G2 - R2 * B2 - G2 * B2)) :
			360.0 - 180.0 / PI * acos((R2 - G2 / 2.0 - B2 / 2.0) / sqrt(R2 * R2 + G2 * G2 + B2 * B2 - R2 * G2 - R2 * B2 - G2 * B2));

		// calculate gradient of each factor and compare with sum of maximum value set by a user and measurement inaccuracy
		if (V1 == 0.0 && V2 == 0.0)
			return 1;
		if (abs(V1 - V2) > max_brightness_gradient * distance)
			return 0;
		if (S1 == 0.0 && S2 == 0.0)
			return 1;
		if (abs(S1 - S2) > max_saturation_gradient * distance + coverage_rate / (255.0 * (V1 + V2) / 2))
			return 0;
		const double HUE_INACCURACY = coverage_rate * (180.0 / (255.0 * (V1 + V2) / 2) + 180.0 / (255.0 * (S1 + S2) / 2));
		if (std::min(abs(H1 - H2), 360.0 - abs(H1 - H2)) > max_hue_gradient * distance + HUE_INACCURACY)
			return 0;
		return 1;
	}

	/*Function finds groups with a number of elements less than minimum given by a user and marks them as ungrouped.*/
	void RemoveTooSmallGroups(Context& context)
	{
		float progress = 0.0;
		// calculate repetition numbers
		std::map<int, int> repetition_numbers;
		for (auto itr = points_graph_structure.begin(), end = points_graph_structure.end(); itr != end; ++itr)
		{
			repetition_numbers[itr->group_number]++;
			context.Feedback().Update(++progress / 2 / points_graph_structure.size());
		}
		// assign number '0' to ungrouped points
		for (auto itr = points_graph_structure.begin(), end = points_graph_structure.end(); itr != end; ++itr)
		{
			if (repetition_numbers[itr->group_number] < min_points_group)
				itr->group_number = 0;
			context.Feedback().Update(++progress / 2 / points_graph_structure.size());
		}
	}
	
	/*Function adds ungrouped points surrounded by neighbours assigned to groups to the locally most popular group.*/
	void FillHolesInClusters(Context& context)
	{
		float progress = 0.0;
		// set initial states
		for (auto itr = points_graph_structure.begin(), end = points_graph_structure.end(); itr != end; ++itr)
		{
			if (itr->group_number != 0 || itr->is_deleted)
			{
				itr->is_calculated = 1;
				context.Feedback().Update(++progress / points_graph_structure.size());
			}
			else
				itr->is_calculated = 0;
		}
		bool any_change = 0;
		do
		{
			any_change = 0;
			std::vector<int> new_group_numbers;
			new_group_numbers.reserve(points_graph_structure.size());
			for (auto itr = points_graph_structure.begin(), end = points_graph_structure.end(); itr != end; ++itr)
			{
				if (itr->is_calculated)
				{
					new_group_numbers.push_back(0);
					continue;
				}
				// calculate repetition number of each group number in the neighbourhood of the point
				std::map<int, int> repetition_number;
				for (auto n_itr = itr->neighbours.begin(), n_end = itr->neighbours.end(); n_itr != n_end; ++n_itr)
				{
					if ((*n_itr)->group_number != 0)
						repetition_number[(*n_itr)->group_number]++;
				}
				// find most popular group in the neighbourhood
				auto max_group = std::max_element(repetition_number.begin(), repetition_number.end(), [](auto left, auto right) {
					return left.second < right.second;
					});
				// assign point to this group if it exist
				if (max_group != repetition_number.end() && repetition_number[0] < join_coefficient * itr->neighbours.size())
				{
					new_group_numbers.push_back(max_group->first);
					itr->is_calculated = 1;
					any_change = 1;
					context.Feedback().Update(++progress / points_graph_structure.size());
				}
				else
					new_group_numbers.push_back(0);
			}
			// write values to the structure at the end of each loop to guarantee isotropic addition of points
			for (int i = 0; i < points_graph_structure.size(); i++)
			{
				if (new_group_numbers[i] != 0)
					points_graph_structure[i].group_number = new_group_numbers[i];
			}
			// end the algorithm if there is no more points to change their group
		} while (any_change);
	}

	/*Function merge clusters with a relatively high number of points on common border.*/
	void MergeClusters(Context& context)
	{
		float progress = 0.0;
		std::vector<std::pair<int, std::vector<Point*>>> clusters;
		// make a list of clusters 
		for (auto itr = points_graph_structure.begin(), end = points_graph_structure.end(); itr != end; ++itr)
		{
			if (itr->group_number == 0)
				continue;

			auto cluster = std::find_if(clusters.begin(), clusters.end(), [itr](auto cluster) {
				return cluster.first == itr->group_number;
				});
			if (cluster == clusters.end())
				clusters.push_back(std::pair<int, std::vector<Point*>>(itr->group_number, std::vector<Point*>({ &*itr })));
			else
				cluster->second.push_back(&*itr);
		}
		// sort the clusters by number of points
		std::sort(clusters.begin(), clusters.end(), [](std::pair<int, std::vector<Point*>> a, std::pair<int, std::vector<Point*>> b)
			{
				return a.second.size() < b.second.size();
			});
		// get a smaller cluster
		for (int i = 0; i < clusters.size(); i++)
		{
			int max_neighbours = 0;
			// get a bigger cluster
			for (int j = i + 1; j < clusters.size(); j++)
			{
				const int BIGGER_GROUP = clusters[j].first;
				// find number of point on the border beetwen them
				int border_points = std::count_if(clusters[i].second.begin(), clusters[i].second.end(), [BIGGER_GROUP](auto point) {
					return std::any_of(point->neighbours.begin(), point->neighbours.end(), [BIGGER_GROUP](auto neigbour) {
						return neigbour->group_number == BIGGER_GROUP;
						});
					});
				// add sizes of all merged part into the smaller cluster
				int size = clusters[i].second.size();
				for (int k = 0; k < i; k++)
				{
					if (clusters[k].first == clusters[i].first)
						size += clusters[k].second.size();
				}
				// compare border size with number of points in the smaller cluster circumference
				if (border_points < sqrt(4 * PI * size) * max_distance / merge_coefficient || border_points < max_neighbours)
					continue;

				// merge the smaller cluster and all clusters merged into them into the bigger cluster
				for (auto itr = clusters[i].second.begin(), end = clusters[i].second.end(); itr != end; ++itr)
				{
					(*itr)->group_number = BIGGER_GROUP;
				}
				for (int k = 0; k < i; k++)
				{
					if (clusters[k].first != clusters[i].first)
						continue;

					for (auto itr = clusters[k].second.begin(), end = clusters[k].second.end(); itr != end; ++itr)
					{
						(*itr)->group_number = BIGGER_GROUP;
					}
					clusters[k].first = BIGGER_GROUP;
				}
				clusters[i].first = BIGGER_GROUP;
				max_neighbours = border_points;
			}
			context.Feedback().Update(++progress / clusters.size());
		}
	}
	
	/*Function changes group numbers to subsequent numbers with '0' assigned to ungrouped points and write them to a given vector of layer values.*/
	void TidyUpGroupNumbers(Context& context, std::vector<StoredReal>& cluster_numbers)
	{
		float progress = 0.0;
		std::map<int, int> repetition_numebers;
		// get the size of each cluster
		for (auto itr = points_graph_structure.begin(), end = points_graph_structure.end(); itr != end; ++itr)
		{
			// ignore ungrouped points
			if (itr->group_number != 0)
				repetition_numebers[itr->group_number]++;
			context.Feedback().Update(++progress / 2 / points_graph_structure.size());
		}
		// sort clusters by number of points
		std::vector<std::map<int, int>::iterator> iters;
		for (auto itr = repetition_numebers.begin(), end = repetition_numebers.end(); itr != end; ++itr)
			iters.push_back(itr);
		std::sort(iters.begin(), iters.end(), [](std::map<int, int>::iterator a, std::map<int, int>::iterator b)
			{
				return a->second > b->second;
			});
		// assign a sequential number to each cluster from the largest to the smallest
		for (auto itr = points_graph_structure.begin(), end = points_graph_structure.end(); itr != end; ++itr)
		{
			for (int j = 0; j < iters.size(); j++)
			{
				if ((iters[j])->first == itr->group_number)
				{
					// ignore number '0' assigned to ungrouped points
					itr->group_number = j + 1;
					break;
				}
			}
			cluster_numbers.push_back(itr->group_number);
			context.Feedback().Update(++progress / 2 / points_graph_structure.size());
		}
	}

	/*Function finds the best fitting cluster to the loudspeaker and write to the given layer '1.0' for points of this cluster and '0.0' for others.*/
	void FindTheSpeaker(Context& context, std::vector<StoredReal>& is_speaker)
	{
		float progress = 0.0;
		std::map<int, std::vector<Point*>> clusters;
		// make a list of clusters 
		for (auto itr = points_graph_structure.begin(), end = points_graph_structure.end(); itr != end; ++itr)
		{
			if (clusters.find(itr->group_number) == clusters.end())
				clusters[itr->group_number] = { &*itr };
			else
				clusters[itr->group_number].push_back(&*itr);
		}
		// find the cluster with the least fitting error
		std::map<int, double> errors;
		for (auto itr = clusters.begin(), end = clusters.end(); itr != end; ++itr)
		{
			errors[itr->first] = GetBoxFittingError(context, itr->second);
			// linearize operation progress
			++progress;
			context.Feedback().Update((2.0 - progress / clusters.size()) * progress / clusters.size());
		}
		auto speaker_itr = std::min_element(errors.begin(), errors.end(), [](auto left, auto right) {
			return left.second < right.second;
			});
		for (auto itr = points_graph_structure.begin(), end = points_graph_structure.end(); itr != end; ++itr)
		{
			if (itr->group_number == speaker_itr->first)
				is_speaker.push_back(1.0);
			else
				is_speaker.push_back(0.0);
		}
		//// print the values of all fitting errors
		//for (auto itr = errors.begin(), end = errors.end(); itr != end; ++itr)
		//{
		//	OGX_LINE.Msg(Info, std::to_wstring(itr->first) + L". " + std::to_wstring(itr->second));
		//}
	}

	/*Function calculates a box fitted to the borders of given points cluster. Returns sum of angles and dimensions relative errors and relative RMSE for this fitting.*/
	double GetBoxFittingError(Context& context, std::vector<Point*>& cluster)
	{
		// find the longest diagonal
		double max_dist = 0.0, local_dist = 0.0;
		Clouds::Point3D point_A, point_B, point_C, point_D;
		for (int i = 0; i < cluster.size(); i++)
		{
			for (int j = i + 1; j < cluster.size(); j++)
			{
				local_dist = Math::CalcPointToPointDistance3D(cluster[i]->xyz.cast<double>(), cluster[j]->xyz.cast<double>());
				if (local_dist > max_dist)
				{
					max_dist = local_dist;
					point_A = cluster[i]->xyz;
					point_B = cluster[j]->xyz;
				}
			}
		}
		const Clouds::Point3D CENTER_POINT = (point_A + point_B) / 2;
		const auto MAIN_DIAGONAL = Math::CalcLine3D(point_A.cast<double>(), point_B.cast<double>());
		// find the farthest vertex from the diagonal
		max_dist = 0.0;
		for (int k = 0; k < cluster.size(); k++)
		{
			local_dist = Math::CalcPointToLineDistance3D(cluster[k]->xyz.cast<double>(), MAIN_DIAGONAL);
			if (local_dist > max_dist)
			{
				max_dist = local_dist;
				point_C = cluster[k]->xyz;
			}
		}
		const auto SECTION_PLANE = Math::CalcPlane3D(point_A.cast<double>(), point_B.cast<double>(), point_C.cast<double>());
		// find the farthest vertex from the cross-section plane including main diagonal
		max_dist = 0.0;
		for (int l = 0; l < cluster.size(); l++)
		{
			local_dist = Math::CalcPointToPlaneDistance3D(cluster[l]->xyz.cast<double>(), SECTION_PLANE.cast<double>());
			local_dist += Math::CalcPointToPointDistance3D(cluster[l]->xyz.cast<double>(), CENTER_POINT.cast<double>());
			if (local_dist > max_dist)
			{
				max_dist = local_dist;
				point_D = cluster[l]->xyz;
			}
		}
		// calculate all vertex of the fitted box
		Clouds::Point3D point_000, point_001, point_010, point_100, point_011, point_101, point_110, point_111;
		point_000 = point_A;
		point_111 = point_B;
		if (Math::CalcPointToPointDistance3D(point_000.cast<double>(), point_D.cast<double>()) <
			Math::CalcPointToPointDistance3D(point_111.cast<double>(), point_D.cast<double>()))
		{
			// D on plane 0XX
			if (Math::CalcPointToPointDistance3D(point_C.cast<double>(), point_D.cast<double>()) <
				Math::CalcPointToPointDistance3D(point_C.cast<double>(), point_000.cast<double>()))
			{
				point_011 = point_C;
				point_100 = 2 * CENTER_POINT - point_C;
			}
			else
			{
				point_100 = point_C;
				point_011 = 2 * CENTER_POINT - point_C;
			}
			if (Math::CalcPointToPointDistance3D(point_D.cast<double>(), point_000.cast<double>()) <
				Math::CalcPointToPointDistance3D(point_D.cast<double>(), point_011.cast<double>()))
			{
				point_001 = point_D;
				point_010 = point_000 + point_011 - point_D;
			}
			else
			{
				point_010 = point_D;
				point_001 = point_000 + point_011 - point_D;
			}
			point_101 = 2 * CENTER_POINT - point_010;
			point_110 = 2 * CENTER_POINT - point_001;
		}
		else
		{
			// D on plane 1XX
			if (Math::CalcPointToPointDistance3D(point_C.cast<double>(), point_D.cast<double>()) <
				Math::CalcPointToPointDistance3D(point_C.cast<double>(), point_111.cast<double>()))
			{
				point_100 = point_C;
				point_011 = 2 * CENTER_POINT - point_C;
			}
			else
			{
				point_011 = point_C;
				point_100 = 2 * CENTER_POINT - point_C;
			}
			if (Math::CalcPointToPointDistance3D(point_D.cast<double>(), point_111.cast<double>()) <
				Math::CalcPointToPointDistance3D(point_D.cast<double>(), point_100.cast<double>()))
			{
				point_110 = point_D;
				point_101 = point_111 + point_100 - point_D;
			}
			else
			{
				point_101 = point_D;
				point_110 = point_111 + point_100 - point_D;
			}
			point_010 = 2 * CENTER_POINT - point_101;
			point_001 = 2 * CENTER_POINT - point_110;
		}
		// calculate all walls of the box
		const auto PLANE_0XX = Math::CalcPlane3D(point_000.cast<double>(), point_010.cast<double>(), point_001.cast<double>());
		const auto PLANE_1XX = Math::CalcPlane3D(point_111.cast<double>(), point_110.cast<double>(), point_101.cast<double>());
		const auto PLANE_X0X = Math::CalcPlane3D(point_000.cast<double>(), point_001.cast<double>(), point_100.cast<double>());
		const auto PLANE_X1X = Math::CalcPlane3D(point_111.cast<double>(), point_011.cast<double>(), point_110.cast<double>());
		const auto PLANE_XX0 = Math::CalcPlane3D(point_000.cast<double>(), point_010.cast<double>(), point_100.cast<double>());
		const auto PLANE_XX1 = Math::CalcPlane3D(point_111.cast<double>(), point_011.cast<double>(), point_101.cast<double>());
		// calculate RMS error of the box fitting
		double RMSE = 0.0;
		for (int i = 0; i < cluster.size(); i++)
		{
			local_dist = std::min({
				Math::CalcPointToPlaneDistance3D(cluster[i]->xyz.cast<double>(), PLANE_0XX.cast<double>()),
				Math::CalcPointToPlaneDistance3D(cluster[i]->xyz.cast<double>(), PLANE_1XX.cast<double>()),
				Math::CalcPointToPlaneDistance3D(cluster[i]->xyz.cast<double>(), PLANE_X0X.cast<double>()),
				Math::CalcPointToPlaneDistance3D(cluster[i]->xyz.cast<double>(), PLANE_X1X.cast<double>()),
				Math::CalcPointToPlaneDistance3D(cluster[i]->xyz.cast<double>(), PLANE_XX0.cast<double>()),
				Math::CalcPointToPlaneDistance3D(cluster[i]->xyz.cast<double>(), PLANE_XX1.cast<double>())
				});
			RMSE += pow(local_dist, 2);
		}
		RMSE = sqrt(RMSE / cluster.size());
		max_dist = Math::CalcPointToLineDistance3D(point_C.cast<double>(), MAIN_DIAGONAL);
		// calculate angle errors
		double relative_angle_error = 0.0;
		const auto X_VEC = Math::Vector3D(point_100.cast<double>() - point_000.cast<double>());
		const auto Y_VEC = Math::Vector3D(point_010.cast<double>() - point_000.cast<double>());
		const auto Z_VEC = Math::Vector3D(point_001.cast<double>() - point_000.cast<double>());
		relative_angle_error += abs(Math::CalcAngleBetweenTwoVectors(X_VEC.cast<double>(), Y_VEC.cast<double>()) - PI / 2) / (PI / 2);
		relative_angle_error += abs(Math::CalcAngleBetweenTwoVectors(Y_VEC.cast<double>(), Z_VEC.cast<double>()) - PI / 2) / (PI / 2);
		relative_angle_error += abs(Math::CalcAngleBetweenTwoVectors(Z_VEC.cast<double>(), X_VEC.cast<double>()) - PI / 2) / (PI / 2);
		// calculate dimension errors
		double relative_demention_error = 0.0;
		relative_demention_error += abs(X_VEC.norm() - speaker_dim_x) / speaker_dim_x;
		relative_demention_error += abs(Y_VEC.norm() - speaker_dim_y) / speaker_dim_y;
		relative_demention_error += abs(Z_VEC.norm() - speaker_dim_z) / speaker_dim_z;
		// return combined relative error
		return (RMSE / max_dist + relative_angle_error + relative_demention_error) / 7.0;
	}
};

OGX_EXPORT_METHOD(LoudspeakerSearching)
