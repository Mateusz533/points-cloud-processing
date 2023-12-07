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

struct LoudspeakerSearching : public Plugin::EasyMethod {
	/*HARDCODED PARAMETERS AND DATA*/
	const double HAUSDORF_THRESHOLD = 50000;	// number of points which change version of Hausdorf algorithm

	/*Structure store data of a 3D point with vector of pointers to its neighbours, which make it an element of a graph structure.*/
	struct Point
	{
		Clouds::Point3D xyz;
		Clouds::Color color;
		std::vector<Point*> neighbours;
		int group_number = 0;
		bool is_calculated = false;
		bool is_deleted = false;

		Point(Clouds::Point3D p_coordinates, Clouds::Color p_color) {
			xyz = p_coordinates;
			color = p_color;
		}
		bool isReady() const {
			return is_calculated || is_deleted;
		}
	};

	enum class PointState : uint8_t {
		DELETED,
		CALCULATED,
		NOT_CALCULATED,
	};

	class FeedbackUpdater {
	public:
		FeedbackUpdater() = delete;
		FeedbackUpdater(Context& context, unsigned int step_number, bool decreasing_growth = false) :
			context(context),
			step_number(step_number),
			decreasing_growth(decreasing_growth)
		{}

		void count() {
			const double last_ratio = 1.0 * counter++ / step_number;
			const double ratio = 1.0 * counter / step_number;
			const double percentage = decreasing_growth ? 100 * (1.0 - pow(1.0 - ratio, 10.0)) : 100 * ratio;
			const double last_percentage = decreasing_growth ? 100 * (1.0 - pow(1.0 - last_ratio, 10.0)) : 100 * last_ratio;
			if (floor(percentage) > floor(last_percentage))
				context.Feedback().Update(0.01 * percentage);
		}
	private:
		Context& context;
		unsigned int counter = 0;
		const unsigned int step_number;
		const bool decreasing_growth;
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
	virtual void DefineParameters(ParameterBank& bank) 	{
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

	virtual void Run(Context& context) {
		// get cloud data
		Clouds::ICloud* original_cloud = GetValidCloud(context, m_node_id);
		OGX_LINE.Msg(Info, L"All data is correct. Processing has started.");

		// reduce the number of points
		ResourceID new_node_id;
		Clouds::ICloud* cloud = CreateRandomlyReducedCloud(context, original_cloud, new_node_id, points_to_remove_ratio);
		OGX_LINE.Msg(Info, L"Reduced cloud has been done.");

		// create points structure from created cloud
		BuildPointsGraphStructure(context, new_node_id, points_graph_structure, max_distance);
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
		RemoveNoise(context, points_graph_structure, min_neighbours, neighbours_to_smooth);

		// reduce noise in colors using median filter
		MedianFilterColors(context, points_graph_structure, median_filter_size);
		OGX_LINE.Msg(Info, L"Noise filtration has been done.");

		// visualize the cloud smoothing
		CopyGraphStructuctureToCloud(points_graph_structure, range);

		// make segmentation using Hausdorf algorithm
		Hausdorf(context, points_graph_structure, HAUSDORF_THRESHOLD);
		OGX_LINE.Msg(Info, L"Segmentation by Hausdorf algorithm has been done.");

		// fill holes with ungrouped points in clusters
		FillHolesInClusters(context, points_graph_structure);

		// remove too small groups
		RemoveTooSmallGroups(context, points_graph_structure);
		OGX_LINE.Msg(Info, L"Removing too small groups has been done.");

		// fill holes made by removing small groups
		FillHolesInClusters(context, points_graph_structure);
		OGX_LINE.Msg(Info, L"Filling holes in clusters has been done.");

		// merge clusters
		MergeClusters(context, points_graph_structure);
		OGX_LINE.Msg(Info, L"Merging clusters has been done.");

		// change group numbers to subsequent numbers with '0' assigned to ungrouped points
		TidyUpGroupNumbers(context, segment_numbers, points_graph_structure);
		OGX_LINE.Msg(Info, L"The ordering of group numbers has been done.");

		// assign all cluster numbers to created layer
		range.SetLayerVals(segment_numbers, *segments_layer);
		OGX_LINE.Msg(Info, L"Whole segmentation process has been finished.");

		// create next layer
		const String NEW_LAYER_NAME = L"Loudspeaker";
		const auto new_layers = cloud->FindLayers(NEW_LAYER_NAME);
		auto speaker_layer = new_layers.empty() ? cloud->CreateLayer(NEW_LAYER_NAME, 0.0) : new_layers[0];
		std::vector<StoredReal> is_speaker;
		is_speaker.reserve(points_graph_structure.size());

		// find the best fitted cluster to the speaker
		FindTheSpeaker(context, is_speaker, points_graph_structure);
		OGX_LINE.Msg(Info, L"Loudspeaker identification has been done.");

		// mark found cluster in created layer
		range.SetLayerVals(is_speaker, *speaker_layer);
		OGX_LINE.Msg(Info, L"Processing has been finished. No errors.");

		//// check time
		//auto start = std::chrono::high_resolution_clock::now();
		//auto stop = std::chrono::high_resolution_clock::now();
		//auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		//auto time = duration.count();
	}

	/*Function reads data of all points from given node, check their correctness and returns a cloud of these points.*/
	Clouds::ICloud* GetValidCloud(Context& context, ResourceID node_id)	{
		auto node = context.m_project->TransTreeFindNode(node_id);
		if (!node)
			ReportError(L"Invalid node!");
		
		auto element = node->GetElement();
		if (!element)
			ReportError(L"Invalid element!");
		
		auto cloud = element->GetData<Clouds::ICloud>();
		if (!cloud)
			ReportError(L"Invalid cloud!");
		
		if (points_to_remove_ratio < 0.0 || points_to_remove_ratio >= 1.0)
			ReportError(L"Invalid points to remove ratio!");
		
		if (max_distance <= 0.0)
			ReportError(L"Invalid maximum distance to neighbour!");
		
		if (min_neighbours < 0)
			ReportError(L"Invalid minimum neigbours number!");
		
		if (min_points_group < 0)
			ReportError(L"Invalid minimum points cluster!");
		
		if (neighbours_to_smooth < 5)
			ReportError(L"Invalid number of neighbours to smooth!");
		
		if (median_filter_size < 2)
			ReportError(L"Invalid median filter size!");
		
		if (max_brightness_gradient < 0.0)
			ReportError(L"Invalid maximum brightness gradient!");
		
		if (max_saturation_gradient < 0.0)
			ReportError(L"Invalid maximum saturation gradient!");
		
		if (max_hue_gradient < 0.0)
			ReportError(L"Invalid maximum hue gradient!");
		
		if (coverage_rate < 0.0)
			ReportError(L"Invalid coverage rate!");
		
		if (join_coefficient < 0.0 || join_coefficient > 1.0)
			ReportError(L"Invalid join coefficient!");
		
		if (merge_coefficient < 0.0)
			ReportError(L"Invalid merge coefficient!");
		
		if (speaker_dim_x <= 0.0 || speaker_dim_x < speaker_dim_y || speaker_dim_x < speaker_dim_z)
			ReportError(L"Invalid speaker dimention X!");
		
		if (speaker_dim_y <= 0.0 || speaker_dim_y < speaker_dim_z)
			ReportError(L"Invalid speaker dimention Y!");
		
		if (speaker_dim_z <= 0.0)
			ReportError(L"Invalid speaker dimention Z!");
		
		return cloud;
	}

	/*Function makes a new node with given ID, assign to them the cloud with random points from a given cloud with ratio given by a user.*/
	Clouds::ICloud* CreateRandomlyReducedCloud(Context& context, Clouds::ICloud*& cloud, ResourceID& node_id, float points_to_remove_ratio) const {
		if (points_to_remove_ratio == 0.0)
			return cloud;

		auto reduced_node = context.m_project->TransTreeFindNode(m_node_id)->CreateChild();
		reduced_node->Rename(L"Reduced cloud");
		auto reduced_element = context.Project().ElementCreate<Clouds::ICloud>();
		reduced_node->SetElement(reduced_element);
		auto reduced_cloud = reduced_element->GetData<Clouds::ICloud>();

		Clouds::PointsRange range;
		cloud->GetAccess().GetAllPoints(range);
		const int range_size = range.size();
		std::vector<Clouds::Point3D> coordinates;
		range.GetXYZ(coordinates);
		std::vector<Clouds::Color> colors;
		range.GetColors(colors);
		std::vector<Clouds::Point3D> reduced_coordinates;
		std::vector<Clouds::Color> reduced_colors;

		std::vector<std::pair<Clouds::Point3D, Clouds::Color>> points;
		points.reserve(range_size);
		std::vector<std::pair<Clouds::Point3D, Clouds::Color>> reduced_points;
		std::transform(coordinates.begin(), coordinates.end(), colors.begin(), std::back_inserter(points),
			[](Clouds::Point3D xyz, Clouds::Color rgb) {
				return std::pair<Clouds::Point3D, Clouds::Color>(xyz, rgb);
			});

		FeedbackUpdater feedback_updater(context, range_size);
		std::mt19937 random_numbers_engine;
		std::uniform_real_distribution<float> distribution(0.0, 1.0);
		reduced_points.reserve(range_size * (1 - points_to_remove_ratio));
		std::copy_if(points.begin(), points.end(), std::back_inserter(reduced_points), [&](const auto& point) {
			feedback_updater.count();
			return distribution(random_numbers_engine) > points_to_remove_ratio;
			});
		for (const std::pair<Clouds::Point3D, Clouds::Color>& point : reduced_points) {
			reduced_coordinates.push_back(point.first);
			reduced_colors.push_back(point.second);
		}

		Clouds::PointsRange reduced_range;
		reduced_cloud->GetAccess().AllocPoints(reduced_coordinates.size(), &reduced_range);
		reduced_range.SetXYZ(reduced_coordinates);
		reduced_range.SetColors(reduced_colors);
		node_id = reduced_node->GetID();

		return reduced_cloud;
	}

	/*Function finds neighbours for all the points from the cloud in the given node and makes the lists of graph elements assigned to the vector in this structure.*/
	void BuildPointsGraphStructure(Context& context, ResourceID child_node_id, std::vector<Point>& graph, float max_distance)	{
		const auto cloud = context.m_project->TransTreeFindNode(child_node_id)->GetElement()->GetData<Clouds::ICloud>();
		Clouds::PointsRange range;
		cloud->GetAccess().GetAllPoints(range);
		const int range_size = range.size();
		std::vector<Clouds::Point3D> coordinates;
		std::vector<Clouds::Color> colors;
		range.GetXYZ(coordinates);
		range.GetColors(colors);

		graph.reserve(range_size);
		std::transform(coordinates.begin(), coordinates.end(), colors.begin(), std::back_inserter(graph),
			[](Clouds::Point3D& xyz, Clouds::Color& color) {
				return Point(xyz, color);
			});

		std::vector<Point*> sorted_by_x_ptrs;
		sorted_by_x_ptrs.reserve(range_size);
		std::transform(graph.begin(), graph.end(), std::back_inserter(sorted_by_x_ptrs), [](Point& point) {
			return &point;
			});
		std::sort(sorted_by_x_ptrs.begin(), sorted_by_x_ptrs.end(), [](Point* a, Point* b) {
			return a->xyz.x() < b->xyz.x();
			});

		FeedbackUpdater feedback_updater(context, range_size);
		for (auto itr = sorted_by_x_ptrs.begin(), end = sorted_by_x_ptrs.end(); itr != end; ++itr) {
			const auto& base_point_ptr = *itr;
			const auto base_point_xyz = base_point_ptr->xyz.cast<double>();
			const float x_max = max_distance + base_point_ptr->xyz.x();

			for (auto n_itr = itr + 1; n_itr != end; ++n_itr) {
				const auto& neighbour_ptr = *n_itr;
				if (neighbour_ptr->xyz.x() > x_max)
					break;

				if (Math::CalcPointToPointDistance3D(base_point_xyz, neighbour_ptr->xyz.cast<double>()) <= max_distance) {
					base_point_ptr->neighbours.push_back(neighbour_ptr);
					neighbour_ptr->neighbours.push_back(base_point_ptr);
				}
			}

			std::sort(base_point_ptr->neighbours.begin(), base_point_ptr->neighbours.end(), [&base_point_ptr](Point* a, Point* b) {
				return Math::CalcPointToPointDistance3D(a->xyz.cast<double>(), base_point_ptr->xyz.cast<double>())
					< Math::CalcPointToPointDistance3D(b->xyz.cast<double>(), base_point_ptr->xyz.cast<double>());
				});
			feedback_updater.count();
		}
	}
	
	/*Function copies all points' data from the graph structure to the given range of the cloud.*/
	void CopyGraphStructuctureToCloud(std::vector<Point>& graph, Clouds::PointsRange& range) {
		std::vector<Clouds::Point3D> coordinates;
		std::vector<Clouds::Color> colors;
		coordinates.reserve(graph.size());
		colors.reserve(graph.size());
		for (const auto& point : graph) {
			coordinates.push_back(point.xyz);
			colors.push_back(point.color);
		}
		range.SetXYZ(coordinates);
		range.SetColors(colors);
	}

	/*Function marks lonely points as 'deleted' and smooths the cloud by projecting all points in the structure to a local fitted sphere or plane.*/
	void RemoveNoise(Context& context, std::vector<Point>& graph, int min_neighbours, int neighbours_to_smooth)	{
		const int number_of_points = graph.size();
		FeedbackUpdater feedback_updater(context, number_of_points);
		for (auto& point : graph) {
			if (point.neighbours.size() < min_neighbours)
				point.is_deleted = true;
		}

		std::vector<Clouds::Point3D> new_coordinates;
		new_coordinates.reserve(number_of_points);
		std::transform(graph.begin(), graph.end(), std::back_inserter(new_coordinates), [&](Point& point) {
			const auto& neighbours = point.neighbours;

			if (neighbours.size() < neighbours_to_smooth || point.is_deleted)
				return point.xyz;

			auto base_point = point.xyz.cast<double>();
			std::vector<decltype(base_point)> n_xyz = { base_point };
			n_xyz.reserve(neighbours_to_smooth + 1);
			std::transform(neighbours.begin(), neighbours.begin() + neighbours_to_smooth, std::back_inserter(n_xyz), [](Point* n_itr) {
				return n_itr->xyz.cast<double>();
				});
			const auto sphere = Math::CalcBestSphere3D(n_xyz.begin(), n_xyz.end());
			const auto plane = Math::CalcBestPlane3D(n_xyz.begin(), n_xyz.end());

			const double sq_sum_sphere = std::accumulate(n_xyz.begin(), n_xyz.end(), 0.0, [&sphere](double init, auto& pt) {
				const auto point_onto_sphere = sphere.Project(pt).cast<double>();
				return init + pow(Math::CalcPointToPointDistance3D(point_onto_sphere, pt), 2);
				});
			const double sq_sum_plane = std::accumulate(n_xyz.begin(), n_xyz.end(), 0.0, [&plane](double init, auto& pt) {
				const auto point_onto_plane = Math::ProjectPointOntoPlane(plane, pt).cast<double>();
				return init + pow(Math::CalcPointToPointDistance3D(point_onto_plane, pt), 2);
				});

			const auto xyz = point.xyz.cast<double>();
			const auto new_xyz = (sq_sum_sphere < sq_sum_plane) ? sphere.Project(xyz) : Math::ProjectPointOntoPlane(plane, xyz);

			feedback_updater.count();
			return Clouds::Point3D(new_xyz.x(), new_xyz.y(), new_xyz.z());
			});

		std::transform(graph.begin(), graph.end(), new_coordinates.begin(), graph.begin(), [](Point& point, Clouds::Point3D& xyz) {
			point.xyz = xyz;
			return point;
			});

		for (auto& point : graph) {
			std::sort(point.neighbours.begin(), point.neighbours.end(), [point](Point* a, Point* b)
				{
					return Math::CalcPointToPointDistance3D(a->xyz.cast<double>(), point.xyz.cast<double>())
						< Math::CalcPointToPointDistance3D(b->xyz.cast<double>(), point.xyz.cast<double>());
				});
		};
	}
	
	/*Function filters all points using the median colors of neighbours.*/
	void MedianFilterColors(Context& context, std::vector<Point>& graph, int median_filter_size) {
		const int number_of_points = graph.size();
		FeedbackUpdater feedback_updater(context, number_of_points);
		std::vector<Clouds::Color> colors;
		colors.reserve(number_of_points);

		std::transform(graph.begin(), graph.end(), colors.begin(), [&](const auto& point) {
			if (point.is_deleted || point.neighbours.size() < median_filter_size) {
				feedback_updater.count();
				return point.color;
			}

			typedef std::pair<double, double> WeightedColor;
			std::vector<WeightedColor> red, green, blue;
			const auto begin = point.neighbours.begin();
			const auto point_xyz = point.xyz.cast<double>();
			const auto half_weight = 0.5 * std::accumulate(begin, begin + median_filter_size, 0.0, [&](double init, Point* neighbour) {
				const double weight = pow(Math::CalcPointToPointDistance3D(point_xyz, neighbour->xyz.cast<double>()), -2.0);
				red.push_back(WeightedColor(neighbour->color.x(), weight));
				green.push_back(WeightedColor(neighbour->color.y(), weight));
				blue.push_back(WeightedColor(neighbour->color.z(), weight));
				return init + weight;
				});

			std::vector<std::vector<WeightedColor>*> rgb = { &red,&green,&blue };
			std::vector<double> median_color;
			std::transform(rgb.begin(), rgb.end(), std::back_inserter(median_color), [&](std::vector<WeightedColor>* channel) {
				std::sort(channel->begin(), channel->end(), [](WeightedColor& left, WeightedColor& right) {
					return left.first > right.first;
					});

				double current_weight = 0.0;
				for (const auto& color : *channel) {
					current_weight += color.second;
					if (current_weight > half_weight)
						return color.first;
				}
				});

			feedback_updater.count();
			return Clouds::Color(median_color[0], median_color[1], median_color[2], 0);
			});

		for (int i = 0; i < number_of_points; ++i)
			graph[i].color = colors[i];
	}

	/*Function groups the points into clusters choosing one of the implemented recursive Hausdorf algorithm depending on the number of points.*/
	void Hausdorf(Context& context, std::vector<Point>& graph, const double hausdorf_threshold) {
		// set initial values
		for (auto& point : graph) {
			point.is_calculated = false;
			point.group_number = 0;
		}
		FeedbackUpdater feedback_updater(context, graph.size());
		int group_counter = 0;

		// realize Hausdorf algorithm for each point
		for (auto& point : graph) {
			if (point.group_number != 0 || point.is_deleted)
				continue;
			++group_counter;

			// check if recursion with this number of points results stack overflow
			if (graph.size() < hausdorf_threshold)
				DFSHausdorfGroup(point, group_counter, feedback_updater);
			else
				MixSearchedHausdorfGroup(point, group_counter, feedback_updater);
		}
	}
	
	/*Function groups points into a cluster with given number, starting from a given point using recursive Hausdorf's algorithm adapted to the graph structure of data. It should not be used for large clouds due to stack overflow.*/
	void DFSHausdorfGroup(Point& start_point, const int& group_counter, FeedbackUpdater& feedback_updater) {
		// relize the algorithm for the point
		start_point.group_number = group_counter;
		GroupNeighbours(start_point, start_point.xyz, feedback_updater);

		// repeate for the neighbours recursively
		for (auto p_neighbour : start_point.neighbours) {
			if (p_neighbour->isReady() || p_neighbour->group_number != start_point.group_number)
				continue;

			DFSHausdorfGroup(*p_neighbour, group_counter, feedback_updater);
		}
	}
	
	/*Function groups points into a cluster with given number, starting from a given point using recursive Hausdorf's algorithm with 2 steps in one iteration of the function. This modification is less clear but prevents stack overflow.*/
	void MixSearchedHausdorfGroup(Point& start_point, const int& group_counter, FeedbackUpdater& feedback_updater) {
		// relize the algorithm for the point
		start_point.group_number = group_counter;
		GroupNeighbours(start_point, start_point.xyz, feedback_updater);

		// repeate for the neighbours in a loop
		std::vector<Point*> neighbours;
		for (auto p_neighbour : start_point.neighbours) {
			if (p_neighbour->isReady() || p_neighbour->group_number != start_point.group_number)
				continue;

			GroupNeighbours(*p_neighbour, start_point.xyz, feedback_updater);
			neighbours.push_back(p_neighbour);
		}

		// repeate for the neighbours of neighbours recursively
		for (auto p_point : neighbours) {
			for (auto p_neighbour : p_point->neighbours) {
				if (p_neighbour->isReady() || p_neighbour->group_number != p_point->group_number)
					continue;

				MixSearchedHausdorfGroup(*p_neighbour, group_counter, feedback_updater);
			}
		}
	}

	void GroupNeighbours(Point& point, Clouds::Point3D base_xyz, FeedbackUpdater& feedback_updater) {
		const auto xyz = base_xyz.cast<double>();
		for (auto p_neighbour : point.neighbours) {
			if (p_neighbour->isReady())
				continue;

			const double distance = Math::CalcPointToPointDistance3D(xyz, p_neighbour->xyz.cast<double>());
			if (CheckColorGradient(point.color, p_neighbour->color, distance))
				p_neighbour->group_number = point.group_number;
		}
		point.is_calculated = true;
		feedback_updater.count();
	}
	
	/*Function returns 'true' if the color gradient of given data converted to HSV are within tolerance and 'false' otherwise.*/	
	bool CheckColorGradient(Clouds::Color color_1, Clouds::Color color_2, double distance) {
		const double R1 = color_1.x();
		const double G1 = color_1.y();
		const double B1 = color_1.z();
		const double R2 = color_2.x();
		const double G2 = color_2.y();
		const double B2 = color_2.z();

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
			return true;
		if (abs(V1 - V2) > max_brightness_gradient * distance)
			return false;
		if (S1 == 0.0 && S2 == 0.0)
			return true;
		if (abs(S1 - S2) > max_saturation_gradient * distance + coverage_rate / (255.0 * (V1 + V2) / 2))
			return false;
		const double HUE_INACCURACY = coverage_rate * (180.0 / (255.0 * (V1 + V2) / 2) + 180.0 / (255.0 * (S1 + S2) / 2));
		if (std::min(abs(H1 - H2), 360.0 - abs(H1 - H2)) > max_hue_gradient * distance + HUE_INACCURACY)
			return false;
		return true;
	}

	/*Function finds groups with a number of elements less than minimum given by a user and marks them as ungrouped.*/
	void RemoveTooSmallGroups(Context& context, std::vector<Point>& graph) {
		FeedbackUpdater feedback_updater(context, 2 * graph.size());
		std::map<int, int> repetition_numbers;
		
		for (const auto& point : graph) {
			++repetition_numbers[point.group_number];
			feedback_updater.count();
		}
		
		for (auto& point: graph) {
			if (repetition_numbers[point.group_number] < min_points_group)
				point.group_number = 0;
			feedback_updater.count();
		}
	}
	
	/*Function adds ungrouped points surrounded by neighbours assigned to groups to the locally most popular group.*/
	void FillHolesInClusters(Context& context, std::vector<Point>& graph) {
		FeedbackUpdater feedback_updater(context, graph.size());

		// set initial states
		for (auto& point : graph) {
			if (point.group_number != 0 || point.is_deleted) {
				point.is_calculated = true;
				feedback_updater.count();
			}
			else
				point.is_calculated = false;
		}

		bool any_change = true;
		while (any_change) {
			any_change = false;
			std::vector<int> new_group_numbers;
			new_group_numbers.reserve(graph.size());

			for (auto& point : graph) {
				if (point.is_calculated) {
					new_group_numbers.push_back(0);
					continue;
				}

				std::map<int, int> repetition_number;
				for (const auto& p_neighbour : point.neighbours) {
					if (p_neighbour->group_number != 0)
						++repetition_number[p_neighbour->group_number];
				}

				auto max_group = std::max_element(repetition_number.begin(), repetition_number.end(), [](auto& a, auto& b) {
					return a.second < b.second;
					});

				// assign point to most popular group if it exist
				if (max_group != repetition_number.end() && repetition_number[0] < join_coefficient * point.neighbours.size()) {
					new_group_numbers.push_back(max_group->first);
					point.is_calculated = true;
					any_change = true;
					feedback_updater.count();
				}
				else
					new_group_numbers.push_back(0);
			}

			// write values to the structure at the end of each loop to guarantee isotropic addition of points
			for (int i = 0; i < graph.size(); i++) {
				if (new_group_numbers[i] != 0)
					graph[i].group_number = new_group_numbers[i];
			}
		}
	}

	/*Function merge clusters with a relatively high number of points on common border.*/
	void MergeClusters(Context& context, std::vector<Point>& graph) {
		typedef std::pair<int, std::vector<Point*>> PointCluster;
		std::vector<PointCluster> clusters;

		for (auto& point : graph) {
			if (point.group_number == 0)
				continue;

			auto cluster = std::find_if(clusters.begin(), clusters.end(), [&point](PointCluster& cluster) {
				return cluster.first == point.group_number;
				});

			if (cluster == clusters.end())
				clusters.push_back(PointCluster(point.group_number, std::vector<Point*>({ &point })));
			else
				cluster->second.push_back(&point);
		}

		std::sort(clusters.begin(), clusters.end(), [](PointCluster& a, PointCluster& b) {
			return a.second.size() < b.second.size();
			});

		FeedbackUpdater feedback_updater(context, clusters.size());

		for (auto itr_smaller = clusters.begin(), end = clusters.end(); itr_smaller != end; ++itr_smaller) {
			PointCluster& smaller_cluster = *itr_smaller;
			int max_neighbours = 0;

			std::for_each(itr_smaller + 1, clusters.end(), [&](PointCluster& bigger_cluster) {
				const int& bigger_group = bigger_cluster.first;

				const int border_points = std::count_if(smaller_cluster.second.begin(), smaller_cluster.second.end(),
					[&](Point* point) {
						return std::any_of(point->neighbours.begin(), point->neighbours.end(), [&](Point* neigbour) {
							return neigbour->group_number == bigger_group;
							});
					});

				// add sizes of all parts merged into the smaller cluster
				const int size = std::accumulate(clusters.begin(), itr_smaller, smaller_cluster.second.size(),
					[&](size_t init, PointCluster& cluster) {
						return init + (cluster.first == smaller_cluster.first) ? cluster.second.size() : 0;
					});

				// compare border size with number of points in the smaller cluster circumference
				if (border_points < sqrt(4 * PI * size) * max_distance / merge_coefficient || border_points < max_neighbours)
					return;

				// merge the smaller cluster and all clusters merged into it, into the bigger cluster
				std::for_each(clusters.begin(), itr_smaller + 1, [&](PointCluster& cluster) {
					if (cluster.first != smaller_cluster.first)
						return;

					for (auto p_point : cluster.second)
						p_point->group_number = bigger_group;

					cluster.first = bigger_group;
					});

				max_neighbours = border_points;
				});

			feedback_updater.count();
		}
	}
	
	/*Function changes group numbers to subsequent numbers with '0' assigned to ungrouped points and write them to a given vector of layer values.*/
	void TidyUpGroupNumbers(Context& context, std::vector<StoredReal>& cluster_numbers, std::vector<Point>& graph) {
		FeedbackUpdater feedback_updater(context, 2 * graph.size());

		std::map<int, int> repetition_numebers;
		for (auto& point : graph) {
			if (point.group_number != 0)
				++repetition_numebers[point.group_number];
			feedback_updater.count();
		}

		std::vector<std::pair<const int, int>*> ordered_groups;
		for (auto& group_size : repetition_numebers)
			ordered_groups.push_back(&group_size);

		std::sort(ordered_groups.begin(), ordered_groups.end(), [](std::pair<const int, int>*& a, std::pair<const int, int>*& b) {
			return a->second > b->second;
			});

		std::map<int, int> new_group_numbers;
		for (int j = 0; j < ordered_groups.size(); ++j) {
			// ignore number '0' assigned to ungrouped points
			new_group_numbers[ordered_groups[j]->first] = j + 1;
		}

		for (auto& point : graph) {
			point.group_number = new_group_numbers[point.group_number];
			cluster_numbers.push_back(point.group_number);
			feedback_updater.count();
		}
	}

	/*Function finds the best fitting cluster to the loudspeaker and write to the given layer '1.0' for points of this cluster and '0.0' for others.*/
	void FindTheSpeaker(Context& context, std::vector<StoredReal>& is_speaker, std::vector<Point>& graph) {
		std::map<int, std::vector<Point*>> clusters;

		for (auto& point : graph) {
			if (clusters.find(point.group_number) == clusters.end())
				clusters[point.group_number] = { &point };
			else
				clusters[point.group_number].push_back(&point);
		}

		FeedbackUpdater feedback_updater(context, clusters.size(), true);

		std::map<int, double> fitting_errors;
		for (auto& cluster : clusters) {
			// ignore ungrouped points
			fitting_errors[cluster.first] = GetBoxFittingError(context, cluster.second);
			feedback_updater.count();
		}

		auto speaker_itr = std::min_element(fitting_errors.begin(), fitting_errors.end(), [](auto& left, auto& right) {
			return left.second < right.second;
			});

		for (auto& point : graph)
			is_speaker.push_back(point.group_number == speaker_itr->first ? 1.0 : 0.0);

		//// print the values of all fitting errors
		//for (auto& itr : fitting_errors)
		//	OGX_LINE.Msg(Info, std::to_wstring(itr.first) + L". " + std::to_wstring(itr.second));
	}

	/*Function calculates a box fitted to the borders of given points cluster. Returns sum of angles and dimensions relative errors and relative RMSE for this fitting.*/
	double GetBoxFittingError(Context& context, const std::vector<Point*>& cluster) {
		Math::Point3D first_diagonal_end, second_diagonal_end;
		double max_distance_between_points = 0.0;
		for (int i = 0; i < cluster.size(); ++i) {
			const auto& first_point = cluster[i]->xyz.cast<double>();

			for (int j = i + 1; j < cluster.size(); ++j) {
				const auto& second_point = cluster[j]->xyz.cast<double>();
				const double local_distance = Math::CalcPointToPointDistance3D(first_point, second_point);

				if (local_distance > max_distance_between_points) {
					max_distance_between_points = local_distance;
					first_diagonal_end = first_point;
					second_diagonal_end = second_point;
				}
			}
		}
		const Math::Point3D center_point = (first_diagonal_end + second_diagonal_end) / 2;
		const Math::Line3D longest_diagonal = Math::CalcLine3D(first_diagonal_end, second_diagonal_end);
		
		Math::Point3D farthest_from_diagonal;
		double max_distance_from_diagonal = 0.0;
		for (const auto p_point : cluster) {
			const auto& point = p_point->xyz.cast<double>();
			const double local_distance = Math::CalcPointToLineDistance3D(point, longest_diagonal);

			if (local_distance > max_distance_from_diagonal) {
				max_distance_from_diagonal = local_distance;
				farthest_from_diagonal = point;
			}
		}
		const Math::Plane3D section_plane = Math::CalcPlane3D(first_diagonal_end, second_diagonal_end, farthest_from_diagonal);
		
		Math::Point3D farthest_from_section_plane;
		double max_distance_from_section_plane_and_center = 0.0;
		for (const auto p_point : cluster) {
			const auto& point = p_point->xyz.cast<double>();
			double local_distance = Math::CalcPointToPlaneDistance3D(point, section_plane);
			local_distance += Math::CalcPointToPointDistance3D(point, center_point);

			if (local_distance > max_distance_from_section_plane_and_center) {
				max_distance_from_section_plane_and_center = local_distance;
				farthest_from_section_plane = point;
			}
		}

		const bool is_first_end_on_0XX = Math::CalcPointToPointDistance3D(first_diagonal_end, farthest_from_section_plane)
			< Math::CalcPointToPointDistance3D(second_diagonal_end, farthest_from_section_plane);

		const Math::Point3D& vertex_000 = is_first_end_on_0XX ? first_diagonal_end : second_diagonal_end;
		const Math::Point3D& vertex_111 = is_first_end_on_0XX ? second_diagonal_end : first_diagonal_end;

		const bool is_far_diag_on_0XX = Math::CalcPointToPointDistance3D(farthest_from_diagonal, farthest_from_section_plane)
			< Math::CalcPointToPointDistance3D(farthest_from_diagonal, vertex_000);

		const Math::Point3D vertex_011 = is_far_diag_on_0XX ? farthest_from_diagonal : 2 * center_point - farthest_from_diagonal;
		const Math::Point3D vertex_100 = is_far_diag_on_0XX ? 2 * center_point - farthest_from_diagonal : farthest_from_diagonal;
		
		const bool is_far_sec_pl_001 = Math::CalcPointToPointDistance3D(farthest_from_section_plane, vertex_000)
			< Math::CalcPointToPointDistance3D(farthest_from_section_plane, vertex_011);

		const Math::Point3D vertex_001 = is_far_sec_pl_001 ? farthest_from_section_plane : vertex_000 + vertex_011 - farthest_from_section_plane;
		const Math::Point3D vertex_010 = is_far_sec_pl_001 ? vertex_000 + vertex_011 - farthest_from_section_plane : farthest_from_section_plane;
	
		const Math::Point3D vertex_101 = 2 * center_point - vertex_010;
		const Math::Point3D vertex_110 = 2 * center_point - vertex_001;

		const Math::Plane3D wall_0XX = Math::CalcPlane3D(vertex_000, vertex_010, vertex_001);
		const Math::Plane3D wall_1XX = Math::CalcPlane3D(vertex_111, vertex_110, vertex_101);
		const Math::Plane3D wall_X0X = Math::CalcPlane3D(vertex_000, vertex_001, vertex_100);
		const Math::Plane3D wall_X1X = Math::CalcPlane3D(vertex_111, vertex_011, vertex_110);
		const Math::Plane3D wall_XX0 = Math::CalcPlane3D(vertex_000, vertex_010, vertex_100);
		const Math::Plane3D wall_XX1 = Math::CalcPlane3D(vertex_111, vertex_011, vertex_101);
		
		const double box_fitting_square_error = std::accumulate(cluster.begin(), cluster.end(), 0.0,
			[&](double init, Point* p_point) {
				const auto& point = p_point->xyz.cast<double>();
				const auto min_point_wall_distance = std::min({
					Math::CalcPointToPlaneDistance3D(point, wall_0XX),
					Math::CalcPointToPlaneDistance3D(point, wall_1XX),
					Math::CalcPointToPlaneDistance3D(point, wall_X0X),
					Math::CalcPointToPlaneDistance3D(point, wall_X1X),
					Math::CalcPointToPlaneDistance3D(point, wall_XX0),
					Math::CalcPointToPlaneDistance3D(point, wall_XX1)
					});
				return init + pow(min_point_wall_distance, 2);
			});
		const double rmse = sqrt(box_fitting_square_error / cluster.size());
		const double max_distance = Math::CalcPointToLineDistance3D(farthest_from_diagonal, longest_diagonal);

		double relative_angle_error = 0.0;
		const Math::Vector3D x_vec(vertex_100 - vertex_000);
		const Math::Vector3D y_vec(vertex_010 - vertex_000);
		const Math::Vector3D z_vec(vertex_001 - vertex_000);
		relative_angle_error += abs(Math::CalcAngleBetweenTwoVectors(x_vec, y_vec) - PI / 2) / (PI / 2);
		relative_angle_error += abs(Math::CalcAngleBetweenTwoVectors(y_vec, z_vec) - PI / 2) / (PI / 2);
		relative_angle_error += abs(Math::CalcAngleBetweenTwoVectors(z_vec, x_vec) - PI / 2) / (PI / 2);
		
		double relative_demention_error = 0.0;
		relative_demention_error += abs(x_vec.norm() - speaker_dim_x) / speaker_dim_x;
		relative_demention_error += abs(y_vec.norm() - speaker_dim_y) / speaker_dim_y;
		relative_demention_error += abs(z_vec.norm() - speaker_dim_z) / speaker_dim_z;		

		const double combined_relative_error = (rmse / max_distance + relative_angle_error + relative_demention_error) / 7.0;
		//// print the values of all fitting errors
		//OGX_LINE.Msg(Info, std::to_wstring(x_vec.norm()) + L", " + std::to_wstring(y_vec.norm()) + L", "
		//	+ std::to_wstring(z_vec.norm()) + L", Angle error: " + std::to_wstring(relative_angle_error) + L", RMSE: "
		//	+ std::to_wstring(rmse / max_distance));
		
		return combined_relative_error;
	}
};

OGX_EXPORT_METHOD(LoudspeakerSearching)
