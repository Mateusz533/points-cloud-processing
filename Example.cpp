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
	std::vector<Point> mPointsGraphStructure;

	// parameters
	ResourceID mNodeId;
	float mMaxDistance = 0.020;
	float mPointsToRemoveRatio = 0.9;
	int mMinNeighbours = 10;
	int mMinPointsGroup = 100;
	int mNeighboursToSmooth = 10;
	int mMedianFlterSize = 10;
	float mMaxBrightnessGradient = 7.5;
	float mMaxSaturationGradient = 2.5;
	float mMaxHueGradient = 180;
	float mCoverageRate = 0.4;					// additional tolerance of hue and saturation differences for mean HSV=(*,1,1), 10 times smaller for HSV=(*,10,10) etc.
	float mJoinCoefficient = 0.8;				// maximum ratio of ungrouped point among the neighbours to join the ungrouped point to any of the clusters
	float mMergeCoefficient = 0.02;				// the higher it is the more clusters will be merged
	float mSpeakerDimX = 0.18;
	float mSpeakerDimY = 0.16;
	float mSpeakerDimZ = 0.12;

	// constructor
	LoudspeakerSearching() : EasyMethod(L"Mateusz Frejlich", L"The program segments the cloud and recognizes the loudspeaker inside it.\nThe default values of all parameters are best fitted to the tested cloud.")
	{
	}

	// add input/output parameters
	virtual void DefineParameters(ParameterBank& bank) 	{
		bank.Add(L"Node ID", mNodeId).AsNode();
		bank.Add(L"Points to remove ratio", mPointsToRemoveRatio).Min(0.0).Max(1.0);
		bank.Add(L"Max. distance to neighbour in a cluster [u]", mMaxDistance).Describe(L"In distance unit used in the cloud.");
		bank.Add(L"Min. number of neighbours", mMinNeighbours).Describe(L"Points below this value will be deleted.");
		bank.Add(L"Min. number of points in a cluster", mMinPointsGroup).Describe(L"Clusters below this value will be canceled.");
		bank.Add(L"Number of neigbours to smooth", mNeighboursToSmooth).Describe(L"Number of neighbours used to smoothing each point by fitting a surface. If there are less neighbours, smoothing will not be done.");
		bank.Add(L"Median filter size", mMedianFlterSize).Describe(L"Number of neighbours used to smoothing the color of each point by median filter. If there are less neighbours, smoothing will not be done.");
		bank.Add(L"Max. hue gradient in a cluster [\260/u]", mMaxHueGradient).Describe(L"One degree per distance unit.");
		bank.Add(L"Max. saturation gradient in a cluster [1/u]", mMaxSaturationGradient).Describe(L"Max. scope value per distance unit.");
		bank.Add(L"Max. brightness gradient in a cluster [1/u]", mMaxBrightnessGradient).Describe(L"Max. scope value per distance unit.");
		bank.Add(L"Coverage rate", mCoverageRate).Describe(L"The higher it is, the higher tolerance based on measurement inaccuracy of hue and saturation will be considered.");
		bank.Add(L"Join coefficient", mJoinCoefficient).Describe(L"The higher it is, the more ungrouped points will be joined to any cluster.");
		bank.Add(L"Merge coefficient", mMergeCoefficient).Describe(L"The higher it is, the more clusters will be merged.");
		bank.Add(L"Loudspeaker dimension X [u]", mSpeakerDimX).Describe(L"The largest dimension");
		bank.Add(L"Loudspeaker dimension Y [u]", mSpeakerDimY).Describe(L"Medium dimention");
		bank.Add(L"Loudspeaker dimension Z [u]", mSpeakerDimZ).Describe(L"The smallest dimension");
	}

	virtual void Run(Context& context) {
		// get cloud data
		Clouds::ICloud* pOriginalCloud = GetValidCloud(context, mNodeId);
		OGX_LINE.Msg(Info, L"All data is correct. Processing has started.");

		// reduce the number of points
		ResourceID newNodeId;
		Clouds::ICloud* pCloud = CreateRandomlyReducedCloud(context, pOriginalCloud, newNodeId, mPointsToRemoveRatio);
		OGX_LINE.Msg(Info, L"Reduced cloud has been done.");

		// create points structure from created cloud
		BuildPointsGraphStructure(context, newNodeId, mPointsGraphStructure, mMaxDistance);
		OGX_LINE.Msg(Info, L"Points structure has been done.");

		// create a new layer
		const String layerName = L"Segmentation";
		auto layers = pCloud->FindLayers(layerName);
		auto* pSegmentsLayer = layers.empty() ? pCloud->CreateLayer(layerName, 0.0) : layers[0];

		// reserve memory for results vector
		Clouds::PointsRange range;
		pCloud->GetAccess().GetAllPoints(range);
		std::vector<StoredReal> segmentNumbers;
		segmentNumbers.reserve(mPointsGraphStructure.size());

		// remove lonely points and smooth the cloud
		RemoveNoise(context, mPointsGraphStructure, mMinNeighbours, mNeighboursToSmooth);

		// reduce noise in colors using median filter
		MedianFilterColors(context, mPointsGraphStructure, mMedianFlterSize);
		OGX_LINE.Msg(Info, L"Noise filtration has been done.");

		// visualize the cloud smoothing
		CopyGraphStructuctureToCloud(mPointsGraphStructure, range);

		// make segmentation using Hausdorf algorithm
		GroupHausdorf(context, mPointsGraphStructure, HAUSDORF_THRESHOLD);
		OGX_LINE.Msg(Info, L"Segmentation by Hausdorf algorithm has been done.");

		// fill holes with ungrouped points in clusters
		FillHolesInClusters(context, mPointsGraphStructure);

		// remove too small groups
		RemoveTooSmallGroups(context, mPointsGraphStructure);
		OGX_LINE.Msg(Info, L"Removing too small groups has been done.");

		// fill holes made by removing small groups
		FillHolesInClusters(context, mPointsGraphStructure);
		OGX_LINE.Msg(Info, L"Filling holes in clusters has been done.");

		// merge clusters
		MergeClusters(context, mPointsGraphStructure);
		OGX_LINE.Msg(Info, L"Merging clusters has been done.");

		// change group numbers to subsequent numbers with '0' assigned to ungrouped points
		TidyUpGroupNumbers(context, segmentNumbers, mPointsGraphStructure);
		OGX_LINE.Msg(Info, L"The ordering of group numbers has been done.");

		// assign all cluster numbers to created layer
		range.SetLayerVals(segmentNumbers, *pSegmentsLayer);
		OGX_LINE.Msg(Info, L"Whole segmentation process has been finished.");

		// create next layer
		const String newLayerName = L"Loudspeaker";
		const auto newLayers = pCloud->FindLayers(newLayerName);
		auto* pSpeakerLayer = newLayers.empty() ? pCloud->CreateLayer(newLayerName, 0.0) : newLayers[0];
		std::vector<StoredReal> isSpeaker;
		isSpeaker.reserve(mPointsGraphStructure.size());

		// find the best fitted cluster to the speaker
		FindTheSpeaker(context, isSpeaker, mPointsGraphStructure);
		OGX_LINE.Msg(Info, L"Loudspeaker identification has been done.");

		// mark found cluster in created layer
		range.SetLayerVals(isSpeaker, *pSpeakerLayer);
		OGX_LINE.Msg(Info, L"Processing has been finished. No errors.");

		//// check time
		//auto start = std::chrono::high_resolution_clock::now();
		//auto stop = std::chrono::high_resolution_clock::now();
		//auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		//auto time = duration.count();
	}

	/*Function reads data of all points from given node, check their correctness and returns a cloud of these points.*/
	Clouds::ICloud* GetValidCloud(Context& rContext, ResourceID nodeId)	{
		auto* pNode = rContext.m_project->TransTreeFindNode(nodeId);
		if (!pNode)
			ReportError(L"Invalid node!");
		
		auto* pElement = pNode->GetElement();
		if (!pElement)
			ReportError(L"Invalid element!");
		
		auto* pCloud = pElement->GetData<Clouds::ICloud>();
		if (!pCloud)
			ReportError(L"Invalid cloud!");
		
		if (mPointsToRemoveRatio < 0.0 || mPointsToRemoveRatio >= 1.0)
			ReportError(L"Invalid points to remove ratio!");
		
		if (mMaxDistance <= 0.0)
			ReportError(L"Invalid maximum distance to neighbour!");
		
		if (mMinNeighbours < 0)
			ReportError(L"Invalid minimum neigbours number!");
		
		if (mMinPointsGroup < 0)
			ReportError(L"Invalid minimum points cluster!");
		
		if (mNeighboursToSmooth < 5)
			ReportError(L"Invalid number of neighbours to smooth!");
		
		if (mMedianFlterSize < 2)
			ReportError(L"Invalid median filter size!");
		
		if (mMaxBrightnessGradient < 0.0)
			ReportError(L"Invalid maximum brightness gradient!");
		
		if (mMaxSaturationGradient < 0.0)
			ReportError(L"Invalid maximum saturation gradient!");
		
		if (mMaxHueGradient < 0.0)
			ReportError(L"Invalid maximum hue gradient!");
		
		if (mCoverageRate < 0.0)
			ReportError(L"Invalid coverage rate!");
		
		if (mJoinCoefficient < 0.0 || mJoinCoefficient > 1.0)
			ReportError(L"Invalid join coefficient!");
		
		if (mMergeCoefficient < 0.0)
			ReportError(L"Invalid merge coefficient!");
		
		if (mSpeakerDimX <= 0.0 || mSpeakerDimX < mSpeakerDimY || mSpeakerDimX < mSpeakerDimZ)
			ReportError(L"Invalid speaker dimention X!");
		
		if (mSpeakerDimY <= 0.0 || mSpeakerDimY < mSpeakerDimZ)
			ReportError(L"Invalid speaker dimention Y!");
		
		if (mSpeakerDimZ <= 0.0)
			ReportError(L"Invalid speaker dimention Z!");
		
		return pCloud;
	}

	/*Function makes a new node with given ID, assign to them the cloud with random points from a given cloud with ratio given by a user.*/
	Clouds::ICloud* CreateRandomlyReducedCloud(Context& rContext, Clouds::ICloud*& rCloud, ResourceID& rNodeId, float pointsToRemoveRatio) const {
		if (pointsToRemoveRatio == 0.0)
			return rCloud;

		auto* pReducedNode = rContext.m_project->TransTreeFindNode(mNodeId)->CreateChild();
		pReducedNode->Rename(L"Reduced cloud");
		auto* pReducedElement = rContext.Project().ElementCreate<Clouds::ICloud>();
		pReducedNode->SetElement(pReducedElement);
		auto* pReducedCloud = pReducedElement->GetData<Clouds::ICloud>();

		Clouds::PointsRange range;
		rCloud->GetAccess().GetAllPoints(range);
		const int rangeSize = range.size();
		std::vector<Clouds::Point3D> coordinates;
		range.GetXYZ(coordinates);
		std::vector<Clouds::Color> colors;
		range.GetColors(colors);
		std::vector<Clouds::Point3D> reducedCoordinates;
		std::vector<Clouds::Color> reducedColors;

		std::vector<std::pair<Clouds::Point3D, Clouds::Color>> points;
		points.reserve(rangeSize);
		std::vector<std::pair<Clouds::Point3D, Clouds::Color>> reducedPoints;
		std::transform(coordinates.begin(), coordinates.end(), colors.begin(), std::back_inserter(points),
			[](Clouds::Point3D xyz, Clouds::Color rgb) {
				return std::pair<Clouds::Point3D, Clouds::Color>(xyz, rgb);
			});

		FeedbackUpdater feedbackUpdater(rContext, rangeSize);
		std::mt19937 randomNumbersEngine;
		std::uniform_real_distribution<float> distribution(0.0, 1.0);
		reducedPoints.reserve(rangeSize * (1 - pointsToRemoveRatio));
		std::copy_if(points.begin(), points.end(), std::back_inserter(reducedPoints), [&](const auto& rPoint) {
			feedbackUpdater.count();
			return distribution(randomNumbersEngine) > pointsToRemoveRatio;
			});
		for (const std::pair<Clouds::Point3D, Clouds::Color>& rPoint : reducedPoints) {
			reducedCoordinates.push_back(rPoint.first);
			reducedColors.push_back(rPoint.second);
		}

		Clouds::PointsRange reducedRange;
		pReducedCloud->GetAccess().AllocPoints(reducedCoordinates.size(), &reducedRange);
		reducedRange.SetXYZ(reducedCoordinates);
		reducedRange.SetColors(reducedColors);
		rNodeId = pReducedNode->GetID();

		return pReducedCloud;
	}

	/*Function finds neighbours for all the points from the cloud in the given node and makes the lists of graph elements assigned to the vector in this structure.*/
	void BuildPointsGraphStructure(Context& rContext, ResourceID childNodeId, std::vector<Point>& rGraph, float maxDistance) {
		auto* pCloud = rContext.m_project->TransTreeFindNode(childNodeId)->GetElement()->GetData<Clouds::ICloud>();
		Clouds::PointsRange range;
		pCloud->GetAccess().GetAllPoints(range);
		const int rangeSize = range.size();
		std::vector<Clouds::Point3D> coordinates;
		std::vector<Clouds::Color> colors;
		range.GetXYZ(coordinates);
		range.GetColors(colors);

		rGraph.reserve(rangeSize);
		std::transform(coordinates.begin(), coordinates.end(), colors.begin(), std::back_inserter(rGraph),
			[](Clouds::Point3D& rXYZ, Clouds::Color& rColor) {
				return Point(rXYZ, rColor);
			});

		std::vector<Point*> ptrsSortedByX;
		ptrsSortedByX.reserve(rangeSize);
		std::transform(rGraph.begin(), rGraph.end(), std::back_inserter(ptrsSortedByX), [](Point& rPoint) {
			return &rPoint;
			});
		std::sort(ptrsSortedByX.begin(), ptrsSortedByX.end(), [](Point* pA, Point* pB) {
			return pA->xyz.x() < pB->xyz.x();
			});

		FeedbackUpdater feedbackUpdater(rContext, rangeSize);
		for (auto itr = ptrsSortedByX.begin(), end = ptrsSortedByX.end(); itr != end; ++itr) {
			const auto& rBasePointPtr = *itr;
			const auto basePointXYZ = rBasePointPtr->xyz.cast<double>();
			const float xMax = maxDistance + rBasePointPtr->xyz.x();

			for (auto neighbourItr = itr + 1; neighbourItr != end; ++neighbourItr) {
				const auto& rNeighbourPtr = *neighbourItr;
				if (rNeighbourPtr->xyz.x() > xMax)
					break;

				if (Math::CalcPointToPointDistance3D(basePointXYZ, rNeighbourPtr->xyz.cast<double>()) <= maxDistance) {
					rBasePointPtr->neighbours.push_back(rNeighbourPtr);
					rNeighbourPtr->neighbours.push_back(rBasePointPtr);
				}
			}

			const auto& rBasePointXYZ = rBasePointPtr->xyz.cast<double>();
			std::sort(rBasePointPtr->neighbours.begin(), rBasePointPtr->neighbours.end(), [&](Point* pA, Point* pB) {
				return Math::CalcPointToPointDistance3D(pA->xyz.cast<double>(), rBasePointXYZ)
					< Math::CalcPointToPointDistance3D(pB->xyz.cast<double>(), rBasePointXYZ);
				});
			feedbackUpdater.count();
		}
	}
	
	/*Function copies all points' data from the graph structure to the given range of the cloud.*/
	void CopyGraphStructuctureToCloud(std::vector<Point>& rGraph, Clouds::PointsRange& rRange) {
		std::vector<Clouds::Point3D> coordinates;
		std::vector<Clouds::Color> colors;
		coordinates.reserve(rGraph.size());
		colors.reserve(rGraph.size());
		for (const auto& rPoint : rGraph) {
			coordinates.push_back(rPoint.xyz);
			colors.push_back(rPoint.color);
		}
		rRange.SetXYZ(coordinates);
		rRange.SetColors(colors);
	}

	/*Function marks lonely points as 'deleted' and smooths the cloud by projecting all points in the structure to a local fitted sphere or plane.*/
	void RemoveNoise(Context& rContext, std::vector<Point>& rGraph, int minNeighbours, int neighboursToSmooth)	{
		const int numberOfPoints = rGraph.size();
		FeedbackUpdater feedbackUpdater(rContext, numberOfPoints);
		for (auto& rPoint : rGraph) {
			if (rPoint.neighbours.size() < minNeighbours)
				rPoint.is_deleted = true;
		}

		std::vector<Clouds::Point3D> newCoordinates;
		newCoordinates.reserve(numberOfPoints);
		std::transform(rGraph.begin(), rGraph.end(), std::back_inserter(newCoordinates), [&](Point& rPoint) {
			const auto& rNeighbours = rPoint.neighbours;

			if (rNeighbours.size() < neighboursToSmooth || rPoint.is_deleted)
				return rPoint.xyz;

			auto basePoint = rPoint.xyz.cast<double>();
			std::vector<decltype(basePoint)> nXYZ = { basePoint };
			nXYZ.reserve(neighboursToSmooth + 1);
			std::transform(rNeighbours.begin(), rNeighbours.begin() + neighboursToSmooth, std::back_inserter(nXYZ), [](Point* nItr) {
				return nItr->xyz.cast<double>();
				});
			const auto sphere = Math::CalcBestSphere3D(nXYZ.begin(), nXYZ.end());
			const auto plane = Math::CalcBestPlane3D(nXYZ.begin(), nXYZ.end());

			const double sqSumSphere = std::accumulate(nXYZ.begin(), nXYZ.end(), 0.0, [&sphere](double init, auto& rPoint) {
				const auto pointOntoSphere = sphere.Project(rPoint).cast<double>();
				return init + pow(Math::CalcPointToPointDistance3D(pointOntoSphere, rPoint), 2);
				});
			const double sqSumPlane = std::accumulate(nXYZ.begin(), nXYZ.end(), 0.0, [&plane](double init, auto& rPoint) {
				const auto point_onto_plane = Math::ProjectPointOntoPlane(plane, rPoint).cast<double>();
				return init + pow(Math::CalcPointToPointDistance3D(point_onto_plane, rPoint), 2);
				});

			const auto xyz = rPoint.xyz.cast<double>();
			const auto newXYZ = (sqSumSphere < sqSumPlane) ? sphere.Project(xyz) : Math::ProjectPointOntoPlane(plane, xyz);

			feedbackUpdater.count();
			return Clouds::Point3D(newXYZ.x(), newXYZ.y(), newXYZ.z());
			});

		std::transform(rGraph.begin(), rGraph.end(), newCoordinates.begin(), rGraph.begin(), [](Point& rPoint, Clouds::Point3D& rXYZ) {
			rPoint.xyz = rXYZ;
			return rPoint;
			});

		for (auto& rPoint : rGraph) {
			std::sort(rPoint.neighbours.begin(), rPoint.neighbours.end(), [rPoint](Point* a, Point* b) {
				return Math::CalcPointToPointDistance3D(a->xyz.cast<double>(), rPoint.xyz.cast<double>())
					< Math::CalcPointToPointDistance3D(b->xyz.cast<double>(), rPoint.xyz.cast<double>());
				});
		}
	}
	
	/*Function filters all points using the median colors of neighbours.*/
	void MedianFilterColors(Context& rContext, std::vector<Point>& rGraph, int medianFilterSize) {
		const int numberOfPoints = rGraph.size();
		FeedbackUpdater feedbackUpdater(rContext, numberOfPoints);
		std::vector<Clouds::Color> colors;
		colors.reserve(numberOfPoints);

		std::transform(rGraph.begin(), rGraph.end(), colors.begin(), [&](const auto& rPoint) {
			if (rPoint.is_deleted || rPoint.neighbours.size() < medianFilterSize) {
				feedbackUpdater.count();
				return rPoint.color;
			}

			typedef std::pair<double, double> WeightedColor;
			std::vector<WeightedColor> red, green, blue;
			const auto begin = rPoint.neighbours.begin();
			const auto pointXYZ = rPoint.xyz.cast<double>();
			const auto halfWeight = 0.5 * std::accumulate(begin, begin + medianFilterSize, 0.0, [&](double init, Point* pNeighbour) {
				const double weight = pow(Math::CalcPointToPointDistance3D(pointXYZ, pNeighbour->xyz.cast<double>()), -2.0);
				red.push_back(WeightedColor(pNeighbour->color.x(), weight));
				green.push_back(WeightedColor(pNeighbour->color.y(), weight));
				blue.push_back(WeightedColor(pNeighbour->color.z(), weight));
				return init + weight;
				});

			std::vector<std::vector<WeightedColor>*> rgb = { &red,&green,&blue };
			std::vector<double> medianColor;
			std::transform(rgb.begin(), rgb.end(), std::back_inserter(medianColor), [&](std::vector<WeightedColor>* pChannel) {
				std::sort(pChannel->begin(), pChannel->end(), [](WeightedColor& rLeft, WeightedColor& rRight) {
					return rLeft.first > rRight.first;
					});

				double currentWeight = 0.0;
				for (const auto& rColor : *pChannel) {
					currentWeight += rColor.second;
					if (currentWeight > halfWeight)
						return rColor.first;
				}
				});

			feedbackUpdater.count();
			return Clouds::Color(medianColor[0], medianColor[1], medianColor[2], 0);
			});

		for (int i = 0; i < numberOfPoints; ++i)
			rGraph[i].color = colors[i];
	}

	/*Function groups the points into clusters choosing one of the implemented recursive Hausdorf algorithm depending on the number of points.*/
	void GroupHausdorf(Context& rContext, std::vector<Point>& rGraph, const double hausdorfThreshold) {
		// set initial values
		for (auto& rPoint : rGraph) {
			rPoint.is_calculated = false;
			rPoint.group_number = 0;
		}
		FeedbackUpdater feedbackUpdater(rContext, rGraph.size());
		int groupCounter = 0;

		// realize Hausdorf algorithm for each point
		for (auto& rPoint : rGraph) {
			if (rPoint.group_number != 0 || rPoint.is_deleted)
				continue;
			++groupCounter;

			// check if recursion with this number of points results stack overflow
			if (rGraph.size() < hausdorfThreshold)
				GroupHausdorfDFS(rPoint, groupCounter, feedbackUpdater);
			else
				GroupHausdorfMixSearch(rPoint, groupCounter, feedbackUpdater);
		}
	}
	
	/*Function groups points into a cluster with given number, starting from a given point using recursive Hausdorf's algorithm adapted to the graph structure of data. It should not be used for large clouds due to stack overflow.*/
	void GroupHausdorfDFS(Point& rStartPoint, const int& rGroupCounter, FeedbackUpdater& rFeedbackUpdater) {
		// relize the algorithm for the point
		rStartPoint.group_number = rGroupCounter;
		GroupNeighbours(rStartPoint, rStartPoint.xyz, rFeedbackUpdater);

		// repeate for the neighbours recursively
		for (auto pNeighbour : rStartPoint.neighbours) {
			if (pNeighbour->isReady() || pNeighbour->group_number != rStartPoint.group_number)
				continue;

			GroupHausdorfDFS(*pNeighbour, rGroupCounter, rFeedbackUpdater);
		}
	}
	
	/*Function groups points into a cluster with given number, starting from a given point using recursive Hausdorf's algorithm with 2 steps in one iteration of the function. This modification is less clear but prevents stack overflow.*/
	void GroupHausdorfMixSearch(Point& rStartPoint, const int& rGroupCounter, FeedbackUpdater& rFeedbackUpdater) {
		// relize the algorithm for the point
		rStartPoint.group_number = rGroupCounter;
		GroupNeighbours(rStartPoint, rStartPoint.xyz, rFeedbackUpdater);

		// repeate for the neighbours in a loop
		std::vector<Point*> neighbours;
		for (auto pNeighbour : rStartPoint.neighbours) {
			if (pNeighbour->isReady() || pNeighbour->group_number != rStartPoint.group_number)
				continue;

			GroupNeighbours(*pNeighbour, rStartPoint.xyz, rFeedbackUpdater);
			neighbours.push_back(pNeighbour);
		}

		// repeate for the neighbours of neighbours recursively
		for (auto pPoint : neighbours) {
			for (auto pNeighbour : pPoint->neighbours) {
				if (pNeighbour->isReady() || pNeighbour->group_number != pPoint->group_number)
					continue;

				GroupHausdorfMixSearch(*pNeighbour, rGroupCounter, rFeedbackUpdater);
			}
		}
	}

	void GroupNeighbours(Point& rPoint, Clouds::Point3D baseXYZ, FeedbackUpdater& rFeedbackUpdater) const {
		const auto xyz = baseXYZ.cast<double>();
		for (auto pNeighbour : rPoint.neighbours) {
			if (pNeighbour->isReady())
				continue;

			const double distance = Math::CalcPointToPointDistance3D(xyz, pNeighbour->xyz.cast<double>());
			if (CheckColorGradient(rPoint.color, pNeighbour->color, distance))
				pNeighbour->group_number = rPoint.group_number;
		}
		rPoint.is_calculated = true;
		rFeedbackUpdater.count();
	}
	
	/*Function returns 'true' if the color gradient of given data converted to HSV are within tolerance and 'false' otherwise.*/	
	bool CheckColorGradient(Clouds::Color color1, Clouds::Color color2, double distance) const {
		const double R1 = color1.x();
		const double G1 = color1.y();
		const double B1 = color1.z();
		const double R2 = color2.x();
		const double G2 = color2.y();
		const double B2 = color2.z();

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
		if (abs(V1 - V2) > mMaxBrightnessGradient * distance)
			return false;
		if (S1 == 0.0 && S2 == 0.0)
			return true;
		if (abs(S1 - S2) > mMaxSaturationGradient * distance + mCoverageRate / (255.0 * (V1 + V2) / 2))
			return false;
		const double hueInaccuracy = mCoverageRate * (180.0 / (255.0 * (V1 + V2) / 2) + 180.0 / (255.0 * (S1 + S2) / 2));
		if (std::min(abs(H1 - H2), 360.0 - abs(H1 - H2)) > mMaxHueGradient * distance + hueInaccuracy)
			return false;
		return true;
	}

	/*Function finds groups with a number of elements less than minimum given by a user and marks them as ungrouped.*/
	void RemoveTooSmallGroups(Context& rContext, std::vector<Point>& rGraph) const {
		FeedbackUpdater feedbackUpdater(rContext, 2 * rGraph.size());
		std::map<int, int> repetitionNumbers;
		
		for (const auto& rPoint : rGraph) {
			++repetitionNumbers[rPoint.group_number];
			feedbackUpdater.count();
		}
		
		for (auto& rPoint: rGraph) {
			if (repetitionNumbers[rPoint.group_number] < mMinPointsGroup)
				rPoint.group_number = 0;
			feedbackUpdater.count();
		}
	}
	
	/*Function adds ungrouped points surrounded by neighbours assigned to groups to the locally most popular group.*/
	void FillHolesInClusters(Context& rContext, std::vector<Point>& rGraph) const {
		FeedbackUpdater feedbackUpdater(rContext, rGraph.size());

		// set initial states
		for (auto& rPoint : rGraph) {
			if (rPoint.group_number != 0 || rPoint.is_deleted) {
				rPoint.is_calculated = true;
				feedbackUpdater.count();
			}
			else
				rPoint.is_calculated = false;
		}

		bool anyChange = true;
		while (anyChange) {
			anyChange = false;
			std::vector<int> newGroupNumbers;
			newGroupNumbers.reserve(rGraph.size());

			for (auto& rPoint : rGraph) {
				if (rPoint.is_calculated) {
					newGroupNumbers.push_back(0);
					continue;
				}

				std::map<int, int> repetitionNumber;
				for (const auto& pNeighbour : rPoint.neighbours) {
					if (pNeighbour->group_number != 0)
						++repetitionNumber[pNeighbour->group_number];
				}

				auto maxGroup = std::max_element(repetitionNumber.begin(), repetitionNumber.end(), [](auto& a, auto& b) {
					return a.second < b.second;
					});

				// assign point to most popular group if it exist
				if (maxGroup != repetitionNumber.end() && repetitionNumber[0] < mJoinCoefficient * rPoint.neighbours.size()) {
					newGroupNumbers.push_back(maxGroup->first);
					rPoint.is_calculated = true;
					anyChange = true;
					feedbackUpdater.count();
				}
				else
					newGroupNumbers.push_back(0);
			}

			// write values to the structure at the end of each loop to guarantee isotropic addition of points
			for (int i = 0; i < rGraph.size(); ++i) {
				if (newGroupNumbers[i] != 0)
					rGraph[i].group_number = newGroupNumbers[i];
			}
		}
	}

	/*Function merge clusters with a relatively high number of points on common border.*/
	void MergeClusters(Context& rContext, std::vector<Point>& rGraph) {
		typedef std::pair<int, std::vector<Point*>> PointCluster;
		std::vector<PointCluster> clusters;

		for (auto& rPoint : rGraph) {
			if (rPoint.group_number == 0)
				continue;

			auto clusterItr = std::find_if(clusters.begin(), clusters.end(), [&rPoint](PointCluster& rCluster) {
				return rCluster.first == rPoint.group_number;
				});

			if (clusterItr == clusters.end())
				clusters.push_back(PointCluster(rPoint.group_number, std::vector<Point*>({ &rPoint })));
			else
				clusterItr->second.push_back(&rPoint);
		}

		std::sort(clusters.begin(), clusters.end(), [](PointCluster& a, PointCluster& b) {
			return a.second.size() < b.second.size();
			});

		FeedbackUpdater feedbackUpdater(rContext, clusters.size());

		for (auto itrSmaller = clusters.begin(), end = clusters.end(); itrSmaller != end; ++itrSmaller) {
			PointCluster& rSmallerCluster = *itrSmaller;
			int maxNeighbours = 0;

			std::for_each(itrSmaller + 1, clusters.end(), [&](PointCluster& bigger_cluster) {
				const int& rBiggerGroup = bigger_cluster.first;

				const int borderPoints = std::count_if(rSmallerCluster.second.begin(), rSmallerCluster.second.end(),
					[&](Point* pPoint) {
						return std::any_of(pPoint->neighbours.begin(), pPoint->neighbours.end(), [&](Point* pNeigbour) {
							return pNeigbour->group_number == rBiggerGroup;
							});
					});

				// add sizes of all parts merged into the smaller cluster
				const int size = std::accumulate(clusters.begin(), itrSmaller, rSmallerCluster.second.size(),
					[&](size_t init, PointCluster& cluster) {
						return init + (cluster.first == rSmallerCluster.first) ? cluster.second.size() : 0;
					});

				// compare border size with number of points in the smaller cluster circumference
				if (borderPoints < sqrt(4 * PI * size) * mMaxDistance / mMergeCoefficient || borderPoints < maxNeighbours)
					return;

				// merge the smaller cluster and all clusters merged into it, into the bigger cluster
				std::for_each(clusters.begin(), itrSmaller + 1, [&](PointCluster& rCluster) {
					if (rCluster.first != rSmallerCluster.first)
						return;

					for (auto pPoint : rCluster.second)
						pPoint->group_number = rBiggerGroup;

					rCluster.first = rBiggerGroup;
					});

				maxNeighbours = borderPoints;
				});

			feedbackUpdater.count();
		}
	}
	
	/*Function changes group numbers to subsequent numbers with '0' assigned to ungrouped points and write them to a given vector of layer values.*/
	void TidyUpGroupNumbers(Context& rContext, std::vector<StoredReal>& rClusterNumbers, std::vector<Point>& rGraph) {
		FeedbackUpdater feedbackUpdater(rContext, 2 * rGraph.size());

		std::map<int, int> repetitionNumbers;
		for (auto& rPoint : rGraph) {
			if (rPoint.group_number != 0)
				++repetitionNumbers[rPoint.group_number];
			feedbackUpdater.count();
		}

		std::vector<std::pair<const int, int>*> orderedGroups;
		for (auto& rGroupSize : repetitionNumbers)
			orderedGroups.push_back(&rGroupSize);

		std::sort(orderedGroups.begin(), orderedGroups.end(), [](std::pair<const int, int>*& a, std::pair<const int, int>*& b) {
			return a->second > b->second;
			});

		std::map<int, int> newGroupNumbers;
		for (int j = 0; j < orderedGroups.size(); ++j) {
			// ignore number '0' assigned to ungrouped points
			newGroupNumbers[orderedGroups[j]->first] = j + 1;
		}

		for (auto& rPoint : rGraph) {
			rPoint.group_number = newGroupNumbers[rPoint.group_number];
			rClusterNumbers.push_back(rPoint.group_number);
			feedbackUpdater.count();
		}
	}

	/*Function finds the best fitting cluster to the loudspeaker and write to the given layer '1.0' for points of this cluster and '0.0' for others.*/
	void FindTheSpeaker(Context& rContext, std::vector<StoredReal>& rIsSpeaker, std::vector<Point>& rGraph) {
		std::map<int, std::vector<Point*>> clusters;

		for (auto& rPoint : rGraph) {
			if (clusters.find(rPoint.group_number) == clusters.end())
				clusters[rPoint.group_number] = { &rPoint };
			else
				clusters[rPoint.group_number].push_back(&rPoint);
		}

		FeedbackUpdater feedbackUpdater(rContext, clusters.size(), true);

		std::map<int, double> fittingErrors;
		for (auto& rCluster : clusters) {
			// ignore ungrouped points
			fittingErrors[rCluster.first] = GetBoxFittingError(rContext, rCluster.second);
			feedbackUpdater.count();
		}

		auto speakerItr = std::min_element(fittingErrors.begin(), fittingErrors.end(), [](auto& rLeft, auto& rRight) {
			return rLeft.second < rRight.second;
			});

		for (auto& rPoint : rGraph)
			rIsSpeaker.push_back(rPoint.group_number == speakerItr->first ? 1.0 : 0.0);

		//// print the values of all fitting errors
		//for (auto& rError : fittingErrors)
		//	OGX_LINE.Msg(Info, std::to_wstring(rError.first) + L". " + std::to_wstring(rError.second));
	}

	/*Function calculates a box fitted to the borders of given points cluster. Returns sum of angles and dimensions relative errors and relative RMSE for this fitting.*/
	double GetBoxFittingError(Context& rContext, const std::vector<Point*>& rCluster) {
		Math::Point3D firstDiagonalEnd, secondDiagonalEnd;
		double maxDistanceBetweenPoints = 0.0;
		for (int i = 0; i < rCluster.size(); ++i) {
			const auto& rFirstPoint = rCluster[i]->xyz.cast<double>();

			for (int j = i + 1; j < rCluster.size(); ++j) {
				const auto& rSecondPoint = rCluster[j]->xyz.cast<double>();
				const double localDistance = Math::CalcPointToPointDistance3D(rFirstPoint, rSecondPoint);

				if (localDistance > maxDistanceBetweenPoints) {
					maxDistanceBetweenPoints = localDistance;
					firstDiagonalEnd = rFirstPoint;
					secondDiagonalEnd = rSecondPoint;
				}
			}
		}
		const Math::Point3D centerPoint = (firstDiagonalEnd + secondDiagonalEnd) / 2;
		const Math::Line3D longestDiagonal = Math::CalcLine3D(firstDiagonalEnd, secondDiagonalEnd);
		
		Math::Point3D farthestFromDiagonal;
		double maxDistanceFromDiagonal = 0.0;
		for (const auto pPoint : rCluster) {
			const auto& rPoint = pPoint->xyz.cast<double>();
			const double localDistance = Math::CalcPointToLineDistance3D(rPoint, longestDiagonal);

			if (localDistance > maxDistanceFromDiagonal) {
				maxDistanceFromDiagonal = localDistance;
				farthestFromDiagonal = rPoint;
			}
		}
		const Math::Plane3D sectionPlane = Math::CalcPlane3D(firstDiagonalEnd, secondDiagonalEnd, farthestFromDiagonal);
		
		Math::Point3D farthestFromSectionPlane;
		double maxDistanceFromSectionPlaneAndCenter = 0.0;
		for (const auto pPoint : rCluster) {
			const auto& rPoint = pPoint->xyz.cast<double>();
			double localDistance = Math::CalcPointToPlaneDistance3D(rPoint, sectionPlane);
			localDistance += Math::CalcPointToPointDistance3D(rPoint, centerPoint);

			if (localDistance > maxDistanceFromSectionPlaneAndCenter) {
				maxDistanceFromSectionPlaneAndCenter = localDistance;
				farthestFromSectionPlane = rPoint;
			}
		}

		const bool isFirstEndOn0YZ = Math::CalcPointToPointDistance3D(firstDiagonalEnd, farthestFromSectionPlane)
			< Math::CalcPointToPointDistance3D(secondDiagonalEnd, farthestFromSectionPlane);

		const Math::Point3D& vertex000 = isFirstEndOn0YZ ? firstDiagonalEnd : secondDiagonalEnd;
		const Math::Point3D& vertex111 = isFirstEndOn0YZ ? secondDiagonalEnd : firstDiagonalEnd;

		const bool isFarDiagOn0YZ = Math::CalcPointToPointDistance3D(farthestFromDiagonal, farthestFromSectionPlane)
			< Math::CalcPointToPointDistance3D(farthestFromDiagonal, vertex000);

		const Math::Point3D vertex011 = isFarDiagOn0YZ ? farthestFromDiagonal : 2 * centerPoint - farthestFromDiagonal;
		const Math::Point3D vertex100 = isFarDiagOn0YZ ? 2 * centerPoint - farthestFromDiagonal : farthestFromDiagonal;
		
		const bool isFarSecPlane001 = Math::CalcPointToPointDistance3D(farthestFromSectionPlane, vertex000)
			< Math::CalcPointToPointDistance3D(farthestFromSectionPlane, vertex011);

		const Math::Point3D vertex001 = isFarSecPlane001 ? farthestFromSectionPlane : vertex000 + vertex011 - farthestFromSectionPlane;
		const Math::Point3D vertex010 = isFarSecPlane001 ? vertex000 + vertex011 - farthestFromSectionPlane : farthestFromSectionPlane;
	
		const Math::Point3D vertex101 = 2 * centerPoint - vertex010;
		const Math::Point3D vertex110 = 2 * centerPoint - vertex001;

		const Math::Plane3D wall0YZ = Math::CalcPlane3D(vertex000, vertex010, vertex001);
		const Math::Plane3D wall1YZ = Math::CalcPlane3D(vertex111, vertex110, vertex101);
		const Math::Plane3D wallX0Z = Math::CalcPlane3D(vertex000, vertex001, vertex100);
		const Math::Plane3D wallX1Z = Math::CalcPlane3D(vertex111, vertex011, vertex110);
		const Math::Plane3D wallXY0 = Math::CalcPlane3D(vertex000, vertex010, vertex100);
		const Math::Plane3D wallXY1 = Math::CalcPlane3D(vertex111, vertex011, vertex101);
		
		const double boxFittingSquareError = std::accumulate(rCluster.begin(), rCluster.end(), 0.0,
			[&](double init, Point* pPoint) {
				const auto& rPoint = pPoint->xyz.cast<double>();
				const auto minPointWallDistance = std::min({
					Math::CalcPointToPlaneDistance3D(rPoint, wall0YZ),
					Math::CalcPointToPlaneDistance3D(rPoint, wall1YZ),
					Math::CalcPointToPlaneDistance3D(rPoint, wallX0Z),
					Math::CalcPointToPlaneDistance3D(rPoint, wallX1Z),
					Math::CalcPointToPlaneDistance3D(rPoint, wallXY0),
					Math::CalcPointToPlaneDistance3D(rPoint, wallXY1)
					});
				return init + pow(minPointWallDistance, 2);
			});
		const double rmse = sqrt(boxFittingSquareError / rCluster.size());
		const double maxDistance = Math::CalcPointToLineDistance3D(farthestFromDiagonal, longestDiagonal);

		double relativeAngleError = 0.0;
		const Math::Vector3D xVec(vertex100 - vertex000);
		const Math::Vector3D yVec(vertex010 - vertex000);
		const Math::Vector3D zVec(vertex001 - vertex000);
		relativeAngleError += abs(Math::CalcAngleBetweenTwoVectors(xVec, yVec) - PI / 2) / (PI / 2);
		relativeAngleError += abs(Math::CalcAngleBetweenTwoVectors(yVec, zVec) - PI / 2) / (PI / 2);
		relativeAngleError += abs(Math::CalcAngleBetweenTwoVectors(zVec, xVec) - PI / 2) / (PI / 2);
		
		double relativeDimentionError = 0.0;
		relativeDimentionError += abs(xVec.norm() - mSpeakerDimX) / mSpeakerDimX;
		relativeDimentionError += abs(yVec.norm() - mSpeakerDimY) / mSpeakerDimY;
		relativeDimentionError += abs(zVec.norm() - mSpeakerDimZ) / mSpeakerDimZ;

		const double combinedRelativeError = (rmse / maxDistance + relativeAngleError + relativeDimentionError) / 7.0;
		//// print the values of all fitting errors
		//OGX_LINE.Msg(Info, std::to_wstring(xVec.norm()) + L", " + std::to_wstring(yVvec.norm()) + L", "
		//	+ std::to_wstring(zVec.norm()) + L", Angle error: " + std::to_wstring(relativeAngleError) + L", RMSE: "
		//	+ std::to_wstring(rmse / maxDistance));
		
		return combinedRelativeError;
	}
};

OGX_EXPORT_METHOD(LoudspeakerSearching)
