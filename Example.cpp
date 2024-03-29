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
#include <queue>
#include <stack>

using namespace ogx;
using namespace Data;

struct LoudspeakerSearching : public Plugin::EasyMethod {
	/*HARDCODED PARAMETERS AND DATA*/
	const double HAUSDORF_THRESHOLD = 50000;	// number of points which change version of Hausdorf algorithm

	/**
	* Store data of a 3D point with vector of pointers to its neighbours, which make it an element of a graph
	* structure.\*
	*/
	struct Point {
		Clouds::Point3D xyz;
		Clouds::Color color;
		std::vector<Point*> neighbours;
		unsigned int groupNumber = 0;

		Point(Clouds::Point3D coordinates, Clouds::Color color) : xyz(coordinates), color(color) {}
		bool isReady() const {
			return mCalculated || mDeleted;
		}
		bool isCalculated() const {
			return mCalculated;
		}
		bool isDeleted() const {
			return mDeleted;
		}
		void setCalculated() {
			mCalculated = true;
		}
		void setUncalculated() {
			mCalculated = false;
		}
		void setDeleted() {
			mDeleted = true;
		}
		bool isUngrouped() const {
			return groupNumber == 0;
		}

	private:
		bool mCalculated = false;
		bool mDeleted = false;
	};

	class FeedbackUpdater {
	public:
		FeedbackUpdater() = delete;
		FeedbackUpdater(Context& context, unsigned int stepNumber, bool decreasingGrowth = false) :
			mContext(context),
			mStepNumber(stepNumber),
			mDecreasingGrowth(decreasingGrowth)
		{}

		void count() {
			update(1);
		}
		void update(int numberOfCounts) {
			const double last_ratio = 1.0 * mCounter / mStepNumber;
			mCounter += numberOfCounts;
			const double ratio = 1.0 * mCounter / mStepNumber;
			const double percentage = mDecreasingGrowth ? 100 * (1.0 - pow(1.0 - ratio, 10.0)) : 100 * ratio;
			const double last_percentage = mDecreasingGrowth ? 100 * (1.0 - pow(1.0 - last_ratio, 10.0)) : 100 * last_ratio;
			if (floor(percentage) > floor(last_percentage))
				mContext.Feedback().Update(0.01 * percentage);			
		}

	private:
		Context& mContext;
		unsigned int mCounter = 0;
		const unsigned int mStepNumber;
		const bool mDecreasingGrowth;
	};

	/**
	* Vector of graph elements with points parameters and pointers to their neighbours inside a spherical search
	* kernel. It is used to save the time for searching neighbours in each operation. It might be made one time for
	* many layers with diffrent parameters such as:
	* max color gradients in HSV
	* min points in a group
	* min neighbous of a point
	* dimensions of the loudspeaker
	* if only search kernel radius is not	changed and operations are done for the same cloud.\*
	*/
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
	/// additional tolerance of hue and saturation differences for mean HSV=(*,1,1),
	/// 10 times smaller for HSV=(*,10,10) etc.\*
	float mCoverageRate = 0.4;
	/// maximum ratio of ungrouped point among the neighbours to join the ungrouped point to any of the clusters
	float mJoinCoefficient = 0.8;
	/// the higher it is the more clusters will be merged
	float mMergeCoefficient = 0.02;
	float mSpeakerDimX = 0.18;
	float mSpeakerDimY = 0.16;
	float mSpeakerDimZ = 0.12;

	// constructor
	LoudspeakerSearching() : EasyMethod(L"Mateusz Frejlich", L"The program segments the cloud and recognizes the \
		loudspeaker inside it.\nThe default values of all parameters are best fitted to the tested cloud.")
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
		Clouds::ICloud* pOriginalCloud = GetValidCloud(context, mNodeId);
		OGX_LINE.Msg(Info, L"All data is correct. Processing has started.");

		ResourceID newNodeId;
		Clouds::ICloud* pCloud = CreateRandomlyReducedCloud(context, pOriginalCloud, newNodeId, mPointsToRemoveRatio);
		OGX_LINE.Msg(Info, L"Reduced cloud has been done.");

		BuildPointsGraphStructure(context, newNodeId, mPointsGraphStructure, mMaxDistance);
		OGX_LINE.Msg(Info, L"Points structure has been done.");

		RemoveNoise(context, mPointsGraphStructure, mMinNeighbours, mNeighboursToSmooth);
		FilterColorsByMedian(context, mPointsGraphStructure, mMedianFlterSize);
		OGX_LINE.Msg(Info, L"Noise filtration has been done.");

		// visualize the cloud smoothing
		Clouds::PointsRange range;
		pCloud->GetAccess().GetAllPoints(range);
		CopyGraphStructuctureToCloud(mPointsGraphStructure, range);

		GroupHausdorf(context, mPointsGraphStructure, HAUSDORF_THRESHOLD);
		OGX_LINE.Msg(Info, L"Segmentation by Hausdorf algorithm has been done.");

		FillHolesInClusters(context, mPointsGraphStructure);
		RemoveTooSmallGroups(context, mPointsGraphStructure);
		OGX_LINE.Msg(Info, L"Removing too small groups has been done.");

		FillHolesInClusters(context, mPointsGraphStructure);
		OGX_LINE.Msg(Info, L"Filling holes in clusters has been done.");

		MergeClusters(context, mPointsGraphStructure);
		OGX_LINE.Msg(Info, L"Merging clusters has been done.");

		TidyUpGroupNumbers(context, mPointsGraphStructure);
		OGX_LINE.Msg(Info, L"The ordering of group numbers has been done.");

		CreateGroupNumbersLayer(pCloud, range, mPointsGraphStructure);
		OGX_LINE.Msg(Info, L"Whole segmentation process has been finished.");

		std::vector<StoredReal> isLoundspeaker;
		isLoundspeaker.reserve(mPointsGraphStructure.size());
		FindTheSpeaker(context, isLoundspeaker, mPointsGraphStructure);
		OGX_LINE.Msg(Info, L"Loudspeaker identification has been done.");

		const String newLayerName = L"Loudspeaker";
		const auto newLayers = pCloud->FindLayers(newLayerName);
		auto* pSpeakerLayer = newLayers.empty() ? pCloud->CreateLayer(newLayerName, 0.0) : newLayers[0];
		range.SetLayerVals(isLoundspeaker, *pSpeakerLayer);
		OGX_LINE.Msg(Info, L"Processing has been finished. No errors.");
	}

	/// Reads data of all points from given node, check their correctness and returns a cloud of these points.
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
	/**
	* Makes a new node with given ID, assign to them the cloud with random points from a given cloud with ratio given
	* by a user.\*
	*/
	Clouds::ICloud* CreateRandomlyReducedCloud(Context& rContext, Clouds::ICloud*& rCloud, ResourceID& rNodeId,
		float pointsToRemoveRatio) const {
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
	/**
	* Finds neighbours for all the points from the cloud in the given node and makes the lists of graph elements
	* assigned to the vector in this structure.\*
	*/
	void BuildPointsGraphStructure(Context& rContext, ResourceID childNodeId, std::vector<Point>& rGraph,
		float maxDistance) {
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
		
		using Cluster = std::vector<Point*>;
		using CubeVector = std::vector<std::vector<std::vector<Cluster>>>;
		class Dimention {
		public:
			Dimention() = delete;
			Dimention(std::vector<Point>& rGraph, double maxDistance, std::function<double(Point&)> dimentionGetter) :
				getDimention(dimentionGetter),
				offset(calculateOffset(rGraph)),
				length(calculateDimention(rGraph)),
				divisions(length / maxDistance),
				indexShifts(calculateIndexShifts(rGraph, maxDistance))
			{}
			int calculateIndex(Point& rPoint) const {
				return std::min(static_cast<int>((getDimention(rPoint) - offset) / length * divisions), divisions - 1);
			}
			void forEachNeighbour(std::function<void(int, int)> fun) {
				for (int baseInd = 0; baseInd < divisions; ++baseInd) {
					for (const int& shift : indexShifts) {
						const int neighbourInd = baseInd + shift;
						if (neighbourInd < 0 || neighbourInd >= divisions)
							break;

						fun(baseInd, neighbourInd);
					}
				}
			}
			size_t size() const {
				return divisions;
			}

		private:
			double calculateOffset(std::vector<Point>& rGraph) const {
				return getDimention(*std::min_element(rGraph.begin(), rGraph.end(), [&](Point& left, Point& right) {
					return getDimention(left) < getDimention(right);
					}));
			}
			double calculateDimention(std::vector<Point>& rGraph) const {
				return getDimention(*std::max_element(rGraph.begin(), rGraph.end(), [&](Point& left, Point& right) {
					return getDimention(left) < getDimention(right);
					})) - offset;
			}
			std::vector<int> calculateIndexShifts(std::vector<Point>& rGraph, double maxDistance) const {
				std::vector<int> shifts({ 0 });
				for (int i = 1; (i - 1) * length / divisions < maxDistance; ++i) {
					shifts.push_back(i);
					shifts.push_back(-i);
				}
				return shifts;
			}
			const std::function<double(Point&)> getDimention;
			const double offset;
			const double length;
			const int divisions;
			const std::vector<int> indexShifts;
		};

		Dimention x(rGraph, maxDistance, [](Point& rPoint)->double {return rPoint.xyz.x(); });
		Dimention y(rGraph, maxDistance, [](Point& rPoint)->double {return rPoint.xyz.y(); });
		Dimention z(rGraph, maxDistance, [](Point& rPoint)->double {return rPoint.xyz.z(); });

		CubeVector clusters(x.size(), std::vector<std::vector<Cluster>>(y.size(), std::vector<Cluster>(z.size())));

		for (auto& rPoint : rGraph) {
			const int indX = x.calculateIndex(rPoint);
			const int indY = y.calculateIndex(rPoint);
			const int indZ = z.calculateIndex(rPoint);
			auto& cluster = clusters[indX][indY][indZ];
			cluster.push_back(&rPoint);
		}
		
		FeedbackUpdater feedbackUpdater(rContext, 2 * rangeSize);

		x.forEachNeighbour([&](int baseX, int neighbourX) {
			y.forEachNeighbour([&](int baseY, int neighbourY) {
				z.forEachNeighbour([&](int baseZ, int neighbourZ) {
					for (auto* pPoint : clusters[baseX][baseY][baseZ]) {
						const auto& first = pPoint->xyz.cast<double>();

						for (auto* pNeighbour : clusters[neighbourX][neighbourY][neighbourZ]) {
							const auto& second = pNeighbour->xyz.cast<double>();

							if (pPoint == pNeighbour)
								feedbackUpdater.count();
							else if (Math::CalcPointToPointDistance3D(first, second) <= maxDistance)
								pPoint->neighbours.push_back(pNeighbour);
						}
					}
					});
				});
			});

		for (auto& rPoint : rGraph) {
			const auto rPointXYZ = rPoint.xyz.cast<double>();
			std::sort(rPoint.neighbours.begin(), rPoint.neighbours.end(), [&](Point* pA, Point* pB) {
				return Math::CalcPointToPointDistance3D(pA->xyz.cast<double>(), rPointXYZ)
					< Math::CalcPointToPointDistance3D(pB->xyz.cast<double>(), rPointXYZ);
				});
			feedbackUpdater.count();
		}
	}
	/// Copies all points' data from the graph structure to the given range of the cloud.
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
	/**
	* Marks lonely points as 'deleted' and smooths the cloud by projecting all points in the structure to a local
	* fitted sphere or plane.\*
	*/
	void RemoveNoise(Context& rContext, std::vector<Point>& rGraph, int minNeighbours, int neighboursToSmooth) {
		const int numberOfPoints = rGraph.size();
		FeedbackUpdater feedbackUpdater(rContext, numberOfPoints * 2);
		for (auto& rPoint : rGraph) {
			if (rPoint.neighbours.size() < minNeighbours)
				rPoint.setDeleted();
		}

		std::vector<Clouds::Point3D> newCoordinates;
		newCoordinates.reserve(numberOfPoints);
		std::transform(rGraph.begin(), rGraph.end(), std::back_inserter(newCoordinates), [&](Point& rPoint) {
			const auto& rNeighbours = rPoint.neighbours;

			if (rNeighbours.size() < neighboursToSmooth || rPoint.isDeleted())
				return rPoint.xyz;

			auto basePoint = rPoint.xyz.cast<double>();
			std::vector<decltype(basePoint)> nXYZ = { basePoint };
			nXYZ.reserve(neighboursToSmooth + 1);
			std::transform(rNeighbours.begin(), rNeighbours.begin() + neighboursToSmooth, std::back_inserter(nXYZ),
				[](Point* nItr) {
					return nItr->xyz.cast<double>();
				});
			const auto sphere = Math::CalcBestSphere3D(nXYZ.begin(), nXYZ.end());
			const auto plane = Math::CalcBestPlane3D(nXYZ.begin(), nXYZ.end());

			const double sqSumSphere = std::accumulate(nXYZ.begin(), nXYZ.end(), 0.0, [&](double init, auto& rPoint) {
				const auto pointOntoSphere = sphere.Project(rPoint).cast<double>();
				return init + pow(Math::CalcPointToPointDistance3D(pointOntoSphere, rPoint), 2);
				});
			const double sqSumPlane = std::accumulate(nXYZ.begin(), nXYZ.end(), 0.0, [&](double init, auto& rPoint) {
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
			feedbackUpdater.count();
		}
	}
	/// Filters all points using the median colors of neighbours.
	void FilterColorsByMedian(Context& rContext, std::vector<Point>& rGraph, int medianFilterSize) {
		const int numberOfPoints = rGraph.size();
		FeedbackUpdater feedbackUpdater(rContext, numberOfPoints);
		std::vector<Clouds::Color> colors;
		colors.reserve(numberOfPoints);

		std::transform(rGraph.begin(), rGraph.end(), colors.begin(), [&](const Point& rPoint) {
			if (rPoint.isDeleted() || rPoint.neighbours.size() < medianFilterSize) {
				feedbackUpdater.count();
				return rPoint.color;
			}

			using WeightedColor = std::pair<double, double>;
			std::vector<WeightedColor> red, green, blue;
			const auto begin = rPoint.neighbours.begin();
			const auto pointXYZ = rPoint.xyz.cast<double>();
			const auto halfWeight = 0.5 * std::accumulate(begin, begin + medianFilterSize, 0.0,
				[&](double init, Point* pNeighbour) {
					const auto& rNeighborXYZ = pNeighbour->xyz.cast<double>();
					const double weight = pow(Math::CalcPointToPointDistance3D(pointXYZ, rNeighborXYZ), -2.0);
					red.push_back(WeightedColor(pNeighbour->color.x(), weight));
					green.push_back(WeightedColor(pNeighbour->color.y(), weight));
					blue.push_back(WeightedColor(pNeighbour->color.z(), weight));
					return init + weight;
				});

			std::vector<std::vector<WeightedColor>*> rgb = { &red,&green,&blue };
			std::vector<double> medianColor;
			std::transform(rgb.begin(), rgb.end(), std::back_inserter(medianColor),
				[&](std::vector<WeightedColor>* pChannel) {
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
	/**
	* Groups the points into clusters choosing one of the implemented Hausdorf algorithm depending on the number
	* of points.\*
	*/
	void GroupHausdorf(Context& rContext, std::vector<Point>& rGraph, const double hausdorfThreshold) {
		// set initial values
		for (auto& rPoint : rGraph) {
			rPoint.setUncalculated();
			rPoint.groupNumber = 0;
		}
		FeedbackUpdater feedbackUpdater(rContext, rGraph.size());
		int groupCounter = 0;

		// realize Hausdorf algorithm for each point
		for (auto& rPoint : rGraph) {
			if (!rPoint.isUngrouped() || rPoint.isDeleted())
				continue;
			++groupCounter;

			// check if recursion with this number of points results stack overflow
			if (rGraph.size() < hausdorfThreshold)
				GroupHausdorfDFS(rPoint, groupCounter, feedbackUpdater);
			else
				GroupHausdorfMixSearch(rPoint, groupCounter, feedbackUpdater);
		}
	}
	/**
	* Groups points into a cluster with given number using Hausdorf's algorithm, starting from a given point with depth
	* first searching of the graph data structure.\*
	*/
	void GroupHausdorfDFS(Point& rStartPoint, const int& rGroupCounter, FeedbackUpdater& rFeedbackUpdater) {
		std::stack<Point*> neighbours({ &rStartPoint });
		rStartPoint.groupNumber = rGroupCounter;

		while (!neighbours.empty()) {
			Point* pBasePoint = neighbours.top();
			neighbours.pop();
			if (pBasePoint->isReady() || pBasePoint->groupNumber != rGroupCounter)
				continue;

			for (auto pNeighbour : pBasePoint->neighbours)
				neighbours.push(pNeighbour);

			GroupNeighbours(pBasePoint, pBasePoint->xyz, rFeedbackUpdater);
		}
	}
	/**
	* Groups points into a cluster with given number using Hausdorf's algorithm, starting from a given point with
	* breadth first searching of the graph data structure.\*
	*/
	void GroupHausdorfBFS(Point& rStartPoint, const int& rGroupCounter, FeedbackUpdater& rFeedbackUpdater) {
		std::queue<Point*> neighbours({ &rStartPoint });
		rStartPoint.groupNumber = rGroupCounter;

		while (!neighbours.empty()) {
			Point* pBasePoint = neighbours.front();
			neighbours.pop();
			if (pBasePoint->isReady() || pBasePoint->groupNumber != rGroupCounter)
				continue;

			for (auto pNeighbour : pBasePoint->neighbours)
				neighbours.push(pNeighbour);

			GroupNeighbours(pBasePoint, pBasePoint->xyz, rFeedbackUpdater);
		}
	}
	/**
	* Groups points into a cluster with given number, starting from a given point using recursive Hausdorf's algorithm
	* with 2 steps in one iteration of the function. This modification is less clear but prevents stack overflow.\*
	*/
	void GroupHausdorfMixSearch(Point& rStartPoint, const int& rGroupCounter, FeedbackUpdater& rFeedbackUpdater) {
		// relize the algorithm for the point
		rStartPoint.groupNumber = rGroupCounter;
		GroupNeighbours(&rStartPoint, rStartPoint.xyz, rFeedbackUpdater);

		// repeate for the neighbours in a loop
		std::vector<Point*> neighbours;
		for (auto pNeighbour : rStartPoint.neighbours) {
			if (pNeighbour->isReady() || pNeighbour->groupNumber != rStartPoint.groupNumber)
				continue;

			GroupNeighbours(pNeighbour, rStartPoint.xyz, rFeedbackUpdater);
			neighbours.push_back(pNeighbour);
		}

		// repeate for the neighbours of neighbours recursively
		for (auto pPoint : neighbours) {
			for (auto pNeighbour : pPoint->neighbours) {
				if (pNeighbour->isReady() || pNeighbour->groupNumber != pPoint->groupNumber)
					continue;

				GroupHausdorfMixSearch(*pNeighbour, rGroupCounter, rFeedbackUpdater);
			}
		}
	}
	/// Assigns group numbers to the neighbours of given point according to Hausdorf algorithm.
	void GroupNeighbours(Point* pPoint, Clouds::Point3D baseXYZ, FeedbackUpdater& rFeedbackUpdater) const {
		const auto xyz = baseXYZ.cast<double>();
		for (auto pNeighbour : pPoint->neighbours) {
			if (pNeighbour->isReady())
				continue;

			const double distance = Math::CalcPointToPointDistance3D(xyz, pNeighbour->xyz.cast<double>());
			if (CheckColorGradient(pPoint->color, pNeighbour->color, distance))
				pNeighbour->groupNumber = pPoint->groupNumber;
		}
		pPoint->setCalculated();
		rFeedbackUpdater.count();
	}
	/// \return 'true' if the color gradient of given data converted to HSV are within tolerance and 'false' otherwise.
	bool CheckColorGradient(Clouds::Color color1, Clouds::Color color2, double distance) const {
		auto rgb2hsv = [](const Clouds::Color& color, double& rHue, double& rSaturation, double& rValue) {
			const double R = color.x();
			const double G = color.y();
			const double B = color.z();
			rValue = std::max({ R, G, B }) / 255.0;
			rSaturation = rValue == 0.0 ? 0.0 : 1.0 - std::min({ R, G, B }) / 255.0 / rValue;
			double aux = 180.0 / PI * acos((R - (G + B) / 2.0) / sqrt(R * R + G * G + B * B - R * G - R * B - G * B));
			rHue = G >= B ? aux : 360.0 - aux;
			};

		double V1, S1, H1, V2, S2, H2;
		rgb2hsv(color1, H1, S1, V1);
		rgb2hsv(color2, H2, S2, V2);

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
	/// Finds groups with a number of elements less than minimum given by a user and marks them as ungrouped.
	void RemoveTooSmallGroups(Context& rContext, std::vector<Point>& rGraph) const {
		FeedbackUpdater feedbackUpdater(rContext, 2 * rGraph.size());
		std::map<int, int> repetitionNumbers;
		
		for (const auto& rPoint : rGraph) {
			++repetitionNumbers[rPoint.groupNumber];
			feedbackUpdater.count();
		}
		
		for (auto& rPoint: rGraph) {
			if (repetitionNumbers[rPoint.groupNumber] < mMinPointsGroup)
				rPoint.groupNumber = 0;
			feedbackUpdater.count();
		}
	}
	/// Adds ungrouped points surrounded by neighbours assigned to groups to the locally most popular group.
	void FillHolesInClusters(Context& rContext, std::vector<Point>& rGraph) const {
		FeedbackUpdater feedbackUpdater(rContext, rGraph.size());

		// set initial states
		for (auto& rPoint : rGraph) {
			if (!rPoint.isUngrouped()) {
				rPoint.setCalculated();
				feedbackUpdater.count();
			}
			else
				rPoint.setUncalculated();
		}

		bool anyChange = true;
		while (anyChange) {
			anyChange = false;
			std::vector<int> newGroupNumbers;
			newGroupNumbers.reserve(rGraph.size());

			std::transform(rGraph.begin(), rGraph.end(), std::back_inserter(newGroupNumbers), [&](Point& rPoint) {
				if (rPoint.isReady())
					return 0;

				std::map<int, int> repetitionNumber;
				for (const auto& pNeighbour : rPoint.neighbours) {
					if (!pNeighbour->isUngrouped())
						++repetitionNumber[pNeighbour->groupNumber];
				}

				auto maxGroup = std::max_element(repetitionNumber.begin(), repetitionNumber.end(), [](auto& rA, auto& rB) {
					return rA.second < rB.second;
					});

				if (maxGroup == repetitionNumber.end())
					return 0;

				if (repetitionNumber[0] >= mJoinCoefficient * rPoint.neighbours.size())
					return 0;

				rPoint.setCalculated();
				feedbackUpdater.count();
				anyChange = true;
				return maxGroup->first;
				});

			// write values to the structure to guarantee isotropic addition of points
			for (int i = 0; i < rGraph.size(); ++i) {
				if (newGroupNumbers[i] != 0)
					rGraph[i].groupNumber = newGroupNumbers[i];
			}
		}
	}
	/// Merges clusters with a relatively high number of points on common border.
	void MergeClusters(Context& rContext, std::vector<Point>& rGraph) {
		using PointCluster = std::pair<int, std::vector<Point*>>;
		std::vector<PointCluster> clusters;

		for (auto& rPoint : rGraph) {
			if (rPoint.isUngrouped())
				continue;

			auto clusterItr = std::find_if(clusters.begin(), clusters.end(), [&rPoint](PointCluster& rCluster) {
				return rCluster.first == rPoint.groupNumber;
				});

			if (clusterItr == clusters.end())
				clusters.push_back(PointCluster(rPoint.groupNumber, std::vector<Point*>({ &rPoint })));
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
							return pNeigbour->groupNumber == rBiggerGroup;
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
						pPoint->groupNumber = rBiggerGroup;

					rCluster.first = rBiggerGroup;
					});

				maxNeighbours = borderPoints;
				});

			feedbackUpdater.count();
		}
	}
	/// Changes group numbers to subsequent numbers with '0' assigned to ungrouped points.
	void TidyUpGroupNumbers(Context& rContext, std::vector<Point>& rGraph) {
		FeedbackUpdater feedbackUpdater(rContext, 4 * rGraph.size());

		std::set<int> groupNumbers;
		for (auto& rPoint : rGraph) {
			groupNumbers.insert(rPoint.groupNumber);
			feedbackUpdater.count();
		}

		std::map<int, int> newGroupNumbers;
		int groupCounter = 0;
		for (const int& rNumber : groupNumbers) 
			newGroupNumbers[rNumber] = groupCounter++;

		for (auto& rPoint : rGraph) {
			rPoint.groupNumber = newGroupNumbers[rPoint.groupNumber];
			feedbackUpdater.count();
		}

		struct GroupNode {
			int number;
			std::vector<GroupNode*> neighbours;
			GroupNode(int number) : number(number) {}
		};

		std::map<int, std::set<int>> neighbourGroups;
		for (auto& rPoint : rGraph) {
			for (auto* pNeighbour : rPoint.neighbours)
				neighbourGroups[rPoint.groupNumber].insert(pNeighbour->groupNumber);

			feedbackUpdater.count();
		}

		std::vector<GroupNode> groupGraph;
		for (int i = 0; i < groupCounter; ++i)
			groupGraph.push_back(GroupNode(i));

		for (auto& rPair : neighbourGroups) {
			for (auto& rNeighbourNumber : rPair.second) {
				// ignore number '0' assigned to ungrouped points
				if (rPair.first != rNeighbourNumber && rNeighbourNumber != 0)
					groupGraph[rPair.first].neighbours.push_back(&groupGraph[rNeighbourNumber]);
			}
		}

		const auto& calculateGroupMatchQuality = [](std::vector<GroupNode>& rGroupGraph) -> int {
			return std::accumulate(rGroupGraph.begin() + 1, rGroupGraph.end(), 0, [](int init, GroupNode& rNode) {
				const auto& rNeighbours = rNode.neighbours;
				return init + std::accumulate(rNeighbours.begin(), rNeighbours.end(), 0, [&](int init, GroupNode* pNeighbour) {
					return init + std::abs(pNeighbour->number - rNode.number);
					});
				});
			};

		bool anyChange = true;
		while (anyChange) {
			anyChange = false;
			for (int i = 1; i < groupGraph.size(); ++i) {
				for (int j = i + 1; j < groupGraph.size(); ++j) {
					const int matchQualityBefore = calculateGroupMatchQuality(groupGraph);
					std::swap(groupGraph[i].number, groupGraph[j].number);
					const int matchQualityAfter = calculateGroupMatchQuality(groupGraph);

					if (matchQualityBefore >= matchQualityAfter)
						std::swap(groupGraph[i].number, groupGraph[j].number);
					else
						anyChange = true;
				}
			}
		}

		for (auto& rPoint : rGraph) {
			rPoint.groupNumber = groupGraph[rPoint.groupNumber].number;
			feedbackUpdater.count();
		}
	}
	/// Creates new layer and writes as values the group numbers of all points.
	void CreateGroupNumbersLayer(Clouds::ICloud* pCloud, Clouds::PointsRange& range, std::vector<Point>& rGraph) {
		const String layerName = L"Segmentation";
		auto layers = pCloud->FindLayers(layerName);
		auto* pSegmentsLayer = layers.empty() ? pCloud->CreateLayer(layerName, 0.0) : layers[0];
		std::vector<StoredReal> segmentNumbers;

		segmentNumbers.reserve(rGraph.size());
		for (auto& rPoint : rGraph)
			segmentNumbers.push_back(rPoint.groupNumber);

		range.SetLayerVals(segmentNumbers, *pSegmentsLayer);
	}
	/**
	* Finds the best fitting cluster to the loudspeaker and write to the given layer '1.0' for points of this cluster
	* and '0.0' for others.\*
	*/
	void FindTheSpeaker(Context& rContext, std::vector<StoredReal>& rIsSpeaker, std::vector<Point>& rGraph) {
		std::map<int, std::vector<Point*>> clusters;

		for (auto& rPoint : rGraph) {
			if (clusters.find(rPoint.groupNumber) == clusters.end())
				clusters[rPoint.groupNumber] = { &rPoint };
			else
				clusters[rPoint.groupNumber].push_back(&rPoint);
		}

		FeedbackUpdater feedbackUpdater(rContext, rGraph.size());

		std::map<int, double> fittingErrors;
		for (auto& rCluster : clusters) {
			// ignore ungrouped points
			if (rCluster.first != 0)
				fittingErrors[rCluster.first] = CalculateBoxFittingError(rCluster.second);

			feedbackUpdater.update(rCluster.second.size());
		}

		auto speakerItr = std::min_element(fittingErrors.begin(), fittingErrors.end(), [](auto& rLeft, auto& rRight) {
			return rLeft.second < rRight.second;
			});

		for (auto& rPoint : rGraph)
			rIsSpeaker.push_back(rPoint.groupNumber == speakerItr->first ? 1.0 : 0.0);

		//// print the values of all fitting errors
		//for (auto& rError : fittingErrors)
		//	OGX_LINE.Msg(Info, std::to_wstring(rError.first) + L". " + std::to_wstring(rError.second));
	}
	/**
	* Calculates a box fitted to the borders of given points cluster.
	* \return Sum of angles and dimensions relative errors and relative RMSE for this fitting.
	*/
	double CalculateBoxFittingError(const std::vector<Point*>& rCluster) {
		const auto maxDiagonalEnds = CalculateMaxDiagonalEnds(rCluster);
		const Math::Point3D centerPoint = (maxDiagonalEnds.first + maxDiagonalEnds.second) / 2;

		const Math::Line3D longestDiagonal = Math::CalcLine3D(maxDiagonalEnds.first, maxDiagonalEnds.second);
		const Math::Point3D farthestFromDiagonal = CalculateFarthestPointFromLine(rCluster, longestDiagonal);

		auto sectionPlane = Math::CalcPlane3D(maxDiagonalEnds.first, maxDiagonalEnds.second, farthestFromDiagonal);
		auto farthestSectionPlaneCorner = CalculateFarthestPointFromPlaneAndPoint(rCluster, sectionPlane, centerPoint);

		std::vector<Math::Point3D> mainVertecies({
			maxDiagonalEnds.first,
			maxDiagonalEnds.second,
			farthestFromDiagonal,
			2 * centerPoint - farthestFromDiagonal
			});
		std::sort(mainVertecies.begin(), mainVertecies.end(), [&](Math::Point3D& left, Math::Point3D& right) {
			return Math::CalcPointToPointDistance3D(left, farthestSectionPlaneCorner)
				< Math::CalcPointToPointDistance3D(right, farthestSectionPlaneCorner);
			});

		const Math::Point3D& vertex000 = farthestSectionPlaneCorner;
		const Math::Point3D& vertex001 = mainVertecies[0];
		const Math::Point3D& vertex010 = mainVertecies[1];
		const Math::Point3D& vertex101 = mainVertecies[2];
		const Math::Point3D& vertex110 = mainVertecies[3];
		const Math::Point3D vertex111 = 2 * centerPoint - vertex000;

		const std::vector<Math::Plane3D> walls({
			Math::CalcPlane3D(vertex000, vertex010, vertex001),
			Math::CalcPlane3D(vertex111, vertex110, vertex101),
			Math::CalcPlane3D(vertex000, vertex001, vertex101),
			Math::CalcPlane3D(vertex111, vertex010, vertex110),
			Math::CalcPlane3D(vertex000, vertex010, vertex110),
			Math::CalcPlane3D(vertex111, vertex001, vertex101),
			});

		const double rmse = CalculateWallsFittingRMSE(rCluster, walls);
		const double maxDistance = Math::CalcPointToLineDistance3D(farthestFromDiagonal, longestDiagonal);

		double relativeAngleError = 0.0;
		const Math::Vector3D xVec(vertex110 - vertex010);
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

	std::pair<Math::Point3D, Math::Point3D> CalculateMaxDiagonalEnds(const std::vector<Point*>& rCluster) {
		const auto xMinMax = std::minmax_element(rCluster.begin(), rCluster.end(), [](Point* pLeft, Point* pRight) {
			return pLeft->xyz.x() < pRight->xyz.x();
			});
		const auto yMinMax = std::minmax_element(rCluster.begin(), rCluster.end(), [](Point* pLeft, Point* pRight) {
			return pLeft->xyz.y() < pRight->xyz.y();
			});
		const auto zMinMax = std::minmax_element(rCluster.begin(), rCluster.end(), [](Point* pLeft, Point* pRight) {
			return pLeft->xyz.z() < pRight->xyz.z();
			});

		std::vector<Point*> borderPoints = {
			*xMinMax.first, *xMinMax.second, *yMinMax.first, *yMinMax.second, *zMinMax.first, *zMinMax.second
		};

		Math::Point3D scopeCenter = Math::Point3D(
			((*xMinMax.first)->xyz.x() + (*xMinMax.second)->xyz.x()) / 2,
			((*yMinMax.first)->xyz.y() + (*yMinMax.second)->xyz.y()) / 2,
			((*zMinMax.first)->xyz.z() + (*zMinMax.second)->xyz.z()) / 2
		);

		auto& rFarthestToCenter = **std::max_element(rCluster.begin(), rCluster.end(), [&](Point* pA, Point* pB) {
			return Math::CalcPointToPointDistance3D(pA->xyz.cast<double>(), scopeCenter)
				< Math::CalcPointToPointDistance3D(pB->xyz.cast<double>(), scopeCenter);
			});
		auto maxDistanceToCenter = Math::CalcPointToPointDistance3D(rFarthestToCenter.xyz.cast<double>(), scopeCenter);

		Math::Point3D firstDiagonalEnd, secondDiagonalEnd;
		double maxDistanceBetweenPoints = 0.0;
		for (int i = 0; i < borderPoints.size(); ++i) {
			const auto& rFirstPoint = borderPoints[i]->xyz.cast<double>();

			for (int j = i + 1; j < borderPoints.size(); ++j) {
				const auto& rSecondPoint = borderPoints[j]->xyz.cast<double>();
				const double localDistance = Math::CalcPointToPointDistance3D(rFirstPoint, rSecondPoint);

				if (localDistance > maxDistanceBetweenPoints) {
					maxDistanceBetweenPoints = localDistance;
					firstDiagonalEnd = rFirstPoint;
					secondDiagonalEnd = rSecondPoint;
				}
			}
		}

		std::vector<Point*> outsideSphere;
		outsideSphere.reserve(rCluster.size());
		std::copy_if(rCluster.begin(), rCluster.end(), std::back_inserter(outsideSphere), [&](Point* pPoint) {
				const auto& rPointXYZ = pPoint->xyz.cast<double>();
				const double distanceToCenter = Math::CalcPointToPointDistance3D(rPointXYZ, scopeCenter);
				return (distanceToCenter + maxDistanceToCenter >= maxDistanceBetweenPoints);
			});

		for (int i = 0; i < outsideSphere.size(); ++i) {
			const auto& rFirstPoint = outsideSphere[i]->xyz.cast<double>();

			for (int j = i + 1; j < outsideSphere.size(); ++j) {
				const auto& rSecondPoint = outsideSphere[j]->xyz.cast<double>();
				const double localDistance = Math::CalcPointToPointDistance3D(rFirstPoint, rSecondPoint);

				if (localDistance > maxDistanceBetweenPoints) {
					maxDistanceBetweenPoints = localDistance;
					firstDiagonalEnd = rFirstPoint;
					secondDiagonalEnd = rSecondPoint;
				}
			}
		}

		return std::pair<Math::Point3D, Math::Point3D>(firstDiagonalEnd, secondDiagonalEnd);
	}
	
	Math::Point3D CalculateFarthestPointFromLine(const std::vector<Point*>& rCluster, const Math::Line3D& rLine) {
		Math::Point3D farthestFromLine;
		double maxDistanceFromLine = 0.0;

		for (const auto pPoint : rCluster) {
			const auto& rPoint = pPoint->xyz.cast<double>();
			const double localDistance = Math::CalcPointToLineDistance3D(rPoint, rLine);

			if (localDistance > maxDistanceFromLine) {
				maxDistanceFromLine = localDistance;
				farthestFromLine = rPoint;
			}
		}

		return farthestFromLine;
	}
	
	Math::Point3D CalculateFarthestPointFromPlaneAndPoint(const std::vector<Point*>& rCluster,
		const Math::Plane3D& rPlane, const Math::Point3D& rPoint) {
		Math::Point3D farthestFromPlaneAndPoint;
		double maxDistanceFromPlaneAndPoint = 0.0;

		for (const auto pNextPoint : rCluster) {
			const auto& rNextPoint = pNextPoint->xyz.cast<double>();
			double localDistance = Math::CalcPointToPlaneDistance3D(rNextPoint, rPlane);
			localDistance += Math::CalcPointToPointDistance3D(rNextPoint, rPoint);

			if (localDistance > maxDistanceFromPlaneAndPoint) {
				maxDistanceFromPlaneAndPoint = localDistance;
				farthestFromPlaneAndPoint = rNextPoint;
			}
		}

		return farthestFromPlaneAndPoint;
	}

	double CalculateWallsFittingRMSE(const std::vector<Point*>& rCluster, const std::vector<Math::Plane3D>& rWalls) {
		const double wallsFittingSquareError = std::accumulate(rCluster.begin(), rCluster.end(), 0.0,
			[&](double init, Point* pPoint) {
				const auto& rPoint = pPoint->xyz.cast<double>();
				std::vector<double> distances;
				std::transform(rWalls.begin(), rWalls.end(), std::back_inserter(distances), [&](Math::Plane3D rWall) {
					return Math::CalcPointToPlaneDistance3D(rPoint, rWall);
					});
				const double& minPointWallDistance = *std::min_element(distances.begin(), distances.end());
				return init + pow(minPointWallDistance, 2);
			});
		return sqrt(wallsFittingSquareError / rCluster.size());
	}
};

OGX_EXPORT_METHOD(LoudspeakerSearching)
