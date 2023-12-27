#pragma once

#ifndef _CLUSTERING_H_
#define _CLUSTERING_H_

#include "raster.h"
#include "feature.h"
#include <functional>
#include <array>

namespace GD {
	using label = size_t;

	struct random_cent {
		Array_features operator()(size_t k, std::vector<Array_features> const&);
	};

	struct kpp_cent {
		Array_features operator()(size_t k, std::vector<Array_features> const&);
	};

	namespace cluster {
		class Kmeans {
		private: // Basic Member
			size_t k;
			size_t maxiter;
			std::function<double(Array_features const&, Array_features const&)> distMetricfunc;

			// Results Member
			std::vector<Array_features> centroids;
			std::vector<std::pair<label, float>> clusterAss;
		private: // Flag About Setting
			bool DistanceMetricSeted = false;
			bool MaxIterSeted = false;

		private: // Basic Function
			void calculateNearestCenter_(std::vector<Array_features> const& dataMat, std::vector<Array_features> const&, std::vector<std::pair<label, float>>&, bool& changed, size_t, size_t);
			void calculateNearestCenter(std::vector<Array_features> const& dataMat, std::vector<Array_features> const&, std::vector<std::pair<label, float>>&, bool& changed, size_t);
			void updateCenter_(std::vector<Array_features>const& dataMat, std::vector<Array_features>& newCentroids, std::vector<int>& clusterCount, size_t, size_t);
			void updateCenter(std::vector<Array_features>const& dataMat, size_t workers);

		protected:
			std::function<std::vector<Array_features>(size_t, std::vector<Array_features> const&)> cents_creater;
		public:
			Kmeans(size_t);

			template<typename metric>
			void setDistanceMetric(metric const& distMetric);

			void setMaxIter(size_t);
			void run(std::vector<Array_features> const&, size_t);

			std::vector<Array_features>&& move_centers();
			std::vector<Array_features> const& get_centers() const;
			std::vector<std::pair<label, float>>&& move_clusterAss();
			std::vector<std::pair<label, float>> const& get_clusterAss() const;

		};
	}
}

#endif // !_CLUSTERING_H_