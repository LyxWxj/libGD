#include "clustering.h"
#include <thread>
#include <numeric>

using namespace GD;
using namespace cluster;
using namespace std;

Kmeans::Kmeans(size_t k): 
	k(k), DistanceMetricSeted(false), MaxIterSeted(false), maxiter(10) {
}

template<typename metric>
void Kmeans::setDistanceMetric(metric const& distanceMetric) {
	this->distMetricfunc = [distanceMetric](Array_features const& a, Array_features const& b) {
		return distanceMetric(a, b);
	};
	this->DistanceMetricSeted = true;
}

void Kmeans::setMaxIter(size_t maxiter) {
	this->maxiter = maxiter;
	this->MaxIterSeted = true;
}

void Kmeans::calculateNearestCenter_(
	vector<Array_features> const& dataMat,
	vector<Array_features> const& centers,
	vector<pair<label, float>>& clusterAssment,
	bool& changed, size_t begin, size_t end
) {
	for (int i = begin; i < end; ++i) {
		double mindst = numeric_limits<double>::max();
		auto const& vec = dataMat[i];
		label oldcluster = this->clusterAss[i].first;
		label newcluster = -1;
		for (int j = 0; j < this->k; ++j) {
			double dist = this->distMetricfunc(vec, centers[i]);
			if (dist < mindst) {
				mindst = dist;
				clusterAssment[i] = make_pair(j, dist);
				newcluster = j;
			}
		}
		changed = true;
	}
}

void Kmeans::calculateNearestCenter(
	vector<Array_features> const& dataMat,
	vector<Array_features> const& centers,
	vector<pair<label, float>>& clusterAssment,
	bool& changed, size_t workers
) {
	if (workers == 1) {
		return calculateNearestCenter_(dataMat, centers, clusterAssment, changed, 0, dataMat.size());
	}
	size_t length = dataMat.size();
	size_t slice = length / workers;
	vector<thread> tasks;
	auto task = [this, &dataMat, &centers, &clusterAssment, &changed](size_t begin, size_t end) {
		this->calculateNearestCenter_(dataMat, centers, clusterAssment, changed, begin, end);
	};

	for (int i = 0; i < workers; ++i) {
		size_t begin = i * slice, end = (i + 1) * slice;
		tasks.emplace_back(thread(task, begin, end));
	}
	for (int i = 0; i < tasks.size(); ++i) {
		tasks[i].join();
	}
}



void Kmeans::updateCenter_(vector<Array_features>const& dataMat, vector<Array_features>& newCentroids, vector<int>& clusterCount, size_t begin, size_t end) {
	for (int i = begin; i < end; ++i) {
		label cluster = clusterAss[i].first;
		newCentroids[cluster] += dataMat[i];
		clusterCount[cluster] += 1;
	}
}

void Kmeans::updateCenter(vector<Array_features> const& dataMat, size_t workers) {
	vector<Array_features> newCentroids(this->k);
	vector<int> clusterCount(this->k);
	if (workers == 1) { 
		updateCenter_(dataMat, newCentroids, clusterCount, 0, dataMat.size());
		GD::operator/=(newCentroids, clusterCount);
		this->centroids = std::move(newCentroids);
		return;
	}
	size_t length = dataMat.size();
	size_t slice = length / workers;
	vector<thread> tasks;
	auto task = [this, &dataMat, &newCentroids, &clusterCount](size_t begin, size_t end) {
		this->updateCenter_(dataMat, newCentroids, clusterCount, begin, end);
	};
	for (int i = 0; i < workers; ++i) {
		size_t begin = i * slice, end = (i + 1) * slice;
		tasks.emplace_back(thread(task, begin, end));
	}
	for (int i = 0; i < tasks.size(); ++i) {
		tasks[i].join();
	}
}

void Kmeans::run(std::vector<Array_features> const& dataMat, size_t workers = 1) {
	if (this->DistanceMetricSeted == false) {
		cerr << "Distance Measurement Ambiguros!" << endl;
		exit(1);
	} 
	if (this->MaxIterSeted == false) {
		cerr << "MaxIteration Unseted!" << endl;
		exit(1);
	}
	this->centroids = cents_creater(k, dataMat);
	bool clusterChanged = true;
	size_t iters = 0;
	while (clusterChanged || iters < maxiter) {
		clusterChanged = false;
		++iters;
		calculateNearestCenter(dataMat, centroids, clusterAss, clusterChanged, workers);
		if (clusterChanged)
			updateCenter(dataMat, workers);
	}
}

vector<Array_features>&& Kmeans::move_centers() {
	return std::move(this->centroids);
}

vector<Array_features> const& Kmeans::get_centers() const {
	return this->centroids;
}

vector<pair<label, float>>&& Kmeans::move_clusterAss() {
	return std::move(this->clusterAss);
}

vector<pair<label, float>> const& Kmeans::get_clusterAss() const {
	return this->clusterAss;
}




