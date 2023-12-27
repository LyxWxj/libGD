#pragma once
#ifndef __TEMPLATES1_HPP__
#define __TEMPLATES1_HPP__

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include <gdal/gdal_priv.h>
#include <gdal/gdal.h>
#include <vector>

namespace GD {


	template<typename T>
	T sum(std::vector<T> const& v) {
		T sum = 0;
		std::for_each(v.begin(), v.end(), [&](T const& e) {sum += e; });
		return sum;
	}

	template<typename T>
	void swap(T& left, T& right) {
		T temp = left;
		left = right;
		right = temp;
	}

	template<typename T>
	int partition(std::vector<T>&vec, int left, int right) {
		int i = left, j = right + 1;
		T pivot = vec[left];
		while (1) {
			while (vec[++i] < pivot) {
				if (i == right) break;
			}
			while (vec[--j] > pivot) {
				if (j == left) break;
			}
			if (i >= j) break;
			swap(vec[i], vec[j]);
		}
		swap(vec[left], vec[j]);
		return j;
	}

	template<typename T>
	T qSelect(std::vector<T>& vec, size_t k) {
		int low = 0, high = vec.size() - 1;
		while (low <= high) {
			int j = partition(vec, low, high);
			if (j == k) {
				return vec[k];
			}
			else if (j > k) {
				high = j - 1;
			}
			else if (j < k) {
				low = j + 1;
			}
		}
		return vec[k];
	}

	template<typename _T>
	_T* Malloc(size_t width, size_t height, size_t bands) {
		return (_T*)CPLMalloc(width * height * bands * sizeof(_T));
	}

	template<typename T>
	constexpr double Type_Max() {
		return std::numeric_limits<T>::max();
	}

	template<typename T>
	constexpr double Type_Min() {
		return std::numeric_limits<T>::min();
	}

	template<typename ContainerA, typename ContainerB>
	double dot_product(ContainerA const& a, ContainerB const& b) {
		double sum = 0;
		size_t size = a.size() > b.size() ? b.size() : a.size();
		for (size_t i = 0; i < size; ++i) {
			sum += a[i] * b[i];
		}
		return sum;
	}
}

#endif // __TEMPLATES1_HPP__