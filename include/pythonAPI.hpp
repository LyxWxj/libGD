#pragma once
#ifndef __PYTHONAPI_H__
#define __PYTHONAPI_H__


#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <list>
#include <vector>
#include <array>
#include <cassert>
#include <memory>
#include "raster.h"
#include "image_fusion.h"
#include "image_segment.h"

namespace py = pybind11;

namespace GD {
	template<typename T>
	py::array_t<T> _ptr_to_arrays_1d(T* data, py::ssize_t col) {
		auto result = py::array_t<T>(col);
		py::buffer_info buf = result.request();
		T* ptr = (T*)buf.ptr;
		memcpy(ptr, data, col * sizeof(T));
		return result;
	}

	template<typename T>
	py::array_t<T> _ptr_to_arrays_2d(T* data, py::ssize_t row, py::ssize_t col) {
		auto result = _ptr_to_arrays_1d(data, row * col);
		result.resize({ row, col });
		return result;
	}


	template<typename T>
	py::array_t<T> ptr_to_arrays_3d(T* data, py::ssize_t height, py::ssize_t width, py::ssize_t band) {
		py::array_t<T> out(band * height * width);
		out.resize({ height, width, band });
		size_t Area = height * width;
		for (int c = 0; c < band; ++c) {
			for (int i = 0; i < height; ++i) {
				for (int j = 0; j < width; ++j) {
					out.mutable_at(i, j, c) = data[c * Area + i * width + j];
				}
			}
		}
		return out;
	}
 
	template<typename T>
	GD::Matrix<T> array_to_Matrix(py::array_t<T> const& arr) {
		py::buffer_info buf = arr.request();
		T* ptr = (T*)buf.ptr;
		GD::Matrix<T> mat(arr.shape(0), arr.shape(1));
		memcpy(mat.Data.data(), ptr, arr.shape(0) * arr.shape(1) * sizeof(T));
		return mat;
	}

	template<typename T>
	GD::RasterPtr array_to_Raster(py::array_t<T> const& arr) {
		py::buffer_info buf = arr.request();
		GD::Raster<T>* raster = new GD::Raster<T>;
		raster->pData = new T[arr.shape(0) * arr.shape(1) * arr.shape(2)];
		raster->lHeight = arr.shape(0);
		raster->lWidth = arr.shape(1);
		if (std::is_same_v<T, unsigned char>) {
			raster->setType(GDT_Byte);
		}
		else if (std::is_same_v<T, float>) {
			raster->setType(GDT_Float32);
		}
		else if (std::is_same_v<T, double>) {
			raster->setType(GDT_Float64);
		}
		else if (std::is_same_v<T, int>){
			raster->setType(GDT_Int32);
		}
		else if (std::is_same_v<T, short>){
			raster->setType(GDT_Int16);
		}
		else if (std::is_same_v<T, unsigned short>){
			raster->setType(GDT_UInt16);
		}
		else if (std::is_same_v<T, unsigned int>){
			raster->setType(GDT_UInt32);
		}
		else {
			throw std::runtime_error("Unsupported data type.");
		}
		if (arr.ndim() == 2)
			raster->nBands = 1;
		else if (arr.ndim() == 3)
			raster->nBands = arr.shape(2);
		raster->area = arr.shape(0) * arr.shape(1);
		raster->pData = GD::Malloc<T>(raster->width(), raster->height(), raster->nBands);
		size_t height = arr.shape(0), width = arr.shape(1), bands = arr.shape(2);
		for(int i = 0; i < height; ++ i) {
			for(int j = 0; j < width; ++j) {
				for(int c = 0; c < bands; ++c) {
					raster->pData[c * raster->area + i * width + j] = arr.at(i, j, c);
				}
			}
		}
		GD::RasterPtr ret (raster);
		return ret;
	}

	template<typename T>
	py::array_t<T> Raster_to_array(GD::Raster<T>* raster) {
		return ptr_to_arrays_3d(raster->pData, raster->height(), raster->width(), raster->bands());
	}

	template<typename T>
	py::array_t<T> filter(py::array_t<T> const& data, py::array_t<float> kernel) {
		GD::RasterPtr raster = array_to_Raster(data);
		GD::Matrix<float> mat = array_to_Matrix(kernel);
		GD::RasterPtr out = raster->filter(mat);
		return ptr_to_arrays_3d((T*)(out->data()), out->height(), out->width(), out->bands());
	}

	template<typename T>
	py::array_t<T> MedianFilter(py::array_t<T> const& data) {
		GD::RasterPtr raster = array_to_Raster(data);
		GD::RasterPtr out = raster->MedianFilter();
		return ptr_to_arrays_3d((T*)(out->data()), out->height(), out->width(), out->bands());
	}

	template<typename T>
	py::array_t<T> GaussFilter(py::array_t<T> const& data, py::ssize_t kernelSize, float stdx, float stdy) {
		GD::RasterPtr raster = array_to_Raster(data);
		GD::RasterPtr out = raster->GaussBlur(kernelSize, stdx, stdy);
		return ptr_to_arrays_3d((T*)(out->data()), out->height(), out->width(), out->bands());
	}

	template<typename T>
	py::array_t<T> Sobel(py::array_t<T> const& data) {
		GD::RasterPtr raster = array_to_Raster(data);
		GD::RasterPtr out = raster->Sobel();
		return ptr_to_arrays_3d((T*)(out->data()), out->height(), out->width(), out->bands());
	}

	template<typename T>
	py::array_t<T> Laplace(py::array_t<T> const& data) {
		GD::RasterPtr raster = array_to_Raster(data);
		GD::RasterPtr out = raster->Laplace_8();
		return ptr_to_arrays_3d((T*)(out->data()), out->height(), out->width(), out->bands());
	}

	template<typename T>
	py::array_t<T> toGray(py::array_t<T> const& data) {
		if (data.ndim() == 2) return data;
		GD::RasterPtr raster = array_to_Raster(data);
		GD::RasterPtr out = raster->toGray();
		return _ptr_to_arrays_2d((T*)(out->data()), out->height(), out->width());
	}

	template<typename T>
	py::array_t<float> IHSFusing(py::array_t<T> const& Ms, py::array_t<T> const& Pans) {
		GD::RasterPtr rasterMs = array_to_Raster(Ms);
		GD::RasterPtr rasterPans = array_to_Raster(Pans);
		GD::IHSFusion fusion;
		RasterPtr result = fusion.FusionOnCPU(rasterMs, rasterPans); // Raster<float>
		auto res_array = ptr_to_arrays_3d((float*) result->data(), result->height(), result->width(), result->bands());
		return res_array;
	}

	template<typename T>
	py::array_t<float> GIHSFusing(py::array_t<T> const& Ms, py::array_t<T> const& Pans) {
		GD::RasterPtr rasterMs = array_to_Raster(Ms);
		GD::RasterPtr rasterPans = array_to_Raster(Pans);
		GD::GIHSFusion fusion;
		RasterPtr result = fusion.FusionOnCPU(rasterMs, rasterPans); // Raster<float>
		auto res_array = ptr_to_arrays_3d((float*)result->data(), result->height(), result->width(), result->bands());
		return res_array;
	}

	template<typename T>
	py::array_t<float> AIHSFusing(py::array_t<T> const& Ms, py::array_t<T> const& Pans) {
		GD::RasterPtr rasterMs = array_to_Raster(Ms);
		GD::RasterPtr rasterPans = array_to_Raster(Pans);
		GD::AIHSFusion fusion;
		fusion.setUseSobol(true);
		fusion.setUseImportanceSampling(true);
		RasterPtr result = fusion.FusionOnCPU(rasterMs, rasterPans); // Raster<float>
		auto res_array = ptr_to_arrays_3d((float*)result->data(), result->height(), result->width(), result->bands());
		return res_array;
	}

	template<typename T>
	py::array_t<float> IAIHSFusing(py::array_t<T> const& Ms, py::array_t<T> const& Pans) {
		GD::RasterPtr rasterMs = array_to_Raster(Ms);
		GD::RasterPtr rasterPans = array_to_Raster(Pans);
		GD::IAIHSFusionWithAlpha fusion;
		fusion.setUseSobol(true);
		fusion.setUseImportanceSampling(true);
		RasterPtr result = fusion.FusionOnCPU(rasterMs, rasterPans); // Raster<float>
		auto res_array = ptr_to_arrays_3d((float*)result->data(), result->height(), result->width(), result->bands());
		return res_array;
	}

	template<typename T>
	py::array_t<T> HistgramEqualize(py::array_t<T> const& data) {
		GD::RasterPtr raster = array_to_Raster(data);
		GD::RasterPtr out = raster->equalizeHist();
		return ptr_to_arrays_3d((T*)(out->data()), out->height(), out->width(), out->bands());
	}

	float flectance (float p) {
		if ( p > 255 ) return 0.f;
		if (p < 0) return 255.f;
		return 255.f - p;
		};

	template<typename Func, typename T>
	auto flectance_img = [](GD::Raster<T>* img, Func f) {
		// 对img原地修改
		size_t width = img->width(), height = img->height(), bands = img->bands();
		size_t area = width * height * bands;
		for (long long i = 0; i < area; ++i) {
			T value = img->pData[i];
			img->pData[i] = (T)f(value);
		}
	};

	template<typename T>
	py::array_t<T> ReflectImg(py::array_t<T> const& data) {
		GD::RasterPtr raster = array_to_Raster(data);
		GD::Raster<T>* raster_ = dynamic_cast<GD::Raster<T>*>(raster.get());
		flectance_img<decltype(flectance), T>(raster_, flectance);
		auto res_array = ptr_to_arrays_3d((T*)raster_->data(), raster_->height(), raster_->width(), raster_->bands());
		return res_array;
	}	

	template<typename T>
	py::array_t<T> SLICSegment(py::array_t<T> const& data, py::array_t<T> const& gradiant, py::ssize_t segments) {
		GD::RasterPtr raster = array_to_Raster(data);
		GD::RasterPtr raster_gradiant = array_to_Raster(gradiant);
		GD::SLIC slic;
		slic.set_raster(raster);
		slic.setSegments(segments);
		slic.set_Gradient_reflected(raster_gradiant);
		RasterPtr result = slic.run();
		auto res_array = ptr_to_arrays_3d((T*)result->data(), result->height(), result->width(), result->bands());
		return res_array;
	}
}



#endif // !__PYTHONAPI_H__