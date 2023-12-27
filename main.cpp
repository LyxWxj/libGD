#include <iostream>
#include "include\templates.hpp"
#include "include\pythonAPI.hpp"
#include <string>
#include <pybind11/embed.h>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>


std::string GDInformation() {
	return "libGD is An Image Proccessing Python library.\nWritten with Cplusplus By Yuxuan Lou SZ University.\nOnly For Study.";
}

std::string test() {
	return "LibGD is working.";
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

PYBIND11_MODULE(GD, m) {
	m.doc() = "libGD"; // optional module docstring
	m.def("filterf32", &GD::filter<float>, "Image Convolution");
	m.def("filterui8", &GD::filter<unsigned char>, "Image Convolution");
	m.def("Medianfilterf32", &GD::MedianFilter<float>, "Image Median Filter");
	m.def("Medianfilterui8", &GD::MedianFilter<unsigned char>, "Image Median Filter");
	m.def("Gaussianfilterf32", &GD::GaussFilter<float>, "Image Gaussian Filter");
	m.def("Gaussianfilterui8", &GD::GaussFilter<unsigned char>, "Image Gaussian Filter");
	m.def("Sobelf32", &GD::Sobel<float>, "Image Sobel Filter");
	m.def("Sobelui8", &GD::Sobel<unsigned char>, "Image Sobel Filter");
	m.def("Laplacef32", &GD::Laplace<float>, "Image Laplace Filter");
	m.def("Laplaceui8", &GD::Laplace<unsigned char>, "Image Laplace Filter");
	m.def("toGrayf32", &GD::toGray<float>, "Image to Gray");
	m.def("toGrayui8", &GD::toGray<unsigned char>, "Image to Gray");
	m.def("IHSFusingf32", &GD::IHSFusing<float>, "Image IHS Fusing");
	m.def("IHSFusingui8", &GD::IHSFusing<unsigned char>, "Image IHS Fusing");
	m.def("GIHSFusingf32", &GD::GIHSFusing<float>, "Image GIHS Fusing");
	m.def("GIHSFusingui8", &GD::GIHSFusing<unsigned char>, "Image GIHS Fusing");
	m.def("AIHSFusingf32", &GD::AIHSFusing<float>, "Image AIHS Fusing");
	m.def("AIHSFusingui8", &GD::AIHSFusing<unsigned char>, "Image AIHS Fusing");
	m.def("IAIHSFusingf32", &GD::IAIHSFusing<float>, "Image IAIHS Fusing");
	m.def("IAIHSFusingui8", &GD::IAIHSFusing<unsigned char>, "Image IAIHS Fusing");
	m.def("SLICSegmentf32", &GD::SLICSegment<float>, "Image SLIC Segment");
	m.def("SLICSegmentui8", &GD::SLICSegment<unsigned char>, "Image SLIC Segment");
	m.def("GDInformation", &GDInformation, "GDInformation");
	m.def("ReflectImgf32", &GD::ReflectImg<float>, "Reflect Image");
	m.def("ReflectImgui8", &GD::ReflectImg<unsigned char>, "Reflect Image");
	m.def("HistgramEqualizef32", &GD::HistgramEqualize<float>, "Histgram Equalize");
	m.def("HistgramEqualizeui8", &GD::HistgramEqualize<unsigned char>, "Histgram Equalize");
	m.def("test", &test, "test");
}

//float flectance(float p) {
//	if (p > 255) return 0.f;
//	if (p < 0) return 255.f;
//	return 255.f - p;
//};
//
//auto flectance_img (GD::RasterPtr img) {
//	 对img原地修改
//	size_t width = img->width(), height = img->height(), bands = img->bands();
//	size_t area = width * height;
//	for (int c = 0; c < bands; ++c) {
//		for (long long i = 0; i < area; ++i) {
//			auto val = img->get(i / width, i % width, c);
//			val = ::flectance(val);
//			img->set(i / width, i % width, c, val);
//		}
//	}
//};
//
//using namespace GD;
//
//int main() {
//	GDALAllRegister();
//	py::initialize_interpreter();
//	auto Ms = GD::open("output//reBlur.tif");
//	auto Pans = GD::open("output//PAN.tif");
//	auto Msarr = GD::Raster_to_array(dynamic_cast<Raster<unsigned char>*>(Ms.get()));
//	auto Pansarr = GD::Raster_to_array(dynamic_cast<Raster<unsigned char>*>(Pans.get()));
//	auto res = GD::GIHSFusing<unsigned char>(Msarr, Pansarr);
//	auto resRaster = GD::array_to_Raster(res);
//	resRaster->save("output//GIHSFusing.tif");
//	res = GD::AIHSFusing<unsigned char>(Msarr, Pansarr);
//	resRaster = GD::array_to_Raster(res);
//	resRaster->save("output//AIHSFusing.tif");
//	res = GD::IAIHSFusing<unsigned char>(Msarr, Pansarr);
//	resRaster = GD::array_to_Raster(res);
//	resRaster->save("output//IAIHSFusing.tif");
//	return 0;
//}
//
