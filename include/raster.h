#ifndef __RASTER__H_
#define __RASTER__H_

#include <gdal/gdal_priv.h>
#include <gdal/gdal.h>
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
#include <cstring>
#include <initializer_list>
#include <random>
#include <thread>
#include "feature.h"
#include <Eigen/Eigen>
#include <Eigen/Core>
#include "templates.hpp"

#define POS(x, y, channel, Width, Area) (x + y * Width + Area * channel)

namespace GD {

	template<typename _Ty>
	struct Matrix {
	public:
		static_assert(std::is_arithmetic<_Ty>::value, "_Ty must be numeric type");
		// data
		std::vector<_Ty*> Data;
		size_t m_width, m_height;
	public:
		Matrix(size_t width = 1, size_t height = 1);
		Matrix(std::initializer_list<std::initializer_list<_Ty>> const& init_list);
	};

	class RasterBase {
	public:
		using RasterPtr = std::shared_ptr<RasterBase>;
		RasterBase() = default;
		RasterBase(const std::string& filepath);
		RasterBase(const char* filepath);
		virtual Array_features get(size_t i, size_t j) const;
		virtual float get(size_t i, size_t j, size_t channel) const;
		virtual size_t width() const;
		virtual size_t height() const;
		virtual size_t bands() const;
		virtual void* data() const;
		virtual GDALDataType type() const;
		virtual void init(GDALDataset* pDatasetRead,
			size_t lWidth, size_t lHeight, size_t nBands,
			GDALDataType datatype);
		virtual void save(std::string const&, std::string const& DriverName = "GTIFF") const;
		virtual std::vector<std::array<size_t, 256>> histgram() const;
		virtual RasterPtr equalizeHist(std::vector<std::array<size_t, 256>> histgram);
		virtual RasterPtr equalizeHist();
		virtual RasterPtr channel_extraction(size_t channel);
		virtual RasterPtr filter(Matrix<float> const& kernel);
		virtual RasterPtr MedianFilter();
		virtual RasterPtr addSalt(size_t num);
		virtual RasterPtr Sobel();
		virtual RasterPtr Laplace_4();
		virtual RasterPtr Laplace_8();
		virtual RasterPtr GaussBlur(size_t size, float stdx, float stdy);
		virtual RasterPtr addGaussNoise(size_t size, float stdy);
		virtual RasterPtr Laplace_4_Sharpening();
		virtual RasterPtr Laplace_8_Sharpening();
		virtual RasterPtr Sobel_Sharpening();
		virtual RasterPtr toGray();
		virtual ~RasterBase();
		virtual void resize(size_t width, size_t height);
		virtual void set(size_t i, size_t j, Array_features const& arr);
		virtual void set(size_t i, size_t j, Eigen::Vector3f const& arr);
		virtual void set(size_t i, size_t j, size_t channel, double value);
		virtual void add(size_t i, size_t j, size_t channel, double const value);
		virtual void add(size_t i, size_t j, double const value);
		virtual void add(size_t i, size_t j, size_t channel, Array_features const& arr);
		virtual size_t type_size()const;
		virtual size_t Area()const;
		virtual RasterPtr copy() const;
		virtual void setnullptr();
	};

	using RasterPtr = std::shared_ptr<RasterBase>;

	template<typename T>
	class Raster : virtual public RasterBase {
	private:
		GDALDataset* pDatasetRead = nullptr;
		GDALDriver* pDriver = nullptr;
		GDALDataType datatype;
		std::vector<float> conv(long long pos, Matrix<float> const& mat);
		RasterPtr Laplace(Matrix<float> const& kernel);
	public:
		size_t lWidth, lHeight, nBands;
		size_t area;
		T* pData = nullptr;
		Array_features get(size_t i, size_t j) const override;
		float get(size_t i, size_t j, size_t channel)const override;
		size_t width() const override;
		size_t height() const override;
		size_t bands() const override;
		void* data() const override;
		GDALDataType type() const override;
		void init(GDALDataset* pDatasetRead,
			size_t lWidth, size_t lHeight, size_t nBands,
			GDALDataType datatype)override;

		void save(std::string const& filepath, std::string const& DriverName = "GTIFF") const override;

		std::vector<std::array<size_t, 256>> histgram() const override;

		RasterPtr equalizeHist(std::vector<std::array<size_t, 256>> histgram) override;

		RasterPtr equalizeHist()override;

		RasterPtr channel_extraction(size_t channel)override;

		RasterPtr filter(Matrix<float> const& kernel)override;

		RasterPtr MedianFilter()override;

		RasterPtr addSalt(size_t num) override;
		RasterPtr Sobel() override;

		RasterPtr GaussBlur(size_t size, float stdx, float stdy)override;

		RasterPtr addGaussNoise(size_t size, float std) override;
		RasterPtr Laplace_4()override;

		RasterPtr Laplace_8()override;

		RasterPtr Laplace_4_Sharpening() override;
		RasterPtr Laplace_8_Sharpening() override;
		RasterPtr Sobel_Sharpening()override;

		RasterPtr toGray() override;

		RasterPtr copy()const override ;

		void resize(size_t width, size_t height) override;

		void set(size_t i, size_t j, Array_features const& arr) override;
		void set(size_t i, size_t j, Eigen::Vector3f const& arr) override;
		void set(size_t i, size_t j, size_t channel, double value)override;
		void add(size_t i, size_t j, size_t channel, double const value)override;
		void add(size_t i, size_t j, double const value)override;
		void add(size_t i, size_t j, size_t channel, Array_features const& arr)override;
		size_t type_size()const override;
		size_t Area()const override;
		void setnullptr() override;
		void setType(GDALDataType type) {
			this->datatype = type;
		}
	public:
		Raster();

		Raster(size_t width, size_t height, size_t nBands, GDALDataType datatype);

		~Raster();

		T& operator()(size_t i, size_t j, size_t band = 0) {
			return pData[POS(j, i, band, lWidth, area)];
		}

		T* operator()(size_t y) {
			return pData + y * lWidth;
		}
	};

	RasterPtr getByType(GDALDataType datatype);

	RasterPtr open(const std::string& filepath);

	RasterPtr operator-(RasterBase const& lhs, RasterBase const& rhs);

	void HistgramMatch(RasterPtr src, RasterPtr dst);

	void RGB2IHS(Eigen::Vector3f const&, Eigen::Vector3f&);
	void IHS2RGB(Eigen::Vector3f const&, Eigen::Vector3f&);
	
}

#endif // __RASTER__H_