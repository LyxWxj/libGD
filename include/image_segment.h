#pragma once
#ifndef IMAGE_SEGMENT_H
#define IMAGE_SEGMENT_H

#include "raster.h"
#include "sampler.h"
#include "feature.h"
#include <mutex>
#include <unordered_set>
#include <set>

namespace GD {
	using label = size_t;
	class ImageSegment {
	public:
		virtual void set_raster(GD::RasterPtr) = 0;
		virtual void set_Gradient_reflected(GD::RasterPtr) = 0;
		virtual void setSegments(size_t segments) = 0;
		virtual GD::RasterPtr run() = 0;
	};

	class SLIC : public ImageSegment {
	public:
		class cluster {
		public:
			std::mutex mu;
			size_t center_pos;
			std::unordered_set<size_t> pixels;
			cluster(size_t center_pos) :center_pos(center_pos) {}
			cluster(cluster const& other) :center_pos(other.center_pos), pixels(other.pixels) {}
			cluster(cluster&& other) noexcept :center_pos(other.center_pos), pixels(std::move(other.pixels)) {}
			void add_pixel(size_t pixel) {
				mu.lock();
				pixels.insert(pixel);
				mu.unlock();
			}
			void remove_pixel(size_t pixel) {
				mu.lock();
				pixels.erase(pixel);
				mu.unlock();
			}
			void clear() {
				pixels.clear();
			}
			size_t get_center_pos() {
				return center_pos;
			}
		};
	private:
		std::vector<std::pair<int, int>> centroids;
		std::vector<std::pair<GD::label, float>> clusterAss;
	private:
		GD::RasterPtr raster;
		GD::RasterPtr gradient_reflected;
		size_t segments;
		std::vector<size_t> centers;
		GD::RasterPtr class_result;
		std::vector<std::array<float, 3>> features;
		int S;
	private:
		static std::array<float, 3> bgr2xyz(std::array<float, 3>const& bgr);
		static float gamma(float);
		static std::array<float, 3> xyz2lab(std::array<float, 3>const&);
		static std::array<float, 3> bgr2lab(std::array<float, 3>const&);
		void calc_features();
		bool assignment();
		bool update();
	public:
		virtual void set_raster(GD::RasterPtr)override;
		template<typename T>
		void set_raster(Raster<T>*);
		virtual void set_Gradient_reflected(GD::RasterPtr)override;
		template<typename T>
		void set_Gradient_reflected(Raster<T>*);
		virtual void setSegments(size_t segments) override;
		virtual GD::RasterPtr run() override;
	};
	template<typename T>
	inline void SLIC::set_raster(Raster<T>* pRaster) {
		RasterPtr raster(pRaster);
		set_raster(raster);
	}
	template<typename T>
	inline void SLIC::set_Gradient_reflected(Raster<T>* pRaster) {
		RasterPtr raster(pRaster);
		set_Gradient_reflected(raster);
	}
}
#endif // IMAGE_SEGMENT_H