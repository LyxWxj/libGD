#pragma once

#ifndef __SAMPLER__H__
#define __SAMPLER__H__

#include <iostream>
#include <vector>
#include <string>
#include <thread>
#include "rand_Engine.h"
#include "raster.h"


namespace GD {
	class Sampler {
	public:
		Sampler() = default;
		~Sampler() = default;
		virtual std::vector<size_t> sample(size_t num) = 0;
		virtual void set_cdf_y_marginal(std::vector<float> const& cdf_y_marginal) = 0;
		virtual void set_cdf_x_condition(std::shared_ptr<GD::Raster<float>> cdf_x_condition) = 0;
	};

	template<typename RandomEngine>
	class UniformSampler : public Sampler {
		// 均匀采样器
	private:
		std::shared_ptr<RandomEngine> RandomCreator;
		size_t width = 0, height = 0;
	public:
		UniformSampler(size_t width, size_t height);
		~UniformSampler() = default;
		std::vector<size_t> sample(size_t num) override;
		void set_cdf_y_marginal(std::vector<float> const& cdf_y_marginal) {} // 不需要设置y的边缘分布
		void set_cdf_x_condition(std::shared_ptr<GD::Raster<float>> cdf_x_condition) {} // 不需要设置x的条件分布
	};

	class ImportanceSamplerHelper {
		// 重要性采样器辅助类, 用于计算重要性采样的PDF, cdf等
	private:
		GD::RasterPtr src = nullptr;
	public:
		ImportanceSamplerHelper() = default;
		~ImportanceSamplerHelper() = default;
		std::shared_ptr<GD::Raster<float>> calcPDF();
		std::vector<float> pdf2cdf(std::vector<float> const& pdf);
		std::vector<float> pdf2cdf(float* pdf, size_t size);
		void setSrc(GD::RasterPtr src);
		std::shared_ptr<GD::Raster<float>>
			calc_pdf_x_condition(GD::Raster<float> const& pdf, std::vector<float> const& pdf_y_marginal);
		std::tuple<std::vector<float>, std::shared_ptr<GD::Raster<float>>> calc_cdf_y_marginal_And_cdf_x_condition();
	};

	template<typename RandomEngine>
	class ImportanceSampler : public Sampler {
		// 重要性采样器
		std::shared_ptr<RandomEngine> RandomCreator;
		size_t width = 0, height = 0; // 采样的区间
	private:
		std::vector<float> const* cdf_y_marginal = nullptr;
		std::shared_ptr<GD::Raster<float>> cdf_x_condition = nullptr;
		size_t sample_from_cdf(std::vector<float> const& array, int dim);
		size_t sample_from_cdf(float* array, size_t size, int dim);
	private:
		size_t lower_bound(std::vector<float> const& arr, float val);
		size_t lower_bound(float* arr, size_t size, float val);
	public:
		ImportanceSampler(size_t width, size_t height);
		~ImportanceSampler() = default;
		std::vector<size_t> sample(size_t num) override;
		void set_cdf_y_marginal(std::vector<float> const& cdf_y_marginal);
		void set_cdf_x_condition(std::shared_ptr<GD::Raster<float>> cdf_x_condition);
	};

	using UnifSampler = UniformSampler<GD::UniformGenerator>;
	using SobolSampler = UniformSampler<GD::LdsGenerator>;
	using ImptcSampler = ImportanceSampler<GD::UniformGenerator>;
	using ImptcSamplerWithSobol = ImportanceSampler<GD::LdsGenerator>;
}

#endif // !__SAMPLER__H__