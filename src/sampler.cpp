#include "sampler.h"

template<typename RandomEngine>
GD::UniformSampler<RandomEngine>::UniformSampler(size_t width, size_t height) :width(width), height(height) {
	if (std::is_same<RandomEngine, GD::UniformGenerator>::value) {
		this->RandomCreator = std::make_shared<RandomEngine>(width, height);
	}
	else if (std::is_same<RandomEngine, GD::LdsGenerator>::value) {
		this->RandomCreator = std::make_shared<RandomEngine>(width * height, 2);
	}
}

template<typename RandomEngine>
std::vector<size_t> GD::UniformSampler<RandomEngine>::sample(size_t nums) {
	std::vector<size_t> samples(nums);
	for (size_t i = 0; i < nums; ++i) {
		float x = RandomCreator->Get(0) * width;
		float y = RandomCreator->Get(1) * height;
		samples[i] = (size_t)(y * width + x);
	}
	return samples;
}

std::shared_ptr<GD::Raster<float>> GD::ImportanceSamplerHelper::calcPDF() {
	// 从一个GD::RasterPtr中计算PDF (归一化)
	if (src == nullptr) {
		std::cerr << "Error: src is nullptr" << std::endl;
		return nullptr;
	}
	size_t width = src->width(), height = src->height();
	std::shared_ptr<GD::Raster<float>> pdf = std::make_shared<GD::Raster<float>>(width, height, 1, GDT_Float32);
	float sum = 0;
	for (size_t i = 0; i < height; ++i) {
		float s = 0;
		for (long long j = 0; j < width; ++j) {
			float value = src->get(i, j)[0];
			s += value;
			(*pdf)(i, j) = value;
		}
		sum += s;
	}
//#pragma omp parallel for num_threads(10)
	for (long long i = 0; i < height; ++i) {
		for (size_t j = 0; j < width; ++j) {
			(*pdf)(i, j) /= sum;
		}
	}
	return pdf;
}

std::vector<float> GD::ImportanceSamplerHelper::pdf2cdf(std::vector<float> const& pdf) {
	// 从PDF计算CDF
	size_t size = pdf.size();
	std::vector<float> cdf(size, 0);
	cdf[0] = pdf[0];
	for (size_t i = 1; i < size; ++i) {
		cdf[i] = cdf[i - 1] + pdf[i];
	}
	return cdf;
}

std::vector<float> GD::ImportanceSamplerHelper::pdf2cdf(float* pdf, size_t size) {
	std::vector<float> cdf(size, 0);
	cdf[0] = pdf[0];
	for (size_t i = 1; i < size; ++i) {
		cdf[i] = cdf[i - 1] + pdf[i];
	}
	return cdf;
}

std::shared_ptr<GD::Raster<float>> GD::ImportanceSamplerHelper::
calc_pdf_x_condition(GD::Raster<float> const& pdf, std::vector<float> const& pdf_y_marginal) {
	std::shared_ptr<GD::Raster<float>> res = std::make_shared<GD::Raster<float>>(pdf.width(), pdf.height(), 1, GDT_Float32);
	return res;
}

void GD::ImportanceSamplerHelper::setSrc(GD::RasterPtr src) {
	this->src = src;
}

void calc_pdf_y_marginal(std::shared_ptr<GD::Raster<float>> pdf, std::vector<float>& pdf_y_marginal, size_t width, size_t height) {
#pragma omp parallel for num_threads(10)
	for (long long i = 0; i < height; ++i) {
		auto p = (*pdf)(i);
		for (long long j = 0; j < width; ++j) {
			pdf_y_marginal[i] += p[j];
		}
	}
}

std::tuple<std::vector<float>, std::shared_ptr<GD::Raster<float>>> GD::ImportanceSamplerHelper::
calc_cdf_y_marginal_And_cdf_x_condition() {

	size_t height = src->height(), width = src->width();
	// 归一化图像
	auto p_pdf = calcPDF();
	// 计算y边缘概率密度
	auto pdf_y_marginal = std::vector<float>(height, 0);
	calc_pdf_y_marginal(p_pdf, pdf_y_marginal, width, height);
	// 计算y边缘概率分布
	auto cdf_y_marginal = pdf2cdf(pdf_y_marginal);
	auto cdf_x_condition = std::make_shared<GD::Raster<float>>(width, height, 1, GDT_Float32);
	// 计算x条件概率密度
	// p_pdf[i][j] 表示在y = i 的条件下, 取 x = j 的条件概率

	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			(*p_pdf)(i, j) /= pdf_y_marginal[i];
		}
	}
	// 计算x条件概率分布
	for (long long i = 0; i < height; ++i) {
		auto p = (*p_pdf)(i);
		auto cdf = pdf2cdf(p, width);
		for (long long j = 0; j < width; ++j) {
			cdf_x_condition->set(i, j, 0, cdf[j]);
		}
	}
	return { cdf_y_marginal, cdf_x_condition };
}

template<typename RandomEngine>
GD::ImportanceSampler<RandomEngine>::ImportanceSampler(size_t width, size_t height) :width(width), height(height) {
	this->RandomCreator = std::make_shared<RandomEngine>(width * height, 2);
}

template<typename RandomEngine>
std::vector<size_t> GD::ImportanceSampler<RandomEngine>::sample(size_t nums) {
	std::vector<size_t> samples(nums);
	float* p = (float*)cdf_x_condition->data();
	for (size_t i = 0; i < nums; ++i) {
		auto y = sample_from_cdf(*cdf_y_marginal, 0);
		auto x = sample_from_cdf(p + y * width, width, 1);
		samples[i] = y * width + x;
	}
	return samples;
}

template<typename RandomEngine>
size_t GD::ImportanceSampler<RandomEngine>::lower_bound(std::vector<float> const& arr, float val) {
	return std::lower_bound(arr.begin(), arr.end(), val) - arr.begin();
}

template<typename RandomEngine>
size_t GD::ImportanceSampler<RandomEngine>::lower_bound(float* arr, size_t size, float val) {
	return std::lower_bound(arr, arr + size, val) - arr;
}

template<typename RandomEngine>
void GD::ImportanceSampler<RandomEngine>::set_cdf_y_marginal(std::vector<float> const& cdf_y_marginal) {
	this->cdf_y_marginal = &cdf_y_marginal;
}

template<typename RandomEngine>
void GD::ImportanceSampler<RandomEngine>::set_cdf_x_condition(std::shared_ptr<GD::Raster<float>> cdf_x_condition) {
	this->cdf_x_condition = cdf_x_condition;
}

template<typename RandomEngine>
size_t GD::ImportanceSampler<RandomEngine>::sample_from_cdf(std::vector<float> const& cdfArray, int dim) {
	float val = this->RandomCreator->Get(dim);
	return lower_bound(cdfArray, val);
}

template<typename RandomEngine>
size_t GD::ImportanceSampler<RandomEngine>::sample_from_cdf(float* cdfArray, size_t size, int dim) {
	float val = this->RandomCreator->Get(dim);
	return lower_bound(cdfArray, size, val);
}

void templateTest() { // 不需要调用, 只提示编译器创建模板的实例
	GD::UnifSampler A = GD::UnifSampler(0,0);
	GD::SobolSampler B = GD::SobolSampler(0, 0);
	GD::ImptcSampler C = GD::ImptcSampler(0, 0);
	GD::ImptcSamplerWithSobol D = GD::ImptcSamplerWithSobol(0, 0);
}