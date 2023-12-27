#include "image_segment.h"
#include <omp.h>
const float param_13 = 1.0f / 3.0f;
const float param_16116 = 16.0f / 116.0f;
const float Xn = 0.950456f;
const float Yn = 1.0f;
const float Zn = 1.088754f;
const float m = 55.0f;

inline float dc (std::array<float, 3>const& labj, std::array<float, 3>const& labi) {
    float d0 = labj[0] - labi[0];
    float d1 = labj[1] - labi[1];
    float d2 = labj[2] - labi[2];
    return d0 * d0 + d1 * d1 + d2 * d2;
};

inline float ds (std::array<unsigned short, 2>const& xyj, std::array<unsigned short, 2>const& xyi) {
    int d0 = xyj[0] - xyi[0];
    int d1 = xyj[1] - xyi[1];
    return d0 * d0 + d1 * d1;
};


void GD::SLIC::setSegments(size_t segs) {
	segments = segs;
	S = 2 * std::sqrt(raster->Area() / segments);
	centroids.resize(segments);
	clusterAss = std::vector<std::pair<GD::label, float>>(raster->width() * raster->height(), std::pair<GD::label, float>(0xf0000000, 1e9));
}

void GD::SLIC::set_Gradient_reflected(GD::RasterPtr Gradient) {
	gradient_reflected = Gradient;
}

void GD::SLIC::set_raster(GD::RasterPtr Raster) {
	raster = Raster;
}

std::array<float,3> GD::SLIC::bgr2xyz(std::array<float,3>const& bgr) {
	float b = bgr[0] / 255.0f, g = bgr[1] / 255.0f, r = bgr[2] / 255.0f;
	float X = 0.412453f * r + 0.357580f * g + 0.180423f * b;
	float Y = 0.212671f * r + 0.715160f * g + 0.072169f * b;
	float Z = 0.019334f * r + 0.119193f * g + 0.950227f * b;
	return std::array<float, 3>{ X, Y, Z };
}

float GD::SLIC::gamma(float x) {
	return x > 0.04045 ? powf((x + 0.055f) / 1.055f, 2.4f) : (x / 12.92);
}

std::array<float, 3> GD::SLIC::xyz2lab(std::array<float, 3> const& xyz)
{
	float X = xyz[0] / Xn, Y = xyz[1] / Yn, Z = xyz[2] / Zn;
	float fX, fY, fZ;
	if (Y > 0.008856f)
		fY = pow(Y, param_13);
	else
		fY = 7.787f * Y + param_16116;

	if (X > 0.008856f)
		fX = pow(X, param_13);
	else
		fX = 7.787f * X + param_16116;

	if (Z > 0.008856)
		fZ = pow(Z, param_13);
	else
		fZ = 7.787f * Z + param_16116;
	float L = 116.0f * fY - 16.0f;
	L = L > 0 ? L : 0;
	float a = 500.0f * (fX - fY);
	float b = 200.0f * (fY - fZ);
	return std::array<float, 3>{ L, a, b };
}

std::array<float, 3> GD::SLIC::bgr2lab(std::array<float, 3> const& bgr) {
	return xyz2lab(bgr2xyz(bgr));
}

void GD::SLIC::calc_features() {
	int width = raster->width(), height = raster->height();
	size_t total = width * height;
	this->features.resize(raster->width() * raster->height());
	for (size_t p = 0; p < total; ++p) {
		int i = p / width, j = p % width;
		auto bgr = raster->get(i, j);
		auto lab = bgr2lab({bgr[0], bgr[1], bgr[2]});
		features[p] = std::array<float, 3>{ lab[0], lab[1], lab[2]};
	}
	for (size_t p = 0; p < this->centers.size(); ++p) {
        size_t pos = centers[p];
        int i = pos / width, j = pos % width;
        centroids[p] = { i, j };
    }
}

bool task_innerHloop(int c,int h, int j, int S, size_t width,
	std::vector<std::array<float, 3>> const& features,
	std::array<float, 3>const & labi,
	std::array<unsigned short, 2>const & xyi,
	std::vector<std::pair<GD::label, float>>& clusterAss
) {
	bool res = false;
	for (int w = j - S; w < j + S; ++w) {
		if (w < 0) {
			continue;
		}else if (w == width) {
			break;
		}
		auto const& labj = features[h * width + w];
		std::array<unsigned short, 2> xyj = { (unsigned short)h, (unsigned short)w };
		float Dc = dc(labj, labi) / (m * m);
		float Ds = ds(xyj, xyi) / (float)(S * S);
		float D = Dc + Ds;
		auto& thisClus = clusterAss[h * width + w];
		if (D < thisClus.second) {
			int old = thisClus.first;
			thisClus.first = c;
			thisClus.second = D;
			res = true;
		}
	}
	return res;
}

bool GD::SLIC::assignment() {
	bool clusterChanged = false;
	size_t width = raster->width(), height = raster->height();
	for (int c = 0; c < centroids.size(); ++ c) {
		auto& center = centroids[c];
		int i = center.first, j = center.second;
		auto const& labi = features[i * width + j];
        auto xyi = std::array<unsigned short, 2>{ (unsigned short)i, (unsigned short)j };
		for (int h = i - S; h < i + S; ++h) {
			if (h < 0) {
				continue;
			} else if (h == height) {
				break;
			}
			clusterChanged = task_innerHloop(c, h, j, S, width, features, labi, xyi, clusterAss);
		}
	}
	return clusterChanged;
}

bool GD::SLIC::update() {
	std::vector<std::vector<std::pair<size_t, float>>> poses_weights(segments);
	size_t Area = raster->Area();
	int width = raster->width(), height = raster->height();
	for (size_t i = 0; i < Area; ++i) {
		int h = i / width, w = i % width;
		auto const& thisClus = clusterAss[i];
		int c = thisClus.first;
		float d = thisClus.second;
		poses_weights[c].push_back({ i, d });
	}
	for (int i = 0; i < segments; ++i) {
		auto& center = centroids[i];
		long long resi = 0, resj = 0;
		double weights_sum = 0;
		size_t size = poses_weights[i].size();
		for (int j = 0; j < size; ++j) {
			auto pos = poses_weights[i][j].first;
			auto weight = poses_weights[i][j].second;
			weights_sum += weight;
			resi += pos / width * weight;
			resj += (pos % width) * weight;
			// resi += weights[i][j] * dists[i][j] / width;
			// resj += weights[i][j] * (dists[i][j] % width);
		}
		if (weights_sum == 0) {
			continue;
		}
		resi /= weights_sum;
		resj /= weights_sum;
		center.first = resi;
		center.second = resj;	
	}
	return true;
}

auto GD::SLIC::run()->RasterPtr {
	if (gradient_reflected == nullptr || raster == nullptr) {
		throw std::runtime_error("GD::SLIC::run() : gradient or raster is null");
	}
	if (gradient_reflected->width() != raster->width() || gradient_reflected->height() != raster->height()) {
		throw std::runtime_error("GD::SLIC::run() : gradient and raster have different sizes");
	}
	if (segments == 0) {
		throw std::runtime_error("GD::SLIC::run() : segments is 0");
	}
	int width = raster->width(), height = raster->height();
	GD::ImptcSamplerWithSobol sampler(width, height);
	GD::ImportanceSamplerHelper helper;
	helper.setSrc(gradient_reflected);
	auto [cdf_y_marginal, cdf_x_condition] = helper.calc_cdf_y_marginal_And_cdf_x_condition();
	sampler.set_cdf_x_condition(cdf_x_condition);
	sampler.set_cdf_y_marginal(cdf_y_marginal);
	this->centers = sampler.sample(segments);
	std::sort(centers.begin(), centers.end());
	class_result = getByType(GDT_Byte);
	class_result->init(nullptr, width, height, 3, GDT_Byte);
	calc_features();
	bool clusterChanged = true;
	while (clusterChanged) {
		clusterChanged = assignment();
		if (clusterChanged)
		update();
	}
	// Ϊÿ��cluster����һ����ɫ
	std::vector<std::array<unsigned char, 3>> colors;
	for (int i = 0; i < segments; i++) {
		colors.push_back(
			std::array<unsigned char, 3>{ (unsigned char)(rand() % 255), (unsigned char)(rand() % 255), (unsigned char)(rand() % 255)});
	}
	// ��ÿ�����ص����ɫ����Ϊ������cluster����ɫ
	size_t total = width * height;
	for (int c = 0; c < 3; ++c) {
		for (size_t i = 0; i < total; ++i) {
			class_result->set(i / width, i % width, c, colors[clusterAss[i].first][c]);
		}
	}
	return class_result;
}