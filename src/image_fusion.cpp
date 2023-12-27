#include "image_fusion.h"
#include "raster.h"
#include "feature.h"
#include "templates.hpp"
#include <algorithm>
#include <ceres/ceres.h>
#include "addresidual.h"

const int N = 12;
const std::array<double, 4> ts = { 2.5, 3.5, 2.0, 4.0 };

GD::Matrix<float> SobelKenel({
        {0, -1, 0},
        {-1, 0, 1},
        {0, 1, 0}
    });

//#define CHECKRasterPtrBand(src1, src2) {if (src1->bands() != 3 && src2->bands() != 1) { \
//                    std::cerr << "Error: src1 bands > 3 src2 bands != 1" << std::endl; \
//                    return nullptr; \
//                } \
//                if (src1->bands() != 1 && src2->bands() > 3) { \
//                    std::cerr << "Error: src1 bands != 1 src2 bands > 3" << std::endl; \
//                    return nullptr; \
//                } \
//                if (src1->type() != src2->type()) { \
//                    std::cerr << "Error: src1 type != src2 type" << std::endl; \
//                    return nullptr; \
//                }}



GD::RasterPtr GD::IHSFusion::FusionOnCPU(GD::RasterPtr src1, GD::RasterPtr src2, size_t workers) {
    //基于IHS变换的图像融合基类
    // CHECKRasterPtrBand(src1, src2);
    if (src1->bands() == 1) swap(src1, src2);
    // src1 is RGB, src2 is PAN
    size_t width = src2->width(), height = src2->height();
    size_t bands = src1->bands();
    src1->resize(width, height);
    GD::RasterPtr dst = GD::getByType(GDT_Float32);
    dst->init(nullptr, width, height, bands, GDT_Float32);

// #pragma omp parallel for num_threads(workers)
    for (long long i = 0; i < height; ++i) {
        for (long long j = 0; j < width; ++j) {
            auto arr = src1->get(i, j);
            float b = arr[0], g = arr[1], r = arr[2];
            auto pan = src2->get(i, j)[0];
            Eigen::Vector3f rgb(r, g, b);
            Eigen::Vector3f ihs; GD::RGB2IHS(rgb, ihs);
            ihs[0] = pan;
            GD::IHS2RGB(ihs, rgb);
            dst->set(i, j, std::vector<float>{rgb[2], rgb[1], rgb[0]});
        }
    }
    return dst;
}

GD::RasterPtr GD::GIHSFusion::FusionOnCPU(GD::RasterPtr src1, GD::RasterPtr src2, size_t workers) {
    //GIHS/FIHS 光谱失真
    // CHECKRasterPtrBand(src1, src2);
    if (src1->bands() == 1) swap(src1, src2);
    // src1 is RGB, src2 is PAN
    size_t width = src2->width(), height = src2->height();
    size_t bands = src1->bands();
    src1->resize(width, height);
    GD::RasterPtr dst = GD::getByType(GDT_Float32);
    dst->init(nullptr, width, height, bands, GDT_Float32);

// #pragma omp parallel for num_threads(workers)
    for (long long i = 0; i < width; ++i) {
        for (long long j = 0; j < height; ++j) {
            auto arr = src1->get(i, j);
            float I = GD::sum(arr) / arr.size();
            auto pan = src2->get(i, j)[0];
            float delta = pan - I;
            for (auto& val : arr) {
                val += delta;
            }
            dst->set(i, j, arr);
        }
    }
    return dst;
}

template<typename Array_T>
static void Normalize(Array_T& arr) {
    float sum = 0;
    for (auto& val : arr) {
        sum += val;
    }
    for (auto& val : arr) {
        val /= sum;
    }
}

auto  flectance = [](float p) {
    return 1.f/(1.f + p);
};

template<typename Func>
auto flectance_img = [](GD::RasterPtr img, Func f) {
    // 对img原地修改
    size_t width = img->width(), height = img->height(), bands = img->bands();
    size_t area = width * height;
    for (int c = 0; c < bands; ++c) {
        for (long long i = 0; i < area; ++i) {
            auto val = img->get(i / width, i % width, c);
            val = f(val);
            img->set(i / width, i % width, c, val);
        }
    }
};

template<typename T>
bool GD::AIHSFusion::CostFunctorB::operator()(const T* const alpha, T* residual) const {
    T temp = T(0);
    const size_t size = rgbs.size();
    const size_t channels = rgbs[0].size();
    for (int i = 0; i < size; ++i) {
        auto const& rgb = rgbs[i];
        const float Pan = Pans[i];
        T I{ 0.0 };
        for (int channel = 0; channel < channels; ++channel) {
            I += alpha[channel] * (T)rgb[channel];
        }
        temp += ceres::pow(I - (T)Pan, 2);
    }
    *residual = temp/(T)size;
    return true;
}

template<size_t size>
void GD::AIHSFusion::optimize_alpha(std::vector<size_t> const& points, GD::RasterPtr src, GD::RasterPtr PAN, std::array<double, size>& alpha) {
    // 初始化Problem
    ceres::Problem problem;
    std::vector<GD::Array_features> rgbs;
    std::vector<float> Pans;
    for (auto pos : points) {
        pos %= PAN->Area();
        int i = pos / PAN->width(), j = pos % PAN->width();
        auto I = src->get(i, j);
        float P = PAN->get(i, j, 0);
        rgbs.emplace_back(std::move(I));
        Pans.emplace_back(std::move(P));
    }
    AddResidual<CostFunctorB>(&problem, alpha, rgbs, Pans);
    for (int i = 0; i < alpha.size(); ++i) {
        problem.SetParameterLowerBound(alpha.data(), i, 0.0);
        problem.SetParameterUpperBound(alpha.data(), i, 1.0);
    }
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
}

template<size_t size>
void processPixel_AIHS(size_t i, size_t j, GD::RasterPtr dPan, GD::RasterPtr dMs, GD::RasterPtr PAN, GD::RasterPtr Ms, std::array<double, size>const& Alphas, std::array<double, 4> const& ts, GD::RasterPtr dst) {
    auto h = [&](double x) {
        return exp(-1e-9 / (1e-10 + x));
        };
    float dp = dPan->get(i, j, 0);
    float P = PAN->get(i, j, 0);
    auto I = Ms->get(i, j);
    float P_I = P - GD::dot_product(Alphas, I);
    float add = h(dp) * P_I;
    for (int m = 0; m < size; ++m) {
        I[m] += add;
    }
    dst->set(i, j, I);
}

template<size_t size>
void processRowAIHS(int i, size_t width, GD::RasterPtr dPan,
    GD::RasterPtr dMs, GD::RasterPtr PAN,
    GD::RasterPtr Ms, 
    std::array<double, size>const& Alphas, 
    std::array<double, 4>const& ts, GD::RasterPtr dst) {
    for (size_t j = 0; j < width; ++j) {
        processPixel_AIHS(i, j, dPan, dMs, PAN, Ms, Alphas, ts, dst);
    }
}

void t1(GD::RasterPtr Ms, GD::RasterPtr& dMs) {
    dMs = Ms->Sobel();
    flectance_img<decltype(flectance)>(dMs, flectance);
}

void t2(GD::RasterPtr Pan, GD::RasterPtr* dPan) {
    *dPan = Pan->Sobel();;
}



GD::RasterPtr GD::AIHSFusion::FusionOnCPU(GD::RasterPtr src1, GD::RasterPtr src2, const size_t workers) {
    auto Ms = src1, PAN = src2;
    RasterPtr dPan, dMs;
    dMs = Ms->Sobel();
    flectance_img<decltype(flectance)>(dMs, flectance);
    std::thread t_(t2, PAN, &dPan);
    size_t width = PAN->width(), height = PAN->height();
    size_t bands = Ms->bands();
    src1->resize(width, height);
    GD::Sampler* sampler = nullptr;
    if (this->useimportanceSampling) {
        if (this->useSobol) {
            GD::ImptcSamplerWithSobol* Imsampler = new ImptcSamplerWithSobol(width, height);
            ImportanceSamplerHelper helper;
            helper.setSrc(dMs);
            auto [cdf_y_marginal, cdf_x_condition] = helper.calc_cdf_y_marginal_And_cdf_x_condition();
            Imsampler->set_cdf_x_condition(cdf_x_condition);
            Imsampler->set_cdf_y_marginal(cdf_y_marginal);
            sampler = Imsampler;
        }
        else {
            GD::ImptcSampler* Imsampler = new ImptcSampler(width, height);
            ImportanceSamplerHelper helper;
            helper.setSrc(dMs);
            auto [cdf_y_marginal, cdf_x_condition] = helper.calc_cdf_y_marginal_And_cdf_x_condition();
            Imsampler->set_cdf_x_condition(cdf_x_condition);
            Imsampler->set_cdf_y_marginal(cdf_y_marginal);
            sampler = Imsampler;
        }
    }
    else {
        if (this->useSobol) {
            sampler = new SobolSampler(width, height);
        }
        else {
            sampler = new UnifSampler(width, height);
        }
    }
    std::array<double, 3> Alphas;
    for (int i = 0; i < bands; ++i) {
        Alphas[i] = 1. / bands;
    }
    int sampleNum = N * sqrt(width * height);
    std::vector<size_t> samplePoses = sampler->sample(sampleNum);
    std::sort(samplePoses.begin(), samplePoses.end());
    optimize_alpha(samplePoses, Ms, PAN, Alphas);
    GD::RasterPtr dst = GD::getByType(GDT_Float32);
    dst->init(nullptr, width, height, bands, GDT_Float32);
    t_.join();
    std::vector<std::thread> threads;

    for (int i = 0; i < height;) {
        int upper = height > i + 50 ? i + 50 : height;
        for (int j = i; j < upper; ++j) {
            threads.emplace_back(std::thread(processRowAIHS<3>, j, width, dPan, dMs, PAN, Ms, Alphas, ts, dst));
        }

        for (int j = i; j < upper; ++j) {
            if (threads[j].joinable()) {
                threads[j].join();
            }
        }
        i += 50;
    }
    
    std::cout << "归一化的参数:";
    Normalize(Alphas);
    for (int i = 0; i < bands; ++i) {
        std::cout << Alphas[i] << "\t";
    }
    std::cout << "\n";
    delete sampler;
    for (int i = 0; i < threads.size(); ++i) {
        if (threads[i].joinable()) {
			threads[i].join();
		}
    }
    return dst;
}

void GD::AIHSFusion::setRecord(std::ofstream& os) {
    this->os = &os;
}

void GD::AIHSFusion::setUseImportanceSampling(bool use) {
    this->useimportanceSampling = use;
}

void GD::AIHSFusion::setUseSobol(bool use) {
    this->useSobol = use;
}

template<size_t size>
void processPixel_IAIHS(size_t i, size_t j, GD::RasterPtr dPan, GD::RasterPtr dMs, GD::RasterPtr PAN, GD::RasterPtr Ms, std::array<double, size>const& Alphas, std::array<double, 4>const& ts, GD::RasterPtr dst) {
    auto h = [&](double x) {
        return exp(-1e-9 / (1e-10 + x));
    };
    float dp = dPan->get(i, j, 0);
    auto dm = dMs->get(i, j);
    float P = PAN->get(i, j, 0);
    auto I = Ms->get(i, j);
    float P_I = P - GD::dot_product(Alphas, I);
    for (int m = 0; m < size; ++m) {
        float wt = 1 - Alphas[m] * h(dp) + Alphas[m] * h(dm[m]);
        I[m] += wt * P_I;
    }
    dst->set(i, j, I);
}

template<size_t size>
void processRowIAIHS(int i, size_t width, GD::RasterPtr dPan, GD::RasterPtr dMs, GD::RasterPtr PAN, GD::RasterPtr Ms, std::array<double, size>const& Alphas, std::array<double, 4>const& ts, GD::RasterPtr dst) {
    for (size_t j = 0; j < width; ++j) {
        processPixel_IAIHS(i, j, dPan, dMs, PAN, Ms, Alphas, ts, dst);
    }
}



GD::RasterPtr GD::IAIHSFusionWithAlpha::FusionOnCPU(GD::RasterPtr src1, GD::RasterPtr src2, size_t workers) {
    if (src1->bands() == 1) swap(src1, src2);
    // src1 is RGB, src2 is PAN
    auto Ms = src1, PAN = src2;
    RasterPtr dPan, dMs_f, dMs;
    dMs = Ms->Sobel();
    dMs_f = dMs->copy();
    flectance_img<decltype(flectance)>(dMs_f, flectance);
    auto t12 = [&]() {
        dPan = PAN->Sobel();
        };
    std::thread t_(t12);
    size_t width = PAN->width(), height = PAN->height();
    size_t bands = Ms->bands();
    src1->resize(width, height);

    GD::Sampler* sampler = nullptr;
    if (this->useimportanceSampling) {
        if (this->useSobol) {
            GD::ImptcSamplerWithSobol* Imsampler = new ImptcSamplerWithSobol(width, height);
            ImportanceSamplerHelper helper;
            helper.setSrc(dMs_f);
            auto [cdf_y_marginal, cdf_x_condition] = helper.calc_cdf_y_marginal_And_cdf_x_condition();
            Imsampler->set_cdf_x_condition(cdf_x_condition);
            Imsampler->set_cdf_y_marginal(cdf_y_marginal);
            sampler = Imsampler;
        }
        else {
            GD::ImptcSampler* Imsampler = new ImptcSampler(width, height);
            ImportanceSamplerHelper helper;
            helper.setSrc(dMs_f);
            auto [cdf_y_marginal, cdf_x_condition] = helper.calc_cdf_y_marginal_And_cdf_x_condition();
            Imsampler->set_cdf_x_condition(cdf_x_condition);
            Imsampler->set_cdf_y_marginal(cdf_y_marginal);
            sampler = Imsampler;
        }
    }
    else {
        if (this->useSobol) {
            sampler = new SobolSampler(width, height);
        }
        else {
            sampler = new UnifSampler(width, height);
        }
    }
    std::array<double, 3> Alphas;
    for (int i = 0; i < bands; ++i) {
        Alphas[i] = 1. / bands;
    }
    int sampleNum = N * sqrt(width * height);
    std::vector<size_t> samplePoses = sampler->sample(sampleNum);

    optimize_alpha(samplePoses, Ms, PAN, Alphas);

    GD::RasterPtr dst = GD::getByType(GDT_Float32);
    dst->init(nullptr, width, height, bands, GDT_Float32);
    t_.join();
    std::vector<std::thread> threads;


    for (int i = 0; i < height;) {
        int upper = height > i + 50 ? i + 50 : height;
        for (int j = i; j < upper; ++j) {
            threads.emplace_back(std::thread(processRowIAIHS<3>, j, width, dPan, dMs, PAN, Ms, Alphas, ts, dst));
        }

        for (int j = i; j < upper; ++j) {
            if (threads[j].joinable()) {
                threads[j].join();
            }
        }
        i += 50;
    }
    std::cout << "归一化的参数:";
    Normalize(Alphas);
    for (int i = 0; i < bands; ++i) {
        std::cout << Alphas[i] << "\t";
    }
    std::cout << "\n";
    delete sampler;
    for (int i = 0; i < threads.size(); ++i) {
        if (threads[i].joinable()) {
            threads[i].join();
        }
    }
    return dst;
}

