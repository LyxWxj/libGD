#pragma once

#ifndef __IMG_FUSION_H__
#define __IMG_FUSION_H__

#include "raster.h"
#include "feature.h"
#include <Eigen/Core>
#include <Eigen/Eigen>
#include "rand_Engine.h"
#include "sampler.h"
#include <ceres/ceres.h>

namespace GD {

	template<size_t N_residual, size_t N_params>
	struct GetCostFunction {

	};

	class IHSFusion {
		//基于IHS变换的图像融合基类
	public:
		IHSFusion() = default;
		virtual ~IHSFusion() = default;
		virtual GD::RasterPtr FusionOnCPU(RasterPtr src1, RasterPtr src2, size_t workers = 10);
	};

	class GIHSFusion : public IHSFusion {
	//GIHS/FIHS 光谱失真
	public:
				GIHSFusion() = default;
		virtual ~GIHSFusion() = default;
		virtual GD::RasterPtr FusionOnCPU(RasterPtr src1, RasterPtr src2, size_t workers = 10) override;
	};

	class AIHSFusion : public IHSFusion {
	protected:
		struct CostFunctorA {
			CostFunctorA(GD::Array_features rgb, float Pan)  {
				this->rgb = rgb;
				this->Pan = Pan;
			}
			template<typename T>
			bool operator()(const T* const alpha, T* residual) const {
				T temp = T(0);
				for (int i = 0; i < 3; ++i) {
					temp += ceres::abs(alpha[i] * (T)rgb[i] - T(Pan));
				}
				*residual = temp;
				return true;
			}
		private:
			GD::Array_features rgb;
			float Pan;
		};

		struct CostFunctorB {
			CostFunctorB(std::vector<GD::Array_features> const& rgbs, std::vector<float> const& Pan) {
				this->rgbs = rgbs;
				this->Pans = Pan;
			}
			CostFunctorB(std::vector<GD::Array_features>&& rgbs, std::vector<float>&& Pan) {
				this->rgbs = rgbs;
				this->Pans = Pan;
			}
			template<typename T>
			bool operator()(const T* const alpha, T* residual) const;
		private:
			std::vector<GD::Array_features> rgbs;
			std::vector<float> Pans;
		};
		// Adaptive IHS Fusion
		// 可以克服GIHS/FIHS的光谱失真问题 基于PAN尽可能逼近强度分量 可减少光谱失真的假设
	protected:
		bool useimportanceSampling = false; // 是否使用重要性采样
		bool useSobol = false; // 是否使用Sobol低差异序列采样
		double gamma = 1.f / 512.f;
		double lambda = 2 * 1e-8;
		std::ofstream* os = nullptr;
	protected:
		template<size_t size>
		void optimize_alpha(std::vector<size_t> const& points, GD::RasterPtr src, GD::RasterPtr PAN, std::array<double, size>& alpha);

	public:
		AIHSFusion() = default;
		virtual ~AIHSFusion() = default;
		virtual GD::RasterPtr FusionOnCPU(RasterPtr src1, RasterPtr src2, size_t workers = 10) override;
		virtual void setRecord(std::ofstream & os);
		virtual void setUseImportanceSampling(bool use);
		virtual void setUseSobol(bool use);
	};

	class IAIHSFusionWithAlpha : public AIHSFusion {
		 // 在IAIHS的基础上提出了全新的边缘检测算子W_i(t)
	public:
				IAIHSFusionWithAlpha() = default;
		virtual ~IAIHSFusionWithAlpha() = default;
		virtual GD::RasterPtr FusionOnCPU(RasterPtr src1, RasterPtr src2, size_t workers = 1) override;
	};
}

#endif // __IMG_FUSION_H__