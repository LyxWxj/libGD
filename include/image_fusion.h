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
		//����IHS�任��ͼ���ںϻ���
	public:
		IHSFusion() = default;
		virtual ~IHSFusion() = default;
		virtual GD::RasterPtr FusionOnCPU(RasterPtr src1, RasterPtr src2, size_t workers = 10);
	};

	class GIHSFusion : public IHSFusion {
	//GIHS/FIHS ����ʧ��
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
		// ���Կ˷�GIHS/FIHS�Ĺ���ʧ������ ����PAN�����ܱƽ�ǿ�ȷ��� �ɼ��ٹ���ʧ��ļ���
	protected:
		bool useimportanceSampling = false; // �Ƿ�ʹ����Ҫ�Բ���
		bool useSobol = false; // �Ƿ�ʹ��Sobol�Ͳ������в���
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
		 // ��IAIHS�Ļ����������ȫ�µı�Ե�������W_i(t)
	public:
				IAIHSFusionWithAlpha() = default;
		virtual ~IAIHSFusionWithAlpha() = default;
		virtual GD::RasterPtr FusionOnCPU(RasterPtr src1, RasterPtr src2, size_t workers = 1) override;
	};
}

#endif // __IMG_FUSION_H__