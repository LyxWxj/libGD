#pragma once

#ifndef _ADD_RESIDUAL_H
#define _ADD_RESIDUAL_H

#include "image_fusion.h"

//template<typename CostFunctor, int... Is>
//auto CreateCostFunction(std::vector<GD::Array_features> const& rgbs, GD::Array_features const& Pans, std::integer_sequence<int, Is...>) {
//	return new ceres::AutoDiffCostFunction<CostFunctor, 1, (Is + 1)...>(new CostFunctor(rgbs, Pans));
//}
//
//template<typename CostFunctor, size_t... Is>
//auto CreateCostFunctions(std::vector<GD::Array_features> const& rgbs, GD::Array_features const& Pans, std::index_sequence<Is...>) {
//	return std::array<ceres::CostFunction*, sizeof...(Is)>
//	{
//		new ceres::
//			AutoDiffCostFunction<CostFunctor, 1, (Is + 1)>(new CostFunctor(rgbs, Pans))...
//	};
//}
//
//template<typename CostFunctor>
//auto CreateCostFunctions(std::vector<GD::Array_features> const& rgbs, GD::Array_features const& Pans) {
//	return CreateCostFunctions<CostFunctor>(rgbs, Pans, std::make_index_sequence<100>());
//}

template<typename CostFunctor, size_t size>
void AddResidual(
	ceres::Problem* problem,
	std::array<double, size>& alphas,
	std::vector<GD::Array_features> const& rgbs,
	GD::Array_features const& Pans) {
	problem->AddResidualBlock(
		new ceres::AutoDiffCostFunction<CostFunctor, 1, size>(new CostFunctor(rgbs, Pans)),
		nullptr, alphas.data()
	);
}

#endif