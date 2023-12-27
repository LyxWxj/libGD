#include "feature.h"


using namespace GD;

Array_features GD::operator+(Array_features const& lhs, Array_features const& rhs) {
	Array_features res;
	for (int i = 0; i < lhs.size(); ++i) {
		res.push_back(lhs[i] + rhs[i]);
	}
	return res;
}

Array_features& GD::operator+=(Array_features& lhs, Array_features const& rhs) {
	for (int i = 0; i < lhs.size(); ++i) {
		lhs[i] += rhs[i];
	}
	return lhs;
}

Array_features& GD::operator/=(Array_features& lhs, float rhs) {
	for (int i = 0; i < lhs.size(); ++i) {
		lhs[i] /= rhs;
	}
	return lhs;
}

Array_features GD::operator/(Array_features const& lhs, Array_features const& rhs) {
	Array_features res;
	for (int i = 0; i < lhs.size(); ++i) {
		res.push_back(lhs[i] / rhs[i]);
	}
	return res;
}

Array_features GD::operator*(Array_features const& lhs, float rhs)
{
	Array_features res;
	for (int i = 0; i < lhs.size(); ++i) {
		res.push_back(lhs[i] * rhs);
	}
	return res;
}

Array_features& GD::operator*=(Array_features& lhs, float rhs)
{
	for (int i = 0; i < lhs.size(); ++i) {
		lhs[i] *= rhs;
	}
	return lhs;
}

Array_features& GD::operator/=(Array_features& lhs, Array_features const& rhs) {
	for (int i = 0; i < lhs.size(); ++i) {
		lhs[i] /= rhs[i];
	}
	return lhs;
}

std::vector<Array_features>& GD::operator/= (std::vector<Array_features>& lhs, Array_features const& rhs) {
	for (int i = 0; i < lhs.size(); ++i) {
		lhs[i] /= rhs;
	}
	return lhs;
}

Array_features GD::operator/(Array_features const& lhs, float rhs) {
	Array_features res;
	for (int i = 0; i < lhs.size(); ++i) {
		res.push_back(lhs[i] / rhs);
	}
	return res;
}

float GD::metrics::dist_eucl::operator() (Array_features const& lhs, Array_features const& rhs) {
	float res{ 0 };
	for (int i = 0; i < lhs.size(); ++i) {
		float diff = lhs[i] - rhs[i];
		res += diff * diff;
	}
	return res;
}

float GD::metrics::dist_Manht::operator() (Array_features const& lhs, Array_features const& rhs) {
	float res{ 0 };
	for (int i = 0; i < lhs.size(); ++i) {
		res += abs(lhs[i] - rhs[i]);
	}
	return res;
}

