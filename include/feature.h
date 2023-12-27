#pragma once
#ifndef __FEATURE__H_
#define __FEATURE__H_

#include <vector>
#include <functional>

namespace GD {

	using Array_features = std::vector<float>;

	Array_features operator+(Array_features const& lhs, Array_features const& rhs);

	Array_features& operator+=(Array_features& lhs, Array_features const& rhs);

	Array_features operator/(Array_features const& lhs, float rhs);

	Array_features& operator/=(Array_features& lhs, float rhs);

	Array_features operator/(Array_features const& lhs, Array_features const& rhs);

	Array_features operator*(Array_features const& lhs, float rhs);

	Array_features& operator*=(Array_features& lhs, float rhs);

	template<typename T, typename Alloc>
	Array_features operator/(Array_features const& lhs, std::vector<T, Alloc> const& rhs) {
		Array_features res(lhs.size());
		for (size_t i = 0; i < lhs.size(); ++i) {
			res[i] = lhs[i] / rhs[i];
		}
		return res;
	}

	template<typename T, typename Alloc>
	Array_features& operator/=(Array_features& lhs, std::vector<T, Alloc> const& rhs) {
	for (size_t i = 0; i < lhs.size(); ++i) {
			lhs[i] /= rhs[i];
		}
		return lhs;
	}

	Array_features& operator/=(Array_features& lhs, Array_features const& rhs);

	std::vector<Array_features>& operator/= (std::vector<Array_features>& lhs, Array_features const& rhs);

	template<typename T, typename Alloc>
	std::vector<Array_features>& operator/= (std::vector<Array_features>& lhs, std::vector<T, Alloc> const& rhs) {
		for (auto& feature : lhs) {
			feature /= rhs;
		}
		return lhs;
	}

	namespace metrics {
		struct dist_eucl {
			float operator() (Array_features const&, Array_features const&);
		};

		struct dist_Manht {
			float operator() (Array_features const&, Array_features const&);
		};
	}
}



#endif // !__FEATURE_H_
