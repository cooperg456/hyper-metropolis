#pragma once

#include <vector>



namespace helpers {

	template <typename T>
	const std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>> & mat) {
		std::vector<std::vector<T>> trans;

		size_t max_row = 0;
		for (const auto & row : mat) {
			max_row = std::max(max_row, row.size());
		}

		trans.resize(max_row);

		for (size_t i = 0; i < mat.size(); i++) {
			for (size_t j = 0; j < mat[i].size(); j++) {
				trans[j].push_back(mat[i][j]);
			}
		}

		return trans;
	}

}