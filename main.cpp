#include <iostream>
#include "Matrix.hpp"
#include <queue>

struct State {
    size_t num_of_moves;
    SquareMatrix<int> matrix;
    bool operator<(const State& that) const {
        return num_of_moves < that.num_of_moves;
    }
};

enum type {
    main_diag,
    another_diag,
};

State calculateRes(const SquareMatrix<int>& matrix);
SquareMatrix<int> getNewMatrix(SquareMatrix<int> matrix, size_t i, size_t j, size_t n, type transpose_type);

int main() {
    size_t n;
    std::cin >> n;
    SquareMatrix<int> matrix(n);
    matrix.read();
    auto res = calculateRes(matrix);
    std::cout << res.num_of_moves << '\n' << res.matrix.getDeterminant() << '\n';
    res.matrix.print();
}

State calculateRes(const SquareMatrix<int>& matrix) {
    std::priority_queue<State> states;
    states.push({0, matrix});
    auto res = states.top();
    size_t k_moves = 0;
    for (size_t k = 0; k < (size_t) std::pow(matrix.getSize(), 4); ++k) {
        auto top_state = states.top();
        states.pop();
        if (top_state.matrix.getDeterminant() > res.matrix.getDeterminant()) {
            res = top_state;
        }
        auto curr_det = top_state.matrix.getDeterminant();
        for (size_t i = 0; i < matrix.getSize(); ++i) {
            for (size_t j = 0; j < matrix.getSize(); ++j) {
                for (size_t n = 0; n < std::min(matrix.getSize() - i, matrix.getSize() - j); ++n) {
                    auto first_new_matrix = getNewMatrix(matrix, i, j, n, main_diag);
                    auto second_new_matrix = getNewMatrix(matrix, i, j, n, another_diag);
                    if (first_new_matrix.getDeterminant() >= curr_det) {
                        states.push({top_state.num_of_moves + 1, first_new_matrix});
                    }
                    if (second_new_matrix.getDeterminant() >= curr_det) {
                        states.push({top_state.num_of_moves + 1, second_new_matrix});
                    }
                }
            }
        }
    }
    return res;
}

SquareMatrix<int> getNewMatrix(SquareMatrix<int> matrix, size_t i, size_t j, size_t n, type transpose_type) {
    switch (transpose_type) {
        case main_diag: {
            for (size_t f = i; f < i + n; ++f) {
                for (size_t s = j; s < f; ++s) {
                    std::swap(matrix((int64_t) f, (int64_t) s), matrix((int64_t) (j + f - i), (int64_t) (i + s - j)));
                }
            }
            return matrix;
        }
        default: {
            for (size_t f = i; f < i + n; ++f) {
                for (size_t s = j; s < f; ++s) {
                    std::swap(matrix((int64_t) f, (int64_t) s), matrix((int64_t) (n + i + j - 1 - f), (int64_t) (n + j + i - 1 - s)));
                }
            }
            return matrix;
        }
    }
}