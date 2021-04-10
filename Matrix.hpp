#ifndef KISS_MATRIX_HPP
#define KISS_MATRIX_HPP

#include <cmath>

class MatrixWrongSizeError : std::logic_error {
  public:
    MatrixWrongSizeError() : std::logic_error("Sizes are not equal.") {}
};

class MatrixIndexError : std::length_error {
  public:
    MatrixIndexError() : std::length_error("No such index.") {}
};

class MatrixIsDegenerateError : std::logic_error {
  public:
    MatrixIsDegenerateError() : std::logic_error("Matrix is degenerate.") {}
};

template <typename T> T getZero() {
    return T(0);
}

template <typename T> T getOne() {
    return T(1);
}

template <typename T>
class Matrix {
  private:
    size_t rowsCnt;
    size_t colsCnt;
  protected:
    T **array;
  public:
    Matrix<T>(size_t, size_t);
    Matrix(const Matrix<T> &);
    ~Matrix();
    Matrix<T> &operator=(const Matrix<T> &);
    size_t getRowsNumber() const;
    size_t getColumnsNumber() const;
    template <typename M>
    friend Matrix<M> &operator+=(Matrix<M> &, const Matrix<M> &);
    template <typename M>
    friend Matrix<M> &operator*=(Matrix<M> &, const M &);
    template <typename M>
    friend Matrix<M> operator*(const Matrix<M> &, const Matrix<M> &);
    Matrix<T> &transpose();
    Matrix<T> getTransposed() const;
    T operator()(int64_t, int64_t) const;
    T &operator()(int64_t, int64_t);
    template <typename M>
    friend std::ostream &operator<<(std::ostream &, const Matrix<M> &);
    template <typename M>
    friend std::istream &operator>>(std::istream &, Matrix<M> &);
    void elTransOfRows(size_t, size_t, const T &, bool);
};
template<typename T>
Matrix<T> operator*(const Matrix<T> &, const T &);
template<typename T>
Matrix<T> operator+(const Matrix<T> &, const Matrix<T> &);
template<typename T>
Matrix<T> &operator-=(Matrix<T> &, const Matrix<T> &);
template<typename T>
Matrix<T> operator*(const T &, const Matrix<T> &);
template<typename T>
Matrix<T> operator-(const Matrix<T> &, const Matrix<T> &);
template<typename T>
Matrix<T> &operator*=(Matrix<T> &, const Matrix<T> &);
template<typename T>
Matrix<T> &operator*=(Matrix<T> &, const int64_t &);
template<typename T>
Matrix<T> operator*(const Matrix<T> &, const int64_t &);
template<typename T>
Matrix<T> operator*(const int64_t &, const Matrix<T> &);

template <typename T>
class SquareMatrix : public Matrix<T> {
  private:
    size_t size;
  public:
    explicit SquareMatrix(size_t);
    SquareMatrix(const SquareMatrix<T> &);
    ~SquareMatrix() = default;
    size_t getSize() const;
    SquareMatrix<T> &operator=(const SquareMatrix<T> &);
    SquareMatrix<T> &transpose();
    SquareMatrix<T> getTransposed() const;
    SquareMatrix<T> getMinor(size_t, size_t) const;
    T getDeterminant() const;
    SquareMatrix<T> &invert();
    SquareMatrix<T> getInverse() const;
    T getTrace() const;
    void read();
    void print() const;
    int getMaxDeterminant() const;

};

template<typename T>
Matrix<T>::Matrix(size_t rows, size_t cols) {
    array = new T *[rows];
    for (size_t i = 0; i < rows; ++i) {
        array[i] = new T[cols];
        for (size_t j = 0; j < cols; ++j) {
            array[i][j] = getZero<T>();
        }
    }
    this->rowsCnt = rows;
    this->colsCnt = cols;
}

template<typename T>
Matrix<T>::~Matrix() {
    for (size_t i = 0; i < this->rowsCnt; ++i) {
        delete[] this->array[i];
    }
    delete[] this->array;
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T> &that) : Matrix<T>(that.rowsCnt, that.colsCnt) {
    for (size_t i = 0; i < this->rowsCnt; ++i) {
        for (size_t j = 0; j < this->colsCnt; ++j) {
            this->array[i][j] = that.array[i][j];
        }
    }
}

template<typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &that) {
    if (&that == this) return *this;
    this->~Matrix<T>();
    this->rowsCnt = that.rowsCnt;
    this->colsCnt = that.colsCnt;
    this->array = new T *[this->rowsCnt];
    for (size_t i = 0; i < rowsCnt; ++i) {
        this->array[i] = new T[colsCnt];
        for (size_t j = 0; j < this->colsCnt; ++j) {
            this->array[i][j] = that.array[i][j];
        }
    }
    return *this;
}

template<typename T>
size_t Matrix<T>::getRowsNumber() const {
    return this->rowsCnt;
}

template<typename T>
size_t Matrix<T>::getColumnsNumber() const {
    return this->colsCnt;
}

template<typename T>
Matrix<T> &operator+=(Matrix<T> &this_m, const Matrix<T> &m) {
    if ((this_m.colsCnt != m.colsCnt) || (this_m.rowsCnt != m.rowsCnt)) throw MatrixWrongSizeError();
    for (size_t i = 0; i < this_m.rowsCnt; ++i) {
        for (size_t j = 0; j < this_m.colsCnt; ++j) {
            this_m.array[i][j] += m.array[i][j];
        }
    }
    return this_m;
}

template<typename T>
Matrix<T> operator+(const Matrix<T> &m1, const Matrix<T> &m2) {
    Matrix<T> res = m1;
    res += m2;
    return res;
}

template<typename T>
Matrix<T> &operator-=(Matrix<T> &this_m, const Matrix<T> &m) {
    this_m += m * T(-1);
    return this_m;
}

template<typename T>
Matrix<T> operator*(const Matrix<T> &m, const T &num){
    Matrix<T> res = m;
    res *= num;
    return res;
}
template<typename T>
Matrix<T> operator*(const T &num, const Matrix<T> &m) {
    return m * num;
}

template<typename T>
Matrix<T> operator-(const Matrix<T> &m1, const Matrix<T> &m2) {
    return m1 + m2 * T(-1);
}

template<typename T>
Matrix<T> operator*(const Matrix<T> &m1, const Matrix<T> &m2) {
    if (m1.colsCnt != m2.rowsCnt) throw MatrixWrongSizeError();
    Matrix<T> res(m1.rowsCnt, m2.colsCnt);
    size_t l = m1.colsCnt;
    for (size_t i = 0; i < res.rowsCnt; ++i) {
        for (size_t j = 0; j < res.colsCnt; ++j) {
            T sum = getZero<T>();
            for (size_t k = 0; k < l; ++k) {
                sum += m1.array[i][k] * m2.array[k][j];
            }
            res.array[i][j] = sum;
        }
    }
    return res;
}

template<typename T>
Matrix<T> &Matrix<T>::transpose() {
    Matrix<T> temp(this->colsCnt, this->rowsCnt);
    for (size_t i = 0; i < this->rowsCnt; ++i) {
        for (size_t j = 0; j < this->colsCnt; ++j) {
            temp.array[j][i] = this->array[i][j];
        }
    }
    (*this) = temp;
    return *this;
}

template<typename T>
Matrix<T> Matrix<T>::getTransposed() const {
    Matrix res = *this;
    res.transpose();
    return res;
}

template<typename T>
T Matrix<T>::operator()(int64_t i, int64_t j) const{
    if (i >= this->rowsCnt || j >= this->colsCnt || i < 0 || j < 0) throw MatrixIndexError();
    return this->array[i][j];
}

template<typename T>
T &Matrix<T>::operator()(int64_t i, int64_t j) {
    if (i >= this->rowsCnt || j >= this->colsCnt || i < 0 || j < 0) throw MatrixIndexError();
    return this->array[i][j];
}

template<typename M>
std::ostream &operator<<(std::ostream &os, const Matrix<M> &m) {
    for (size_t i = 0; i < m.rowsCnt; ++i) {
        for (size_t j = 0; j < m.colsCnt; ++j) {
            os << m.array[i][j] << ' ';
        }
        os << '\n';
    }
    return os;
}

template<typename M>
std::istream &operator>>(std::istream &is, Matrix<M> &m) {
    for (size_t i = 0; i < m.rowsCnt; ++i) {
        for (size_t j = 0; j < m.colsCnt; ++j) {
            is >> m.array[i][j];
        }
    }
    return is;
}

template<typename T>
Matrix<T> &operator*=(Matrix<T> &m, const T &num) {
    for (size_t i = 0; i < m.rowsCnt; ++i) {
        for (size_t j = 0; j < m.colsCnt; ++j) {
            m.array[i][j] *= static_cast<T>(num);
        }
    }
    return m;
}

template<typename T>
void Matrix<T>::elTransOfRows(size_t row1, size_t row2, const T &num, bool swap) {
    T *buf = new T[this->colsCnt];
    for (size_t i = 0; i < this->colsCnt; ++i) {
        buf[i] = (*this)(row2, i) * num;
    }
    if (swap) {
        for (size_t i = 0; i < this->colsCnt; ++i) {
            (*this)(row2, i) = (*this)(row1, i);
            (*this)(row1, i) = buf[i];
        }
    } else {
        for (size_t i = 0; i < this->colsCnt; ++i) {
            (*this)(row1, i) += buf[i];
        }
    }
}

template<typename M>
Matrix<M> &operator*=(Matrix<M> &this_m, const Matrix<M> &m) {
    this_m = this_m * m;
    return this_m;
}

template<typename T>
Matrix<T> &operator*=(Matrix<T> &m, const int64_t &num) {
    return m *= T(num);
}
template<typename T>
Matrix<T> operator*(const Matrix<T> &m, const int64_t &num) {
    return m * T(num);
}
template<typename T>
Matrix<T> operator*(const int64_t &num, const Matrix<T> &m) {
    return m * num;
}

template<typename T>
SquareMatrix<T>::SquareMatrix(size_t size) : Matrix<T>::Matrix(size, size), size(size) {}

template<typename T>
SquareMatrix<T>::SquareMatrix(const SquareMatrix<T> &m) : Matrix<T>::Matrix(m), size(m.size) {}

template<typename T>
size_t SquareMatrix<T>::getSize() const {
    return this->size;
}

template<typename T>
SquareMatrix<T> &operator*=(SquareMatrix<T> &this_sm, const SquareMatrix<T> &sm) {
    Matrix<T> &this_m = this_sm;
    const Matrix<T> &m = sm;
    this_m *= m;
    return this_sm;
}

template<typename M>
SquareMatrix<M> &operator*=(SquareMatrix<M> &this_sm, const M &num) {
    Matrix<M> &this_m = this_sm;
    this_m *= num;
    return this_sm;
}

template<typename T>
SquareMatrix<T> operator*(const SquareMatrix<T> &sm, const T &num) {
    SquareMatrix<T> res = sm;
    res *= num;
    return res;
}

template<typename T>
SquareMatrix<T> operator*(const T &num, const SquareMatrix<T> &sm) {
    return sm * num;
}

template<typename T>
SquareMatrix<T> operator*(const SquareMatrix<T> &m1, const SquareMatrix<T> &m2) {
    SquareMatrix<T> res(m1);
    res *= m2;
    return res;
}

template<typename T>
SquareMatrix<T> &SquareMatrix<T>::operator=(const SquareMatrix<T> &that) {
    this->Matrix<T>::operator=(that);
    this->size = that.size;
    return *this;
}

template<typename T>
SquareMatrix<T> &SquareMatrix<T>::transpose() {
    this->Matrix<T>::transpose();
    return *this;
}

template<typename T>
SquareMatrix<T> SquareMatrix<T>::getTransposed() const {
    SquareMatrix<T> res = *this;
    return res.transpose();
}

template<typename T>
SquareMatrix<T> SquareMatrix<T>::getMinor(size_t i, size_t j) const {
    SquareMatrix<T> res(this->size - 1);
    for (size_t k = 0; k < i; ++k) {
        for (size_t l = 0; l < j; ++l) {
            res(k, l) = (*this)(k, l);
        }
    }
    for (size_t k = i + 1; k < this->size; ++k) {
        for (size_t l = 0; l < j; ++l) {
            res(k - 1, l) = (*this)(k, l);
        }
    }
    for (size_t k = 0; k < i; ++k) {
        for (size_t l = j + 1; l < this->size; ++l) {
            res(k, l - 1) = (*this)(k, l);
        }
    }
    for (size_t k = i + 1; k < this->size; ++k) {
        for (size_t l = j + 1; l < this->size; ++l) {
            res(k - 1, l - 1) = (*this)(k, l);
        }
    }
    return res;
}

template<typename T>
T SquareMatrix<T>::getDeterminant() const {
    T det = getZero<T>();
    if (this->size == 1) {
        return (*this)(0, 0);
    }
    for (int i = 0; i < this->size; ++i) {
        det += (*this)(0, i) * static_cast<T>(pow(-1, i)) * (this->getMinor(0, i)).getDeterminant();
    }
    return det;
}

template<typename T>
SquareMatrix<T> Identity(size_t s) {
    SquareMatrix<T> res(s);
    for (size_t i = 0; i < s; ++i) {
        res(i, i) = getOne<T>();
    }
    return res;
}

template<typename T>
SquareMatrix<T> &SquareMatrix<T>::invert() {
    SquareMatrix<T> id = Identity<T>(this->size);
    for (size_t i = 0; i < this->size; ++i) {
        if ((*this)(i, i) == getZero<T>()) {
            size_t j = i;
            while (j < this->size && (*this)(j, i) == getZero<T>()) {
                ++j;
            }
            if (j == this->size) throw MatrixIsDegenerateError();
            this->elTransOfRows(i, j, getOne<T>(), true);
            id.elTransOfRows(i, j, getOne<T>(), true);
        }
        id.elTransOfRows(i, i, getOne<T>() / (*this)(i, i), true);
        this->elTransOfRows(i, i, getOne<T>() / (*this)(i, i), true);
        for (size_t j = 0; j < i; ++j) {
            id.elTransOfRows(j, i, static_cast<T>(-1) * (*this)(j, i), false);
            this->elTransOfRows(j, i, static_cast<T>(-1) * (*this)(j, i), false);
        }
        for (size_t j = i + 1; j < this->size; ++j) {
            id.elTransOfRows(j, i, static_cast<T>(-1) * (*this)(j, i), false);
            this->elTransOfRows(j, i, static_cast<T>(-1) * (*this)(j, i), false);
        }
    }
    *this = id;
    return *this;
}

template<typename T>
SquareMatrix<T> SquareMatrix<T>::getInverse() const {
    SquareMatrix<T> res = *this;
    return res.invert();
}

template<typename T>
T SquareMatrix<T>::getTrace() const {
    T trace = getZero<T>();
    for (size_t i = 0; i < this->size; ++i) {
        trace += (*this)(i, i);
    }
    return trace;
}

template<typename T>
void SquareMatrix<T>::read() {
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            std::cin >> Matrix<T>::array[i][j];
        }
    }
}

template<typename T>
int SquareMatrix<T>::getMaxDeterminant() const {
    int upper_bound = 1;
    int num_of_zeros = 0;
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            if (Matrix<T>::array[i][j]) {
                upper_bound *= Matrix<T>::array[i][j];
            } else {
                ++num_of_zeros;
            }
        }
    }
    if (num_of_zeros > size * (size - 1)) {
        return 0;
    }
    upper_bound = std::abs(upper_bound);
    int lower_bound = 0;

}

template<typename T>
void SquareMatrix<T>::print() const {
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            std::cout << Matrix<T>::array[i][j] << " ";
        }
        std::cout << "\n";
    }
}

template<typename T>
SquareMatrix<T> &operator+=(SquareMatrix<T> &this_sm, const SquareMatrix<T> &sm) {
    Matrix<T> &this_m = this_sm;
    const Matrix<T> &m = sm;
    this_m += m;
    return this_sm;
}

template<typename T>
SquareMatrix<T> &operator-=(SquareMatrix<T> &this_sm, const SquareMatrix<T> &sm) {
    this_sm += T(-1) * sm;
    return this_sm;
}

template<typename T>
SquareMatrix<T> operator+(const SquareMatrix<T> &sm1, const SquareMatrix<T> &sm2) {
    SquareMatrix<T> res = sm1;
    res += sm2;
    return res;
}

template<typename T>
SquareMatrix<T> operator-(const SquareMatrix<T> &m1, const SquareMatrix<T> &m2) {
    SquareMatrix<T> res = m1;
    res -= m2;
    return res;
}

template<typename T>
SquareMatrix<T> &operator*=(SquareMatrix<T> &sm, const int64_t &num) {
    return sm *= T(num);
}
template<typename T>
SquareMatrix<T> operator*(const SquareMatrix<T> &sm, const int64_t &num) {
    return sm * T(num);
}
template<typename T>
SquareMatrix<T> operator*(const int64_t &num, const SquareMatrix<T> &sm) {
    return sm * num;
}

#endif //KISS_MATRIX_HPP
