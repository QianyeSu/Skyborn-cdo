#ifndef MATRIX_VIEW_H
#define MATRIX_VIEW_H

// Modified code from https://github.com/pwwiur/Matrix

#include <stddef.h>

template <typename T>
class MatrixView
{
private:
  T *m_array{};
  size_t m_numRows{};
  size_t m_numColumns{};

  class VectorView
  {
  public:
    VectorView(T *arr, size_t columns, size_t x) : px_arr{ arr }, px_columns{ columns }, px_x{ x } {}

    // clang-format off
    T &operator[](size_t y) const { return px_arr[y + px_columns * px_x]; }
    // clang-format on

  private:
    T *px_arr{};
    size_t px_columns{};
    size_t px_x{};
  };

public:
  MatrixView(T *array, size_t rows, size_t columns) : m_array{ array }, m_numRows{ rows }, m_numColumns{ columns }
  {
    if (rows == 0 || columns == 0) throw "Matrix initialization arguments rows and columns must be greater than zero.";
  }

  // clang-format off
  const VectorView operator[](size_t x) const { return VectorView(m_array, m_numColumns, x); }
  VectorView operator[](size_t x) { return VectorView(m_array, m_numColumns, x); }
  // clang-format on
};

#endif  // MATRIX_VIEW_H
