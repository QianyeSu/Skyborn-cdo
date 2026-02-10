/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef VECTOR3D_H
#define VECTOR3D_H

#define _USE_MATH_DEFINES

#include <cstdio>
#include <cmath>

// clang-format off
class  // Vector3d
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
Vector3d
// clang-format on
{
private:
  double X{ 0.0 };
  double Y{ 0.0 };
  double Z{ 0.0 };

public:
  Vector3d(Vector3d const &) = default;
  Vector3d(Vector3d &&) = default;
  Vector3d &operator=(Vector3d const &) = default;
  Vector3d &operator=(Vector3d &&) = default;

  Vector3d() noexcept {}
  constexpr explicit Vector3d(double const &x, double const &y, double const &z) : X(x), Y(y), Z(z) {}

  [[nodiscard]] constexpr double &
  operator[](size_t index) noexcept
  {
    switch (index)
      {
      case 0: return X;
      case 1: return Y;
      case 2: return Z;
      }

    return X;
  }

  /*
  [[nodiscard]] constexpr const double &
  operator[](size_t index) noexcept
  {
    switch (index)
      {
      case 0: return X;
      case 1: return Y;
      case 2: return Z;
      }

    return X;
  }
  */
  [[nodiscard]] friend constexpr auto
  operator+(Vector3d const &lhs, Vector3d const &rhs) noexcept
  {
    return Vector3d{ lhs.X + rhs.X, lhs.Y + rhs.Y, lhs.Z + rhs.Z };
  }

  [[nodiscard]] friend constexpr auto
  operator-(Vector3d const &lhs, Vector3d const &rhs) noexcept
  {
    return Vector3d{ lhs.X - rhs.X, lhs.Y - rhs.Y, lhs.Z - rhs.Z };
  }

  // clang-format off
  [[nodiscard]] constexpr Vector3d operator+() const { return *this; }
  [[nodiscard]] constexpr Vector3d operator-() const { return Vector3d(-X, -Y, -Z); }
  // [[nodiscard]] constexpr Vector3d operator-() const { return { -X, -Y, -Z }; }
  // clang-format on

  // Calculate the cross/outer/vector product
  [[nodiscard]] friend constexpr auto
  operator%(Vector3d const &lhs, Vector3d const &rhs) noexcept
  {
    return Vector3d{ lhs.Y * rhs.Z - lhs.Z * rhs.Y, lhs.Z * rhs.X - lhs.X * rhs.Z, lhs.X * rhs.Y - lhs.Y * rhs.X };
  }

  // Division by scalars
  friend constexpr Vector3d
  operator/(Vector3d const &vec, const double &scalar) noexcept
  {
    return Vector3d(vec.X / scalar, vec.Y / scalar, vec.Z / scalar);
  }

  constexpr Vector3d &
  operator/=(double scalar) noexcept
  {
    X /= scalar;
    Y /= scalar;
    Z /= scalar;
    return *this;
  }

  // Calculate the dot/inner/scalar  product
  constexpr double
  operator*(Vector3d const &other) const noexcept
  {
    return (X * other.X) + (Y * other.Y) + (Z * other.Z);
  }

  double
  magnitude() const noexcept
  {
    return std::hypot(X, Y, Z);
  }

  Vector3d
  normalised() const noexcept
  {
    return Vector3d(*this) / magnitude();
  }

  void
  d_normalize()
  {
    auto dnorm = std::sqrt(*this * *this);
    *this = *this / dnorm;
  }

  double
  longitude() const noexcept
  {
    return std::atan2(Y, X);
  }

  double
  latitude() const noexcept
  {
    return M_PI_2 - std::acos(Z);
  }
};

static inline Vector3d
vector_product(Vector3d const &v0, Vector3d const &v1, Vector3d const &v2)
{
  // e1, e2: edges of the underlying planar triangle: v1-v0 ands v2-v0, respectively
  auto e1 = v1 - v0;
  auto e2 = v2 - v0;
  auto cu = e1 % e2;
  if ((cu * v0) < 0.0) cu = -cu;
  cu.d_normalize();
  return cu;
}

static inline Vector3d
circum_center_mean(Vector3d const &v0, Vector3d const &v1, Vector3d const &v2)
{
  /*
    v0, v1, v2: the coordinates of the three triangle vertices (_dmo,nit vectors) in
    counter clockwise order center: the coordinates of the circumcenter unless co-linear
  */
  // cu0, cu1, cu2: vector product of center:  e1 x e2

  auto cu0 = vector_product(v0, v1, v2);
  auto cu1 = vector_product(v1, v2, v0);
  auto cu2 = vector_product(v2, v0, v1);

  auto center = cu0 + cu1 + cu2;
  center.d_normalize();
  return center;
}

#endif
