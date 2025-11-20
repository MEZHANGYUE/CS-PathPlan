/**
 * @file
 * @brief Defines the Vec2 class.
 */

#pragma once

#include <cmath>
#include <sstream>
#include <string>

namespace math_util
{

template <typename T>
constexpr inline T kGeometryEps();

template <>
constexpr inline double kGeometryEps()
{
  return 1e-10;
}

template <>
constexpr inline float kGeometryEps()
{
  return 1e-8f;
}

/**
 * @class Vec2
 * @brief Implements a templated class of 2-dimensional vectors.
 */
template <typename T>
class Vec2
{
public:
  //! Constructor which takes x- and y-coordinates.
  constexpr Vec2(const T x, const T y) noexcept : x_(x), y_(y) {}

  //! Constructor from PointENU object.
  // explicit Vec2(const common::PointENU & p) noexcept
  // : x_(static_cast<T>(p.x())), y_(static_cast<T>(p.y()))
  // {
  // }

  //! Constructor returning the zero vector.
  constexpr Vec2() noexcept : Vec2(static_cast<T>(0.0), static_cast<T>(0.0)) {}

  //! Creates a unit-vector with a given angle to the positive x semi-axis
  static Vec2 CreateUnitVec2d(const T angle) { return Vec2(std::cos(angle), std::sin(angle)); }

  //! Getter for x component
  inline T x() const { return x_; }

  //! Getter for y component
  inline T y() const { return y_; }

  //! Setter for x component
  inline void set_x(const T x) { x_ = x; }

  //! Setter for y component
  inline void set_y(const T y) { y_ = y; }

  //! Gets the length of the vector
  inline T Length() const { return std::sqrt(x_ * x_ + y_ * y_); }

  //! Gets the squared length of the vector
  inline T LengthSquare() const { return x_ * x_ + y_ * y_; }

  //! Gets the angle between the vector and the positive x semi-axis
  inline T Angle() const { return std::atan2(y_, x_); }

  //! Returns the unit vector that is co-linear with this vector
  void Normalize()
  {
    const T l = Length();
    if (l > kGeometryEps<T>()) {
      x_ /= l;
      y_ /= l;
    }
  }

  //! Returns the distance to the given vector
  T DistanceTo(const Vec2 & other) const
  {
    const T dx = x_ - other.x_;
    const T dy = y_ - other.y_;
    return std::sqrt(dx * dx + dy * dy);
  }

  //! Returns the squared distance to the given vector
  T DistanceSquareTo(const Vec2 & other) const
  {
    const T dx = x_ - other.x_;
    const T dy = y_ - other.y_;
    return dx * dx + dy * dy;
  }

  //! Returns the "cross" product between these two Vec2 (non-standard).
  T CrossProd(const Vec2 & other) const { return x_ * other.y_ - y_ * other.x_; }

  //! Returns the inner product between these two Vec2.
  T InnerProd(const Vec2 & other) const { return x_ * other.x_ + y_ * other.y_; }

  //! rotate the vector by angle.
  Vec2 rotate(const T angle) const
  {
    return Vec2(
      x_ * std::cos(angle) - y_ * std::sin(angle), x_ * std::sin(angle) + y_ * std::cos(angle));
  }

  //! rotate the vector itself by angle.
  void SelfRotate(const T angle)
  {
    const T tmp_x = x_;
    x_ = x_ * std::cos(angle) - y_ * std::sin(angle);
    y_ = tmp_x * std::sin(angle) + y_ * std::cos(angle);
  }

  //! transform the vector to the new frame defined by the origin of
  // the new frame presented in the current frame and the angle of rotation
  // from the current frame to the new frame
  Vec2 TransformToNewFrame(const Vec2 & origin, const T angle) const
  {
    Vec2 tmp(x_ - origin.x(), y_ - origin.y());
    tmp.SelfRotate(-angle);
    return tmp;
  }

  //! rotate the vector around the given origin point, clockwise
  Vec2 RotateAroundPoint(const Vec2 & origin, const T angle) const
  {
    Vec2 tmp = TransformToNewFrame(origin, angle);
    return Vec2(tmp.x() + origin.x(), tmp.y() + origin.y());
  }

  //! transform the vector itself to the new frame defined by the origin of
  // the new frame presented in the current frame and the angle of rotation
  // from the current frame to the new frame
  void SelfTransform(const Vec2 & origin, const T angle)
  {
    x_ -= origin.x();
    y_ -= origin.y();
    SelfRotate(-angle);
  }

  //! Sums two Vec2
  Vec2 operator+(const Vec2 & other) const { return Vec2(x_ + other.x(), y_ + other.y()); }

  //! Subtracts two Vec2
  Vec2 operator-(const Vec2 & other) const { return Vec2(x_ - other.x(), y_ - other.y()); }

  //! Multiplies Vec2 by a scalar
  Vec2 operator*(const T ratio) const { return Vec2(x_ * ratio, y_ * ratio); }

  //! Divides Vec2 by a scalar
  Vec2 operator/(const T ratio) const { return Vec2(x_ / ratio, y_ / ratio); }

  //! Sums another Vec2 to the current one
  Vec2 & operator+=(const Vec2 & other)
  {
    x_ += other.x();
    y_ += other.y();
    return *this;
  }

  //! Subtracts another Vec2 to the current one
  Vec2 & operator-=(const Vec2 & other)
  {
    x_ -= other.x();
    y_ -= other.y();
    return *this;
  }

  //! Multiplies this Vec2 by a scalar
  Vec2 & operator*=(const T ratio)
  {
    x_ *= ratio;
    y_ *= ratio;
    return *this;
  }

  //! Divides this Vec2 by a scalar
  Vec2 & operator/=(const T ratio)
  {
    x_ /= ratio;
    y_ /= ratio;
    return *this;
  }

  //! Compares two Vec2
  inline bool operator==(const Vec2 & other) const
  {
    return (
      std::abs(x_ - other.x()) < kGeometryEps<T>() && std::abs(y_ - other.y()) < kGeometryEps<T>());
  }
  inline bool operator!=(const Vec2 & other) const { return !operator==(other); }

  // Negates Vec2
  Vec2 operator-() const { return Vec2(-x_, -y_); }

  //! Return a PointENU representation of this Vec2
  // PointENU ToPointENU() const
  // {
  //   PointENU p;
  //   p.set_x(x_);
  //   p.set_y(y_);
  //   return p;
  // }

  //! Returns a human-readable string representing this object
  std::string DebugString() const
  {
    std::ostringstream buf;
    buf << "vec2d ( x = " << x_ << "  y = " << y_ << " )";
    return buf.str();
  }

protected:
  T x_ = static_cast<T>(0.0);
  T y_ = static_cast<T>(0.0);
};

//! Multiplies the given Vec2 by a given scalar
template <typename TScaler, typename TVec>
Vec2<TVec> operator*(const TScaler ratio, const Vec2<TVec> & vec)
{
  return vec * static_cast<TVec>(ratio);
}

}  // namespace math_util
