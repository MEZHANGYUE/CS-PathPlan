/**
 * @file
 * @brief Define the LineSegment2 class.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <utility>

#include "math_utils.hpp"
#include "vec2.hpp"

namespace math_util
{

/**
 * @class LineSegment2
 * @brief Line segment in 2-D.
 */
template <class T>
class LineSegment2
{
public:
  /**
   * @brief Empty constructor.
   */
  LineSegment2() { unit_direction_ = Vec2<T>(static_cast<T>(1.0), static_cast<T>(0.0)); }

  /**
   * @brief Constructor with start point and end point.
   * @param start The start point of the line segment.
   * @param end The end point of the line segment.
   */
  LineSegment2(const Vec2<T> & start, const Vec2<T> & end) : start_(start), end_(end)
  {
    const T dx = end_.x() - start_.x();
    const T dy = end_.y() - start_.y();
    length_ = std::sqrt(dx * dx + dy * dy);
    unit_direction_ =
      (length_ <= kGeometryEps<T>() ? Vec2<T>() : Vec2<T>(dx / length_, dy / length_));
    heading_ = unit_direction_.Angle();
  }

  LineSegment2(const T start_x, const T start_y, const T end_x, const T end_y)
  : LineSegment2(Vec2<T>(start_x, start_y), Vec2<T>(end_x, end_y))
  {
  }

  /**
   * @brief Get the start point.
   * @return The start point of the line segment.
   */
  const Vec2<T> & start() const { return start_; }

  /**
   * @brief Get the end point.
   * @return The end point of the line segment.
   */
  const Vec2<T> & end() const { return end_; }

  /**
   * @brief Get the unit direction from the start point to the end point.
   * @return The unit direction of the line segment.
   */
  const Vec2<T> & unit_direction() const { return unit_direction_; }

  /**
   * @brief Get the unit normal to the unit_direction above.
   * The normal is in the CCW direction, i.e. if the unit direction
   * is the x-axis, then the normal will be the y-axis.
   * @return The CCW unit normal of the line segment.
   */
  const Vec2<T> ortho_normal() const { return unit_direction_.rotate(M_PI / 2); }

  /**
   * @brief Get the center of the line segment.
   * @return The center of the line segment.
   */
  Vec2<T> center() const { return (start_ + end_) / static_cast<T>(2.0); }

  /**
   * @brief Get the heading of the line segment.
   * @return The heading, which is the angle between unit direction and x-axis.
   */
  inline T heading() const { return heading_; }

  /**
   * @brief Get the cosine of the heading.
   * @return The cosine of the heading.
   */
  inline T cos_heading() const { return unit_direction_.x(); }

  /**
   * @brief Get the sine of the heading.
   * @return The sine of the heading.
   */
  inline T sin_heading() const { return unit_direction_.y(); }

  /**
   * @brief Get the length of the line segment.
   * @return The length of the line segment.
   */
  inline T length() const { return length_; }

  /**
   * @brief Get the square of length of the line segment.
   * @return The square of length of the line segment.
   */
  inline T length_sqr() const { return length_ * length_; }

  /**
   * @brief Compute the shortest distance from a point on the line segment
   *        to a point in 2-D.
   * @param point The point to compute the distance to.
   * @return The shortest distance from points on the line segment to point.
   */
  T DistanceTo(const Vec2<T> & point) const
  {
    if (length_ <= kGeometryEps<T>()) {
      return point.DistanceTo(start_);
    }
    const T x0 = point.x() - start_.x();
    const T y0 = point.y() - start_.y();
    const T proj = x0 * unit_direction_.x() + y0 * unit_direction_.y();
    if (proj <= static_cast<T>(0.0)) {
      return std::sqrt(x0 * x0 + y0 * y0);
    }
    if (proj >= length_) {
      return point.DistanceTo(end_);
    }
    return std::abs(x0 * unit_direction_.y() - y0 * unit_direction_.x());
  }

  /**
   * @brief Compute the shortest distance from a point on the line segment
   *        to a point in 2-D, and get the nearest point on the line segment.
   * @param point The point to compute the distance to.
   * @param nearest_pt The nearest point on the line segment
   *        to the input point.
   * @return The shortest distance from points on the line segment
   *         to the input point.
   */
  T DistanceTo(const Vec2<T> & point, Vec2<T> * const nearest_pt) const
  {
    CHECK_NOTNULL(nearest_pt);
    if (length_ <= kGeometryEps<T>()) {
      *nearest_pt = start_;
      return point.DistanceTo(start_);
    }
    const T x0 = point.x() - start_.x();
    const T y0 = point.y() - start_.y();
    const T proj = x0 * unit_direction_.x() + y0 * unit_direction_.y();
    if (proj < static_cast<T>(0.0)) {
      *nearest_pt = start_;
      return std::sqrt(x0 * x0 + y0 * y0);
    }
    if (proj > length_) {
      *nearest_pt = end_;
      return point.DistanceTo(end_);
    }
    *nearest_pt = start_ + unit_direction_ * proj;
    return std::abs(x0 * unit_direction_.y() - y0 * unit_direction_.x());
  }

  /**
   * @brief Compute the square of the shortest distance from a point
   *        on the line segment to a point in 2-D.
   * @param point The point to compute the squared of the distance to.
   * @return The square of the shortest distance from points
   *         on the line segment to the input point.
   */
  T DistanceSquareTo(const Vec2<T> & point) const
  {
    if (length_ <= kGeometryEps<T>()) {
      return point.DistanceSquareTo(start_);
    }
    const T x0 = point.x() - start_.x();
    const T y0 = point.y() - start_.y();
    const T proj = x0 * unit_direction_.x() + y0 * unit_direction_.y();
    if (proj <= static_cast<T>(0.0)) {
      return Square(x0) + Square(y0);
    }
    if (proj >= length_) {
      return point.DistanceSquareTo(end_);
    }
    return Square(x0 * unit_direction_.y() - y0 * unit_direction_.x());
  }

  /**
   * @brief Compute the square of the shortest distance from a point
   *        on the line segment to a point in 2-D,
   *        and get the nearest point on the line segment.
   * @param point The point to compute the squared of the distance to.
   * @param nearest_pt The nearest point on the line segment
   *        to the input point.
   * @return The shortest distance from points on the line segment
   *         to the input point.
   */
  T DistanceSquareTo(const Vec2<T> & point, Vec2<T> * const nearest_pt) const
  {
    CHECK_NOTNULL(nearest_pt);
    if (length_ <= kGeometryEps<T>()) {
      *nearest_pt = start_;
      return point.DistanceSquareTo(start_);
    }
    const T x0 = point.x() - start_.x();
    const T y0 = point.y() - start_.y();
    const T proj = x0 * unit_direction_.x() + y0 * unit_direction_.y();
    if (proj <= static_cast<T>(0.0)) {
      *nearest_pt = start_;
      return Square(x0) + Square(y0);
    }
    if (proj >= length_) {
      *nearest_pt = end_;
      return point.DistanceSquareTo(end_);
    }
    *nearest_pt = start_ + unit_direction_ * proj;
    return Square(x0 * unit_direction_.y() - y0 * unit_direction_.x());
  }

  /**
   * @brief Check if a point is within the line segment.
   * @param point The point to check if it is within the line segment.
   * @return Whether the input point is within the line segment or not.
   */
  bool IsPointIn(const Vec2<T> & point) const
  {
    if (length_ <= kGeometryEps<T>()) {
      return std::abs(point.x() - start_.x()) <= kGeometryEps<T>() &&
             std::abs(point.y() - start_.y()) <= kGeometryEps<T>();
    }
    constexpr T kCrossProdEpsilon = 1e-5;  // use a slightly larger value for stability
    const T prod = CrossProd(point, start_, end_);
    if (std::abs(prod) > kCrossProdEpsilon) {
      return false;
    }
    return IsWithin(point.x(), start_.x(), end_.x()) && IsWithin(point.y(), start_.y(), end_.y());
  }

  /**
   * @brief Check if the line segment has an intersect
   *        with another line segment in 2-D.
   * @param other_segment The line segment to check if it has an intersect.
   * @return Whether the line segment has an intersect
   *         with the input other_segment.
   */
  bool HasIntersect(const LineSegment2 & other_segment) const
  {
    Vec2<T> point;
    return GetIntersect(other_segment, &point);
  }

  /**
   * @brief Compute the intersect with another line segment in 2-D if any.
   * @param other_segment The line segment to compute the intersect.
   * @param point the computed intersect between the line segment and
   *        the input other_segment.
   * @return Whether the line segment has an intersect
   *         with the input other_segment.
   */
  bool GetIntersect(const LineSegment2 & other_segment, Vec2<T> * const point) const
  {
    // CHECK_NOTNULL(point);
    if (IsPointIn(other_segment.start())) {
      *point = other_segment.start();
      return true;
    }
    if (IsPointIn(other_segment.end())) {
      *point = other_segment.end();
      return true;
    }
    if (other_segment.IsPointIn(start_)) {
      *point = start_;
      return true;
    }
    if (other_segment.IsPointIn(end_)) {
      *point = end_;
      return true;
    }
    if (length_ <= kGeometryEps<T>() || other_segment.length() <= kGeometryEps<T>()) {
      return false;
    }
    const T cc1 = CrossProd(start_, end_, other_segment.start());
    const T cc2 = CrossProd(start_, end_, other_segment.end());
    if (cc1 * cc2 >= -kGeometryEps<T>()) {
      return false;
    }
    const T cc3 = CrossProd(other_segment.start(), other_segment.end(), start_);
    const T cc4 = CrossProd(other_segment.start(), other_segment.end(), end_);
    if (cc3 * cc4 >= -kGeometryEps<T>()) {
      return false;
    }
    const T ratio = cc4 / (cc4 - cc3);
    *point = Vec2<T>(
      (start_.x() - end_.x()) * ratio + end_.x(), (start_.y() - end_.y()) * ratio + end_.y());
    return true;
  }

  /**
   * @brief Compute the projection of a vector onto the line segment.
   * @param point The end of the vector (starting from the start point of the
   *        line segment) to compute the projection onto the line segment.
   * @return The projection of the vector, which is from the start point of
   *         the line segment to the input point, onto the unit direction.
   */
  T ProjectOntoUnit(const Vec2<T> & point) const
  {
    return unit_direction_.InnerProd(point - start_);
  }

  /**
   * @brief Compute the cross product of a vector onto the line segment.
   * @param point The end of the vector (starting from the start point of the
   *        line segment) to compute the cross product onto the line segment.
   * @return The cross product of the unit direction and
   *         the vector, which is from the start point of
   *         the line segment to the input point.
   */
  T ProductOntoUnit(const Vec2<T> & point) const
  {
    return unit_direction_.CrossProd(point - start_);
  }

  /**
   * @brief Compute perpendicular foot of a point in 2-D on the straight line
   *        expanded from the line segment.
   * @param point The point to compute the perpendicular foot from.
   * @param foot_point The computed perpendicular foot from the input point to
   *        the straight line expanded from the line segment.
   * @return The distance from the input point to the perpendicular foot.
   */
  T GetPerpendicularFoot(const Vec2<T> & point, Vec2<T> * const foot_point) const
  {
    CHECK_NOTNULL(foot_point);
    if (length_ <= kGeometryEps<T>()) {
      *foot_point = start_;
      return point.DistanceTo(start_);
    }
    const T x0 = point.x() - start_.x();
    const T y0 = point.y() - start_.y();
    const T proj = x0 * unit_direction_.x() + y0 * unit_direction_.y();
    *foot_point = start_ + unit_direction_ * proj;
    return std::abs(x0 * unit_direction_.y() - y0 * unit_direction_.x());
  }

  /**
   * @brief Translate the line segment by the given vector.
   * @param delta Vector to translate by
   */
  void Translate(const Vec2<T> & delta)
  {
    start_ += delta;
    end_ += delta;
  }

  /**
   * @brief Rotate the line segment itself by the given vector.
   * @param theta Angle to rotate by
   */
  void RotateAboutOrigin(const T theta)
  {
    Rotate(Vec2<T>(static_cast<T>(0.0), static_cast<T>(0.0)), theta);
  }

  /**
   * @brief Rotate the line segment itself counter clock-wise
   * around the given center by the given angle.
   * @param center Center of rotation. theta Angle to rotate by.
   */
  void Rotate(const Vec2<T> & center, const T theta)
  {
    Vec2<T> v_center_to_start = start_ - center;
    v_center_to_start.SelfRotate(theta);
    start_ = center + v_center_to_start;

    Vec2<T> v_center_to_end = end_ - center;
    v_center_to_end.SelfRotate(theta);
    end_ = center + v_center_to_end;

    const T dx = end_.x() - start_.x();
    const T dy = end_.y() - start_.y();

    unit_direction_ =
      (length_ <= kGeometryEps<T>() ? Vec2<T>() : Vec2<T>(dx / length_, dy / length_));
    heading_ = unit_direction_.Angle();
  }

  /**
   * @brief: checks whether an other segment crosses this segment from our left
   * hand side.
   * @param other
   * @param intersection
   * @return True if other crosses this from the left hand side.
   */
  bool CrossFromLeft(const LineSegment2<T> & other, Vec2<T> * const intersection) const
  {
    if (!GetIntersect(other, intersection)) {
      return false;
    }
    return CrossProd(start_, end_, other.start()) >= 0;
  }

  /**
   * @brief: checks whether an other segment crosses this segment from our right
   * hand side.
   * @param other
   * @param intersection
   * @return True if other crosses this from the right hand side.
   */
  bool CrossFromRight(const LineSegment2<T> & other, Vec2<T> * const intersection) const
  {
    if (!GetIntersect(other, intersection)) {
      return false;
    }
    return CrossProd(start_, end_, other.start()) <= 0;
  }

  /**
   * @brief Get the debug string including the essential information.
   * @return Information of the line segment for debugging.
   */
  std::string DebugString() const
  {
    std::ostringstream buf;
    buf << "segment2d ( start = " << start_.DebugString() << "  end = " << end_.DebugString()
        << " )";
    return buf.str();
  }

  /**
   * @brief Get the shortest distance to the other line segment
   * @param other
   * @return The shortest distance to the other line segment
   */
  T DistanceTo(const LineSegment2<T> & other) const
  {
    if (HasIntersect(other)) {
      return static_cast<T>(0.0);
    }
    T dist1 = DistanceTo(other.start());
    T dist2 = DistanceTo(other.end());
    T dist3 = other.DistanceTo(start());
    T dist4 = other.DistanceTo(end());
    return std::min<T>(dist1, std::min(dist2, std::min(dist3, dist4)));
  }

private:
  static bool IsWithin(T val, T bound1, T bound2)
  {
    if (bound1 > bound2) {
      std::swap(bound1, bound2);
    }
    return val >= bound1 - kGeometryEps<T>() && val <= bound2 + kGeometryEps<T>();
  }

private:
  Vec2<T> start_;
  Vec2<T> end_;
  Vec2<T> unit_direction_;
  T heading_ = static_cast<T>(0.0);
  T length_ = static_cast<T>(0.0);
};

}  // namespace math_util
