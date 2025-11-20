#pragma once

#include <iostream>
#include <cmath>
#include <functional>
#include <cassert>
#include <utility>
#include "vec2.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923  // �� / 2
#endif

namespace math_util
{

using namespace std;

inline float wrapTo2Pi(float a)
{
  bool was_neg = a < 0;
  a = fmod(a, (float)(2.0 * M_PI));
  if (was_neg) a += (float)(2.0 * M_PI);
  return a;
}

//        inline float wrapToPi(float a) {
//            return wrapTo2Pi(a + (float) (M_PI)) - (float) (M_PI);
//        }
inline float wrapToPi(float a)
{
  a = std::fmod(a + M_PI, 2.0 * M_PI);
  if (a < 0.0) {
    a += (2.0 * M_PI);
  }
  return a - M_PI;
}

inline float degToRad(float a) { return wrapTo2Pi(a) * 0.0174532925f; }

inline float radToDeg(float a) { return wrapTo2Pi(a) * 57.2957795f; }

inline std::pair<double, double> cartesian2Polar(double x, double y)
{
  double r = std::sqrt(x * x + y * y);
  double theta = std::atan2(y, x);
  return std::make_pair(r, theta);
}

class QuadPoly
{
public:
  volatile float a, b, c;
  static const char * TAG;

  QuadPoly(float a, float b, float c) { setPoly(a, b, c); }

  ~QuadPoly(){};

  void setPoly(float a_, float b_, float c_)
  {
    a = a_;
    b = b_;
    c = c_;
  }

  float max() { return (4 * a * c - b * b) / (2 * a); }

  float rootPos()
  {
    //        LOGD("QuadPoly::rootPos");
    float s = getSqrItem();
    //!! Important fix for clang!
    //!! if a is very small but not zero, O2 solver will always return 0.
    if (fabsf(a) <= 0.0001f) {
      if (b != 0) {
        float solve_x = -c / b;
        //                LOGD("O1 Solve: %f", solve_x);
        return solve_x;
      }
      //            LOGD("NAN 0");
      return NAN;
    }

    if (std::isnan(s) == 0) {
      float solve_x = (float)((-b + sqrt(s)) / (2 * a));
      //            LOGD("O2 Solve: %f, [a:%f, b:%f, c:%f, sqrt(s):%f]",
      //                 solve_x, a, b, c, sqrt(s));
      return solve_x;
    }
    //        LOGD("NAN 1");
    return NAN;
  }

  float getSqrItem()
  {
    float tmp = (b * b - 4 * a * c);

    if (tmp < 0) {
      //            std::cout << "No real root.\n";
      //            toSrting();
      return NAN;
    }
    return tmp;
  }

  float eval(float x) { return a * x * x + b * x + c; }

  void toSrting() { std::cout << " a = " << a << "; b = " << b << "; c = " << c << std::endl; }
};

inline float FastInvSqrtF32(float x)
{
  long i;
  float x2, y;
  const float threehalfs = 1.5F;

  x2 = x * 0.5F;
  y = x;
  i = *(long *)&y;            // evil floating point bit level hacking
  i = 0x5f3759df - (i >> 1);  // what the fuck?
  y = *(float *)&i;
  y = y * (threehalfs - (x2 * y * y));  // 1st iteration
  //	y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed

  return y;
}

inline float FastSqrtF32(float x) { return 1.0f / FastInvSqrtF32(x); }

inline float NumDiffO1(std::function<float(float)> & f, float point, float inc)
{
  return (f(point + inc) - f(point)) / inc;
}

inline float NumDiffO2(std::function<float(float)> & f, float point, float inc)
{
  float g1 = (f(point + inc) - f(point)) / inc;
  float g2 = (f(point) - f(point - inc)) / inc;

  return (g2 - g1) / inc;
}

#define MIN(a, b) ((a) > (b) ? (b) : (a))
#define MAX(a, b) ((a) < (b) ? (b) : (a))
#define ARRAY_LEN(arr) (sizeof(arr) / sizeof(arr[0]))
#define SQR_SUM(a, b) ((a) * (a) + (b) * (b))

inline float HuberLoss(float x, float thres = 0.4)
{
  if (std::fabs(x) < thres) {
    return x * x;
  }
  return (std::fabs(x) - thres) * 2.0f * thres + (thres * thres);
}

typedef struct
{
  float range[2];
  int n;
} Linspace;

inline std::vector<float> linspace(Linspace & input)
{
  assert(input.range[1] >= input.range[0]);

  std::vector<float> result;

  float step = 0;

  if (input.n > 1) {
    step = (input.range[1] - input.range[0]) / (float)(input.n - 1);
  }

  float x = input.range[0];

  for (int i = 0; i < input.n; i++, x += step) {
    result.emplace_back(x);
  }

  return result;
}

/**
 * @brief Cross product between two 2-D vectors from the common start point,
 *        and end at two other points.
 * @param start_point The common start point of two vectors in 2-D.
 * @param end_point_1 The end point of the first vector.
 * @param end_point_2 The end point of the second vector.
 *
 * @return The cross product result.
 */
template <class T>
inline T CrossProd(
  const Vec2<T> & start_point, const Vec2<T> & end_point_1, const Vec2<T> & end_point_2)
{
  return (end_point_1 - start_point).CrossProd(end_point_2 - start_point);
}

/**
 * @brief Inner product between two 2-D vectors from the common start point,
 *        and end at two other points.
 * @param start_point The common start point of two vectors in 2-D.
 * @param end_point_1 The end point of the first vector.
 * @param end_point_2 The end point of the second vector.
 *
 * @return The inner product result.
 */
template <class T>
inline T InnerProd(
  const Vec2<T> & start_point, const Vec2<T> & end_point_1, const Vec2<T> & end_point_2)
{
  return (end_point_1 - start_point).InnerProd(end_point_2 - start_point);
}

/**
 * @brief Cross product between two vectors.
 *        One vector is formed by 1st and 2nd parameters of the function.
 *        The other vector is formed by 3rd and 4th parameters of the function.
 * @param x0 The x coordinate of the first vector.
 * @param y0 The y coordinate of the first vector.
 * @param x1 The x coordinate of the second vector.
 * @param y1 The y coordinate of the second vector.
 *
 * @return The cross product result.
 */
template <class T>
inline T CrossProd(const T x0, const T y0, const T x1, const T y1)
{
  return x0 * y1 - x1 * y0;
}

/**
 * @brief Inner product between two vectors.
 *        One vector is formed by 1st and 2nd parameters of the function.
 *        The other vector is formed by 3rd and 4th parameters of the function.
 * @param x0 The x coordinate of the first vector.
 * @param y0 The y coordinate of the first vector.
 * @param x1 The x coordinate of the second vector.
 * @param y1 The y coordinate of the second vector.
 *
 * @return The inner product result.
 */
template <class T>
inline T InnerProd(const T x0, const T y0, const T x1, const T y1)
{
  return x0 * x1 + y0 * y1;
}

/**
 * @brief Wrap angle to [0, 2 * PI).
 * @param angle the original value of the angle.
 * @return The wrapped value of the angle.
 */
double WrapAngle(const double angle);

/**
 * @brief Normalize angle to [-PI, PI).
 * @param angle the original value of the angle.
 * @return The normalized value of the angle.
 */
double NormalizeAngle(const double angle);

/**
 * @brief Normalize angle to [0, PI_2].
 * @param angle the original value of the angle.
 * @return The normalized value of the angle.
 */
double NormalizeToAcuteAngle(const double angle);

/**
 * @brief Normalize angle to [-180, 180).
 * @param angle the original value of the angle in degree.
 * @return The normalized value of the angle in degree.
 */
double NormalizeAngleInDeg(const double angle);

/**
 * @brief Calculate the difference between angle from and to
 * @param from the start angle
 * @param from the end angle
 * @return The difference between from and to. The range is between [-PI, PI).
 */
double AngleDiff(const double from, const double to);

/**
 * @brief Get a random integer between two integer values by a random seed.
 * @param s The lower bound of the random integer.
 * @param t The upper bound of the random integer.
 * @param random_seed The random seed.
 * @return A random integer between s and t based on the input random_seed.
 */
int RandomInt(const int s, const int t, unsigned int rand_seed = 1);

/**
 * @brief Get a random double between two integer values by a random seed.
 * @param s The lower bound of the random double.
 * @param t The upper bound of the random double.
 * @param random_seed The random seed.
 * @return A random double between s and t based on the input random_seed.
 */
double RandomDouble(const double s, const double t, unsigned int rand_seed = 1);

/**
 * @brief Compute squared value.
 * @param value The target value to get its squared value.
 * @return Squared value of the input value.
 */
template <typename T>
inline T Square(const T value)
{
  return value * value;
}

/**
 * @brief Clamp a value between two bounds.
 *        If the value goes beyond the bounds, return one of the bounds,
 *        otherwise, return the original value.
 * @param value The original value to be clamped.
 * @param bound1 One bound to clamp the value.
 * @param bound2 The other bound to clamp the value.
 * @return The clamped value.
 */
template <typename T>
T Clamp(const T value, T bound1, T bound2)
{
  if (bound1 > bound2) {
    std::swap(bound1, bound2);
  }

  if (value < bound1) {
    return bound1;
  } else if (value > bound2) {
    return bound2;
  }
  return value;
}

template <typename T>
T Clamp(const T value, T bound1, T bound2, bool * const is_clamped)
{
  if (bound1 > bound2) {
    std::swap(bound1, bound2);
  }

  if (value < bound1) {
    *is_clamped = true;
    return bound1;
  } else if (value > bound2) {
    *is_clamped = true;
    return bound2;
  }
  *is_clamped = false;
  return value;
}

// Gaussian
double Gaussian(const double u, const double std, const double x);

// Sigmoid
double Sigmoid(const double x);

// Rotate a 2d vector counter-clockwise by theta
// Eigen::Vector2d RotateVector2d(const Eigen::Vector2d & v_in, const double theta);

inline std::pair<double, double> RFUToFLU(const double x, const double y)
{
  return std::make_pair(y, -x);
}

inline std::pair<double, double> FLUToRFU(const double x, const double y)
{
  return std::make_pair(-y, x);
}

inline void L2Norm(int feat_dim, float * feat_data)
{
  if (feat_dim == 0) {
    return;
  }
  // feature normalization
  float l2norm = 0.0f;
  for (int i = 0; i < feat_dim; ++i) {
    l2norm += feat_data[i] * feat_data[i];
  }
  if (l2norm == 0) {
    float val = 1.f / std::sqrt(static_cast<float>(feat_dim));
    for (int i = 0; i < feat_dim; ++i) {
      feat_data[i] = val;
    }
  } else {
    l2norm = std::sqrt(l2norm);
    for (int i = 0; i < feat_dim; ++i) {
      feat_data[i] /= l2norm;
    }
  }
}

inline double L2Norm(const std::vector<double> & data)
{
  if (data.size() == 0) {
    return 0.0;
  }
  double l2norm = 0.0f;
  for (size_t i = 0; i < data.size(); ++i) {
    l2norm += data[i] * data[i];
  }
  return std::sqrt(l2norm);
}

// Cartesian coordinates to Polar coordinates
std::pair<double, double> Cartesian2Polar(double x, double y);

/**
 * @brief Compute barycentric coordinates of `p` inside `a`, `b` and `c`.
 *        The coordinates are allowed to be negative, which is useful for
 *        extrapolation.
 * @return `true` if barycentric coordiante is well-defined, otherwise `false`
 */
template <class T>
bool ComputeBarycentricCoordinate2(
  const Vec2<T> & p, const Vec2<T> & a, const Vec2<T> & b, const Vec2<T> & c, T * u, T * v, T * w)
{
  const auto v0 = b - a;
  const auto v1 = c - a;
  const auto v2 = p - a;
  const T denom = CrossProd(v0.x(), v0.y(), v1.x(), v1.y());
  constexpr T kEps = 1e-7;
  if (std::abs(denom) <= kEps) {
    // This means that a, b and c are co-linear
    return false;
  }
  *v = CrossProd(v2.x(), v2.y(), v1.x(), v1.y()) / denom;
  *w = CrossProd(v0.x(), v0.y(), v2.x(), v2.y()) / denom;
  *u = static_cast<T>(1.0) - *v - *w;
  return true;
}

}  // namespace math_util
