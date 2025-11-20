#include "linear_interpolation.hpp"
#include <algorithm>
#include <cmath>

namespace math_util
{

namespace
{

inline double NormalizeAngle(const double angle)
{
  double a = std::fmod(angle + M_PI, 2.0 * M_PI);
  if (a < 0.0) {
    a += (2.0 * M_PI);
  }
  return a - M_PI;
}

constexpr double kEpsilon = 1e-10;

}  // namespace

double slerp(const double a0, const double a1, const double ratio)
{
  return slerp(a0, 0.0, a1, 1.0, ratio);
}

double slerp(const double a0, const double t0, const double a1, const double t1, const double t)
{
  if (std::abs(t1 - t0) <= kEpsilon) {
    std::cout << "input time difference is too small" << std::endl;
    return NormalizeAngle(a0);
  }
  const double a0_n = NormalizeAngle(a0);
  const double a1_n = NormalizeAngle(a1);
  double d = a1_n - a0_n;
  if (d > M_PI) {
    d = d - 2 * M_PI;
  } else if (d < -M_PI) {
    d = d + 2 * M_PI;
  }

  const double r = (t - t0) / (t1 - t0);
  const double a = a0_n + d * r;
  return NormalizeAngle(a);
}

}