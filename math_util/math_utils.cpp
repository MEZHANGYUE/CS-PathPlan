#include "math_utils.hpp"

#include <cmath>
#include <utility>
#include <random>
namespace math_util
{

double WrapAngle(const double angle)
{
  const double new_angle = std::fmod(angle, M_PI * 2.0);
  return new_angle < 0 ? new_angle + M_PI * 2.0 : new_angle;
}

double NormalizeAngle(const double angle)
{
  double a = std::fmod(angle + M_PI, 2.0 * M_PI);
  if (a < 0.0) {
    a += (2.0 * M_PI);
  }
  return a - M_PI;
}

double NormalizeToAcuteAngle(const double angle)
{
  double a = std::fmod(angle + M_PI_2, M_PI);
  if (a < 0.0) {
    a += M_PI;
  }
  return std::fabs(a - M_PI_2);
}

double NormalizeAngleInDeg(const double angle)
{
  double a = angle;
  while (a < -180.0) {
    a += 360.0;
  }
  while (a >= 180.0) {
    a -= 360.0;
  }
  return a;
}

double AngleDiff(const double from, const double to) { return NormalizeAngle(to - from); }

#ifdef __linux__
int RandomInt(const int s, const int t, unsigned int rand_seed)
{
  if (s >= t) {
    return s;
  }
  return s + rand_r(&rand_seed) % (t - s + 1);
}

double RandomDouble(const double s, const double t, unsigned int rand_seed)
{
  return s + (t - s) / 16383.0 * (rand_r(&rand_seed) & 16383);
}
#else
// ��ƽ̨�̰߳�ȫ�����������
int RandomInt(const int s, const int t, unsigned int rand_seed)
{
    if (s >= t) return s;
    std::mt19937 gen(rand_seed);
    std::uniform_int_distribution<int> dist(s, t);
    return dist(gen);
}

// ��ƽ̨�̰߳�ȫ�������������
double RandomDouble(const double s, const double t, unsigned int rand_seed)
{
    std::mt19937 gen(rand_seed);
    std::uniform_real_distribution<double> dist(s, t);
    return dist(gen);
}
#endif
// Gaussian
double Gaussian(const double u, const double std, const double x)
{
  return (1.0 / std::sqrt(2 * M_PI * std * std)) * std::exp(-(x - u) * (x - u) / (2 * std * std));
}

// Sigmoid
double Sigmoid(const double x) { return 1.0 / (1.0 + std::exp(-x)); }

// Eigen::Vector2d RotateVector2d(const Eigen::Vector2d & v_in, const double theta)
// {
//   const double cos_theta = std::cos(theta);
//   const double sin_theta = std::sin(theta);

//   auto x = cos_theta * v_in.x() - sin_theta * v_in.y();
//   auto y = sin_theta * v_in.x() + cos_theta * v_in.y();

//   return {x, y};
// }

std::pair<double, double> Cartesian2Polar(double x, double y)
{
  double r = std::sqrt(x * x + y * y);
  double theta = std::atan2(y, x);
  return std::make_pair(r, theta);
}

}  // namespace math_util
