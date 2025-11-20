#include "bezier.hpp"

namespace math_util
{

Bezier::Bezier()
{
  is_init_ = false;
  result_path_.clear();
}

void Bezier::SetConfig(const BezierConfig & config)
{
  config_ = config;
  return;
}

void Bezier::Init(const C_Point & p1, const C_Point & p2, const double path_resolution)
{
  start_pt_ = p1;
  end_pt_ = p2;
  result_path_.clear();
  is_init_ = true;
  path_resolution_ = path_resolution;
  return;
}

int Bezier::GeneratePath()
{
  if (!is_init_) return -1;
  // calculate the ctrl point
  C_Point p0, p1, p2, p3;
  p0 = start_pt_;
  p3 = end_pt_;
  double dis_p0_p3 = std::hypot(p0.x - p3.x, p0.y - p3.y);
  if (dis_p0_p3 < 1e-1) return -1;

  p1.x = p0.x + std::cos(p0.heading) * dis_p0_p3 * 1.0 / 3.0;
  p1.y = p0.y + std::sin(p0.heading) * dis_p0_p3 * 1.0 / 3.0;
  p2.x = p3.x - std::cos(p3.heading) * dis_p0_p3 * 1.0 / 3.0;
  p2.y = p3.y - std::sin(p3.heading) * dis_p0_p3 * 1.0 / 3.0;

  double dis = std::hypot(p2.x - p1.x, p2.y - p1.y) + dis_p0_p3 * 2.0 / 3.0;
  double resolution = path_resolution_ / dis;

  C_Point tmp;
  double t, it;
  for (t = 0.0; t <= 1.0; t += resolution) {
    it = 1.0 - t;
    tmp.x = it * it * it * p0.x + 3 * it * it * t * p1.x + 3 * it * t * t * p2.x + t * t * t * p3.x;
    tmp.y = it * it * it * p0.y + 3 * it * it * t * p1.y + 3 * it * t * t * p2.y + t * t * t * p3.y;
    result_path_.push_back(tmp);
  }
  return 0;
}

bool Bezier::GetResult()
{
  if (!is_init_) return false;
  if (result_path_.empty()) return false;
  return true;
}

}  // namespace math_util