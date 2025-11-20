#ifndef SRC_BEZIER_H
#define SRC_BEZIER_H

#include <iostream>
#include <vector>
#include <cmath>

namespace math_util
{

struct C_Point
{
  double x;        ///< x_value
  double y;        ///< y value
  double heading;  ///< rad
  C_Point()
  {
    x = 0.0;
    y = 0.0;
    heading = 0.0;
  }
  C_Point(const double & _x, const double & _y, const double & _heading)
  {
    x = _x;
    y = _y;
    heading = _heading;
  }
  C_Point(const double & _x, const double & _y)
  {
    x = _x;
    y = _y;
    heading = 0.0;
  }
  void Set(const double & _x, const double & _y, const double & _heading)
  {
    x = _x;
    y = _y;
    heading = _heading;
  }
  void SetX(const double & _x) { x = _x; }
  void SetY(const double & _y) { y = _y; }
  void SetHeading(const double & _heading) { heading = _heading; }

  double GetX() { return x; }
  double GetY() { return y; }
  double GetHeading() { return heading; }

  void operator=(const C_Point & p)
  {
    x = p.x;
    y = p.y;
    heading = p.heading;
  }
  C_Point operator-(const C_Point & p)
  {
    C_Point tmp;
    tmp.x = x - p.x;
    tmp.y = y - p.y;
    tmp.heading = heading - p.heading;
    return tmp;
  }
  C_Point operator+(const C_Point & p)
  {
    C_Point tmp;
    tmp.x = x + p.x;
    tmp.y = y + p.y;
    tmp.heading = heading + p.heading;
    return tmp;
  }
};

struct BezierConfig
{
  double min_radius;

  BezierConfig() { min_radius = 1.0; }
};

class Bezier
{
public:
  Bezier();
  void SetConfig(const BezierConfig & config);
  void Init(const C_Point & p1, const C_Point & p2, const double path_resolution);
  int GeneratePath();
  bool GetResult();
  std::vector<C_Point> GetResultPath() { return result_path_; }

private:
  BezierConfig config_;
  bool is_init_;
  C_Point start_pt_;
  C_Point end_pt_;
  std::vector<C_Point> result_path_;
  double path_resolution_;
};

}  // namespace math_util

#endif