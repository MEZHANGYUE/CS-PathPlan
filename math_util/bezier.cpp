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

  double k = 1.0 / 3.0;
  double min_r = config_.min_radius;
  
  // Iteratively adjust k to satisfy min turning radius if needed
  // We check curvature at t=0, 0.5, 1.0
  // Increasing k generally increases radius of curvature at endpoints for typical turns
  for (int iter = 0; iter < 10; ++iter) {
      p1.x = p0.x + std::cos(p0.heading) * dis_p0_p3 * k;
      p1.y = p0.y + std::sin(p0.heading) * dis_p0_p3 * k;
      p1.z = p0.z + (p3.z - p0.z) * 1.0 / 3.0; // Z interpolation usually linear-ish

      p2.x = p3.x - std::cos(p3.heading) * dis_p0_p3 * k;
      p2.y = p3.y - std::sin(p3.heading) * dis_p0_p3 * k;
      p2.z = p0.z + (p3.z - p0.z) * 2.0 / 3.0;

      if (min_r <= 1.0) break; // No constraint or very small

      // Check curvature
      bool satisfied = true;
      for (double t : {0.0, 0.5, 1.0}) {
          double it = 1.0 - t;
          // First derivative
          double dx = 3*it*it*(p1.x - p0.x) + 6*it*t*(p2.x - p1.x) + 3*t*t*(p3.x - p2.x);
          double dy = 3*it*it*(p1.y - p0.y) + 6*it*t*(p2.y - p1.y) + 3*t*t*(p3.y - p2.y);
          double dz = 3*it*it*(p1.z - p0.z) + 6*it*t*(p2.z - p1.z) + 3*t*t*(p3.z - p2.z);
          
          // Second derivative
          double ddx = 6*it*(p2.x - 2*p1.x + p0.x) + 6*t*(p3.x - 2*p2.x + p1.x);
          double ddy = 6*it*(p2.y - 2*p1.y + p0.y) + 6*t*(p3.y - 2*p2.y + p1.y);
          double ddz = 6*it*(p2.z - 2*p1.z + p0.z) + 6*t*(p3.z - 2*p2.z + p1.z);

          // Cross product magnitude |v x a|
          double cx = dy*ddz - dz*ddy;
          double cy = dz*ddx - dx*ddz;
          double cz = dx*ddy - dy*ddx;
          double cross_norm = std::sqrt(cx*cx + cy*cy + cz*cz);
          
          double vel_norm = std::sqrt(dx*dx + dy*dy + dz*dz);
          double vel_norm3 = vel_norm * vel_norm * vel_norm;

          if (vel_norm3 > 1e-6) {
              double curvature = cross_norm / vel_norm3;
              if (curvature > 1.0 / min_r) {
                  satisfied = false;
                  break;
              }
          }
      }

      if (satisfied) break;
      
      k += 0.02;
      if (k > 0.45) { // Limit k to avoid loops or cusps
          k = 0.45;
          break; 
      }
  }
  
  // Re-calculate final points with chosen k
  p1.x = p0.x + std::cos(p0.heading) * dis_p0_p3 * k;
  p1.y = p0.y + std::sin(p0.heading) * dis_p0_p3 * k;
  p1.z = p0.z + (p3.z - p0.z) * 1.0 / 3.0;

  p2.x = p3.x - std::cos(p3.heading) * dis_p0_p3 * k;
  p2.y = p3.y - std::sin(p3.heading) * dis_p0_p3 * k;
  p2.z = p0.z + (p3.z - p0.z) * 2.0 / 3.0;

  double dis = std::hypot(p2.x - p1.x, p2.y - p1.y) + dis_p0_p3 * 2.0 / 3.0;
  double resolution = path_resolution_ / dis;

  C_Point tmp;
  double t, it;
  for (t = 0.0; t <= 1.0; t += resolution) {
    it = 1.0 - t;
    tmp.x = it * it * it * p0.x + 3 * it * it * t * p1.x + 3 * it * t * t * p2.x + t * t * t * p3.x;
    tmp.y = it * it * it * p0.y + 3 * it * it * t * p1.y + 3 * it * t * t * p2.y + t * t * t * p3.y;
    tmp.z = it * it * it * p0.z + 3 * it * it * t * p1.z + 3 * it * t * t * p2.z + t * t * t * p3.z;
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

Eigen::MatrixXd Bezier::GenerateTrajectoryMatrix(const Eigen::MatrixXd &Path, const std::string &yaml_path, double sample_distance_override, double v_avg_override)
{
    if (Path.rows() < 2) {
        return Eigen::MatrixXd(0, 3);
    }

    double resolution = 1.0;
    if (sample_distance_override > 0) {
        resolution = sample_distance_override;
    }
    
    // Clear previous state
    result_path_.clear();
    
    std::vector<C_Point> all_points;
    int num_points = Path.rows();

    // Calculate headings
    std::vector<double> headings(num_points);
    for (int i = 0; i < num_points; ++i) {
        double dx, dy;
        if (i == 0) {
            dx = Path(1, 0) - Path(0, 0);
            dy = Path(1, 1) - Path(0, 1);
        } else if (i == num_points - 1) {
            dx = Path(i, 0) - Path(i - 1, 0);
            dy = Path(i, 1) - Path(i - 1, 1);
        } else {
            dx = Path(i + 1, 0) - Path(i - 1, 0);
            dy = Path(i + 1, 1) - Path(i - 1, 1);
        }
        headings[i] = std::atan2(dy, dx);
    }

    // Generate segments
    for (int i = 0; i < num_points - 1; ++i) {
        C_Point p1(Path(i, 0), Path(i, 1), Path(i, 2), headings[i]);
        C_Point p2(Path(i + 1, 0), Path(i + 1, 1), Path(i + 1, 2), headings[i + 1]);

        this->Init(p1, p2, resolution);
        if (this->GeneratePath() == 0) {
            std::vector<C_Point> segment_path = this->GetResultPath();
            // Skip first point for segments after the first one to avoid duplicates
            size_t start_k = (i == 0) ? 0 : 1;
            for (size_t k = start_k; k < segment_path.size(); ++k) {
                all_points.push_back(segment_path[k]);
            }
        } else {
            // Fallback: add end point if generation fails
             C_Point pt(Path(i+1, 0), Path(i+1, 1), Path(i+1, 2), 0);
             all_points.push_back(pt);
        }
    }

    // Convert to Eigen::MatrixXd
    Eigen::MatrixXd result(all_points.size(), 3);
    for (size_t i = 0; i < all_points.size(); ++i) {
        result(i, 0) = all_points[i].x;
        result(i, 1) = all_points[i].y;
        result(i, 2) = all_points[i].z;
    }

    return result;
}

}  // namespace math_util