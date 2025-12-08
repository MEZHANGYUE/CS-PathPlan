#ifndef MATH_UTIL_ALTITUDE_OPTIMIZER_HPP
#define MATH_UTIL_ALTITUDE_OPTIMIZER_HPP

#include <string>
#include <vector>
#include <Eigen/Dense>

// Simple, dependency-light elevation and costmap utilities plus an altitude optimizer
namespace math_util {

struct GeoTransform {
  double origin_x; // top-left X
  double pixel_w;  // pixel width (x resolution)
  double rot_x;    // rotation
  double origin_y; // top-left Y
  double rot_y;
  double pixel_h;  // pixel height (negative if north-up)
};

// ElevationMap: stores a single-band float elevation grid with geotransform
class ElevationMap {
public:
  ElevationMap();
  ~ElevationMap();

  // load a GeoTIFF using GDAL (requires HAVE_GDAL), or load a PGM/ASCII raster as fallback
  bool loadFromTiff(const std::string &path);
  bool loadFromPGM(const std::string &path);

  bool isValid() const { return valid_; }

  // print summary info about the loaded elevation map (dimensions, geotransform, data stats)
  void printInfo() const;

  int getWidth() const { return width_; }
  int getHeight() const { return height_; }
  double getResolution() const { return std::abs(geo_.pixel_w); }
  GeoTransform getGeoTransform() const { return geo_; }

  // get elevation at world coordinates (x,y) in same CRS as GeoTransform; bilinear interpolation
  bool getElevationAt(double x, double y, double &elev) const;

private:
  // Helper to perform in-code downsampling for large rasters
  // poDS_ptr is passed as void* to avoid exposing GDAL headers in this hpp
  bool performDownsampling(void* poDS_ptr, int full_w, int full_h, uint64_t bytes_needed, uint64_t target_bytes, const std::string& path);

  int width_;
  int height_;
  std::vector<float> data_; // row-major, top-left origin
  GeoTransform geo_;
  bool valid_;
};

// CostMap: simple cost map storing 0..254 costs and 255=no-information
class CostMap {
public:
  CostMap() : width_(0), height_(0), resolution_(0.0), origin_x_(0.0), origin_y_(0.0) {}
  void create(int width, int height, double resolution, double origin_x, double origin_y);
  void setCost(int x, int y, unsigned char c);
  unsigned char getCost(int x, int y) const;
  int getWidth() const { return width_; }
  int getHeight() const { return height_; }
  double getResolution() const { return resolution_; }
  double getOriginX() const { return origin_x_; }
  double getOriginY() const { return origin_y_; }

  // fill from elevation map by linear scaling between elev_min..elev_max -> 0..254
  void fromElevationMap(const ElevationMap &emap, double elev_min, double elev_max);

private:
  int width_, height_;
  double resolution_, origin_x_, origin_y_;
  std::vector<unsigned char> costs_; // row-major, top-left origin
};

// AltitudeOptimizer: optimize waypoint altitudes z given xy coordinates and terrain
// Minimizes quadratic objective: lambda_smooth * ||D2 z||^2 + lambda_follow * ||z - elev||^2 + lambda_avoid * w .* (z - (elev + safety))^2
class AltitudeOptimizer {
public:
  AltitudeOptimizer();
    //无人机球体直径,最大爬升下滑比例,
  struct Params {
    double lambda_smooth = 1.0; // smoothing weight
    double lambda_follow = 10.0; // follow terrain weight (aim for z ~= elev + uav_R)
    double max_climb_rate = 2.0; // max vertical change per horizontal meter (m/m)
    double uav_R = 0.5; // UAV effective radius/height clearance (meters)
  };

  // Optimize heights for input waypoints (each waypoint: [x,y,z])
  // xy coordinates are in same CRS as ElevationMap GeoTransform.
  // Returns optimized z values (same length as waypoints) and true on success
  bool load_elevation_data(const std::string &elevation_file);
  bool optimizeHeights(const std::vector<Eigen::Vector3d> &waypoints,
                       const ElevationMap &emap,
                       const CostMap &cmap,
                       const Params &p,
                       std::vector<double> &out_z) const;
  // access loaded elevation map (if any)
  const ElevationMap & getElevationMap() const { return elevation_map_; }
  
private:
  // store loaded elevation map when load_elevation_data is called
  ElevationMap elevation_map_;
};

} // namespace math_util

#endif // MATH_UTIL_ALTITUDE_OPTIMIZER_HPP
