#ifndef ELEVATION_COST_MAP_H
#define ELEVATION_COST_MAP_H

#include <vector>
#include <string>
#include <limits>
#include <cstdint>

#ifdef HAVE_GDAL
#include <gdal_priv.h>
#include <cpl_conv.h>
#endif

struct GeoTransform {
  double origin_x; // top-left X
  double pixel_w;  // pixel width (x resolution)
  double rot_x;    // rotation
  double origin_y; // top-left Y
  double rot_y;
  double pixel_h;  // pixel height (negative if north-up)
};

class ElevationCostMap {
public:
    ElevationCostMap();
    ~ElevationCostMap();

    // Elevation map
    bool loadElevationData(const std::string &path);
    bool performDownsampling(void* poDS_ptr, int full_w, int full_h, uint64_t bytes_needed, uint64_t target_bytes, const std::string& path);
    void printElevationInfo() const;
    bool getElevationAt(double x, double y, double &elev) const;
    bool isElevationValid() const { return elev_valid_; }

    // Cost map (stores elevation-derived costs)
    void createCostMap(int width, int height, double resolution, double origin_x, double origin_y);
    void setCost(int x, int y, float c);
    float getCost(int x, int y) const;
    bool getCostAt(double x, double y, float &val) const;

private:
    // Elevation data
    int elev_width_ = 0;
    int elev_height_ = 0;
    std::vector<float> elev_data_; // row-major, top-left origin
    GeoTransform elev_geo_;
    bool elev_valid_ = false;

    // Cost map data
    int cost_width_ = 0;
    int cost_height_ = 0;
    double cost_resolution_ = 0.0;
    double cost_origin_x_ = 0.0;
    double cost_origin_y_ = 0.0;
    std::vector<float> cost_data_; // row-major, top-left origin
};

#endif
