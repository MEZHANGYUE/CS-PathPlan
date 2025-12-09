#ifndef UAV_PATH_PLANNING_H
#define UAV_PATH_PLANNING_H
#include "json.hpp"
#include "elog.h"
#include <yaml-cpp/yaml.h>
// #include "config.h"
#include "math_util/coordinate_transform.hpp"
#include <string>
#include <vector>
#include<iostream>
#include "math_util/minimum_snap.hpp"
// #include "math_util/altitude_optimizer.hpp"
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std;
using json = nlohmann::json;
using namespace math_util;

struct GeoTransform {
  double origin_x; // top-left X
  double pixel_w;  // pixel width (x resolution)
  double rot_x;    // rotation
  double origin_y; // top-left Y
  double rot_y;
  double pixel_h;  // pixel height (negative if north-up)
};

struct InputData
{
    double distance_points;
    double leader_speed = 30.0; // m/s, average speed override read from input JSON (formerly V_avg)
    double leader_fly_high;
    int formation_model;
    int formation_using;
    int uav_leader_id;
    std::pair<double, double> height_list;
    std::pair<double, double> ready_high_list;
    std::vector<WGS84Coord> high_zhandou_point_wgs84;
    std::vector<WGS84Coord> leader_midway_point_wgs84;
    std::vector<WGS84Coord> ready_zone;
    WGS84Coord uav_leader_start_point_wgs84 = {0.0, 0.0, 0.0};
    std::pair<int, std::pair<int, int>> uavs_plane_data;
};

struct OutputData
{
    double distance_points;
    double leader_fly_high;
    int formation_model;
    int formation_using;
    int uav_leader_id;
    std::pair<double, double> height_list;
    std::pair<double, double> ready_high_list;
    std::vector<WGS84Coord> high_zhandou_point_wgs84;
    std::vector<WGS84Coord> leader_midway_point_wgs84;
    std::vector<WGS84Coord> ready_zone;
    WGS84Coord uav_leader_start_point_wgs84 = {0.0, 0.0, 0.0};
    std::pair<int, std::pair<int, int>> uavs_plane_data;
};
// WGS84椭球体参数
constexpr double WGS84_A = 6378137.0;           // 长半轴 (米)
constexpr double WGS84_E2 = 0.006694379990141;  // 第一偏心率平方
constexpr double WGS84_E = 0.0818191908426;     // 第一偏心率

// 计算卯酉圈曲率半径
inline double calcN(double lat_rad) {
    double sin_lat = sin(lat_rad);
    return WGS84_A / sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat);
}

// 经纬度坐标结构体
struct WGS84Point {
    double lon; // 经度 (度)
    double lat; // 纬度 (度)
    double alt; // 高度 (米)
};

// 东北天坐标结构体
struct ENUPoint {
    double east;  // 东向 (米)
    double north; // 北向 (米)
    double up;    // 天向 (米)
};

// ECEF坐标结构体
struct ECEFPoint {
    double x;
    double y;
    double z;
};

// 将角度转换为弧度
inline double deg2rad(double deg) {
    return deg * M_PI / 180.0;
}

// 将弧度转换为角度
inline double rad2deg(double rad) {
    return rad * 180.0 / M_PI;
}

// UavPathPlanner: 将 uavPathPlanning 中的功能封装为类
class UavPathPlanner {
public:
    UavPathPlanner();
    ~UavPathPlanner();
    std::vector<ENUPoint> Trajectory_ENU={}; //单独优化高度时使用,非空时才能调用高度优化
    // Minisnap 轨迹生成接口
    std::vector<ENUPoint> Minisnap_3D(std::vector<ENUPoint> origin_waypoints, double distance_, double V_avg_override = -1.0);
    std::vector<ENUPoint> Minisnap_EN(std::vector<ENUPoint> origin_waypoints, double distance_, double V_avg_override = -1.0);

    // 主规划接口
    bool getPlan(json &input_json, json &output_json, bool use3D = true);
    //高度优化接口
    bool runAltitudeOptimization(const std::string &elev_file);
    // 辅助函数
    bool loadData(InputData &input_data, json &input_json);
    bool putWGS84ToJson(json &j, const std::string &key, const std::vector<WGS84Point> &traj);
    json generateFollowerTrajectories(const json &input_json, const InputData &input_data,
                                      const std::vector<ENUPoint> &Trajectory_ENU,
                                      const std::vector<WGS84Point> &Trajectory_WGS84);

    // 坐标变换
    ENUPoint wgs84ToENU(const WGS84Point& target, const WGS84Point& reference);
    WGS84Point enuToWGS84(const ENUPoint& enu, const WGS84Point& reference);
    std::vector<ENUPoint> wgs84ToENU_Batch(const std::vector<WGS84Point>& targets, const WGS84Point& reference);
    std::vector<WGS84Point> enuToWGS84_Batch(const std::vector<ENUPoint>& targets, const WGS84Point& reference);
    WGS84Point ecefToWGS84(const ECEFPoint& ecef);
    std::array<std::array<double, 3>, 3> computeENURotationMatrix(double lat_rad, double lon_rad);
    std::array<std::array<double, 3>, 3> computeENURotationMatrixInverse(double lat_rad, double lon_rad);
    ENUPoint ecefToENU(const ECEFPoint& delta_ecef, double ref_lat_rad, double ref_lon_rad);
    ECEFPoint enuToECEF(const ENUPoint& enu, double ref_lat_rad, double ref_lon_rad);

    // 保存 JSON
    static inline bool saveJsonToFile(const json &j, const std::string &filename) {
        try {
            std::ofstream file(filename);
            if (!file.is_open()) {
                std::cerr << "Error: Cannot open file " << filename << std::endl;
                return false;
            }
            file << j.dump(4);
            file.close();
            std::cout << "Successfully saved JSON to " << filename << std::endl;
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Error saving JSON to file: " << e.what() << std::endl;
            return false;
        }
    }

private:
    // 内部状态
    TrajectoryGeneratorTool generator_;
    WGS84Point origin_;
    // 提取的高度优化器函数仍在 cpp 中实现为私有方法
    // 现在只输入高程文件路径（例如 .tif），函数直接使用类成员 Trajectory_ENU 进行高度优化
    // bool runAltitudeOptimization(const std::string &elev_file);
    // 内部辅助：从 JSON 读取 ENU 路径与 distance
    std::vector<ENUPoint> getENUFromJSON(const json& j, const std::string& key);
    double getDistanceFromJSON(const json& j, const std::string& key);

    // Altitude Optimization Params
    struct AltitudeParams {
        double lambda_smooth = 1.0; // smoothing weight
        double lambda_follow = 10.0; // follow terrain weight (aim for z ~= elev + uav_R)
        double max_climb_rate = 2.0; // max vertical change per horizontal meter (m/m)
        double uav_R = 0.5; // UAV effective radius/height clearance (meters)
    };

    // Elevation Map Data
    int elev_width_ = 0;
    int elev_height_ = 0;
    std::vector<float> elev_data_; // row-major, top-left origin
    GeoTransform elev_geo_;
    bool elev_valid_ = false;

    // Helpers for elevation
    bool performDownsampling(void* poDS_ptr, int full_w, int full_h, uint64_t bytes_needed, uint64_t target_bytes, const std::string& path);
    bool getElevationAt(double x, double y, double &elev) const;
    bool optimizeHeights(const std::vector<Eigen::Vector3d> &waypoints, const AltitudeParams &p, std::vector<double> &out_z);

    // ECEF转经纬度（迭代法）
    ECEFPoint wgs84ToECEF(const WGS84Point& lla);

    // Altitude Optimization methods
    bool loadElevationData(const std::string &path);
    
    // Elevation Map Info
    void printElevationInfo() const;
    bool isElevationValid() const { return elev_valid_; }
};
#endif
