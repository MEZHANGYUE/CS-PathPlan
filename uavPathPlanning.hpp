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
    double min_turning_radius = 0.0; // Minimum turning radius constraint
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
    // Bezier 轨迹生成接口
    std::vector<ENUPoint> Bezier_3D(std::vector<ENUPoint> origin_waypoints, double distance_, double V_avg_override = -1.0, double min_radius = 0.0);

    // 主规划接口
    bool getPlan(json &input_json, json &output_json, bool use3D = true, std::string algorithm = "minimum_snap");
    //高度优化接口
    bool runAltitudeOptimization(const std::string &elev_file, json &output_json, const json &input_json);
    // 辅助函数
    bool loadData(InputData &input_data, json &input_json);
    bool putWGS84ToJson(json &j, const std::string &key, const std::vector<WGS84Point> &traj);
    json generateFollowerTrajectories(const json &input_json, const InputData &input_data,
                                      const std::vector<ENUPoint> &Trajectory_ENU,
                                      const std::vector<WGS84Point> &Trajectory_WGS84);

    // 计算轨迹最小转弯半径
    double calculateMinTurningRadius(const std::vector<ENUPoint>& path);

    // 生成 arc-line-arc（相切圆 - 直线 - 相切圆）路径
    // p0: 起点, heading0: 起点初始朝向（弧度）
    // p1: 终点（中间点），p2: 用于确定 p1 的朝向（p1->p2 的方向为切线方向）
    // radius: 相切圆半径（来自 config.yaml -> path_planning.min_turning_radius）
    // resolution: 采样距离（米）
    std::vector<ENUPoint> generateArcLineArc(const ENUPoint &p0, double heading0,
                                             const ENUPoint &p1, const ENUPoint &p2,
                                             double radius, double resolution = 1.0);

    // 计算第二段过渡轨迹（切圆切入优化）并更新第三段巡逻轨迹
    void computeTransitionAndRotatePatrol(const ENUPoint& p0, double heading0, double minR, double resolution, 
                                          const std::vector<ENUPoint>& Patrol_Path, json& output_json);

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

    // CostMap: simple cost map storing height values (float)
    class CostMap {
    public:
      CostMap() : width_(0), height_(0), resolution_(0.0), origin_x_(0.0), origin_y_(0.0) {}
      void create(int width, int height, double resolution, double origin_x, double origin_y);
      void setCost(int x, int y, float c);
      float getCost(int x, int y) const;
      int getWidth() const { return width_; }
      int getHeight() const { return height_; }
      double getResolution() const { return resolution_; }
      double getOriginX() const { return origin_x_; }
      double getOriginY() const { return origin_y_; }
    
      // get value at world coordinates (x,y)
      bool getValueAt(double x, double y, float &val) const;
    
      // fill from elevation map (copy elevation values)
      // Note: In this merged version, we pass the UavPathPlanner instance or relevant data directly
      // But since CostMap is now nested or we just use UavPathPlanner's elevation data, 
      // we might not need a separate fromElevationMap if we access UavPathPlanner's data.
      // However, to keep logic similar to before, we can let it access UavPathPlanner's elevation data.
      // For simplicity, let's keep CostMap as a helper class or struct within UavPathPlanner or just use UavPathPlanner's methods.
      // The user asked to merge ElevationMap and CostMap classes into uavPathPlanning.
      // ElevationMap logic is already largely merged into UavPathPlanner (elev_data_, etc).
      // CostMap logic (grid of costs) can also be merged.
      
    private:
      int width_, height_;
      double resolution_, origin_x_, origin_y_;
      std::vector<float> costs_; // row-major, top-left origin
    };

private:
    // 内部状态
    TrajectoryGeneratorTool generator_;
    WGS84Point origin_;
    // 提取的高度优化器函数仍在 cpp 中实现为私有方法
    // 现在只输入高程文件路径（例如 .tif），函数直接使用类成员 Trajectory_ENU 进行高度优化
    // bool runAltitudeOptimization(const std::string &elev_file);
    // 内部辅助：从 JSON 读取 ENU 路径与 distance
    std::vector<ENUPoint> getENUFromJSON(const json& j, const std::string& key1, const std::string& key2);
    double getDistanceFromJSON(const json& j, const std::string& key);

    // Altitude Optimization Params
    struct AltitudeParams {
        double lambda_smooth = 1.0; // smoothing weight
        double lambda_follow = 0.0; // follow terrain weight (aim for z ~= elev + uav_R)
        double max_climb_rate = 2.0; // max vertical change per horizontal meter (m/m)
        double uav_R = 2.0; // UAV effective radius/height clearance (meters)
        double safe_distance = 50.0; // Safe distance from terrain (meters)
    };

    // Elevation Map Data
    int elev_width_ = 0;
    int elev_height_ = 0;
    std::vector<float> elev_data_; // row-major, top-left origin
    GeoTransform elev_geo_;
    bool elev_valid_ = false;

    // CostMap Data
    int cost_width_ = 0;
    int cost_height_ = 0;
    double cost_resolution_ = 0.0;
    double cost_origin_x_ = 0.0;
    double cost_origin_y_ = 0.0;
    std::vector<float> cost_data_; // row-major, top-left origin

    // Helpers for elevation
    bool performDownsampling(void* poDS_ptr, int full_w, int full_h, uint64_t bytes_needed, uint64_t target_bytes, const std::string& path);
    bool getElevationAt(double x, double y, double &elev) const;
    bool optimizeHeights(const std::vector<Eigen::Vector3d> &waypoints, const AltitudeParams &p, std::vector<double> &out_z);
    bool optimizeHeightsGlobalSmooth(const std::vector<double> &input_z, const std::vector<Eigen::Vector3d> &waypoints, const AltitudeParams &p, std::vector<double> &out_z);

    // ECEF转经纬度（迭代法）
    ECEFPoint wgs84ToECEF(const WGS84Point& lla);

    // Altitude Optimization methods
    bool loadElevationData(const std::string &path);
    
    // Elevation Map Info
    void printElevationInfo() const;
    bool isElevationValid() const { return elev_valid_; }

    // CostMap methods
    void createCostMap(int width, int height, double resolution, double origin_x, double origin_y);
    void setCost(int x, int y, float c);
    float getCost(int x, int y) const;
    bool getCostAt(double x, double y, float &val) const;
    // void initCostMapFromElevation(); // Removed in favor of local ENU map
    void buildLocalENUCostMap(double margin, double resolution);

    json generateVShapeTrajectories(const json &uavs_ids, const json &uav_starts, 
                                    const ENUPoint &leader_start_enu, const Eigen::Matrix2d &R0,
                                    const std::vector<Eigen::Vector2d> &leader_xy, 
                                    const std::vector<double> &leader_headings,
                                    const std::vector<ENUPoint> &Trajectory_ENU,
                                    double safety_distance);

    json generateLineShapeTrajectories(const json &uavs_ids, const json &uav_starts, 
                                       const ENUPoint &leader_start_enu, const Eigen::Matrix2d &R0,
                                       const std::vector<Eigen::Vector2d> &leader_xy, 
                                       const std::vector<double> &leader_headings,
                                       const std::vector<ENUPoint> &Trajectory_ENU,
                                       double safety_distance);
};
#endif
