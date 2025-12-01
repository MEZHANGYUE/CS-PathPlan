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
using namespace std;
using json = nlohmann::json;
using namespace math_util;

struct InputData
{
    double distance_points;
    double leader_speed = 5.0; // m/s, average speed override read from input JSON (formerly V_avg)
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
extern std::vector<ENUPoint>  Trajectory_ENU;

// Minisnap_3D: origin_waypoints, sampling distance (m), optional average speed override (m/s)
std::vector<ENUPoint>  Minisnap_3D(std::vector<ENUPoint> origin_waypoints, double distance_, double V_avg_override = -1.0);

bool getPlan(json &input_json, json &output_json);

bool loadData(InputData &input_data, json &input_json);
bool putWGS84ToJson(json &j, const std::string &key, const std::vector<WGS84Point> &traj);
// 生成跟随者（followers）的编队轨迹，返回用于写入 output_json 的 plane_array
json generateFollowerTrajectories(const json &input_json, const InputData &input_data,
                                  const std::vector<ENUPoint> &Trajectory_ENU,
                                  const std::vector<WGS84Point> &Trajectory_WGS84);
Eigen::VectorXd minimumSnapTimeAllocation(const Eigen::MatrixXd& route, double V_avg, double min_turn_radius = 200);
Eigen::VectorXd minimumSnapTimeAllocation(const Eigen::MatrixXd& route);
// 主转换函数：WGS84经纬度转东北天坐标
ENUPoint wgs84ToENU(const WGS84Point& target, const WGS84Point& reference);
// 东北天坐标转WGS84经纬度
WGS84Point enuToWGS84(const ENUPoint& enu, const WGS84Point& reference);
// 批量转换
std::vector<ENUPoint> wgs84ToENU_Batch(const std::vector<WGS84Point>& targets, 
                                      const WGS84Point& reference);
std::vector<WGS84Point> enuToWGS84_Batch(const std::vector<ENUPoint>& targets, 
                                      const WGS84Point& reference);
inline bool saveJsonToFile(const json &j, const std::string &filename) {
    try {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open file " << filename << std::endl;
            return false;
        }
        
        // 修正：使用 dump() 方法将 json 对象转换为字符串
        file << j.dump(4);  // 缩进4个空格，美化输出
        file.close();
        
        std::cout << "Successfully saved JSON to " << filename << std::endl;
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error saving JSON to file: " << e.what() << std::endl;
        return false;
    }
}                                     
#endif
