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
#include "math_util/polygon2d.hpp"
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cstdint>
#include <memory>
#include <array>

using namespace std;
using json = nlohmann::json;
using namespace math_util;

class ElevationCostMap;

struct ProhibitedZone {
    std::vector<WGS84Coord> polygon;
    std::pair<double, double> height_range;
};

// 飞行区域通用结构体（用于表示准备区域、战斗区域等）
struct FlightZone {
    int zone_id = 0;                            // 区域ID
    std::string zone_type;                      // 区域类型（如 "ready_zone", "battle_zone"）
    std::vector<WGS84Coord> polygon;            // 区域多边形边界定点集合
    std::pair<double, double> height_range;     // 区域的高度范围 [下限, 上限]
    int link_flag = 0;                          // 联结标志或其他自定义数据
};

// 轨迹行（可用于输入/输出）：
// - input_json["using_midway_lines"] 的一行
// - output_json["using_midway_lines"] 的一行
// 形状为：[uav_id, segment_id, [lon,lat,alt], [lon,lat,alt], ...]
struct UavTrajectoryLine {
    int uav_id = 0;
    int segment_id = 0;
    std::vector<WGS84Coord> points;
};

struct InputData
{
    double distance_points = 0.0;
    double leader_speed = 30.0; // m/s, average speed override read from input JSON (formerly V_avg)
    double min_turning_radius = 0.0; // Minimum turning radius constraint
    double leader_fly_high;
    int formation_model;
    int formation_using;
    int uav_leader_id;
    std::pair<double, double> height_list;
    FlightZone ready_zone;
    std::vector<WGS84Coord> high_zhandou_point_wgs84;
    std::vector<WGS84Coord> leader_midway_point_wgs84;
    std::vector<WGS84Coord> uav_start_point_wgs84;
    std::vector<int> uavs_id;
    std::vector<int> ready_id;
    std::vector<int> uav_leader_ids;
    std::vector<std::array<int, 3>> uavs_plane_data_list; // [uav_id, segment_id, point_idx]
    // 输入历史轨迹（若 input.json 提供 using_midway_lines）
    std::vector<UavTrajectoryLine> using_midway_lines;
    // input.json 可能提供的“使用中的无人机列表”（为空通常表示使用全部）
    std::vector<int> using_uav_list;

    // battle_zone 字段（包含多边形和高度）
    std::vector<FlightZone> battle_zones;
    // battle_zone_list：输入可能提供的辅助字段（当前规划逻辑不依赖它），按 int 列表透传解析
    std::vector<int> battle_zone_list;
    WGS84Coord uav_leader_start_point_wgs84 = {0.0, 0.0, 0.0};
    std::pair<int, std::pair<int, int>> uavs_plane_data;
    bool has_prohibited_zone = false;
    std::vector<ProhibitedZone> prohibited_zones;
    bool has_check_prohibited_zone = false;
    std::vector<ProhibitedZone> check_prohibited_zones;
    std::vector<WGS84Coord> existing_midway_lines; // 新增

    // 可由 input_json 覆盖 config 的参数
    double formation_distance = -1.0;
    double position_misalignment = -1.0;
    double uav_R = -1.0;
    int uav_formation_max_row = 0;

    // 高度优化参数覆盖
    double ao_uav_R = -1.0;
    double ao_safe_distance = -1.0;
    double ao_lambda_follow = -1.0;
    double ao_lambda_smooth = -1.0;
    double ao_max_climb_rate = -1.0;
};

struct OutputData
{
    // 这些字段均直接对应 uavPathPlanning.cpp 中构建的 output_json 内容。
    // 注意：OutputData 本身不包含 json 类型；如需序列化/反序列化，请在 cpp 中做转换。

    // output_json["abnormal_uav_plane"]: [int, ...]
    std::vector<int> abnormal_uav_plane;

    // output_json["using_uav_list"]: [int, ...]
    std::vector<int> using_uav_list;

    // output_json["ready_id"]: [int, ...]
    std::vector<int> ready_id;

    // output_json["midway_point_num"]: [int, ...]
    std::vector<int> midway_point_num;

    // output_json["leader_show_points"]: [[lon,lat,alt], ...]
    std::vector<WGS84Coord> leader_show_points;

    // output_json["uav_leader_plane{1,2,3}"]: [[lon,lat,alt], ...]
    std::vector<WGS84Coord> uav_leader_plane1;
    std::vector<WGS84Coord> uav_leader_plane2;
    std::vector<WGS84Coord> uav_leader_plane3;

    // output_json["uav_plane{1,2,3}"]: [[uav_id, [lon,lat,alt], ...], ...]
    // 这里用 UavTrajectoryLine 统一表达（segment_id 分别固定为 1/2/3）。
    std::vector<UavTrajectoryLine> uav_plane1;
    std::vector<UavTrajectoryLine> uav_plane2;
    std::vector<UavTrajectoryLine> uav_plane3;

    // output_json["using_midway_lines"]: [[uav_id, segment_id, [lon,lat,alt]...], ...]
    std::vector<UavTrajectoryLine> using_midway_lines;
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
    struct PlannerConfig {
        struct AltitudeOptimization {
            bool enabled = false;
            std::string elevation_file;
            double lambda_smooth = 1.0;
            double lambda_follow = 0.0;
            double max_climb_rate = 2.0;
            double uav_R = 2.0;
            double safe_distance = 50.0;
        } altitude_optimization;

        struct PathPlanning {
            double position_misalignment = 0.0;
            double min_turning_radius = 0.0;
            double patrol_width = 0.0;
            std::string patrol_mode = "BOW";
            double patrol_region_shrink_distance = 0.0; // 巡逻区域内缩距离（米），0 表示不内缩
            double formation_distance = 50.0;
            int uav_formation_max_row = 8;
            double distance_points = 0.0; // config.yaml -> path_planning.Distance_Points (or distance_points)
            double prohibited_zone_conflict_distance = 50.0; // 禁飞区冲突判定距离阈值（米）
        } path_planning;

        // minimum_snap 轨迹生成参数（config.yaml.minimum_snap）
        MinimumSnapConfig minimum_snap;

        bool loaded = false;
        std::string loaded_from;
        std::string load_error;
    };

    UavPathPlanner();
    ~UavPathPlanner();

    // 统一加载 config.yaml（只读一次，其他地方直接用 config() 取参数）
    bool loadFromYAML(const std::string &config_path = "");
    const PlannerConfig &config() const { return config_; }

    std::vector<ENUPoint> Trajectory_ENU={}; //单独优化高度时使用,非空时才能调用高度优化
    // Minisnap 轨迹生成接口
    std::vector<ENUPoint> Minisnap_3D(std::vector<ENUPoint> origin_waypoints, double distance_, double V_avg_override = -1.0);
    std::vector<ENUPoint> Minisnap_EN(std::vector<ENUPoint> origin_waypoints, double distance_, double V_avg_override = -1.0);
    // Bezier 轨迹生成接口
    std::vector<ENUPoint> Bezier_3D(std::vector<ENUPoint> origin_waypoints, double distance_, double V_avg_override = -1.0, double min_radius = 0.0);

    // 主规划接口
    bool getPlan(json &input_json, json &output_json, bool use3D = true, std::string algorithm = "minimum_snap");
    
    // 准备规划航点（整合中途点与区域边界点）
    std::vector<ENUPoint> preparePlanningWaypoints(int &midwaypoint_num, int &zhandoupoint_num);

    // 检查历史/当前航线与 check_prohibited_zone_wgs84 的冲突并输出 abnormal_uav_plane
    bool check_change(const json &input_json, json &output_json);
    //高度优化接口 以 Trajectory_ENU 为主，优化后还要联动 follower
    bool runAltitudeOptimization(const std::string &elev_file, OutputData &output_data);
    // 辅助函数
    bool loadData(InputData &input_data, json &input_json);
    // 将强类型 OutputData 写入 output_json（字段格式与当前约定保持一致）
    bool outputDataToJson(const OutputData &output_data, json &output_json);
    bool putWGS84ToJson(json &j, const std::string &key, const std::vector<WGS84Point> &traj);
    json generateFollowerTrajectories(const InputData &input_data,
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
                                          const std::vector<ENUPoint>& Patrol_Path, OutputData &output_data);

    // 将轨迹信息追加到 output_json 的 using_midway_lines 字段
    void appendTrajectoryToOutput(json &output_json, int uav_id, int segment_id, const std::vector<WGS84Point> &traj);

    // 避障处理
    std::vector<ENUPoint> avoidProhibitedZones(const std::vector<ENUPoint>& path);

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
    PlannerConfig config_;
    InputData input_data_; // 存储解析后的输入数据，避免重复解析
    OutputData output_data_; // 存储解析后的输出数据，避免重复解析
    // 提取的高度优化器函数仍在 cpp 中实现为私有方法
    // 现在只输入高程文件路径（例如 .tif），函数直接使用类成员 Trajectory_ENU 进行高度优化
    // bool runAltitudeOptimization(const std::string &elev_file);
    // 第三段巡逻轨迹：按 patrol_mode 选择计算方法（巡逻区域点集直接传入）
    std::vector<ENUPoint> computePatrolPathByMode(const std::vector<ENUPoint>& patrol_zone,
                                                  double distance,
                                                  const std::string& patrol_mode,
                                                  const std::vector<ENUPoint>& trajectory_enu);

    // 兼容封装：从 enu_waypoints 尾部截取 patrolpoint_num 个点作为巡逻区域，然后调用 computePatrolPathByMode。
    std::vector<ENUPoint> computePatrolPathFromWaypoints(const std::vector<ENUPoint>& enu_waypoints,
                                                         int patrolpoint_num,
                                                         double distance,
                                                         const std::string& patrol_mode,
                                                         const std::vector<ENUPoint>& trajectory_enu);

    // SINGLE 巡逻路径生成（从巡逻区域多边形点集生成闭合巡逻轨迹）
    std::vector<ENUPoint> gen_single_patrol(const std::vector<ENUPoint> &patrol_zone,
                                            double distance,
                                            const std::vector<ENUPoint> &trajectory_enu);

    // BOW 巡逻路径生成（暂未实现，先留空结构）
    std::vector<ENUPoint> gen_bow_patrol(const std::vector<ENUPoint> &patrol_zone,
                                         double distance,
                                         const std::vector<ENUPoint> &trajectory_enu);

    // CIRCULAR 巡逻路径生成（暂未实现，先留空结构）
    std::vector<ENUPoint> gen_circular_patrol(const std::vector<ENUPoint> &patrol_zone,
                                              double distance,
                                              const std::vector<ENUPoint> &trajectory_enu);

    // 巡逻区域内缩（仅对 east/north 做 2D 内缩；up 保持不变）
    void shrinkPolygon(std::vector<ENUPoint> &polygon, double shrink_meters) const;

    // 编队结束后：为指定无人机选择一个 battle_zone（通常按 uavs_id 下标映射到 battle_zone_wgs84）。
    const FlightZone *selectBattleZoneForUav(int uav_id) const;

    // 检查 battle_zone 在给定高度层(target_up, ENU.up)是否可用：
    // - 高度参数是否有效
    // - 在该高度层是否与禁飞区(prohibited_zones)水平重叠（且高度范围覆盖该层）
    // 函数内部会打印检查结果。
    bool check_battle_zone(int uav_id, const FlightZone &battle_zone, double target_up) const;

    // Altitude Optimization Params
    struct AltitudeParams {
        double lambda_smooth = 1.0; // smoothing weight
        double lambda_follow = 0.0; // follow terrain weight (aim for z ~= elev + uav_R)
        double max_climb_rate = 2.0; // max vertical change per horizontal meter (m/m)
        double uav_R = 2.0; // UAV effective radius/height clearance (meters)
        double safe_distance = 50.0; // Safe distance from terrain (meters)
    };

    // Elevation + Cost map (delegated to ElevationCostMap)
    std::unique_ptr<ElevationCostMap> elev_cost_map_;

    // Helpers for altitude optimization
    AltitudeParams makeAltitudeParams() const;
    bool optimizeSegmentAltitudeENU(std::vector<ENUPoint> &segment_enu);
    // 以 OutputData 为主，优化后直接回写 OutputData（并同步 using_midway_lines）
    bool optimizeAndApplyOutputSegment(OutputData &output_data, const char *key, int segment_id, bool keep_closed_equal_height);
    // 联合优化多段路径高度，并支持特定段（如巡逻段）等高约束
    bool optimizeAndApplyJointSegments(OutputData &output_data, const std::vector<std::string> &keys, const std::vector<int> &segment_ids, int equal_height_segment_idx = -1);
    bool optimizeHeights(const std::vector<Eigen::Vector3d> &waypoints, const AltitudeParams &p, std::vector<double> &out_z);
    bool optimizeHeightsGlobalSmooth(const std::vector<double> &input_z, const std::vector<Eigen::Vector3d> &waypoints, const AltitudeParams &p, std::vector<double> &out_z);

    // ECEF转经纬度（迭代法）
    ECEFPoint wgs84ToECEF(const WGS84Point& lla);

    // Elevation/CostMap methods are provided directly by ElevationCostMap
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

    // formation_model = 3: 一字形编队(竖) / trail
    // uav_formation_max_row: 每一列最多支持的“总行数”（主列包含长机，因此主列最多容纳 max_row-1 架僚机；其余列每列最多 max_row 架僚机）。最小为 1。
    json generateVerticalLineShapeTrajectories(const json &uavs_ids, const json &uav_starts,
                                               const ENUPoint &leader_start_enu, const Eigen::Matrix2d &R0,
                                               const std::vector<Eigen::Vector2d> &leader_xy,
                                               const std::vector<double> &leader_headings,
                                               const std::vector<ENUPoint> &Trajectory_ENU,
                                               int uav_formation_max_row,
                                               double safety_distance);

    // formation_model = 4: 三角形编队 / triangle
    json generateTriangleShapeTrajectories(const json &uavs_ids, const json &uav_starts,
                                           const ENUPoint &leader_start_enu, const Eigen::Matrix2d &R0,
                                           const std::vector<Eigen::Vector2d> &leader_xy,
                                           const std::vector<double> &leader_headings,
                                           const std::vector<ENUPoint> &Trajectory_ENU,
                                           double safety_distance);
};
#endif
