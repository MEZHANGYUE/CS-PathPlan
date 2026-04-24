#include "uavPathPlanning.hpp"
#include "elevation_cost_map.hpp"
#include "math_util/minimum_snap.hpp"
#include "math_util/bezier.hpp"
#include <yaml-cpp/yaml.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <sys/stat.h>
#include <iomanip> // Added for std::setprecision
#include <queue>
#include <map>
#include <set>
#include <unordered_map>
#include <cctype>
#include "math_util/polygon2d.hpp"
#include "algorithms/clipper.hpp"

#ifdef HAVE_GDAL
#include <gdal_priv.h>
#include <cpl_conv.h>
#endif

// (moved into UavPathPlanner::generator_)

namespace {
inline bool fileExists(const std::string &path) {
    struct stat buffer;
    return ::stat(path.c_str(), &buffer) == 0;
}

template <typename T>
inline void yamlAssignIfPresent(const YAML::Node &node, const char *key, T &out) {
    try {
        if (node && node[key]) out = node[key].as<T>();
    } catch (...) {
        // ignore conversion errors
    }
}

inline void yamlAssignIfPresent(const YAML::Node &node, const char *key, std::string &out) {
    try {
        if (node && node[key]) out = node[key].as<std::string>();
    } catch (...) {
        // ignore conversion errors
    }
}

inline void yamlAssignVec3IfPresent(const YAML::Node &node, const char *key, Eigen::Vector3d &out) {
    try {
        if (!node || !node[key] || !node[key].IsSequence() || node[key].size() < 3) return;
        out = Eigen::Vector3d(node[key][0].as<double>(), node[key][1].as<double>(), node[key][2].as<double>());
    } catch (...) {
        // ignore conversion errors
    }
}

struct CheckZoneSpec {
    std::vector<WGS84Point> polygon_wgs;
    double min_h = -std::numeric_limits<double>::infinity();
    double max_h = std::numeric_limits<double>::infinity();
};

inline bool altitudeRangeOverlap(double a_min, double a_max, double b_min, double b_max) {
    return !(a_max < b_min || b_max < a_min);
}

inline bool rangeCovers(double v, double lo, double hi) {
    if (std::isnan(v)) return false;
    if (std::isnan(lo) || std::isnan(hi)) return false;
    if (!std::isfinite(lo) && !std::isfinite(hi)) return true;
    if (!std::isfinite(lo) && std::isfinite(hi)) return v <= hi;
    if (std::isfinite(lo) && !std::isfinite(hi)) return v >= lo;
    if (lo <= hi) return (v >= lo && v <= hi);
    return (v >= hi && v <= lo);
}

inline Polygon2d makePolygon2dFromWgs84(const std::vector<WGS84Coord> &poly_wgs,
                                       const WGS84Point &origin,
                                       UavPathPlanner *planner) {
    std::vector<Vec2d> pts;
    pts.reserve(poly_wgs.size());
    for (const auto &p : poly_wgs) {
        WGS84Point w{p.lon, p.lat, p.alt};
        ENUPoint e = planner->wgs84ToENU(w, origin);
        pts.emplace_back(e.east, e.north);
    }
    return Polygon2d(pts);
}

inline bool polygonsOverlap2D(const Polygon2d &a, const Polygon2d &b) {
    const auto &ap = a.points();
    const auto &bp = b.points();
    if (ap.size() < 3 || bp.size() < 3) return false;

    // 1) any vertex-in-polygon
    for (const auto &p : ap) {
        if (b.IsPointIn(p)) return true;
    }
    for (const auto &p : bp) {
        if (a.IsPointIn(p)) return true;
    }

    // 2) any edge intersection
    std::vector<Vec2d> intersections;
    intersections.reserve(4);
    for (size_t i = 0; i < ap.size(); ++i) {
        const Vec2d &p1 = ap[i];
        const Vec2d &p2 = ap[(i + 1) % ap.size()];
        LineSegment2d seg(p1, p2);
        if (b.Intersections(seg, &intersections) && !intersections.empty()) return true;
    }
    return false;
}

inline bool parseWGS84PointArray(const json &arr, WGS84Point &out) {
    if (!arr.is_array() || arr.size() < 2 || !arr[0].is_number() || !arr[1].is_number()) {
        return false;
    }
    out.lon = arr[0].get<double>();
    out.lat = arr[1].get<double>();
    out.alt = (arr.size() >= 3 && arr[2].is_number()) ? arr[2].get<double>() : 0.0;
    return true;
}

inline bool parseWGS84PointValue(const json &value, WGS84Point &out) {
    if (parseWGS84PointArray(value, out)) {
        return true;
    }

    if (!value.is_object()) {
        return false;
    }

    const char *lon_keys[] = {"lon", "lng", "x", "longitude"};
    const char *lat_keys[] = {"lat", "y", "latitude"};
    const char *alt_keys[] = {"alt", "z", "height", "altitude"};

    bool has_lon = false;
    bool has_lat = false;
    out.alt = 0.0;

    for (const char *key : lon_keys) {
        if (value.contains(key) && value[key].is_number()) {
            out.lon = value[key].get<double>();
            has_lon = true;
            break;
        }
    }
    for (const char *key : lat_keys) {
        if (value.contains(key) && value[key].is_number()) {
            out.lat = value[key].get<double>();
            has_lat = true;
            break;
        }
    }
    for (const char *key : alt_keys) {
        if (value.contains(key) && value[key].is_number()) {
            out.alt = value[key].get<double>();
            break;
        }
    }

    return has_lon && has_lat;
}

inline double degToRad(double deg) {
    return deg * M_PI / 180.0;
}

inline double wgs84DistanceSquaredMeters(const WGS84Point &a, const WGS84Point &b) {
    constexpr double kEarthRadius = 6378137.0;
    const double lat1 = degToRad(a.lat);
    const double lat2 = degToRad(b.lat);
    const double dlat = lat2 - lat1;
    const double dlon = degToRad(b.lon - a.lon);
    const double x = dlon * std::cos((lat1 + lat2) * 0.5) * kEarthRadius;
    const double y = dlat * kEarthRadius;
    const double z = b.alt - a.alt;
    return x * x + y * y + z * z;
}

inline void appendTrajectoryPointsFromJson(const json &segment, std::vector<WGS84Point> &trajectory) {
    if (!segment.is_array()) {
        return;
    }

    for (const auto &pt_json : segment) {
        WGS84Point pt;
        if (!parseWGS84PointValue(pt_json, pt)) {
            continue;
        }
        trajectory.push_back(pt);
    }
}

// 把输入的 leader_midway_point_wgs84（长机中途点）映射到最终生成的航迹点序列上，记录“每个中途点在航迹里的下标”。
// 结构：json array[int]；找不到就填 -1。
inline json buildMidwayPointNum(const InputData &input_data, const json &output_json) {
    json midway_point_num = json::array();
    if (input_data.leader_midway_point_wgs84.empty()) {
        return midway_point_num;
    }

    std::vector<WGS84Point> leader_trajectory;
    appendTrajectoryPointsFromJson(output_json.value("uav_leader_plane1", json::array()), leader_trajectory);
    appendTrajectoryPointsFromJson(output_json.value("uav_leader_plane2", json::array()), leader_trajectory);
    appendTrajectoryPointsFromJson(output_json.value("uav_leader_plane3", json::array()), leader_trajectory);

    for (const auto &midway : input_data.leader_midway_point_wgs84) {
        WGS84Point midway_point{midway.lon, midway.lat, midway.alt};
        if (leader_trajectory.empty()) {
            midway_point_num.push_back(-1);
            continue;
        }

        int best_index = -1;
        double best_dist_sq = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < leader_trajectory.size(); ++i) {
            const double dist_sq = wgs84DistanceSquaredMeters(midway_point, leader_trajectory[i]);
            if (dist_sq < best_dist_sq) {
                best_dist_sq = dist_sq;
                best_index = static_cast<int>(i);
            }
        }
        midway_point_num.push_back(best_index);
    }

    return midway_point_num;
}

inline bool parseHeightRangeArray(const json &arr, double &min_h, double &max_h) {
    if (!arr.is_array() || arr.size() < 2 || !arr[0].is_number() || !arr[1].is_number()) {
        return false;
    }
    min_h = arr[0].get<double>();
    max_h = arr[1].get<double>();
    if (min_h > max_h) std::swap(min_h, max_h);
    return true;
}

inline std::vector<CheckZoneSpec> parseCheckZonesFromInputData(const InputData &input_data) {
    std::vector<CheckZoneSpec> zones;
    zones.reserve(input_data.check_prohibited_zones.size());

    for (const auto &zone : input_data.check_prohibited_zones) {
        CheckZoneSpec z;
        z.min_h = zone.height_range.first;
        z.max_h = zone.height_range.second;
        for (const auto &p : zone.polygon) {
            z.polygon_wgs.push_back({p.lon, p.lat, p.alt});
        }
        if (z.polygon_wgs.size() >= 3) {
            zones.push_back(std::move(z));
        }
    }

    return zones;
}

struct UavProgress {
    int segment_id = 0;
    int point_idx = 0; // 1-based
};

inline std::unordered_map<int, UavProgress> parseUavProgressFromInputData(const InputData &input_data) {
    std::unordered_map<int, UavProgress> progress;
    for (const auto &item : input_data.uavs_plane_data_list) {
        int uav_id = item[0];
        UavProgress cur{item[1], item[2]};

        auto it = progress.find(uav_id);
        if (it == progress.end()) {
            progress[uav_id] = cur;
            continue;
        }

        const UavProgress &old = it->second;
        if (cur.segment_id > old.segment_id ||
            (cur.segment_id == old.segment_id && cur.point_idx > old.point_idx)) {
            it->second = cur;
        }
    }
    return progress;
}

inline double computeTailHeadingRobust(const std::vector<ENUPoint> &path, double fallback_heading = 0.0) {
    if (path.size() < 2) return fallback_heading;

    // 使用末端若干段的加权方向向量，抑制末尾单段抖动/反向小回摆
    const int max_segments = 8;
    Eigen::Vector2d acc(0.0, 0.0);
    int used = 0;

    for (int i = static_cast<int>(path.size()) - 1; i > 0 && used < max_segments; --i) {
        double dx = path[i].east - path[i - 1].east;
        double dy = path[i].north - path[i - 1].north;
        double d = std::hypot(dx, dy);
        if (d < 1e-3) continue;  // 跳过几乎重合点

        // 越靠近末端权重越大
        double w = 1.0 + 0.25 * used;
        acc.x() += w * (dx / d);
        acc.y() += w * (dy / d);
        ++used;
    }

    if (used == 0 || acc.norm() < 1e-9) {
        // 回退到最后一个非零线段
        for (int i = static_cast<int>(path.size()) - 1; i > 0; --i) {
            double dx = path[i].east - path[i - 1].east;
            double dy = path[i].north - path[i - 1].north;
            if (std::hypot(dx, dy) > 1e-3) {
                return std::atan2(dy, dx);
            }
        }
        return fallback_heading;
    }

    return std::atan2(acc.y(), acc.x());
}
} // namespace

const FlightZone *UavPathPlanner::selectBattleZoneForUav(int uav_id) const {
    if (this->input_data_.battle_zones.empty()) return nullptr;

    // Prefer index mapping: uavs_id[i] -> battle_zones[i]
    for (size_t i = 0; i < this->input_data_.uavs_id.size(); ++i) {
        if (this->input_data_.uavs_id[i] == uav_id) {
            if (i < this->input_data_.battle_zones.size()) return &this->input_data_.battle_zones[i];
            break;
        }
    }
    // Fallback: first battle zone
    return &this->input_data_.battle_zones.front();
}

bool UavPathPlanner::check_battle_zone(int uav_id, const FlightZone &battle_zone, double target_up) const {
    bool height_ok = std::isfinite(target_up);
    bool poly_ok = (battle_zone.polygon.size() >= 3);
    bool overlap = false;

    if (!poly_ok) {
        std::cout << "[BattleZoneCheck] uav_id=" << uav_id
                  << " zone_id=" << battle_zone.zone_id
                  << " FAIL: battle_zone polygon invalid (<3 points)." << std::endl;
        return false;
    }

    // Build battle polygon in ENU for 2D overlap checks
    // Note: battle_zone polygon points may have no altitude; we only care about east/north.
    std::vector<WGS84Coord> battle_poly_wgs = battle_zone.polygon;
    Polygon2d battle_poly_enu = makePolygon2dFromWgs84(battle_poly_wgs, this->origin_, const_cast<UavPathPlanner *>(this));

    if (!this->input_data_.prohibited_zones.empty()) {
        for (size_t zi = 0; zi < this->input_data_.prohibited_zones.size(); ++zi) {
            const auto &pz = this->input_data_.prohibited_zones[zi];
            if (pz.polygon.size() < 3) continue;

            const double zmin = pz.height_range.first;
            const double zmax = pz.height_range.second;
            if (height_ok && !rangeCovers(target_up, zmin, zmax)) {
                continue; // this altitude layer does not intersect the prohibited zone
            }

            std::vector<WGS84Coord> pz_poly_wgs = pz.polygon;
            Polygon2d pz_poly_enu = makePolygon2dFromWgs84(pz_poly_wgs, this->origin_, const_cast<UavPathPlanner *>(this));
            if (polygonsOverlap2D(battle_poly_enu, pz_poly_enu)) {
                overlap = true;
                std::cout << "[BattleZoneCheck] uav_id=" << uav_id
                          << " zone_id=" << battle_zone.zone_id
                          << " overlap prohibited_zone idx=" << zi
                          << " at target_up=" << target_up
                          << " (pz.height_range=[" << zmin << "," << zmax << "])" << std::endl;
                break;
            }
        }
    }

    const bool ok = height_ok && !overlap;
    std::cout << "[BattleZoneCheck] uav_id=" << uav_id
              << " zone_id=" << battle_zone.zone_id
              << " target_up=" << target_up
              << " height_ok=" << (height_ok ? 1 : 0)
              << " overlap_prohibited=" << (overlap ? 1 : 0)
              << " => " << (ok ? "OK" : "FAIL") << std::endl;
    return ok;
}

// 统一加载 config.yaml（只读一次，其他地方直接用 this->config_ 的参数）
bool UavPathPlanner::loadFromYAML(const std::string &config_path)
{
    PlannerConfig cfg; // defaults
    std::vector<std::string> candidates;

    if (!config_path.empty()) {
        candidates.push_back(config_path);
    } else {
        candidates = {"config.yaml", "../config.yaml", "../../config.yaml"};
    }

    std::string found;
    for (const auto &p : candidates) {
        if (fileExists(p)) {
            found = p;
            break;
        }
    }

    if (found.empty()) {
        cfg.loaded = false;
        cfg.loaded_from.clear();
        cfg.load_error = "config.yaml not found";
        this->config_ = cfg;
        return false;
    }

    try {
        YAML::Node root = YAML::LoadFile(found);
        cfg.loaded = true;
        cfg.loaded_from = found;
        cfg.load_error.clear();

        if (root["altitude_optimization"]) {
            auto alt = root["altitude_optimization"];
            yamlAssignIfPresent(alt, "enabled", cfg.altitude_optimization.enabled);
            yamlAssignIfPresent(alt, "elevation_file", cfg.altitude_optimization.elevation_file);
            yamlAssignIfPresent(alt, "lambda_smooth", cfg.altitude_optimization.lambda_smooth);
            yamlAssignIfPresent(alt, "lambda_follow", cfg.altitude_optimization.lambda_follow);
            yamlAssignIfPresent(alt, "max_climb_rate", cfg.altitude_optimization.max_climb_rate);
            yamlAssignIfPresent(alt, "uav_R", cfg.altitude_optimization.uav_R);
            yamlAssignIfPresent(alt, "safe_distance", cfg.altitude_optimization.safe_distance);
        }

        if (root["path_planning"]) {
            auto pp = root["path_planning"];
            yamlAssignIfPresent(pp, "minimum_snap_config_file", cfg.path_planning.minimum_snap_config_file);
            yamlAssignIfPresent(pp, "patrol_region_shrink_distance", cfg.path_planning.patrol_region_shrink_distance);
            yamlAssignIfPresent(pp, "position_misalignment", cfg.path_planning.position_misalignment);
            yamlAssignIfPresent(pp, "min_turning_radius", cfg.path_planning.min_turning_radius);
            yamlAssignIfPresent(pp, "patrol_width", cfg.path_planning.patrol_width);
            yamlAssignIfPresent(pp, "patrol_mode", cfg.path_planning.patrol_mode);
            yamlAssignIfPresent(pp, "formation_distance", cfg.path_planning.formation_distance);
            yamlAssignIfPresent(pp, "uav_formation_max_row", cfg.path_planning.uav_formation_max_row);
            yamlAssignIfPresent(pp, "prohibited_zone_conflict_distance", cfg.path_planning.prohibited_zone_conflict_distance);

            // 兼容历史字段
            if (pp["Distance_Points"]) {
                yamlAssignIfPresent(pp, "Distance_Points", cfg.path_planning.distance_points);
            } else {
                yamlAssignIfPresent(pp, "distance_points", cfg.path_planning.distance_points);
            }
        }

        // minimum_snap 参数：仅从 path_planning.minimum_snap_config_file 指向的独立配置文件读取
        if (!cfg.path_planning.minimum_snap_config_file.empty()) {
            std::string resolved = cfg.path_planning.minimum_snap_config_file;
            try {
                if (!fileExists(resolved)) {
                    std::cerr << "Warning: minimum_snap_config_file not found: " << resolved << std::endl;
                } else {
                    YAML::Node ms_root = YAML::LoadFile(resolved);
                    YAML::Node ms = ms_root;
                    // allow wrapper style: minimum_snap: { ... }
                    if (ms_root["minimum_snap"]) ms = ms_root["minimum_snap"];

                    yamlAssignIfPresent(ms, "order", cfg.minimum_snap.order);
                    yamlAssignIfPresent(ms, "path_weight", cfg.minimum_snap.path_weight);
                    yamlAssignIfPresent(ms, "vel_zero_weight", cfg.minimum_snap.vel_zero_weight);
                    yamlAssignIfPresent(ms, "V_avg", cfg.minimum_snap.V_avg);
                    yamlAssignIfPresent(ms, "min_time_s", cfg.minimum_snap.min_time_s);
                    yamlAssignIfPresent(ms, "sample_distance", cfg.minimum_snap.sample_distance);
                    yamlAssignVec3IfPresent(ms, "start_vel", cfg.minimum_snap.start_vel);
                    yamlAssignVec3IfPresent(ms, "end_vel", cfg.minimum_snap.end_vel);
                    yamlAssignVec3IfPresent(ms, "start_acc", cfg.minimum_snap.start_acc);
                    yamlAssignVec3IfPresent(ms, "end_acc", cfg.minimum_snap.end_acc);

                    std::cerr << "Loaded minimum_snap config from " << resolved << std::endl;
                }
            } catch (const std::exception &e) {
                std::cerr << "Warning: failed to load minimum_snap config from " << resolved << ": " << e.what() << std::endl;
            }
        }
        this->config_ = cfg;
        std::cerr << "Loaded config from " << cfg.loaded_from << std::endl;
        return true;
    } catch (const std::exception &e) {
        cfg.loaded = false;
        cfg.loaded_from = found;
        cfg.load_error = e.what();
        this->config_ = cfg;
        std::cerr << "Warning: failed to load config from " << found << ": " << e.what() << std::endl;
        return false;
    }
}

// 经纬度转ECEF坐标
ECEFPoint UavPathPlanner::wgs84ToECEF(const WGS84Point& lla) {
    double lat_rad = deg2rad(lla.lat);
    double lon_rad = deg2rad(lla.lon);
    
    double N = calcN(lat_rad);
    double cos_lat = cos(lat_rad);
    double sin_lat = sin(lat_rad);
    double cos_lon = cos(lon_rad);
    double sin_lon = sin(lon_rad);
    
    ECEFPoint ecef;
    ecef.x = (N + lla.alt) * cos_lat * cos_lon;
    ecef.y = (N + lla.alt) * cos_lat * sin_lon;
    ecef.z = (N * (1 - WGS84_E2) + lla.alt) * sin_lat;
    
    return ecef;
}

// UavPathPlanner 构造/析构
UavPathPlanner::UavPathPlanner() {
    origin_.lon = 0.0;
    origin_.lat = 0.0;
    origin_.alt = 0.0;
    elev_cost_map_ = std::make_unique<ElevationCostMap>();

    // 统一读取 config.yaml（若不存在则使用默认值）
    this->loadFromYAML();
}

UavPathPlanner::~UavPathPlanner() = default;

// ECEF转经纬度（迭代法）
WGS84Point UavPathPlanner::ecefToWGS84(const ECEFPoint& ecef) {
    // 初始估计
    double p = sqrt(ecef.x * ecef.x + ecef.y * ecef.y);
    double theta = atan2(ecef.z * WGS84_A, p * WGS84_A * (1 - WGS84_E2));
    
    double lat_rad = atan2(ecef.z + WGS84_E2 * WGS84_A * (1 - WGS84_E2) * pow(sin(theta), 3) / (1 - WGS84_E2),
                          p - WGS84_E2 * WGS84_A * pow(cos(theta), 3));
    
    // 迭代计算以提高精度
    const int max_iterations = 10;
    const double tolerance = 1e-12;
    
    for (int i = 0; i < max_iterations; ++i) {
        double N = calcN(lat_rad);
        double alt = p / cos(lat_rad) - N;
        
        double lat_new = atan2(ecef.z, p * (1 - WGS84_E2 * N / (N + alt)));
        
        if (fabs(lat_new - lat_rad) < tolerance) {
            lat_rad = lat_new;
            break;
        }
        lat_rad = lat_new;
    }
    
    double lon_rad = atan2(ecef.y, ecef.x);
    
    // 计算最终高度
    double N = calcN(lat_rad);
    double alt;
    if (p < 1e-12) { // 接近极点
        alt = fabs(ecef.z) - WGS84_A * sqrt(1 - WGS84_E2);
    } else {
        alt = p / cos(lat_rad) - N;
    }
    
    WGS84Point lla;
    lla.lat = rad2deg(lat_rad);
    lla.lon = rad2deg(lon_rad);
    lla.alt = alt;
    
    return lla;
}

// 计算东北天坐标系的旋转矩阵
std::array<std::array<double, 3>, 3> UavPathPlanner::computeENURotationMatrix(double lat_rad, double lon_rad) {
    std::array<std::array<double, 3>, 3> R;
    
    double cos_lat = cos(lat_rad);
    double sin_lat = sin(lat_rad);
    double cos_lon = cos(lon_rad);
    double sin_lon = sin(lon_rad);
    
    // East轴单位向量 (东方向)
    R[0][0] = -sin_lon;
    R[0][1] = cos_lon;
    R[0][2] = 0.0;
    
    // North轴单位向量 (北方向)
    R[1][0] = -sin_lat * cos_lon;
    R[1][1] = -sin_lat * sin_lon;
    R[1][2] = cos_lat;
    
    // Up轴单位向量 (天方向)
    R[2][0] = cos_lat * cos_lon;
    R[2][1] = cos_lat * sin_lon;
    R[2][2] = sin_lat;
    
    return R;
}

// 计算东北天坐标系的逆旋转矩阵（转置）
std::array<std::array<double, 3>, 3> UavPathPlanner::computeENURotationMatrixInverse(double lat_rad, double lon_rad) {
    std::array<std::array<double, 3>, 3> R;
    
    double cos_lat = cos(lat_rad);
    double sin_lat = sin(lat_rad);
    double cos_lon = cos(lon_rad);
    double sin_lon = sin(lon_rad);
    
    // 旋转矩阵是正交矩阵，逆矩阵等于转置
    R[0][0] = -sin_lon;
    R[0][1] = -sin_lat * cos_lon;
    R[0][2] = cos_lat * cos_lon;
    
    R[1][0] = cos_lon;
    R[1][1] = -sin_lat * sin_lon;
    R[1][2] = cos_lat * sin_lon;
    
    R[2][0] = 0.0;
    R[2][1] = cos_lat;
    R[2][2] = sin_lat;
    
    return R;
}

// ECEF坐标差转东北天坐标
ENUPoint UavPathPlanner::ecefToENU(const ECEFPoint& delta_ecef, double ref_lat_rad, double ref_lon_rad) {
    auto R = UavPathPlanner::computeENURotationMatrix(ref_lat_rad, ref_lon_rad);
    
    ENUPoint enu;
    enu.east = R[0][0] * delta_ecef.x + R[0][1] * delta_ecef.y + R[0][2] * delta_ecef.z;
    enu.north = R[1][0] * delta_ecef.x + R[1][1] * delta_ecef.y + R[1][2] * delta_ecef.z;
    enu.up = R[2][0] * delta_ecef.x + R[2][1] * delta_ecef.y + R[2][2] * delta_ecef.z;
    
    return enu;
}

// 东北天坐标转ECEF坐标差
ECEFPoint UavPathPlanner::enuToECEF(const ENUPoint& enu, double ref_lat_rad, double ref_lon_rad) {
    auto R_inv = UavPathPlanner::computeENURotationMatrixInverse(ref_lat_rad, ref_lon_rad);
    
    ECEFPoint delta_ecef;
    delta_ecef.x = R_inv[0][0] * enu.east + R_inv[0][1] * enu.north + R_inv[0][2] * enu.up;
    delta_ecef.y = R_inv[1][0] * enu.east + R_inv[1][1] * enu.north + R_inv[1][2] * enu.up;
    delta_ecef.z = R_inv[2][0] * enu.east + R_inv[2][1] * enu.north + R_inv[2][2] * enu.up;
    
    return delta_ecef;
}

// 主转换函数：WGS84经纬度转东北天坐标
ENUPoint UavPathPlanner::wgs84ToENU(const WGS84Point& target, const WGS84Point& reference) {
    // 将参考点和目标点转换为ECEF坐标
    ECEFPoint ref_ecef = UavPathPlanner::wgs84ToECEF(reference);
    ECEFPoint target_ecef = UavPathPlanner::wgs84ToECEF(target);
    
    // 计算ECEF坐标差
    ECEFPoint delta_ecef;
    delta_ecef.x = target_ecef.x - ref_ecef.x;
    delta_ecef.y = target_ecef.y - ref_ecef.y;
    delta_ecef.z = target_ecef.z - ref_ecef.z;
    
    // 转换为东北天坐标
    double ref_lat_rad = deg2rad(reference.lat);
    double ref_lon_rad = deg2rad(reference.lon);
    
    return UavPathPlanner::ecefToENU(delta_ecef, ref_lat_rad, ref_lon_rad);
}

// 东北天坐标转WGS84经纬度
WGS84Point UavPathPlanner::enuToWGS84(const ENUPoint& enu, const WGS84Point& reference) {
    // 将参考点转换为ECEF坐标
    ECEFPoint ref_ecef = UavPathPlanner::wgs84ToECEF(reference);
    
    // 将东北天坐标转换为ECEF坐标差
    double ref_lat_rad = deg2rad(reference.lat);
    double ref_lon_rad = deg2rad(reference.lon);
    ECEFPoint delta_ecef = UavPathPlanner::enuToECEF(enu, ref_lat_rad, ref_lon_rad);
    
    // 计算目标点的ECEF坐标
    ECEFPoint target_ecef;
    target_ecef.x = ref_ecef.x + delta_ecef.x;
    target_ecef.y = ref_ecef.y + delta_ecef.y;
    target_ecef.z = ref_ecef.z + delta_ecef.z;
    
    // 将ECEF坐标转换为经纬度
    return ecefToWGS84(target_ecef);
}
// 批量WGS84转ENU
std::vector<ENUPoint> UavPathPlanner::wgs84ToENU_Batch(const std::vector<WGS84Point>& targets, 
                                      const WGS84Point& reference) {
    std::vector<ENUPoint> results;
    results.reserve(targets.size()); // 预分配空间以提高效率
    
    for (const auto& target : targets) {
        results.push_back(wgs84ToENU(target, reference));
    }
    
    return results;
}

// 批量ENU转WGS84
std::vector<WGS84Point> UavPathPlanner::enuToWGS84_Batch(const std::vector<ENUPoint>& targets, 
                                        const WGS84Point& reference) {
    std::vector<WGS84Point> results;
    results.reserve(targets.size()); // 预分配空间以提高效率
    
    for (const auto& target : targets) {
        results.push_back(enuToWGS84(target, reference));
    }
    
    return results;
}

// 生成 arc-line-arc（相切圆 - 直线 - 相切圆）路径实现,用于任务区域过渡
std::vector<ENUPoint> UavPathPlanner::generateArcLineArc(const ENUPoint &p0, double heading0,
                                                         const ENUPoint &p1, const ENUPoint &p2,
                                                         double radius, double resolution)
{
    std::vector<ENUPoint> path;
    if (radius <= 0.0) {
        // no radius constraint, fallback to straight line
        double dx = p1.east - p0.east;
        double dy = p1.north - p0.north;
        double dist = std::hypot(dx, dy);
        int steps = std::max(1, (int)std::ceil(dist / resolution));
        for (int i = 0; i <= steps; ++i) {
            double t = (double)i / (double)steps;
            ENUPoint q{p0.east + t * dx, p0.north + t * dy, p0.up + t * (p1.up - p0.up)};
            path.push_back(q);
        }
        return path;
    }

    // headings
    double h0 = heading0;
    double h1 = std::atan2(p2.north - p1.north, p2.east - p1.east);

    auto rotate90 = [](double ax, double ay, int sign)->std::pair<double,double>{
        if (sign >= 0) return std::make_pair(-ay, ax); // +90
        else return std::make_pair(ay, -ax); // -90
    };

    // Try combinations of left/right turn directions to find a feasible tangent
    // s = +1 => left-turn (CCW), s = -1 => right-turn (CW)
    bool found = false;
    std::pair<double,double> C1, C2; // centers
    std::pair<double,double> T1, T2; // tangent points
    int best_s0 = 0;
    int best_s1 = 0;
    double best_cost = std::numeric_limits<double>::infinity();

    auto tangent_at = [&](double theta, int sign)->std::pair<double,double>{
        if (sign > 0) return std::make_pair(-sin(theta), cos(theta));
        else return std::make_pair(sin(theta), -cos(theta));
    };

    struct TangentLine { std::pair<double,double> t1, t2; };

    for (int s0 : {1, -1}) {
        auto n0 = rotate90(std::cos(h0), std::sin(h0), s0);
        std::pair<double,double> c1 = {p0.east + radius * n0.first, p0.north + radius * n0.second};

        for (int s1 : {1, -1}) {
            auto n1 = rotate90(std::cos(h1), std::sin(h1), s1);
            std::pair<double,double> c2 = {p1.east + radius * n1.first, p1.north + radius * n1.second};

            double vx = c2.first - c1.first;
            double vy = c2.second - c1.second;
            double d = std::hypot(vx, vy);
            if (d < 1e-6) continue;

            std::vector<TangentLine> candidates;

            if (s0 == s1) {
                // External tangents (equal radius)
                for (int sign : {1, -1}) {
                    auto vperp = rotate90(vx / d, vy / d, sign);
                    candidates.push_back({
                        {c1.first + radius * vperp.first, c1.second + radius * vperp.second},
                        {c2.first + radius * vperp.first, c2.second + radius * vperp.second}
                    });
                }
            } else {
                // Internal tangents (equal radius) exist only when circles are separated enough
                if (d <= 2.0 * radius + 1e-9) continue;
                double phi = std::atan2(vy, vx);
                double alpha = std::acos((2.0 * radius) / d);
                for (int sign : {1, -1}) {
                    double ang = phi + sign * alpha;
                    double ux = std::cos(ang);
                    double uy = std::sin(ang);
                    candidates.push_back({
                        {c1.first + radius * ux, c1.second + radius * uy},
                        {c2.first - radius * ux, c2.second - radius * uy}
                    });
                }
            }

            for (const auto &line : candidates) {
                double lx = line.t2.first - line.t1.first;
                double ly = line.t2.second - line.t1.second;
                double l_len = std::hypot(lx, ly);
                if (l_len < 1e-6) continue;
                double l_dx = lx / l_len;
                double l_dy = ly / l_len;

                // Check alignment at T1
                double theta_t1 = std::atan2(line.t1.second - c1.second, line.t1.first - c1.first);
                auto tan1 = tangent_at(theta_t1, s0);
                if (tan1.first * l_dx + tan1.second * l_dy < 0.99) continue;

                // Check alignment at T2
                double theta_t2 = std::atan2(line.t2.second - c2.second, line.t2.first - c2.first);
                auto tan2 = tangent_at(theta_t2, s1);
                if (tan2.first * l_dx + tan2.second * l_dy < 0.99) continue;

                // Evaluate total cost (arc lengths + line)
                double theta0 = std::atan2(p0.north - c1.second, p0.east - c1.first);
                double delta0 = theta_t1 - theta0;
                while (delta0 <= -M_PI) delta0 += 2 * M_PI;
                while (delta0 > M_PI) delta0 -= 2 * M_PI;
                if (s0 > 0 && delta0 < 0) delta0 += 2 * M_PI;
                if (s0 < 0 && delta0 > 0) delta0 -= 2 * M_PI;

                double theta1 = std::atan2(p1.north - c2.second, p1.east - c2.first);
                double delta1 = theta1 - theta_t2;
                while (delta1 <= -M_PI) delta1 += 2 * M_PI;
                while (delta1 > M_PI) delta1 -= 2 * M_PI;
                if (s1 > 0 && delta1 < 0) delta1 += 2 * M_PI;
                if (s1 < 0 && delta1 > 0) delta1 -= 2 * M_PI;

                double cost = std::fabs(delta0) * radius + l_len + std::fabs(delta1) * radius;
                if (cost < best_cost) {
                    best_cost = cost;
                    found = true;
                    C1 = c1;
                    C2 = c2;
                    T1 = line.t1;
                    T2 = line.t2;
                    best_s0 = s0;
                    best_s1 = s1;
                }
            }
        }
    }

    if (!found) {
        // fallback: straight line
        double dx = p1.east - p0.east;
        double dy = p1.north - p0.north;
        double dist = std::hypot(dx, dy);
        int steps = std::max(1, (int)std::ceil(dist / resolution));
        for (int i = 0; i <= steps; ++i) {
            double t = (double)i / (double)steps;
            ENUPoint q{p0.east + t * dx, p0.north + t * dy, p0.up + t * (p1.up - p0.up)};
            path.push_back(q);
        }
        return path;
    }

    // sample arc from p0 to T1 on circle C1
    double theta0 = atan2(p0.north - C1.second, p0.east - C1.first);
    double theta_t1 = atan2(T1.second - C1.second, T1.first - C1.first);
    double delta0 = theta_t1 - theta0;
    while (delta0 <= -M_PI) delta0 += 2*M_PI;
    while (delta0 > M_PI) delta0 -= 2*M_PI;

    // Fix direction based on best_s0
    if (best_s0 > 0 && delta0 < 0) delta0 += 2*M_PI;
    if (best_s0 < 0 && delta0 > 0) delta0 -= 2*M_PI;

    double arc_len0 = std::fabs(delta0) * radius;
    int steps0 = std::max(1, (int)std::ceil(arc_len0 / resolution));
    for (int i = 0; i <= steps0; ++i) {
        double t = (double)i / (double)steps0;
        double theta = theta0 + delta0 * t;
        ENUPoint q{C1.first + radius * cos(theta), C1.second + radius * sin(theta), p0.up + (p1.up - p0.up) * ( (double)i / steps0 * 0.1 )};
        path.push_back(q);
    }

    // straight from T1 to T2
    double lx = T2.first - T1.first; double ly = T2.second - T1.second;
    double ldist = std::hypot(lx, ly);
    int lsteps = std::max(1, (int)std::ceil(ldist / resolution));
    for (int i = 1; i <= lsteps; ++i) {
        double t = (double)i / (double)lsteps;
        ENUPoint q{T1.first + t * lx, T1.second + t * ly, p0.up + t * (p1.up - p0.up)};
        path.push_back(q);
    }

    // arc from T2 to p1
    double theta_t2 = atan2(T2.second - C2.second, T2.first - C2.first);
    double theta1 = atan2(p1.north - C2.second, p1.east - C2.first);
    double delta1 = theta1 - theta_t2;
    while (delta1 <= -M_PI) delta1 += 2*M_PI;
    while (delta1 > M_PI) delta1 -= 2*M_PI;

    // Fix direction based on best_s1
    if (best_s1 > 0 && delta1 < 0) delta1 += 2*M_PI;
    if (best_s1 < 0 && delta1 > 0) delta1 -= 2*M_PI;

    double arc_len1 = std::fabs(delta1) * radius;
    int steps1 = std::max(1, (int)std::ceil(arc_len1 / resolution));
    for (int i = 1; i <= steps1; ++i) {
        double t = (double)i / (double)steps1;
        double theta = theta_t2 + delta1 * t;
        ENUPoint q{C2.first + radius * cos(theta), C2.second + radius * sin(theta), p1.up};
        path.push_back(q);
    }

    return path;
}

// UavPathPlanner::runAltitudeOptimization: 现在只接受一个高程文件路径参数 (.tif 或 PGM)
UavPathPlanner::AltitudeParams UavPathPlanner::makeAltitudeParams() const
{
    AltitudeParams params;
    params.uav_R = this->config_.altitude_optimization.uav_R;
    params.safe_distance = this->config_.altitude_optimization.safe_distance;
    params.lambda_follow = this->config_.altitude_optimization.lambda_follow;
    params.lambda_smooth = this->config_.altitude_optimization.lambda_smooth;
    params.max_climb_rate = this->config_.altitude_optimization.max_climb_rate;

    if (this->input_data_.ao_uav_R > 0.0) params.uav_R = this->input_data_.ao_uav_R;
    if (this->input_data_.ao_safe_distance > 0.0) params.safe_distance = this->input_data_.ao_safe_distance;
    if (this->input_data_.ao_lambda_follow >= 0.0) params.lambda_follow = this->input_data_.ao_lambda_follow;
    if (this->input_data_.ao_lambda_smooth > 0.0) params.lambda_smooth = this->input_data_.ao_lambda_smooth;
    if (this->input_data_.ao_max_climb_rate > 0.0) params.max_climb_rate = this->input_data_.ao_max_climb_rate;

    return params;
}

bool UavPathPlanner::optimizeSegmentAltitudeENU(std::vector<ENUPoint> &segment_enu)
{
    if (segment_enu.empty()) return false;

    std::vector<Eigen::Vector3d> wpts;
    wpts.reserve(segment_enu.size());
    for (const auto &p : segment_enu) {
        wpts.emplace_back(p.east, p.north, p.up);
    }

    AltitudeParams params = this->makeAltitudeParams();

    std::vector<double> out_z;
    if (!this->optimizeHeights(wpts, params, out_z)) {
        return false;
    }

    for (size_t i = 0; i < out_z.size() && i < segment_enu.size(); ++i) {
        segment_enu[i].up = out_z[i];
    }

    // Second pass: Global Smoothing
    std::vector<double> out_z_smooth;
    AltitudeParams params_smooth = params;
    params_smooth.lambda_smooth *= 10.0;
    params_smooth.max_climb_rate *= 0.5;
    if (this->optimizeHeightsGlobalSmooth(out_z, wpts, params_smooth, out_z_smooth)) {
        for (size_t i = 0; i < out_z_smooth.size() && i < segment_enu.size(); ++i) {
            segment_enu[i].up = out_z_smooth[i];
        }
    }

    return true;
}
// 以 OutputData 为主，优化后直接回写 OutputData（并同步 using_midway_lines）
bool UavPathPlanner::optimizeAndApplyOutputSegment(OutputData &output_data, const char *key, int segment_id, bool keep_closed_equal_height)
{
    if (key == nullptr) return false;

    auto resolveVec = [&](const std::string &k) -> std::vector<WGS84Coord> * {
        if (k == "uav_leader_plane1") return &output_data.uav_leader_plane1;
        if (k == "uav_leader_plane2") return &output_data.uav_leader_plane2;
        if (k == "uav_leader_plane3") return &output_data.uav_leader_plane3;
        return nullptr;
    };

    std::vector<WGS84Coord> *seg_vec = resolveVec(key);
    if (seg_vec == nullptr || seg_vec->empty()) return false;

    std::vector<WGS84Point> seg_wgs;
    seg_wgs.reserve(seg_vec->size());
    for (const auto &pt : *seg_vec) {
        seg_wgs.push_back({pt.lon, pt.lat, pt.alt});
    }
    if (seg_wgs.empty()) return false;

    std::vector<ENUPoint> seg_enu = this->wgs84ToENU_Batch(seg_wgs, this->origin_);
    if (!this->optimizeSegmentAltitudeENU(seg_enu)) return false;

    if (keep_closed_equal_height && seg_enu.size() >= 2) {
        const auto &p0 = seg_enu.front();
        const auto &pn = seg_enu.back();
        if (std::hypot(p0.east - pn.east, p0.north - pn.north) < 1e-6) {
            seg_enu.back().up = seg_enu.front().up;
        }
    }

    seg_wgs = this->enuToWGS84_Batch(seg_enu, this->origin_);
    seg_vec->clear();
    seg_vec->reserve(seg_wgs.size());
    for (const auto &p : seg_wgs) {
        seg_vec->emplace_back(WGS84Coord(p.lon, p.lat, p.alt));
    }

    this->upsertUsingMidwayLine(output_data, this->input_data_.uav_leader_id, segment_id, seg_wgs);

    return true;
}

bool UavPathPlanner::optimizeAndApplyJointSegments(OutputData &output_data, const std::vector<std::string> &keys, const std::vector<int> &segment_ids, int equal_height_segment_idx)
{
    if (keys.empty() || keys.size() != segment_ids.size()) return false;

    auto resolveVec = [&](const std::string &k) -> std::vector<WGS84Coord> * {
        if (k == "uav_leader_plane1") return &output_data.uav_leader_plane1;
        if (k == "uav_leader_plane2") return &output_data.uav_leader_plane2;
        if (k == "uav_leader_plane3") return &output_data.uav_leader_plane3;
        return nullptr;
    };

    std::vector<WGS84Point> joint_wgs;
    std::vector<size_t> segment_end_indices;

    for (const auto &key : keys) {
        std::vector<WGS84Coord> *seg = resolveVec(key);
        if (seg == nullptr || seg->empty()) return false;
        for (const auto &pt : *seg) {
            joint_wgs.push_back({pt.lon, pt.lat, pt.alt});
        }
        segment_end_indices.push_back(joint_wgs.size());
    }

    if (joint_wgs.empty()) return false;

    std::vector<ENUPoint> joint_enu = this->wgs84ToENU_Batch(joint_wgs, this->origin_);
    
    std::vector<Eigen::Vector3d> wpts;
    wpts.reserve(joint_enu.size());
    for (const auto &p : joint_enu) {
        wpts.emplace_back(p.east, p.north, p.up);
    }

    AltitudeParams params = this->makeAltitudeParams();
    std::vector<double> out_z;

    // 如果指定了某一段需要等高（如巡逻段），在优化前处理
    // 目前 optimizeHeights 是通用的 QP，如果需要严格等高，可以在优化后强制拉平，或者在 QP 中加入等式约束。
    // 这里采用简单策略：先联合优化保证连续性，然后对等高段取最大高度拉平，再做一次全局平滑。

    if (!this->optimizeHeights(wpts, params, out_z)) {
        return false;
    }

    // 处理等高约束 (uav_leader_plane3)
    if (equal_height_segment_idx >= 0 && equal_height_segment_idx < (int)segment_end_indices.size()) {
        size_t start_idx = (equal_height_segment_idx == 0) ? 0 : segment_end_indices[equal_height_segment_idx - 1];
        size_t end_idx = segment_end_indices[equal_height_segment_idx];
        if (end_idx > start_idx) {
            double max_h = -1e9;
            for (size_t i = start_idx; i < end_idx; ++i) {
                max_h = std::max(max_h, out_z[i]);
            }
            for (size_t i = start_idx; i < end_idx; ++i) {
                out_z[i] = max_h;
            }
        }
    }

    // 全局平滑
    std::vector<double> out_z_smooth;
    AltitudeParams params_smooth = params;
    params_smooth.lambda_smooth *= 10.0;
    params_smooth.max_climb_rate *= 0.5;
    if (this->optimizeHeightsGlobalSmooth(out_z, wpts, params_smooth, out_z_smooth)) {
        out_z = out_z_smooth;
    }

    // 再次检查等高约束（平滑可能会破坏它，如果是严格等高，平滑后需要再次拉平）
    if (equal_height_segment_idx >= 0 && equal_height_segment_idx < (int)segment_end_indices.size()) {
        size_t start_idx = (equal_height_segment_idx == 0) ? 0 : segment_end_indices[equal_height_segment_idx - 1];
        size_t end_idx = segment_end_indices[equal_height_segment_idx];
        if (end_idx > start_idx) {
            double final_h = out_z[start_idx]; // 取平滑后的第一个点高度作为基准
            for (size_t i = start_idx; i < end_idx; ++i) {
                out_z[i] = final_h;
            }
        }
    }

    // 保证分段衔接点高度连续：如果前一段末点与后一段起点在平面位置上相同/极近，则强制两者高度一致。
    // 典型场景：plane2 末尾切入点 == plane3 起点。
    for (size_t si = 1; si < segment_end_indices.size(); ++si) {
        const size_t boundary = segment_end_indices[si - 1];
        if (boundary == 0 || boundary >= joint_enu.size()) continue;
        const ENUPoint &a = joint_enu[boundary - 1];
        const ENUPoint &b = joint_enu[boundary];
        const double dxy = std::hypot(a.east - b.east, a.north - b.north);
        if (dxy < 0.5) {
            out_z[boundary - 1] = out_z[boundary];
        }
    }

    for (size_t i = 0; i < out_z.size(); ++i) {
        joint_enu[i].up = out_z[i];
    }

    joint_wgs = this->enuToWGS84_Batch(joint_enu, this->origin_);

    // 回写到 OutputData 和 using_midway_lines
    size_t current_offset = 0;
    for (size_t i = 0; i < keys.size(); ++i) {
        size_t next_offset = segment_end_indices[i];
        std::vector<WGS84Point> seg_wgs(joint_wgs.begin() + current_offset, joint_wgs.begin() + next_offset);

        std::vector<WGS84Coord> *dst = resolveVec(keys[i]);
        if (dst) {
            dst->clear();
            dst->reserve(seg_wgs.size());
            for (const auto &p : seg_wgs) {
                dst->emplace_back(WGS84Coord(p.lon, p.lat, p.alt));
            }
        }

        int seg_id = segment_ids[i];
        this->upsertUsingMidwayLine(output_data, this->input_data_.uav_leader_id, seg_id, seg_wgs);
        current_offset = next_offset;
    }

    return true;
}

bool UavPathPlanner::runAltitudeOptimization(const std::string &elev_file, ::OutputData &output_data)
{
    try {
        if (elev_file.empty()) return false;
        if (this->Trajectory_ENU.empty()) {
            std::cerr << "runAltitudeOptimization: Trajectory_ENU is empty, nothing to optimize\n";
            return false;
        }

        if (!elev_cost_map_) {
            elev_cost_map_ = std::make_unique<ElevationCostMap>();
        }
        // Avoid redundant I/O: if already loaded/valid, reuse it.
        if (!elev_cost_map_->isElevationValid()) {
            if (!elev_cost_map_->loadElevationData(elev_file)) {
                std::cerr << "AltitudeOptimizer: failed to load elevation file: " << elev_file << "\n";
                return false;
            }
        }

        // Instead of copying the whole map (which is in WGS84), build a local ENU costmap
        // around the trajectory to solve coordinate system mismatch.
        this->buildLocalENUCostMap(1000.0, 10.0);

        if (!this->optimizeSegmentAltitudeENU(this->Trajectory_ENU)) {
            std::cerr << "AltitudeOptimizer: optimization failed\n";
            return false;
        }

        // Convert updated ENU trajectory to WGS84
        std::vector<WGS84Point> Trajectory_WGS84 = this->enuToWGS84_Batch(this->Trajectory_ENU, this->origin_);

        this->writeLeaderPlane1(output_data, Trajectory_WGS84);
        this->writeFollowerPlane1(output_data, this->Trajectory_ENU, Trajectory_WGS84);
        std::cerr << "AltitudeOptimizer: updated output_data with new trajectories.\n";

        return true;
    } catch (const std::exception &e) {
        std::cerr << "runAltitudeOptimization exception: " << e.what() << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}
//高度优化,根据代价地图高程值
bool UavPathPlanner::optimizeHeights(const std::vector<Eigen::Vector3d> &waypoints,
                                       const AltitudeParams &p,
                                       std::vector<double> &out_z)
{
  size_t n = waypoints.size();
  if (n == 0) return false;
  // Build Hessian H and rhs b for quadratic objective: 0.5*z^T H z - b^T z
  Eigen::SparseMatrix<double> H(n, n);
  std::vector<Eigen::Triplet<double>> trips;
  Eigen::VectorXd b = Eigen::VectorXd::Zero(n);

  // smoothing: discrete second derivative penalty -> use operator L where (L z)_i = z_{i-1} - 2 z_i + z_{i+1}
  // For interior points i=1..n-2, add lambda_smooth * (L^T L) contributions
  if (n >= 3 && p.lambda_smooth > 0.0) {
    for (size_t i = 1; i + 1 < n; ++i) {
      double s = p.lambda_smooth;
      trips.emplace_back(i-1, i-1, s * 1.0);
      trips.emplace_back(i-1, i,   s * -2.0);
      trips.emplace_back(i-1, i+1, s * 1.0);

      trips.emplace_back(i, i-1,   s * -2.0);
      trips.emplace_back(i, i,     s * 4.0);
      trips.emplace_back(i, i+1,   s * -2.0);

      trips.emplace_back(i+1, i-1, s * 1.0);
      trips.emplace_back(i+1, i,   s * -2.0);
      trips.emplace_back(i+1, i+1, s * 1.0);
    }
  }

  // terrain follow weight: lambda_follow * (z - elev)^2 -> contributes lambda_follow * I to H and lambda_follow * elev to b
  for (size_t i = 0; i < n; ++i) {
    double wx = waypoints[i](0);
    double wy = waypoints[i](1);
    double elev = 0.0;
    bool has_elev = false;
    float celev;
    
    // Try CostMap first (which is now in ENU)
    if (elev_cost_map_ && elev_cost_map_->getCostAt(wx, wy, celev)) {
        elev = static_cast<double>(celev);
        has_elev = true;
    } else {
        // Fallback to ElevationMap (need to convert ENU wx,wy to WGS84)
        ENUPoint enu_pt;
        enu_pt.east = wx;
        enu_pt.north = wy;
        enu_pt.up = 0.0;
        WGS84Point wgs_pt = this->enuToWGS84(enu_pt, this->origin_);
        
        if (elev_cost_map_ && elev_cost_map_->getElevationAt(wgs_pt.lon, wgs_pt.lat, elev)) {
            has_elev = true;
        }
    }

    if (has_elev) {
      double s = p.lambda_follow;
      // follow target is terrain + UAV clearance (uav_R)
      // If original height is already safe (higher than elev + safe_distance), use original height as target
      // to avoid pulling the drone down unnecessarily.
      double safe_h = elev + p.safe_distance;
      double target = std::max(waypoints[i](2), safe_h);
      
      trips.emplace_back(i, i, s);
      b(i) += s * target;
    }
  }

  // max climb rate soft-constraint: penalize large vertical differences between consecutive points
  // term for edge (i,i+1): ((z_{i+1}-z_i)/(dist * max_climb_rate))^2
  // which expands to weight * (z_i^2 - 2 z_i z_{i+1} + z_{i+1}^2)
  if (p.max_climb_rate > 0.0) {
    for (size_t i = 0; i + 1 < n; ++i) {
      double x0 = waypoints[i](0), y0 = waypoints[i](1);
      double x1 = waypoints[i+1](0), y1 = waypoints[i+1](1);
      double dist = std::hypot(x1 - x0, y1 - y0);
      if (dist <= 1e-9) continue;
      double denom = (dist * p.max_climb_rate);
      if (denom <= 1e-12) continue;
      double w = 1.0 / (denom * denom); // weight for (z_{i+1}-z_i)^2
      trips.emplace_back(i, i,     w);
      trips.emplace_back(i, i+1,  -w);
      trips.emplace_back(i+1, i,  -w);
      trips.emplace_back(i+1, i+1, w);
    }
  }

  // add tiny diagonal regularization triplets
  for (size_t i = 0; i < n; ++i) {
    trips.emplace_back(i, i, 1e-8);
  }
  // assemble H
  H.setFromTriplets(trips.begin(), trips.end());

  // Solve H z = b
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(H);
  if (solver.info() != Eigen::Success) {
    cerr << "AltitudeOptimizer: solver decomposition failed" << endl;
    return false;
  }
  Eigen::VectorXd z = solver.solve(b);
  if (solver.info() != Eigen::Success) {
    cerr << "AltitudeOptimizer: solver failed" << endl;
    return false;
  }
  out_z.resize(n);
  for (size_t i = 0; i < n; ++i) out_z[i] = z(i);

  // Post-optimization check: ensure z > elev + uav_R
  for (size_t i = 0; i < n; ++i) {
      double wx = waypoints[i](0);
      double wy = waypoints[i](1);
      float celev;
      double min_h = -std::numeric_limits<double>::infinity();
      
      if (elev_cost_map_ && elev_cost_map_->getCostAt(wx, wy, celev)) {
          min_h = static_cast<double>(celev) + p.safe_distance;
      } else {
          double elev;
          // Fallback to ElevationMap (need to convert ENU wx,wy to WGS84)
          ENUPoint enu_pt;
          enu_pt.east = wx;
          enu_pt.north = wy;
          enu_pt.up = 0.0;
          WGS84Point wgs_pt = this->enuToWGS84(enu_pt, this->origin_);

          if (elev_cost_map_ && elev_cost_map_->getElevationAt(wgs_pt.lon, wgs_pt.lat, elev)) {
              min_h = elev + p.safe_distance;
          }
      }
      
      if (min_h > -std::numeric_limits<double>::infinity() && out_z[i] < min_h) {
          out_z[i] = min_h;
      }
  }

  return true;
}
//全局高度优化平滑(二次优化)
bool UavPathPlanner::optimizeHeightsGlobalSmooth(const std::vector<double> &input_z,
                                                 const std::vector<Eigen::Vector3d> &waypoints,
                                                 const AltitudeParams &p,
                                                 std::vector<double> &out_z)
{
    auto start_time = std::chrono::high_resolution_clock::now();
    size_t n = input_z.size();
    if (n == 0 || waypoints.size() != n) return false;

    // Iterative approach to handle inequality constraint z >= input_z
    // We solve unconstrained (w.r.t inequality) first, then add penalties for violations
    
    std::vector<double> current_z = input_z; // Start with input
    std::vector<bool> active_constraints(n, false);
    
    // Max iterations
    int max_iter = 10;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        Eigen::SparseMatrix<double> H(n, n);
        std::vector<Eigen::Triplet<double>> trips;
        Eigen::VectorXd b = Eigen::VectorXd::Zero(n);

        // 1. Smoothing term (minimize 2nd derivative)
        if (n >= 3 && p.lambda_smooth > 0.0) {
            for (size_t i = 1; i + 1 < n; ++i) {
                double s = p.lambda_smooth;
                trips.emplace_back(i-1, i-1, s * 1.0);
                trips.emplace_back(i-1, i,   s * -2.0);
                trips.emplace_back(i-1, i+1, s * 1.0);

                trips.emplace_back(i, i-1,   s * -2.0);
                trips.emplace_back(i, i,     s * 4.0);
                trips.emplace_back(i, i+1,   s * -2.0);

                trips.emplace_back(i+1, i-1, s * 1.0);
                trips.emplace_back(i+1, i,   s * -2.0);
                trips.emplace_back(i+1, i+1, s * 1.0);
            }
        }

        // 2. Climb rate term (minimize 1st derivative)
        if (p.max_climb_rate > 0.0) {
            for (size_t i = 0; i + 1 < n; ++i) {
                double x0 = waypoints[i](0), y0 = waypoints[i](1);
                double x1 = waypoints[i+1](0), y1 = waypoints[i+1](1);
                double dist = std::hypot(x1 - x0, y1 - y0);
                if (dist <= 1e-9) continue;
                double denom = (dist * p.max_climb_rate);
                if (denom <= 1e-12) continue;
                double w = 1.0 / (denom * denom); 
                trips.emplace_back(i, i,     w);
                trips.emplace_back(i, i+1,  -w);
                trips.emplace_back(i+1, i,  -w);
                trips.emplace_back(i+1, i+1, w);
            }
        }

        // 3. Fix start and end points (Hard constraints via penalty or direct elimination)
        // Here using very large penalty to approximate hard constraint
        double fix_weight = 1e10;
        trips.emplace_back(0, 0, fix_weight);
        b(0) += fix_weight * input_z[0];
        
        trips.emplace_back(n-1, n-1, fix_weight);
        b(n-1) += fix_weight * input_z[n-1];

        // 4. Inequality constraints (z >= input_z) via penalty for active set
        double constraint_weight = 1e8;
        for (size_t i = 1; i < n - 1; ++i) {
            if (active_constraints[i]) {
                trips.emplace_back(i, i, constraint_weight);
                b(i) += constraint_weight * input_z[i];
            }
        }

        // Regularization
        for (size_t i = 0; i < n; ++i) trips.emplace_back(i, i, 1e-8);

        H.setFromTriplets(trips.begin(), trips.end());
        
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(H);
        if (solver.info() != Eigen::Success) return false;
        Eigen::VectorXd z = solver.solve(b);
        if (solver.info() != Eigen::Success) return false;

        // Check violations
        bool violation = false;
        for (size_t i = 0; i < n; ++i) {
            current_z[i] = z(i);
            if (current_z[i] < input_z[i] - 1e-3) { // Tolerance
                if (!active_constraints[i]) {
                    active_constraints[i] = true;
                    violation = true;
                }
            }
        }

        if (!violation) break; // Converged
    }
    
    // Final check to ensure strictly >= input_z
    for (size_t i = 0; i < n; ++i) {
        if (current_z[i] < input_z[i]) current_z[i] = input_z[i];
    }
    
    out_z = current_z;
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cerr << "AltitudeOptimizer: second pass (global smoothing) took " << elapsed.count() << "s\n";
    return true;
}
// patrol_mode 支持"CIRCULAR"（回形）和"BOW"（弓形）,"SINGLE"(原始的单圈巡逻)  ,传入的区域 patrol_zone 已经经过收缩处理
std::vector<ENUPoint> UavPathPlanner::gen_single_patrol(const std::vector<ENUPoint> &patrol_zone,
                                                        double distance,
                                                        const std::vector<ENUPoint> &trajectory_enu)
{
    std::vector<ENUPoint> patrol_path;
    if (patrol_zone.size() < 3) {
        std::cerr << "gen_single_patrol failed: patrol_zone.size() < 3 (size=" << patrol_zone.size() << ")" << std::endl;
        return patrol_path;
    }

    std::vector<ENUPoint> patrol_waypoints = patrol_zone;
    patrol_waypoints.push_back(patrol_waypoints[0]); // 闭合巡逻路径

    // 保证闭合处切向连续：P0 -> P1 -> ... -> P0 -> P1
    if (patrol_waypoints.size() > 2) {
        patrol_waypoints.push_back(patrol_waypoints[1]);
    }

    std::vector<ENUPoint> patrol_path_full = Minisnap_3D(patrol_waypoints, distance, this->input_data_.leader_speed);
    if (patrol_path_full.empty()) {
        std::cerr << "gen_single_patrol failed: Minisnap_3D returned empty patrol path. patrol_waypoints.size()="
                  << patrol_waypoints.size() << ", distance=" << distance
                  << ", leader_speed=" << this->input_data_.leader_speed << std::endl;
        return patrol_path;
    }

    if (patrol_waypoints.size() > 2) {
        ENUPoint target_p = patrol_waypoints[patrol_waypoints.size() - 2]; // 倒数第二个航点(P0)
        size_t best_idx = patrol_path_full.size() - 1;
        double min_dist = std::numeric_limits<double>::max();
        size_t search_start = patrol_path_full.size() / 2; // 避免匹配到起点的 P0

        for (size_t i = patrol_path_full.size(); i-- > search_start;) {
            double dx = patrol_path_full[i].east - target_p.east;
            double dy = patrol_path_full[i].north - target_p.north;
            double dz = patrol_path_full[i].up - target_p.up;
            double d = dx * dx + dy * dy + dz * dz;

            if (d < min_dist) {
                min_dist = d;
                best_idx = i;
            }
        }

        if (best_idx < patrol_path_full.size()) {
            patrol_path.assign(patrol_path_full.begin(), patrol_path_full.begin() + best_idx + 1);
        } else {
            patrol_path = patrol_path_full;
        }
    } else {
        patrol_path = patrol_path_full;
    }

    // 第三段高度与第一段（已优化后的）高度保持一致
    if (!trajectory_enu.empty()) {
        const double h_plane1 = trajectory_enu.back().up;
        for (auto &pt : patrol_path) {
            pt.up = h_plane1;
        }
    }

    if (!patrol_path.empty()) {
        patrol_path.push_back(patrol_path[0]); // 闭合巡逻路径
    } else {
        std::cerr << "gen_single_patrol failed: final patrol_path is empty." << std::endl;
    }

    return patrol_path;
}

std::vector<ENUPoint> UavPathPlanner::gen_bow_patrol(const std::vector<ENUPoint> &patrol_zone,
                                                     double distance,
                                                     const std::vector<ENUPoint> &trajectory_enu)
{
    std::vector<ENUPoint> patrol_path;
    if (patrol_zone.size() < 3) {
        std::cerr << "gen_bow_patrol failed: patrol_zone.size() < 3 (size=" << patrol_zone.size() << ")" << std::endl;
        return patrol_path;
    }

    const double patrol_width = this->config_.path_planning.patrol_width;
    if (!(patrol_width > 1e-6)) {
        std::cerr << "gen_bow_patrol failed: invalid patrol_width=" << patrol_width << std::endl;
        return patrol_path;
    }

    const double resolution = (distance > 1e-6) ? distance : 1.0;   //轨迹点密度
    const double keep_up = !trajectory_enu.empty() ? trajectory_enu.back().up : patrol_zone.front().up;  //起始点(如果有的话)的高度

    // Build polygon (already shrunk by caller)
    std::vector<Vec2d> poly_pts;
    poly_pts.reserve(patrol_zone.size());
    for (const auto &p : patrol_zone) {
        poly_pts.emplace_back(p.east, p.north);
    }
    Polygon2d poly(poly_pts);  //二维点向量转换成多边形
    if (poly.points().size() < 3) {
        std::cerr << "gen_bow_patrol failed: Polygon2d has <3 points." << std::endl;
        return patrol_path;
    }

    // 1) Scan along the dominant direction of the ACTUAL polygon boundary.
    // 不再使用最小外接矩形，而是用实际边界最长边作为扫描主方向。
    double scan_heading = 0.0;
    double longest_edge_len = 0.0;
    const auto &poly_boundary = poly.points();
    for (size_t i = 0; i < poly_boundary.size(); ++i) {
        const Vec2d &a = poly_boundary[i];
        const Vec2d &b = poly_boundary[(i + 1) % poly_boundary.size()];
        const double dx = b.x() - a.x();
        const double dy = b.y() - a.y();
        const double edge_len = std::hypot(dx, dy);
        if (edge_len > longest_edge_len + 1e-6) {
            longest_edge_len = edge_len;
            scan_heading = std::atan2(dy, dx);
        }
    }
    if (!(longest_edge_len > 1e-6) || !std::isfinite(scan_heading)) {
        std::cerr << "gen_bow_patrol failed: invalid dominant boundary edge." << std::endl;
        return patrol_path;
    }

    while (scan_heading > M_PI) scan_heading -= 2.0 * M_PI;
    while (scan_heading <= -M_PI) scan_heading += 2.0 * M_PI;

    const Vec2d d = Vec2d::CreateUnitVec2d(scan_heading);
    const Vec2d n(-d.y(), d.x());
    Vec2d origin(0.0, 0.0);
    for (const auto &pt : poly_boundary) {
        origin += pt;
    }
    origin /= static_cast<double>(poly_boundary.size());

    auto toLocal = [&](const Vec2d &p) -> Vec2d {
        Vec2d q = p - origin;
        return Vec2d(q.InnerProd(d), q.InnerProd(n));
    };

    auto toWorld = [&](const Vec2d &pl) -> Vec2d {
        return origin + d * pl.x() + n * pl.y();
    };

    auto appendPoint = [&](const ENUPoint &p) {
        if (!patrol_path.empty()) {
            const ENUPoint &last = patrol_path.back();
            const double dx = p.east - last.east;
            const double dy = p.north - last.north;
            const double dz = p.up - last.up;
            if (dx * dx + dy * dy + dz * dz < 1e-12) return;
        }
        patrol_path.push_back(p);
    };

    auto appendLine = [&](const ENUPoint &a, const ENUPoint &b) {
        const double dx = b.east - a.east;
        const double dy = b.north - a.north;
        const double dz = b.up - a.up;
        const double len = std::hypot(dx, dy);
        const int steps = std::max(1, static_cast<int>(std::ceil(len / resolution)));
        for (int i = 0; i <= steps; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(steps);
            ENUPoint p{a.east + t * dx, a.north + t * dy, a.up + t * dz};
            appendPoint(p);
        }
    };

    auto appendUTurnArcLocal = [&](const Vec2d &p0_l, int dir_sign, const Vec2d &p1_l, const ENUPoint &p0_world_ref) {
        // Rounded U-turn: half circle between scanlines at x = p0_l.x().
        // p0_l and p1_l are expected to differ mainly in y.
        const double x_c = p0_l.x();
        const double y_c = 0.5 * (p0_l.y() + p1_l.y());
        const double r = 0.5 * std::abs(p1_l.y() - p0_l.y());
        if (!(r > 1e-6)) {
            return;
        }

        const double theta0 = std::atan2(p0_l.y() - y_c, p0_l.x() - x_c);
        const double theta1_raw = std::atan2(p1_l.y() - y_c, p1_l.x() - x_c);

        // Decide CW/CCW so that tangent at start matches dir_sign.
        const Vec2d tan_ccw(-std::sin(theta0), std::cos(theta0));
        const bool ccw = (tan_ccw.x() * static_cast<double>(dir_sign) > 0.0);

        double theta1 = theta1_raw;
        double delta = 0.0;
        if (ccw) {
            while (theta1 < theta0) theta1 += 2.0 * M_PI;
            delta = theta1 - theta0;
        } else {
            while (theta1 > theta0) theta1 -= 2.0 * M_PI;
            delta = theta1 - theta0;
        }

        const double arc_len = std::abs(delta) * r;
        const int steps = std::max(1, static_cast<int>(std::ceil(arc_len / resolution)));
        for (int i = 1; i <= steps; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(steps);
            const double theta = theta0 + delta * t;
            Vec2d pl(x_c + r * std::cos(theta), y_c + r * std::sin(theta));
            Vec2d pw = toWorld(pl);
            ENUPoint p{pw.x(), pw.y(), p0_world_ref.up};
            appendPoint(p);
        }
    };

    struct OrientedSeg {
        Vec2d a_w;
        Vec2d b_w;
        double xmin_l = 0.0;
        double xmax_l = 0.0;
        double y_l = 0.0;
    };

    // Determine scanline range in local frame
    double min_x = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    for (const auto &pt : poly.points()) {
        Vec2d pl = toLocal(pt);
        min_x = std::min(min_x, pl.x());
        max_x = std::max(max_x, pl.x());
        min_y = std::min(min_y, pl.y());
        max_y = std::max(max_y, pl.y());
    }
    if (!(std::isfinite(min_x) && std::isfinite(max_x) && std::isfinite(min_y) && std::isfinite(max_y))) {
        std::cerr << "gen_bow_patrol failed: invalid polygon bounds." << std::endl;
        return patrol_path;
    }

    // 如果巡逻区域在“实际边界主方向”的法向投影宽度不足以容纳至少两条巡逻线间距，则直接返回空轨迹。
    const double short_side = max_y - min_y;
    if (!(std::isfinite(short_side)) || short_side < 2.0 * patrol_width - 1e-6) {
        std::cerr << "gen_bow_patrol: short side (" << short_side
                  << ") < 2*patrol_width (" << (2.0 * patrol_width)
                  << "), return empty patrol path." << std::endl;
        return patrol_path;
    }

    const double margin = std::max(patrol_width * 2.0, 10.0);
    const double x0 = min_x - margin;
    const double x1 = max_x + margin;

    // Extra scanline rule (towards max_y side):
    // if remaining margin to top edge + shrink_distance > patrol_width,
    // add one more scanline by expanding the already-shrunk polygon outward.
    const double shrink_dist = this->config_.path_planning.patrol_region_shrink_distance;
    bool need_extra_scanline = false;
    Polygon2d expanded_poly;
    if (shrink_dist > 1e-6) {
        const double k = std::floor((max_y - min_y) / patrol_width);
        const double last_y = min_y + k * patrol_width;
        const double remain = max_y - last_y;
        if (remain + shrink_dist > patrol_width + 1e-6) {
            need_extra_scanline = true;

            auto expandPolygon = [&](const std::vector<ENUPoint> &polygon, double expand_meters) -> std::vector<ENUPoint> {
                std::vector<ENUPoint> out;
                if (!(expand_meters > 0.0) || polygon.size() < 3) return out;

                constexpr double kScale = 1000.0;
                ClipperLib::Path subj;
                subj.reserve(polygon.size());
                for (const auto &pt : polygon) {
                    subj.emplace_back(
                        static_cast<ClipperLib::cInt>(std::llround(pt.east * kScale)),
                        static_cast<ClipperLib::cInt>(std::llround(pt.north * kScale))
                    );
                }

                ClipperLib::ClipperOffset co;
                co.AddPath(subj, ClipperLib::jtMiter, ClipperLib::etClosedPolygon);
                ClipperLib::Paths solution;
                const double offset = expand_meters * kScale;
                co.Execute(solution, offset);
                if (solution.empty()) return out;

                size_t best_idx = 0;
                double best_area = 0.0;
                for (size_t i = 0; i < solution.size(); ++i) {
                    const double a = std::abs(ClipperLib::Area(solution[i]));
                    if (a > best_area) {
                        best_area = a;
                        best_idx = i;
                    }
                }

                out.reserve(solution[best_idx].size());
                for (const auto &ipt : solution[best_idx]) {
                    ENUPoint p;
                    p.east = static_cast<double>(ipt.X) / kScale;
                    p.north = static_cast<double>(ipt.Y) / kScale;
                    p.up = keep_up;
                    out.push_back(p);
                }
                return out;
            };

            std::vector<ENUPoint> expanded_zone = expandPolygon(patrol_zone, shrink_dist);
            if (expanded_zone.size() >= 3) {
                std::vector<Vec2d> epts;
                epts.reserve(expanded_zone.size());
                for (const auto &p : expanded_zone) {
                    epts.emplace_back(p.east, p.north);
                }
                expanded_poly = Polygon2d(epts);
            } else {
                need_extra_scanline = false;
            }
        }
    }

    bool has_prev = false;
    ENUPoint prev_end{};
    Vec2d prev_end_l{};
    int prev_dir_sign = +1;

    // 2) Cover-all: keep all overlap segments per scanline.
    const double scan_y_max = need_extra_scanline ? (max_y + patrol_width + 1e-6) : (max_y + 1e-6);
    for (double y = min_y; y <= scan_y_max; y += patrol_width) {
        Vec2d p_start_w = toWorld(Vec2d(x0, y));
        Vec2d p_end_w = toWorld(Vec2d(x1, y));
        LineSegment2d scan_seg(p_start_w, p_end_w);

        const bool use_expanded = (need_extra_scanline && (y > max_y + 1e-6));
        std::vector<LineSegment2d> overlaps = use_expanded ? expanded_poly.GetAllOverlaps(scan_seg)
                                                           : poly.GetAllOverlaps(scan_seg);
        // For the extra scanline, we only extend outward in the scanline-spacing direction (local y).
        // Do NOT extend the segment length along the scan direction: trim overlaps back to the shrunk
        // polygon's local x-range [min_x, max_x].
        if (use_expanded && !overlaps.empty()) {
            std::vector<LineSegment2d> trimmed;
            trimmed.reserve(overlaps.size());
            for (const auto &seg : overlaps) {
                const Vec2d a = seg.start();
                const Vec2d b = seg.end();
                const Vec2d al = toLocal(a);
                const Vec2d bl = toLocal(b);
                double sx0 = std::min(al.x(), bl.x());
                double sx1 = std::max(al.x(), bl.x());
                const double ix0 = std::max(sx0, min_x);
                const double ix1 = std::min(sx1, max_x);
                if (ix1 - ix0 <= 1e-6) continue;
                Vec2d ta_w = toWorld(Vec2d(ix0, y));
                Vec2d tb_w = toWorld(Vec2d(ix1, y));
                trimmed.emplace_back(ta_w, tb_w);
            }
            overlaps.swap(trimmed);
        }
        if (overlaps.empty()) {
            continue;
        }

        std::vector<OrientedSeg> row;
        row.reserve(overlaps.size());
        for (const auto &seg : overlaps) {
            const Vec2d a = seg.start();
            const Vec2d b = seg.end();
            Vec2d al = toLocal(a);
            Vec2d bl = toLocal(b);
            OrientedSeg os;
            os.a_w = a;
            os.b_w = b;
            os.xmin_l = std::min(al.x(), bl.x());
            os.xmax_l = std::max(al.x(), bl.x());
            os.y_l = y;
            row.push_back(os);
        }

        // Determine snake direction for this row
        const int row_idx = static_cast<int>(std::llround((y - min_y) / patrol_width));
        const bool forward = (row_idx % 2 == 0);
        const int dir_sign = forward ? +1 : -1;

        std::sort(row.begin(), row.end(), [&](const OrientedSeg &s1, const OrientedSeg &s2) {
            if (forward) {
                return s1.xmin_l < s2.xmin_l;
            }
            return s1.xmax_l > s2.xmax_l;
        });

        auto segStartEnd = [&](const OrientedSeg &s) -> std::pair<ENUPoint, ENUPoint> {
            Vec2d al = toLocal(s.a_w);
            Vec2d bl = toLocal(s.b_w);
            Vec2d start_w = s.a_w;
            Vec2d end_w = s.b_w;
            if (dir_sign > 0) {
                if (al.x() > bl.x()) {
                    start_w = s.b_w;
                    end_w = s.a_w;
                }
            } else {
                if (al.x() < bl.x()) {
                    start_w = s.b_w;
                    end_w = s.a_w;
                }
            }
            ENUPoint s0{start_w.x(), start_w.y(), keep_up};
            ENUPoint s1p{end_w.x(), end_w.y(), keep_up};
            return {s0, s1p};
        };

        // 3) Rounded connection between rows (U-turn arc); do not close the curve for now.
        // Connect from previous row end to this row start.
        {
            auto [row_first_start, row_first_end] = segStartEnd(row.front());
            if (has_prev) {
                Vec2d cur_start_l = toLocal(Vec2d(row_first_start.east, row_first_start.north));

                // U-turn arc between scanlines at x = prev_end_l.x(), then optional x-align to current start.
                Vec2d align_end_l(prev_end_l.x(), cur_start_l.y());
                appendUTurnArcLocal(prev_end_l, prev_dir_sign, align_end_l, prev_end);

                Vec2d align_end_w = toWorld(align_end_l);
                ENUPoint align_end{align_end_w.x(), align_end_w.y(), keep_up};
                if (std::hypot(align_end.east - row_first_start.east, align_end.north - row_first_start.north) > 1e-6) {
                    appendLine(align_end, row_first_start);
                }
            } else {
                appendPoint(row_first_start);
            }
        }

        // Traverse all segments in this row, keeping all overlaps.
        for (size_t j = 0; j < row.size(); ++j) {
            auto [s0, s1p] = segStartEnd(row[j]);

            // If we're not exactly at this segment start, connect with a straight line.
            if (!patrol_path.empty()) {
                const ENUPoint &last = patrol_path.back();
                if (std::hypot(last.east - s0.east, last.north - s0.north) > 1e-6) {
                    appendLine(last, s0);
                }
            } else {
                appendPoint(s0);
            }

            appendLine(s0, s1p);
        }

        // Update prev
        if (!patrol_path.empty()) {
            prev_end = patrol_path.back();
            prev_end_l = toLocal(Vec2d(prev_end.east, prev_end.north));
            prev_dir_sign = dir_sign;
            has_prev = true;
        }
    }

    // Ensure all points share the same altitude as the optimized previous segment.
    for (auto &pt : patrol_path) {
        pt.up = keep_up;
    }

    // 最后一个点到第一个点的路径用generateArcLineArc生成，确保切向连续和平滑过渡
    if (patrol_path.size() >= 3) {
        const ENUPoint p0 = patrol_path[patrol_path.size() - 1];
        const ENUPoint p0_prev = patrol_path[patrol_path.size() - 2];
        const ENUPoint p1 = patrol_path.front();
        const ENUPoint p2 = patrol_path[1];

        const double close_dx = p1.east - p0.east;
        const double close_dy = p1.north - p0.north;
        if (std::hypot(close_dx, close_dy) > 1e-3) {
            const double seg_dx = p0.east - p0_prev.east;
            const double seg_dy = p0.north - p0_prev.north;
            double heading0 = 0.0;
            if (std::hypot(seg_dx, seg_dy) > 1e-6) {
                heading0 = std::atan2(seg_dy, seg_dx);
            } else {
                // fallback: use first segment direction
                heading0 = std::atan2(p2.north - p1.north, p2.east - p1.east) + M_PI;
            }

            double radius = this->config_.path_planning.min_turning_radius;
            if (!(radius > 1e-6)) {
                radius = 0.5 * patrol_width;
            }

            ENUPoint sp0{p0.east, p0.north, keep_up};
            ENUPoint sp1{p1.east, p1.north, keep_up};
            ENUPoint sp2{p2.east, p2.north, keep_up};
            std::vector<ENUPoint> close_path = this->generateArcLineArc(sp0, heading0, sp1, sp2, radius, resolution);
            if (!close_path.empty()) {
                // avoid duplicating the starting point
                for (size_t i = 1; i < close_path.size(); ++i) {
                    appendPoint(close_path[i]);
                }
            }
        }
    }

    return patrol_path;
}

std::vector<ENUPoint> UavPathPlanner::gen_circular_patrol(const std::vector<ENUPoint> &patrol_zone,
                                                          double distance,
                                                          const std::vector<ENUPoint> &trajectory_enu)
{
    (void)patrol_zone;
    (void)distance;
    (void)trajectory_enu;
    std::vector<ENUPoint> patrol_path;
    return patrol_path;
}

void UavPathPlanner::shrinkPolygon(std::vector<ENUPoint> &polygon, double shrink_meters) const
{
    if (shrink_meters <= 0.0) return;
    if (polygon.size() < 3) return;

    constexpr double kScale = 1000.0; // 与 algorithms/lloydsVoronoiPartition.h 一致

    ClipperLib::Path subj;
    subj.reserve(polygon.size());
    for (const auto &pt : polygon) {
        subj.emplace_back(
            static_cast<ClipperLib::cInt>(std::llround(pt.east * kScale)),
            static_cast<ClipperLib::cInt>(std::llround(pt.north * kScale))
        );
    }

    ClipperLib::ClipperOffset co;
    co.AddPath(subj, ClipperLib::jtMiter, ClipperLib::etClosedPolygon);

    ClipperLib::Paths solution;
    const double offset = -shrink_meters * kScale; // 向内收缩
    co.Execute(solution, offset);

    if (solution.empty()) {
        std::cerr << "shrinkPolygon: offset produced empty result (shrink_meters="
                  << shrink_meters << ")" << std::endl;
        return;
    }

    // 若产生多个多边形，取面积绝对值最大的一个
    size_t best_idx = 0;
    double best_area = 0.0;
    for (size_t i = 0; i < solution.size(); ++i) {
        const double a = std::abs(ClipperLib::Area(solution[i]));
        if (a > best_area) {
            best_area = a;
            best_idx = i;
        }
    }

    const double keep_up = polygon.empty() ? 0.0 : polygon.front().up;
    std::vector<ENUPoint> out;
    out.reserve(solution[best_idx].size());
    for (const auto &ipt : solution[best_idx]) {
        ENUPoint p;
        p.east = static_cast<double>(ipt.X) / kScale;
        p.north = static_cast<double>(ipt.Y) / kScale;
        p.up = keep_up;
        out.push_back(p);
    }

    if (out.size() < 3) {
        std::cerr << "shrinkPolygon: shrunk polygon has <3 points, keep original (points="
                  << out.size() << ")" << std::endl;
        return;
    }

    polygon.swap(out);
}

std::vector<ENUPoint> UavPathPlanner::computePatrolPathByMode(const std::vector<ENUPoint>& patrol_zone,
                                                              double distance,
                                                              const std::string& patrol_mode,
                                                              const std::vector<ENUPoint>& trajectory_enu)
{
    std::vector<ENUPoint> patrol_path;
    std::vector<ENUPoint> shrunk_zone = patrol_zone;
    if (patrol_zone.size() < 3) {
        std::cerr << "computePatrolPathByMode failed: patrol_zone.size() < 3 (size="
                  << patrol_zone.size() << ")" << std::endl;
        return patrol_path;
    }

    std::string mode = patrol_mode;
    if (mode.empty()) {
        mode = "SINGLE";
    }
    std::transform(mode.begin(), mode.end(), mode.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    shrinkPolygon(shrunk_zone, this->config_.path_planning.patrol_region_shrink_distance); // 先对巡逻区域进行适当收缩，避免生成的路径过于贴近边界
    if (mode == "SINGLE") {
        patrol_path = this->gen_single_patrol(shrunk_zone, distance, trajectory_enu);
    } else if (mode == "BOW") {
        patrol_path = this->gen_bow_patrol(shrunk_zone, distance, trajectory_enu);
        if (patrol_path.empty()) {
            std::cerr << "patrol_mode='BOW' is not implemented yet, fallback to SINGLE." << std::endl;
            patrol_path = this->gen_single_patrol(shrunk_zone, distance, trajectory_enu);
        }
    } else if (mode == "CIRCULAR") {
        patrol_path = this->gen_circular_patrol(shrunk_zone, distance, trajectory_enu);
        if (patrol_path.empty()) {
            std::cerr << "patrol_mode='CIRCULAR' is not implemented yet, fallback to SINGLE." << std::endl;
            patrol_path = this->gen_single_patrol(shrunk_zone, distance, trajectory_enu);
        }
    } else {
        std::cerr << "Unknown patrol_mode='" << mode << "', fallback to SINGLE." << std::endl;
        patrol_path = this->gen_single_patrol(shrunk_zone, distance, trajectory_enu);
    }
    //最后一个点到第一个点的路径用Minisnap_3D生成，确保切向连续和平滑过渡
    // patrol_path.push_back(patrol_path.front()); // 确保路径闭合

    return patrol_path;
}

//检查历史和当前航段是否进入异常区域
bool UavPathPlanner::check_change(const InputData &input_data, OutputData &output_data)
{
    output_data.abnormal_uav_plane.clear();

    if (output_data.using_midway_lines.empty()) {
        std::cout << "check_change: no using_midway_lines in output_data, skipping check.\n";
        return true;
    }

    const auto zones_wgs = parseCheckZonesFromInputData(input_data);
    if (zones_wgs.empty()) {
        std::cout << "check_change: no check zones defined in input_data, skipping check.\n";
        return true;
    }

    // 若 origin_ 仍为默认值，尝试从 using_midway_lines 第一个点初始化一个参考点
    if (std::abs(this->origin_.lon) < 1e-12 && std::abs(this->origin_.lat) < 1e-12) {
        for (const auto &line : output_data.using_midway_lines) {
            if (line.points.empty()) continue;
            const auto &p0 = line.points.front();
            this->origin_.lon = p0.lon;
            this->origin_.lat = p0.lat;
            this->origin_.alt = 0.0;
            break;
        }
    }

    struct ENUZone {
        Polygon2d poly;
        double min_h = -std::numeric_limits<double>::infinity();
        double max_h = std::numeric_limits<double>::infinity();
    };

    std::vector<ENUZone> zones_enu;
    zones_enu.reserve(zones_wgs.size());
    for (const auto &zw : zones_wgs) {
        std::vector<Vec2d> pts;
        pts.reserve(zw.polygon_wgs.size());
        for (const auto &pw : zw.polygon_wgs) {
            ENUPoint ep = this->wgs84ToENU(pw, this->origin_);
            pts.emplace_back(ep.east, ep.north);
        }
        if (pts.size() < 3) continue;
        zones_enu.push_back({Polygon2d(pts), zw.min_h, zw.max_h});
    }
    if (zones_enu.empty()) {
        return true;
    }

    const auto progress_map = parseUavProgressFromInputData(input_data);
    std::set<int> bad_uavs;

    for (const auto &line : output_data.using_midway_lines) {
        const int uav_id = line.uav_id;
        const int segment_id = line.segment_id;

        std::vector<WGS84Point> wpts;
        wpts.reserve(line.points.size());
        for (const auto &pt : line.points) {
            wpts.push_back({pt.lon, pt.lat, pt.alt});
        }
        if (wpts.size() < 2) continue;

        // 基于 uavs_plane_data 剔除“已经飞过/已完成”的航段
        size_t start_idx = 0; // 从哪个点开始检查（检查线段为 [i, i+1]）
        auto pit = progress_map.find(uav_id);
        if (pit != progress_map.end()) {
            const UavProgress &pr = pit->second;
            if (segment_id < pr.segment_id) {
                continue; // 整段已完成，跳过检测
            }
            if (segment_id == pr.segment_id) {
                if (pr.point_idx >= static_cast<int>(wpts.size())) {
                    continue; // 本段也已完成
                }
                if (pr.point_idx > 1) {
                    // 1-based point idx -> 0-based start point
                    start_idx = static_cast<size_t>(pr.point_idx - 1);
                }
            }
        }
        if (start_idx >= wpts.size() - 1) {
            continue;
        }

        std::vector<ENUPoint> epts;
        epts.reserve(wpts.size());
        for (const auto &p : wpts) {
            epts.push_back(this->wgs84ToENU(p, this->origin_));
        }

        bool collided = false;
        for (size_t i = start_idx; i + 1 < epts.size() && !collided; ++i) {
            const auto &a = epts[i];
            const auto &b = epts[i + 1];
            LineSegment2d seg(Vec2d(a.east, a.north), Vec2d(b.east, b.north));
            double seg_min_h = std::min(a.up, b.up);
            double seg_max_h = std::max(a.up, b.up);

            for (const auto &zone : zones_enu) {
                if (!altitudeRangeOverlap(seg_min_h, seg_max_h, zone.min_h, zone.max_h)) {
                    continue;
                }

                // 规则1：点在区域内
                if (zone.poly.IsPointIn(seg.start()) || zone.poly.IsPointIn(seg.end())) {
                    collided = true;
                    break;
                }

                // 规则2：线段与多边形边界相交/重叠/穿越
                if (zone.poly.DistanceTo(seg) <= 1e-8) {
                    collided = true;
                    break;
                }
            }
        }

        if (collided) {
            bad_uavs.insert(uav_id);
        }
    }

    output_data.abnormal_uav_plane.assign(bad_uavs.begin(), bad_uavs.end());
    return true;
}

std::vector<ENUPoint> UavPathPlanner::preparePlanningWaypoints(int &midwaypoint_num, int &zhandoupoint_num)
{
    std::vector<ENUPoint> Enu_waypoint;

    // 记录路径点和区域边界点数量
    midwaypoint_num = static_cast<int>(this->input_data_.leader_midway_point_wgs84.size());
    zhandoupoint_num = static_cast<int>(this->input_data_.high_zhandou_point_wgs84.size());
    std::cout << "leader_midway_point_wgs84 count: " << midwaypoint_num << std::endl;
    std::cout << "high_zhandou_point_wgs84 count: " << zhandoupoint_num << std::endl;

    // Enu_waypoint 包含路径点和区域边界点
    std::vector<WGS84Point> wgs84_points;
    wgs84_points.reserve(this->input_data_.leader_midway_point_wgs84.size() +
                         this->input_data_.high_zhandou_point_wgs84.size());
    double last_alt = 0.0;

    for (const auto &p : this->input_data_.leader_midway_point_wgs84) {
        WGS84Point wp{p.lon, p.lat, p.alt};
        last_alt = wp.alt;
        wgs84_points.push_back(wp);
    }

    std::vector<WGS84Point> add_wgs84_points;
    add_wgs84_points.reserve(this->input_data_.high_zhandou_point_wgs84.size());
    for (const auto &p : this->input_data_.high_zhandou_point_wgs84) {
        // 战斗区域边界点：先使用与中途点一致的基准高度，具体抬高在 plane3 统一按 target_up 赋值。
        add_wgs84_points.push_back({p.lon, p.lat, last_alt});
    }

    if (!wgs84_points.empty() && !add_wgs84_points.empty()) {
        WGS84Point last_pt = wgs84_points.back();
        int min_idx = -1;
        double min_dist = std::numeric_limits<double>::max();
        for (size_t i = 0; i < add_wgs84_points.size(); ++i) {
            ENUPoint enu = this->wgs84ToENU(add_wgs84_points[i], last_pt);
            double dist = enu.east * enu.east + enu.north * enu.north + enu.up * enu.up;
            if (dist < min_dist) {
                min_dist = dist;
                min_idx = static_cast<int>(i);
            }
        }

        if (min_idx != -1) {
            std::rotate(add_wgs84_points.begin(), add_wgs84_points.begin() + min_idx, add_wgs84_points.end());

            if (wgs84_points.size() >= 2 && add_wgs84_points.size() >= 2) {
                WGS84Point prev_pt = wgs84_points[wgs84_points.size() - 2];
                ENUPoint vec_in = this->wgs84ToENU(last_pt, prev_pt);

                ENUPoint vec_next = this->wgs84ToENU(add_wgs84_points[1], last_pt);
                ENUPoint vec_prev = this->wgs84ToENU(add_wgs84_points.back(), last_pt);

                double dot_next = vec_in.east * vec_next.east + vec_in.north * vec_next.north;
                double dot_prev = vec_in.east * vec_prev.east + vec_in.north * vec_prev.north;

                if (dot_prev > dot_next) {
                    std::reverse(add_wgs84_points.begin() + 1, add_wgs84_points.end());
                }
            }
        }
    }

    wgs84_points.insert(wgs84_points.end(), add_wgs84_points.begin(), add_wgs84_points.end());
    if (!wgs84_points.empty()) {
        Enu_waypoint = this->wgs84ToENU_Batch(wgs84_points, this->origin_);
    }

    // Filter Enu_waypoint to remove points too close to the next one
    if (Enu_waypoint.size() > 1) {
        std::vector<ENUPoint> filtered_wpts;
        filtered_wpts.reserve(Enu_waypoint.size());
        double min_dist = 200.0; // Minimum distance in meters for waypoints
        //只过滤路径点的间距,不影响区域边界点
        for (size_t i = 0; i < midwaypoint_num - 1; ++i) {
            const auto& p1 = Enu_waypoint[i];
            const auto& p2 = Enu_waypoint[i+1];
            double dist2d = std::hypot(p1.east - p2.east, p1.north - p2.north);

            if (dist2d > min_dist) {
                filtered_wpts.push_back(p1);
            } else {
                std::cerr << "preparePlanningWaypoints: merging waypoint " << i << " to next (dist=" << dist2d << "m)\n";
            }
        }
        // 将剩余的点（包括最后一个 midway 点和所有的 zhandou 点）添加到结果中
        int start_idx = (midwaypoint_num > 0) ? (midwaypoint_num - 1) : 0;
        for (size_t i = start_idx; i < Enu_waypoint.size(); ++i) {
            filtered_wpts.push_back(Enu_waypoint[i]);
        }

        if (filtered_wpts.size() < static_cast<size_t>(midwaypoint_num + zhandoupoint_num)) {
            std::cerr << "preparePlanningWaypoints: filtered waypoints from " << (midwaypoint_num + zhandoupoint_num)
                      << " to " << filtered_wpts.size() << " points.\n";
            Enu_waypoint = filtered_wpts;
        }
    }

    return Enu_waypoint;
}

namespace {
inline bool containsIdLocal(const std::vector<int> &v, int id) {
    for (int x : v) if (x == id) return true;
    return false;
}

inline void appendUniqueIdLocal(std::vector<int> &v, int id) {
    if (!containsIdLocal(v, id)) v.push_back(id);
}
} // namespace

void UavPathPlanner::upsertUsingMidwayLine(OutputData &output_data, int uav_id, int segment_id,
                                           const std::vector<WGS84Point> &traj) const {
    std::vector<WGS84Coord> pts;
    pts.reserve(traj.size());
    for (const auto &p : traj) {
        pts.emplace_back(WGS84Coord(p.lon, p.lat, p.alt));
    }

    for (auto &line : output_data.using_midway_lines) {
        if (line.uav_id == uav_id && line.segment_id == segment_id) {
            line.points = std::move(pts);
            return;
        }
    }

    UavTrajectoryLine new_line;
    new_line.uav_id = uav_id;
    new_line.segment_id = segment_id;
    new_line.points = std::move(pts);
    output_data.using_midway_lines.push_back(std::move(new_line));
}

void UavPathPlanner::writeLeaderPlane1(OutputData &output_data, const std::vector<WGS84Point> &traj) const {
    output_data.uav_leader_plane1.clear();
    output_data.uav_leader_plane1.reserve(traj.size());
    for (const auto &p : traj) {
        output_data.uav_leader_plane1.emplace_back(WGS84Coord(p.lon, p.lat, p.alt));
    }
    this->upsertUsingMidwayLine(output_data, this->input_data_.uav_leader_id, 1, traj);
}

void UavPathPlanner::writeLeaderSegment(std::vector<WGS84Coord> &dst_segment, OutputData &output_data,
                                        int segment_id, const std::vector<WGS84Point> &traj,
                                        bool sync_using_midway_line) const {
    dst_segment.clear();
    dst_segment.reserve(traj.size());
    for (const auto &p : traj) {
        dst_segment.emplace_back(WGS84Coord(p.lon, p.lat, p.alt));
    }

    if (sync_using_midway_line) {
        this->upsertUsingMidwayLine(output_data, this->input_data_.uav_leader_id, segment_id, traj);
    }
}

void UavPathPlanner::writeFollowerPlane1(OutputData &output_data, const std::vector<ENUPoint> &leader_traj_enu,
                                         const std::vector<WGS84Point> &leader_traj_wgs) {
    if (this->input_data_.formation_using != 1) {
        output_data.uav_plane1.clear();
        return;
    }

    try {
        json plane_array = this->generateFollowerTrajectories(this->input_data_, leader_traj_enu, leader_traj_wgs);

        output_data.uav_plane1.clear();
        if (plane_array.is_array()) {
            for (const auto &follower : plane_array) {
                if (!follower.is_array() || follower.size() < 2 || !follower[0].is_number_integer()) {
                    continue;
                }

                const int uid = follower[0].get<int>();
                UavTrajectoryLine line;
                line.uav_id = uid;
                line.segment_id = 1;

                std::vector<WGS84Point> follower_traj;
                follower_traj.reserve(follower.size() - 1);
                for (size_t i = 1; i < follower.size(); ++i) {
                    WGS84Point p;
                    if (!parseWGS84PointValue(follower[i], p)) {
                        continue;
                    }
                    follower_traj.push_back(p);
                    line.points.emplace_back(WGS84Coord(p.lon, p.lat, p.alt));
                }
                output_data.uav_plane1.push_back(std::move(line));
                if (!follower_traj.empty()) {
                    this->upsertUsingMidwayLine(output_data, uid, 1, follower_traj);
                }
            }
        }
    } catch (const std::exception &e) {
        std::cerr << "调用 generateFollowerTrajectories 出错: " << e.what() << std::endl;
    }
}

bool UavPathPlanner::getFollowerStartWgs84(int uid, WGS84Point &out) const {
    for (size_t i = 0; i < this->input_data_.uavs_id.size(); ++i) {
        if (this->input_data_.uavs_id[i] == uid) {
            if (i < this->input_data_.uav_start_point_wgs84.size()) {
                const auto &p = this->input_data_.uav_start_point_wgs84[i];
                out = {p.lon, p.lat, p.alt};
                return true;
            }
            break;
        }
    }
    return false;
}

void UavPathPlanner::adjustFollowerStartAltitudeIfNeeded(WGS84Point &p, bool formation_enabled) const {
    if (formation_enabled) return;

    // 参考长机起点高度：优先用已生成的 plane1 起点（可能已做高度优化），否则退回输入起点。
    double leader_ref_alt = this->input_data_.uav_leader_start_point_wgs84.alt;
    if (!this->output_data_.uav_leader_plane1.empty() &&
        std::isfinite(this->output_data_.uav_leader_plane1.front().alt)) {
        leader_ref_alt = this->output_data_.uav_leader_plane1.front().alt;
    }
    // 非编队时，输入可能只给了长机起点经纬度（alt=0）。此时使用 leader_midway 的首点高度作为参考高度（仅高度基准，不参与路径）。
    if ((!std::isfinite(leader_ref_alt) || std::abs(leader_ref_alt) < 1e-6) &&
        !this->input_data_.leader_midway_point_wgs84.empty() &&
        std::isfinite(this->input_data_.leader_midway_point_wgs84.front().alt) &&
        this->input_data_.leader_midway_point_wgs84.front().alt > 0.0) {
        leader_ref_alt = this->input_data_.leader_midway_point_wgs84.front().alt;
    }

    // 先保证“至少与长机同高”（仅抬升，不下调）
    if (std::isfinite(leader_ref_alt) && (!std::isfinite(p.alt) || p.alt < leader_ref_alt)) {
        p.alt = leader_ref_alt;
    }

    // 如果高程不可用，就到此为止（用户仍能看到僚机不再贴地）
    if (!this->elev_cost_map_ || !this->elev_cost_map_->isElevationValid()) return;

    AltitudeParams params = this->makeAltitudeParams();
    double min_clearance = params.safe_distance;
    if (!(min_clearance > 0.0) && (params.uav_R > 0.0)) min_clearance = params.uav_R;
    if (!(min_clearance > 0.0)) return;

    // 计算长机参考离地净空 AGL（至少 min_clearance）
    double leader_clearance = min_clearance;
    {
        const auto &ls = this->input_data_.uav_leader_start_point_wgs84;
        double leader_elev = 0.0;
        if (std::isfinite(leader_ref_alt) && this->elev_cost_map_->getElevationAt(ls.lon, ls.lat, leader_elev)) {
            const double c = leader_ref_alt - leader_elev;
            if (std::isfinite(c) && c > leader_clearance) leader_clearance = c;
        }
    }

    // 抬升僚机到“地面 + 长机 AGL”
    double elev = 0.0;
    if (!this->elev_cost_map_->getElevationAt(p.lon, p.lat, elev)) return;
    const double min_alt = elev + leader_clearance;
    if (std::isfinite(min_alt) && (!std::isfinite(p.alt) || p.alt < min_alt)) {
        p.alt = min_alt;
    }
}

UavPathPlanner::AltitudeOptContext UavPathPlanner::prepareAltitudeOptimizationContext(bool ensure_elevation_loaded,
                                                                                      bool print_config_log) {
    AltitudeOptContext ctx;
    ctx.enabled = this->config_.altitude_optimization.enabled;
    ctx.elevation_file = this->config_.altitude_optimization.elevation_file;

    if (ensure_elevation_loaded && ctx.enabled && !ctx.elevation_file.empty()) {
        if (!elev_cost_map_) elev_cost_map_ = std::make_unique<ElevationCostMap>();
        if (!elev_cost_map_->isElevationValid()) {
            if (!elev_cost_map_->loadElevationData(ctx.elevation_file)) {
                std::cerr << "Warning: failed to load elevation file: " << ctx.elevation_file << std::endl;
            }
        }
    }

    if (print_config_log) {
        if (!this->config_.loaded) {
            std::cerr << "Warning: config.yaml not loaded (" << this->config_.load_error << ")" << std::endl;
        } else {
            std::cout << "Loaded config from " << this->config_.loaded_from << std::endl;
            std::cout << "altitude_opt_enabled: " << ctx.enabled << std::endl;
            if (ctx.enabled && ctx.elevation_file.empty()) {
                std::cerr << "Warning: 'elevation_file' parameter missing in config.yaml" << std::endl;
            }
        }
    }

    return ctx;
}

void UavPathPlanner::appendUavSegmentLine(std::vector<UavTrajectoryLine> &dst, int uid, int segment_id,
                                         const std::vector<WGS84Point> &traj) const {
    UavTrajectoryLine line;
    line.uav_id = uid;
    line.segment_id = segment_id;
    line.points.reserve(traj.size());
    for (const auto &p : traj) {
        line.points.emplace_back(WGS84Coord(p.lon, p.lat, p.alt));
    }
    dst.push_back(std::move(line));
}

bool UavPathPlanner::buildTransitionAndRotatePatrol(const ENUPoint& p0, double heading0, double minR, double resolution,
                                                    const std::vector<ENUPoint>& patrol_path,
                                                    std::vector<ENUPoint> &out_transition_path,
                                                    std::vector<ENUPoint> &out_rotated_patrol) const {
    out_transition_path.clear();
    out_rotated_patrol.clear();
    if (patrol_path.empty()) return false;

    if (!(minR > 1e-6)) {
        ENUPoint p1 = patrol_path.front();
        double dx = p1.east - p0.east;
        double dy = p1.north - p0.north;
        double dist = std::hypot(dx, dy);
        int steps = std::max(1, static_cast<int>(std::ceil(dist / resolution)));
        for (int i = 0; i <= steps; ++i) {
            double t = static_cast<double>(i) / steps;
            ENUPoint q{p0.east + t * dx, p0.north + t * dy, p0.up + t * (p1.up - p0.up)};
            out_transition_path.push_back(q);
        }
        out_rotated_patrol = patrol_path;
        return false;
    }

    double best_score = std::numeric_limits<double>::max();
    size_t best_idx = 0;
    double best_arc_len = 0.0;
    double best_line_len = 0.0;
    int best_s = 0;
    double best_cx = 0.0;
    double best_cy = 0.0;
    double best_theta_start = 0.0;
    bool found_any = false;

    for (int s : {1, -1}) {
        double cx = p0.east - s * minR * std::sin(heading0);
        double cy = p0.north + s * minR * std::cos(heading0);
        double theta_start = std::atan2(p0.north - cy, p0.east - cx);

        for (size_t i = 0; i < patrol_path.size(); ++i) {
            const auto &pt = patrol_path[i];
            const auto &next_pt = patrol_path[(i + 1) % patrol_path.size()];
            double patrol_dx = next_pt.east - pt.east;
            double patrol_dy = next_pt.north - pt.north;
            double patrol_len = std::hypot(patrol_dx, patrol_dy);
            if (patrol_len < 1e-3) continue;
            patrol_dx /= patrol_len;
            patrol_dy /= patrol_len;

            double v_cx = pt.east - cx;
            double v_cy = pt.north - cy;
            double dist_cp = std::hypot(v_cx, v_cy);
            if (dist_cp <= minR) continue;

            double alpha = std::atan2(v_cy, v_cx);
            double beta = std::acos(minR / dist_cp);
            double candidates[2] = {alpha + beta, alpha - beta};
            for (double theta : candidates) {
                double tx = cx + minR * std::cos(theta);
                double ty = cy + minR * std::sin(theta);

                double lx = pt.east - tx;
                double ly = pt.north - ty;
                double l_len = std::hypot(lx, ly);
                if (l_len < 1e-3) continue;
                double l_dx = lx / l_len;
                double l_dy = ly / l_len;

                double tan_x = -s * std::sin(theta);
                double tan_y = s * std::cos(theta);
                if (tan_x * l_dx + tan_y * l_dy < 0.99) continue;

                double alignment = l_dx * patrol_dx + l_dy * patrol_dy;
                if (alignment < 0.8) continue;

                double d_theta = theta - theta_start;
                if (s > 0) {
                    while (d_theta <= 0) d_theta += 2 * M_PI;
                    while (d_theta > 2 * M_PI) d_theta -= 2 * M_PI;
                } else {
                    while (d_theta >= 0) d_theta -= 2 * M_PI;
                    while (d_theta < -2 * M_PI) d_theta += 2 * M_PI;
                }
                double arc_len = std::abs(d_theta) * minR;
                double penalty = 1000.0 * (1.0 - alignment);
                double total_cost = arc_len + l_len + penalty;

                if (total_cost < best_score) {
                    best_score = total_cost;
                    best_idx = i;
                    best_arc_len = arc_len;
                    best_line_len = l_len;
                    best_s = s;
                    best_cx = cx;
                    best_cy = cy;
                    best_theta_start = theta_start;
                    found_any = true;
                }
            }
        }
    }

    if (found_any) {
        int steps_arc = std::max(1, static_cast<int>(std::ceil(best_arc_len / resolution)));
        double d_theta_total = (best_s > 0) ? (best_arc_len / minR) : -(best_arc_len / minR);
        for (int i = 0; i <= steps_arc; ++i) {
            double t = static_cast<double>(i) / steps_arc;
            double ang = best_theta_start + d_theta_total * t;
            ENUPoint pt;
            pt.east = best_cx + minR * std::cos(ang);
            pt.north = best_cy + minR * std::sin(ang);
            pt.up = p0.up + (patrol_path[best_idx].up - p0.up) * (t * best_arc_len / (best_arc_len + best_line_len));
            out_transition_path.push_back(pt);
        }

        ENUPoint t_end = out_transition_path.back();
        ENUPoint p_target = patrol_path[best_idx];
        int steps_line = std::max(1, static_cast<int>(std::ceil(best_line_len / resolution)));
        for (int i = 1; i <= steps_line; ++i) {
            double t = static_cast<double>(i) / steps_line;
            ENUPoint pt;
            pt.east = t_end.east + t * (p_target.east - t_end.east);
            pt.north = t_end.north + t * (p_target.north - t_end.north);
            pt.up = t_end.up + t * (p_target.up - t_end.up);
            out_transition_path.push_back(pt);
        }

        out_rotated_patrol.reserve(patrol_path.size() + 1);
        for (size_t i = 0; i < patrol_path.size(); ++i) {
            out_rotated_patrol.push_back(patrol_path[(best_idx + i) % patrol_path.size()]);
        }
        if (!out_rotated_patrol.empty()) {
            out_rotated_patrol.push_back(out_rotated_patrol.front());
        }
        return true;
    }

    std::cerr << "Warning: Failed to find valid tangent transition, falling back to straight line." << std::endl;
    ENUPoint p1 = patrol_path[0];
    double dx = p1.east - p0.east;
    double dy = p1.north - p0.north;
    double dist = std::hypot(dx, dy);
    int steps = std::max(1, static_cast<int>(std::ceil(dist / resolution)));
    for (int i = 0; i <= steps; ++i) {
        double t = static_cast<double>(i) / steps;
        ENUPoint q{p0.east + t * dx, p0.north + t * dy, p0.up + t * (p1.up - p0.up)};
        out_transition_path.push_back(q);
    }
    out_rotated_patrol = patrol_path;
    return false;
}

double UavPathPlanner::computeActualMaxClimbRate(const std::vector<ENUPoint> &path) const {
    double max_rate = 0.0;
    if (path.size() < 2) return max_rate;

    for (size_t i = 1; i < path.size(); ++i) {
        const double dx = path[i].east - path[i - 1].east;
        const double dy = path[i].north - path[i - 1].north;
        const double dist_xy = std::hypot(dx, dy);
        if (dist_xy <= 1e-6) continue;
        const double dz = path[i].up - path[i - 1].up;
        max_rate = std::max(max_rate, std::abs(dz) / dist_xy);
    }
    return max_rate;
}

void UavPathPlanner::enforceTransitionClimbRateAndBorrowPatrolPrefix(std::vector<ENUPoint> &transition_path,
                                                                     std::vector<ENUPoint> &patrol_path,
                                                                     const std::string &log_label) const {
    if (transition_path.empty() || patrol_path.empty()) return;

    const double max_climb_rate = this->makeAltitudeParams().max_climb_rate;
    if (!(max_climb_rate > 0.0)) {
        std::cout << log_label << " actual max climb rate: "
                  << this->computeActualMaxClimbRate(transition_path) << std::endl;
        return;
    }

    std::vector<ENUPoint> patrol_core = patrol_path;
    auto sameXY = [](const ENUPoint &a, const ENUPoint &b) {
        return std::hypot(a.east - b.east, a.north - b.north) <= 1e-6;
    };
    const bool patrol_closed = patrol_core.size() >= 2 && sameXY(patrol_core.front(), patrol_core.back());
    if (patrol_closed) {
        patrol_core.pop_back();
    }
    if (patrol_core.empty()) {
        std::cout << log_label << " actual max climb rate: "
                  << this->computeActualMaxClimbRate(transition_path) << std::endl;
        return;
    }

    const double target_up = patrol_core.front().up;
    auto advanceTowardTarget = [&](double current_up, double dist_xy) {
        const double delta_limit = max_climb_rate * std::max(0.0, dist_xy);
        if (target_up >= current_up) return std::min(target_up, current_up + delta_limit);
        return std::max(target_up, current_up - delta_limit);
    };

    for (size_t i = 1; i < transition_path.size(); ++i) {
        const double dx = transition_path[i].east - transition_path[i - 1].east;
        const double dy = transition_path[i].north - transition_path[i - 1].north;
        const double dist_xy = std::hypot(dx, dy);
        transition_path[i].up = advanceTowardTarget(transition_path[i - 1].up, dist_xy);
    }

    const auto reachedTarget = [&](double up) {
        return std::abs(up - target_up) <= 1e-6;
    };

    if (!reachedTarget(transition_path.back().up)) {
        double loop_length = 0.0;
        for (size_t i = 0; i < patrol_core.size(); ++i) {
            const ENUPoint &a = patrol_core[i];
            const ENUPoint &b = patrol_core[(i + 1) % patrol_core.size()];
            loop_length += std::hypot(b.east - a.east, b.north - a.north);
        }
        if (loop_length <= 1e-6) {
            std::cerr << log_label << " failed to extend plane2: patrol loop length is zero." << std::endl;
        } else {
            const double remaining_h = std::abs(target_up - transition_path.back().up);
            const int max_loops = std::max(1, static_cast<int>(std::ceil(remaining_h / (max_climb_rate * loop_length))) + 1);

            ENUPoint current = transition_path.back();
            int current_idx = 0;
            bool reached = false;

            for (int loop = 0; loop < max_loops && !reached; ++loop) {
                for (size_t step = 0; step < patrol_core.size(); ++step) {
                    const int next_idx = (current_idx + 1) % static_cast<int>(patrol_core.size());
                    const ENUPoint &next_patrol = patrol_core[static_cast<size_t>(next_idx)];
                    const double dx = next_patrol.east - current.east;
                    const double dy = next_patrol.north - current.north;
                    const double dist_xy = std::hypot(dx, dy);
                    if (dist_xy <= 1e-6) {
                        current = {next_patrol.east, next_patrol.north, current.up};
                        current_idx = next_idx;
                        continue;
                    }

                    const double next_up = advanceTowardTarget(current.up, dist_xy);
                    if (!reachedTarget(next_up)) {
                        ENUPoint appended{next_patrol.east, next_patrol.north, next_up};
                        transition_path.push_back(appended);
                        current = appended;
                        current_idx = next_idx;
                        continue;
                    }

                    const double delta_up = std::abs(target_up - current.up);
                    const double step_up = std::abs(next_up - current.up);
                    double t = (step_up > 1e-9) ? (delta_up / step_up) : 1.0;
                    t = std::clamp(t, 0.0, 1.0);

                    ENUPoint split_point;
                    split_point.east = current.east + t * dx;
                    split_point.north = current.north + t * dy;
                    split_point.up = target_up;
                    if (!sameXY(split_point, transition_path.back()) || !reachedTarget(transition_path.back().up)) {
                        transition_path.push_back(split_point);
                    }

                    std::vector<ENUPoint> rebuilt_patrol;
                    rebuilt_patrol.reserve(patrol_core.size() + 2);
                    rebuilt_patrol.push_back(split_point);
                    rebuilt_patrol.push_back({next_patrol.east, next_patrol.north, target_up});
                    for (size_t k = 1; k < patrol_core.size(); ++k) {
                        const size_t idx = (static_cast<size_t>(next_idx) + k) % patrol_core.size();
                        rebuilt_patrol.push_back({patrol_core[idx].east, patrol_core[idx].north, target_up});
                    }
                    rebuilt_patrol.push_back(split_point);
                    patrol_path = std::move(rebuilt_patrol);
                    reached = true;
                    break;
                }
            }

            if (!reached) {
                std::cerr << log_label << " warning: borrowed full patrol prefix loops but still did not reach patrol altitude." << std::endl;
                patrol_path = patrol_core;
                for (auto &pt : patrol_path) pt.up = target_up;
                if (!patrol_path.empty()) patrol_path.push_back(patrol_path.front());
            }
        }
    } else {
        for (auto &pt : patrol_path) {
            pt.up = target_up;
        }
    }

    std::cout << log_label << " actual max climb rate: "
              << this->computeActualMaxClimbRate(transition_path) << std::endl;
}

void UavPathPlanner::generateLeaderPlane2Plane3NonFormation(const WGS84Point &leader_start_wgs, double distance) {
    this->output_data_.uav_leader_plane2.clear();
    this->output_data_.uav_leader_plane3.clear();

    ENUPoint p0 = this->wgs84ToENU(leader_start_wgs, this->origin_);
    const double zhandou_target_up = p0.up + this->input_data_.leader_fly_high;

    std::vector<WGS84Point> leader_battle_wgs;
    leader_battle_wgs.reserve(this->input_data_.high_zhandou_point_wgs84.size());
    for (const auto &pt : this->input_data_.high_zhandou_point_wgs84) {
        leader_battle_wgs.push_back({pt.lon, pt.lat, 0.0});
    }
    std::vector<ENUPoint> leader_battle_enu = this->wgs84ToENU_Batch(leader_battle_wgs, this->origin_);
    for (auto &pt : leader_battle_enu) pt.up = zhandou_target_up;

    std::vector<ENUPoint> ctx_enu;
    ctx_enu.push_back(p0);
    std::vector<ENUPoint> patrol_path = this->computePatrolPathByMode(
        leader_battle_enu,
        distance,
        this->config_.path_planning.patrol_mode,
        ctx_enu);
    for (auto &pt : patrol_path) pt.up = zhandou_target_up;

    if (patrol_path.empty()) {
        std::cerr << "Warning: non-formation leader patrol (plane3) empty; leader plane2/3 not generated." << std::endl;
        return;
    }

    double heading0 = 0.0;
    {
        const ENUPoint &p1 = patrol_path.front();
        const double dx = p1.east - p0.east;
        const double dy = p1.north - p0.north;
        heading0 = (std::hypot(dx, dy) > 1e-6) ? std::atan2(dy, dx) : 0.0;
    }
    const double radius = std::max(0.0, this->input_data_.min_turning_radius);
    const double resolution = (distance > 0.0) ? distance : 300.0;
    std::vector<ENUPoint> transition;
    std::vector<ENUPoint> rotated_patrol;
    this->buildTransitionAndRotatePatrol(p0, heading0, radius, resolution, patrol_path, transition, rotated_patrol);
    if (transition.empty()) return;

    if (!rotated_patrol.empty()) {
        this->enforceTransitionClimbRateAndBorrowPatrolPrefix(transition, rotated_patrol,
                                                              "leader plane2(non-formation)");
    }

    std::vector<WGS84Point> trans_wgs = this->enuToWGS84_Batch(transition, this->origin_);
    std::vector<WGS84Point> patrol_wgs = this->enuToWGS84_Batch(rotated_patrol.empty() ? patrol_path : rotated_patrol, this->origin_);
    this->writeLeaderSegment(this->output_data_.uav_leader_plane3, this->output_data_, 3, patrol_wgs);
    this->writeLeaderSegment(this->output_data_.uav_leader_plane2, this->output_data_, 2, trans_wgs);
}

void UavPathPlanner::generateFollowerPlane1(const std::vector<ENUPoint> &leader_traj_enu,
                                           const std::vector<WGS84Point> &leader_traj_wgs) {
    this->writeFollowerPlane1(this->output_data_, leader_traj_enu, leader_traj_wgs);
}

void UavPathPlanner::generateFollowerPlane2Plane3(bool formation_enabled, double final_heading, double distance,
                                                  std::vector<int> &out_final_ready_ids) {
    this->output_data_.uav_plane2.clear();
    this->output_data_.uav_plane3.clear();

    out_final_ready_ids = this->input_data_.ready_id;
    std::vector<int> battle_ids;

    // 决定每个僚机去 ready_zone 还是 battle_zone
    std::vector<int> candidate_uavs;
    if (formation_enabled) {
        if (!this->output_data_.uav_plane1.empty()) {
            candidate_uavs.reserve(this->output_data_.uav_plane1.size());
            for (const auto &line : this->output_data_.uav_plane1) candidate_uavs.push_back(line.uav_id);
        }
    } else {
        candidate_uavs = this->input_data_.uavs_id;
    }

    if (!candidate_uavs.empty()) {
        for (int uid : candidate_uavs) {

            if (containsIdLocal(this->input_data_.ready_id, uid)) {
                std::cout << "[DestDecide] uav_id=" << uid << " forced to ready_zone (in ready_id)." << std::endl;
                appendUniqueIdLocal(out_final_ready_ids, uid);
                continue;
            }

            const FlightZone *bz = this->selectBattleZoneForUav(uid);
            if (bz == nullptr) {
                std::cout << "[DestDecide] uav_id=" << uid << " no battle_zone_wgs84 provided -> assign to ready_zone." << std::endl;
                appendUniqueIdLocal(out_final_ready_ids, uid);
                continue;
            }

            ENUPoint p0{};
            double heading0 = 0.0;
            std::vector<ENUPoint> ctx_enu;
            if (!this->getFollowerCurrentState(uid, formation_enabled, final_heading, p0, heading0, ctx_enu)) {
                std::cout << "[DestDecide] uav_id=" << uid << " no valid current state -> assign to ready_zone." << std::endl;
                appendUniqueIdLocal(out_final_ready_ids, uid);
                continue;
            }

            const double battle_relative_h = 0.5 * (bz->height_range.first + bz->height_range.second);
            const double battle_target_up = p0.up + battle_relative_h;

            const bool battle_ok = this->check_battle_zone(uid, *bz, battle_target_up);
            if (!battle_ok) {
                std::cout << "[DestDecide] uav_id=" << uid << " battle_zone check failed -> assign to ready_zone." << std::endl;
                appendUniqueIdLocal(out_final_ready_ids, uid);
                continue;
            }

            std::cout << "[DestDecide] uav_id=" << uid << " assigned to battle_zone." << std::endl;
            battle_ids.push_back(uid);
        }
    }

    // 先为 battle_ids 生成 plane2/plane3；失败的再回退到 ready_ids。
    if (!battle_ids.empty()) {
        for (int rid : battle_ids) {
            const FlightZone *bz = this->selectBattleZoneForUav(rid);
            if (bz == nullptr || bz->polygon.size() < 3) {
                appendUniqueIdLocal(out_final_ready_ids, rid);
                continue;
            }

            ENUPoint p0{};
            double heading0 = 0.0;
            std::vector<ENUPoint> ctx_enu;
            if (!this->getFollowerCurrentState(rid, formation_enabled, final_heading, p0, heading0, ctx_enu)) {
                appendUniqueIdLocal(out_final_ready_ids, rid);
                continue;
            }

            const double battle_relative_h = 0.5 * (bz->height_range.first + bz->height_range.second);
            const double battle_target_up = p0.up + battle_relative_h;
            if (!this->check_battle_zone(rid, *bz, battle_target_up)) {
                appendUniqueIdLocal(out_final_ready_ids, rid);
                continue;
            }

            // battle_zone -> ENU（目标高度层为 battle_target_up）
            std::vector<WGS84Point> battle_zone_wgs;
            battle_zone_wgs.reserve(bz->polygon.size());
            for (const auto &pt : bz->polygon) {
                battle_zone_wgs.push_back({pt.lon, pt.lat, battle_target_up});
            }
            if (battle_zone_wgs.size() < 3) {
                appendUniqueIdLocal(out_final_ready_ids, rid);
                continue;
            }
            std::vector<ENUPoint> battle_zone_enu = this->wgs84ToENU_Batch(battle_zone_wgs, this->origin_);

            // 第三段：战斗区域巡逻轨迹
            std::vector<ENUPoint> battle_patrol = this->computePatrolPathByMode(
                battle_zone_enu,
                distance,
                this->config_.path_planning.patrol_mode,
                ctx_enu);

            if (battle_patrol.empty()) {
                std::cerr << "battle_id=" << rid << " failed to generate battle patrol path, fallback to ready_zone." << std::endl;
                appendUniqueIdLocal(out_final_ready_ids, rid);
                continue;
            }
            for (auto &pt : battle_patrol) pt.up = battle_target_up;

            // 非编队/无历史航段时，用指向巡逻起点的方向作为 heading0
            if (!std::isfinite(heading0) || ctx_enu.size() < 2) {
                ENUPoint p1 = battle_patrol.front();
                const double dx = p1.east - p0.east;
                const double dy = p1.north - p0.north;
                if (std::hypot(dx, dy) > 1e-6) {
                    heading0 = std::atan2(dy, dx);
                } else {
                    heading0 = 0.0;
                }
            }

            const double radius = std::max(0.0, this->input_data_.min_turning_radius);
            const double resolution = (distance > 0.0) ? distance : 300.0;
            std::vector<ENUPoint> battle_transition;
            std::vector<ENUPoint> rotated_battle_patrol;
            this->buildTransitionAndRotatePatrol(p0, heading0, radius, resolution,
                                                 battle_patrol, battle_transition, rotated_battle_patrol);
            if (battle_transition.empty()) {
                std::cerr << "battle_id=" << rid << " failed to generate battle transition path, fallback to ready_zone." << std::endl;
                appendUniqueIdLocal(out_final_ready_ids, rid);
                continue;
            }

            if (!rotated_battle_patrol.empty()) {
                this->enforceTransitionClimbRateAndBorrowPatrolPrefix(
                    battle_transition, rotated_battle_patrol,
                    std::string("uav ") + std::to_string(rid) + " battle plane2");
            }

            std::vector<WGS84Point> battle_transition_wgs = this->enuToWGS84_Batch(battle_transition, this->origin_);
            std::vector<WGS84Point> battle_patrol_wgs = this->enuToWGS84_Batch(
                rotated_battle_patrol.empty() ? battle_patrol : rotated_battle_patrol, this->origin_);

            this->upsertUsingMidwayLine(this->output_data_, rid, 2, battle_transition_wgs);
            this->upsertUsingMidwayLine(this->output_data_, rid, 3, battle_patrol_wgs);

            this->appendUavSegmentLine(this->output_data_.uav_plane2, rid, 2, battle_transition_wgs);
            this->appendUavSegmentLine(this->output_data_.uav_plane3, rid, 3, battle_patrol_wgs);
        }
    }

    // 再为 out_final_ready_ids 生成 plane2/plane3
    if (!out_final_ready_ids.empty() && this->input_data_.ready_zone.polygon.size() >= 3) {

        // ready_zone 相对高度：使用 ready_high_list 的平均值
        const double ready_relative_h = 0.5 * (this->input_data_.ready_zone.height_range.first + this->input_data_.ready_zone.height_range.second);

        for (int rid : out_final_ready_ids) {
            // 如果该机已经生成了 battle plane2/3，则跳过 ready 生成
            bool already_has_plane2 = false;
            for (const auto &l2 : this->output_data_.uav_plane2) {
                if (l2.uav_id == rid) { already_has_plane2 = true; break; }
            }
            if (already_has_plane2) continue;

            ENUPoint p0{};
            double heading0 = 0.0;
            std::vector<ENUPoint> ctx_enu;
            if (!this->getFollowerCurrentState(rid, formation_enabled, final_heading, p0, heading0, ctx_enu)) {
                std::cerr << "ready_id=" << rid << " no valid current state, skip ready plane2/3 generation." << std::endl;
                continue;
            }

            const double ready_target_up = p0.up + ready_relative_h;

            std::vector<WGS84Point> ready_zone_wgs;
            ready_zone_wgs.reserve(this->input_data_.ready_zone.polygon.size());
            for (const auto &pt : this->input_data_.ready_zone.polygon) {
                ready_zone_wgs.push_back({pt.lon, pt.lat, ready_target_up});
            }
            if (ready_zone_wgs.size() < 3) {
                std::cerr << "ready_zone invalid (<3 points), skip ready_id=" << rid << std::endl;
                continue;
            }

            std::vector<ENUPoint> ready_zone_enu = this->wgs84ToENU_Batch(ready_zone_wgs, this->origin_);

            std::vector<ENUPoint> ready_patrol = this->computePatrolPathByMode(
                ready_zone_enu,
                distance,
                this->config_.path_planning.patrol_mode,
                ctx_enu);

            if (ready_patrol.empty()) {
                std::cerr << "ready_id=" << rid << " failed to generate ready patrol path." << std::endl;
                continue;
            }

            for (auto &pt : ready_patrol) {
                pt.up = ready_target_up;
            }

            if (!std::isfinite(heading0) || ctx_enu.size() < 2) {
                ENUPoint p1 = ready_patrol.front();
                const double dx = p1.east - p0.east;
                const double dy = p1.north - p0.north;
                if (std::hypot(dx, dy) > 1e-6) {
                    heading0 = std::atan2(dy, dx);
                } else {
                    heading0 = 0.0;
                }
            }

            const double radius = std::max(0.0, this->input_data_.min_turning_radius);
            const double resolution = (distance > 0.0) ? distance : 300.0;
            std::vector<ENUPoint> ready_transition;
            std::vector<ENUPoint> rotated_ready_patrol;
            this->buildTransitionAndRotatePatrol(p0, heading0, radius, resolution,
                                                 ready_patrol, ready_transition, rotated_ready_patrol);

            if (ready_transition.empty()) {
                std::cerr << "ready_id=" << rid << " failed to generate ready transition path." << std::endl;
                continue;
            }

            if (!rotated_ready_patrol.empty()) {
                this->enforceTransitionClimbRateAndBorrowPatrolPrefix(
                    ready_transition, rotated_ready_patrol,
                    std::string("uav ") + std::to_string(rid) + " ready plane2");
            }

            std::vector<WGS84Point> ready_transition_wgs = this->enuToWGS84_Batch(ready_transition, this->origin_);
            std::vector<WGS84Point> ready_patrol_wgs = this->enuToWGS84_Batch(
                rotated_ready_patrol.empty() ? ready_patrol : rotated_ready_patrol, this->origin_);

            this->upsertUsingMidwayLine(this->output_data_, rid, 2, ready_transition_wgs);
            this->upsertUsingMidwayLine(this->output_data_, rid, 3, ready_patrol_wgs);

            this->appendUavSegmentLine(this->output_data_.uav_plane2, rid, 2, ready_transition_wgs);
            this->appendUavSegmentLine(this->output_data_.uav_plane3, rid, 3, ready_patrol_wgs);
        }
    }
}

bool UavPathPlanner::getFollowerCurrentState(int uid, bool formation_enabled, double final_heading,
                                             ENUPoint &p0, double &heading0, std::vector<ENUPoint> &ctx_enu) {
    heading0 = 0.0;
    ctx_enu.clear();

    // 1) 优先使用 plane1（仅在编队模式且 plane1 存在时）
    if (formation_enabled && !this->output_data_.uav_plane1.empty()) {
        const UavTrajectoryLine *follower_line = nullptr;
        for (const auto &line : this->output_data_.uav_plane1) {
            if (line.uav_id == uid) { follower_line = &line; break; }
        }
        if (follower_line != nullptr && follower_line->points.size() >= 2) {
            std::vector<WGS84Point> wgs;
            wgs.reserve(follower_line->points.size());
            for (const auto &pt : follower_line->points) {
                wgs.push_back({pt.lon, pt.lat, pt.alt});
            }
            ctx_enu = this->wgs84ToENU_Batch(wgs, this->origin_);
            if (ctx_enu.size() >= 2) {
                p0 = ctx_enu.back();
                heading0 = computeTailHeadingRobust(ctx_enu, final_heading);
                return true;
            }
        }
    }

    // 2) 回退：使用起点
    WGS84Point start_wgs;
    if (!this->getFollowerStartWgs84(uid, start_wgs)) return false;
    this->adjustFollowerStartAltitudeIfNeeded(start_wgs, formation_enabled);
    p0 = this->wgs84ToENU(start_wgs, this->origin_);
    ctx_enu.push_back(p0);
    heading0 = 0.0;
    return true;
}
//规划函数 用于总调用
bool UavPathPlanner::getPlan(json &input_json, json &output_json, bool use3D, std::string algorithm)
{
    std::vector<ENUPoint> Enu_waypoint;
    int midwaypoint_num = 0;
    int zhandoupoint_num = 0;
    double distance;
    // InputData input_data; // Removed local variable
    if (!loadData(this->input_data_, input_json)) // Use member variable
    {
        std::cerr << "Failed to load intput json data." << std::endl;
    }
    else
    {
        std::cerr << "Successfully load intput json data." << std::endl;
    }

    // 输出结果统一写入 output_data_，最后再一次性转换为 output_json
    this->output_data_ = OutputData{};
    // using_midway_lines 采用“输入历史 + 本次新增”的迭代模式（仅使用 loadData 已解析数据）
    this->output_data_.using_midway_lines = this->input_data_.using_midway_lines;

    const bool formation_enabled = (this->input_data_.formation_using == 1);

    // 长机起点（非编队时用于 origin_ 与长机 plane2/3 起点复用；高度优化逻辑只在非编队时修改）
    WGS84Point leader_start_wgs{this->input_data_.uav_leader_start_point_wgs84.lon,
                                this->input_data_.uav_leader_start_point_wgs84.lat,
                                this->input_data_.uav_leader_start_point_wgs84.alt};

    // ---------- 高度优化开关 / 高程数据加载（非编队也可能依赖地形查询用于起点抬升与后续 plane2/3 优化） ----------
    const auto altitude_ctx = this->prepareAltitudeOptimizationContext(
        /*ensure_elevation_loaded=*/true,
        /*print_config_log=*/false);
    const bool altitude_opt_enabled = altitude_ctx.enabled;
    const std::string elevation_file = altitude_ctx.elevation_file;

    // ---------- 原点选择 ----------
    // 统一使用长机起点作为 origin_（仅经纬度参与 ENU 原点，alt 置 0）
    // 非编队时仍保留现有起点高度优化，用于后续 plane2/3 起点与离地净空约束。
    if (!formation_enabled) {
        // 若输入未提供起点高度，优先用 leader_midway 的首点高度作为参考高度（仅用于高度基准，不参与路径）
        if ((!std::isfinite(leader_start_wgs.alt) || std::abs(leader_start_wgs.alt) < 1e-6) &&
            !this->input_data_.leader_midway_point_wgs84.empty() &&
            std::isfinite(this->input_data_.leader_midway_point_wgs84.front().alt) &&
            this->input_data_.leader_midway_point_wgs84.front().alt > 0.0) {
            leader_start_wgs.alt = this->input_data_.leader_midway_point_wgs84.front().alt;
        }
        // 如果高程可用，保证起点不贴地（至少 safe_distance）
        if (elev_cost_map_ && elev_cost_map_->isElevationValid()) {
            AltitudeParams params = this->makeAltitudeParams();
            double min_clearance = params.safe_distance;
            if (!(min_clearance > 0.0) && (params.uav_R > 0.0)) min_clearance = params.uav_R;
            if (min_clearance > 0.0) {
                double elev = 0.0;
                if (elev_cost_map_->getElevationAt(leader_start_wgs.lon, leader_start_wgs.lat, elev)) {
                    const double min_alt = elev + min_clearance;
                    if (std::isfinite(min_alt) && (!std::isfinite(leader_start_wgs.alt) || leader_start_wgs.alt < min_alt)) {
                        leader_start_wgs.alt = min_alt;
                    }
                }
            }
        }
    }
    this->origin_ = leader_start_wgs;
    this->origin_.alt = 0.0;

    // ---------- 规划点准备（仅编队模式才需要 leader_midway_point_wgs84） ----------
    if (formation_enabled) {
        Enu_waypoint = this->preparePlanningWaypoints(midwaypoint_num, zhandoupoint_num);
    } else {
        midwaypoint_num = 0;
        zhandoupoint_num = static_cast<int>(this->input_data_.high_zhandou_point_wgs84.size());
        std::cout << "formation_using!=1: skip leader_midway_point_wgs84; high_zhandou_point_wgs84 count: "
                  << zhandoupoint_num << std::endl;
    }

    distance = this->input_data_.distance_points;
    if (!(distance > 0.0)) {
        // 最终兜底默认值
        distance = 300.0;
        std::cout << "Distance_Points not provided; using default distance_points=" << distance << std::endl;
    }

    std::vector<WGS84Point> Trajectory_WGS84;

    // ---------- 长机 plane1（仅编队模式才生成/优化） ----------
    if (formation_enabled) {

        // 提取用于轨迹生成的点（排除末尾的 zhandou points）
        std::vector<ENUPoint> planning_waypoints;
        if (Enu_waypoint.size() >= static_cast<size_t>(zhandoupoint_num)) {
            planning_waypoints.assign(Enu_waypoint.begin(), Enu_waypoint.end() - zhandoupoint_num);
        } else {
            planning_waypoints = Enu_waypoint;
        }

        // 避障处理  todo : 所有的轨迹都需要避障处理，而不仅仅是长机 plane1
        if (this->input_data_.has_prohibited_zone) {
            std::cout << "Applying prohibited zone avoidance..." << std::endl;
            planning_waypoints = this->avoidProhibitedZones(planning_waypoints);
        }

        if (algorithm == "minimum_snap") {
            if (use3D) {
                Trajectory_ENU = this->Minisnap_3D(planning_waypoints, distance, this->input_data_.leader_speed);
            } else {
                Trajectory_ENU = this->Minisnap_EN(planning_waypoints, distance, this->input_data_.leader_speed);
            }
        } else if (algorithm == "bspline") {
            std::cerr << "bspline algorithm not implemented yet." << std::endl;
            return false;
        } else if (algorithm == "bezier") {
            Trajectory_ENU = this->Bezier_3D(planning_waypoints, distance, this->input_data_.leader_speed, this->input_data_.min_turning_radius);
        } else {
            std::cerr << "Unknown algorithm: " << algorithm << ". Please use one of: minimum_snap, bspline, bezier" << std::endl;
            return false;
        }

        // ENU -> WGS84，并写入长机 plane1
        Trajectory_WGS84 = this->enuToWGS84_Batch(Trajectory_ENU, this->origin_);
        this->writeLeaderPlane1(this->output_data_, Trajectory_WGS84);
    } else {
        // 非编队：长机 plane1 必须为空
        this->Trajectory_ENU.clear();
        this->output_data_.uav_leader_plane1.clear();
    }

    this->prepareAltitudeOptimizationContext(
        /*ensure_elevation_loaded=*/false,
        /*print_config_log=*/true);

    // 仅在编队模式对长机 plane1 做高度优化（非编队不走 leader_midway/plane1）
    if (formation_enabled) {
        if (altitude_opt_enabled) {
            if (!elevation_file.empty()) {
                std::cout << "Running Altitude Optimization with file: " << elevation_file << std::endl;
                auto start_time = std::chrono::high_resolution_clock::now();
                if (!this->runAltitudeOptimization(elevation_file, this->output_data_)) {
                    std::cerr << "Failed to Altitude Optimization!" << std::endl;
                }
                auto end_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end_time - start_time;
                std::cout << "Altitude optimization time: " << elapsed.count() << "s" << std::endl;
            } else {
                std::cerr << "Altitude optimization enabled but no elevation file specified." << std::endl;
            }
        } else {
            std::cout << "Altitude optimization skipped." << std::endl;
        }
    }

    if (formation_enabled) {
        // 高度优化可能更新了 output_data_.uav_leader_plane1；这里仅同步一份 WGS84Point 供后续生成僚机 plane1 使用。
        Trajectory_WGS84.clear();
        Trajectory_WGS84.reserve(this->output_data_.uav_leader_plane1.size());
        for (const auto &p : this->output_data_.uav_leader_plane1) {
            Trajectory_WGS84.push_back({p.lon, p.lat, p.alt});
        }
    }

    // 计算并输出最小转弯半径（仅对 plane1 轨迹有意义）
    if (formation_enabled) {
        double min_radius = this->calculateMinTurningRadius(Trajectory_ENU);
        if (min_radius > 0) {
            std::cout << "Minimum turning radius: " << min_radius << " m" << std::endl;
        } else {
            std::cout << "Minimum turning radius: N/A (Straight line or insufficient points)" << std::endl;
        }
    }

    // 计算 Trajectory_ENU 最终点的朝向（仅编队会用到）
    double final_heading = 0.0;
    if (formation_enabled) {
        if (Trajectory_ENU.size() >= 2) {
            final_heading = computeTailHeadingRobust(Trajectory_ENU, final_heading);
        }
    }
    // std::cout << "Final heading: " << final_heading << " rad (" << rad2deg(final_heading) << " deg)" << std::endl;

    // 生成并写入僚机 plane1（编队才有）
    this->generateFollowerPlane1(Trajectory_ENU, Trajectory_WGS84);

    std::vector<ENUPoint> Patrol_Path; // 提升作用域以供第二段轨迹计算使用

    // ---------- 非编队：长机直接去自己的战斗区域（high_zhandou_point_wgs84） ----------
    if (!formation_enabled) {
        // 非编队：长机直接去自己的战斗区域（high_zhandou_point_wgs84）
        this->generateLeaderPlane2Plane3NonFormation(leader_start_wgs, distance);
    }

    /*                   第三段段任务区域巡逻轨迹计算（仅编队模式）               */
    if (formation_enabled && zhandoupoint_num != 0)   //存在高战斗区域点
    {
        std::cout<<"开始计算第三段任务区域巡逻轨迹"<<std::endl;
        if (Enu_waypoint.size() < static_cast<size_t>(zhandoupoint_num)) {
            std::cerr << "plane3 leader patrol failed: Enu_waypoint.size()=" << Enu_waypoint.size()
                      << " < zhandoupoint_num=" << zhandoupoint_num << std::endl;
        } else {
            std::vector<ENUPoint> zhandou_zone_enu(Enu_waypoint.end() - zhandoupoint_num, Enu_waypoint.end());

            // 长机战斗区域目标高度：和僚机一致，基于 plane1 末点高度抬高 leader_fly_high
            double base_up = 0.0;
            if (!Trajectory_ENU.empty()) {
                base_up = Trajectory_ENU.back().up;
            } else if (midwaypoint_num > 0 && Enu_waypoint.size() >= static_cast<size_t>(midwaypoint_num)) {
                base_up = Enu_waypoint[static_cast<size_t>(midwaypoint_num - 1)].up;
            }
            const double zhandou_target_up = base_up + this->input_data_.leader_fly_high;
            for (auto &pt : zhandou_zone_enu) {
                pt.up = zhandou_target_up;
            }

            Patrol_Path = this->computePatrolPathByMode(
                zhandou_zone_enu,
                distance,
                this->config_.path_planning.patrol_mode,
                Trajectory_ENU
            );

            // 第三段巡逻高度统一赋值为战斗区域目标高度
            for (auto &pt : Patrol_Path) {
                pt.up = zhandou_target_up;
            }
        }

        if (!Patrol_Path.empty()) {
            std::vector<WGS84Point> Patrol_Path_WGS84 = this->enuToWGS84_Batch(Patrol_Path, this->origin_);
            this->writeLeaderSegment(this->output_data_.uav_leader_plane3, this->output_data_, 3,
                                     Patrol_Path_WGS84, /*sync_using_midway_line=*/Trajectory_ENU.empty());
        } else {
            std::cerr << "Warning: failed to generate patrol path in plane3." << std::endl;
        }
    }


    /*                   第二段任务区域过渡轨迹计算      zhandoupoint_num           */
    if(formation_enabled && zhandoupoint_num != 0 && !Trajectory_ENU.empty() && !Patrol_Path.empty())   //存在高战斗区域点且有前序轨迹
    {
        std::cout<<"开始计算第二段任务区域过渡轨迹 (切圆切入优化)"<<std::endl;
        ENUPoint p0 = Trajectory_ENU.back(); // 起点 ENU
        double heading0 = final_heading; // 起点朝向

        //僚机和长机巡逻区域不一样,所以注释掉
        // // 计算并更新跟随者轨迹
        // json plane_array = this->generateFollowerTrajectories(input_json, input_data, Trajectory_ENU, Trajectory_WGS84);
        // output_json["uav_plane1"] = plane_array;

        // 计算第二段过渡轨迹（切圆切入优化）并更新第三段巡逻轨迹
        computeTransitionAndRotatePatrol(p0, heading0, this->input_data_.min_turning_radius, distance, Patrol_Path, this->output_data_);
    }

    // 最后再对第二段和第三段做高度优化
    if (altitude_opt_enabled && elev_cost_map_ && elev_cost_map_->isElevationValid()) {
        std::vector<std::string> keys = {"uav_leader_plane2", "uav_leader_plane3"};
        std::vector<int> seg_ids = {2, 3};
        // 联合优化第二段和第三段，其中第三段（index 1）需要保持等高点高度相同
        this->optimizeAndApplyJointSegments(this->output_data_, keys, seg_ids, 1);
    }

    // 僚机 plane2/3（战斗区/准备区）生成：内部会完成目的地分配并写入 output_data_.uav_plane2/uav_plane3，
    // 同时返回“最终 ready_id”（包含 battle_zone 不满足而改派的无人机）。
    std::vector<int> final_ready_ids;
    this->generateFollowerPlane2Plane3(formation_enabled, final_heading, distance, final_ready_ids);

    // ready_id
    // 输出时回填“最终 ready_id”（包含 battle_zone 不满足而改派的无人机）
    this->output_data_.ready_id = final_ready_ids;

    // leader_show_points：编队模式显示 leader_midway_point_wgs84 + 战斗区点；非编队只显示战斗区点
    this->output_data_.leader_show_points.clear();
    if (formation_enabled) {
        this->output_data_.leader_show_points.reserve(this->input_data_.leader_midway_point_wgs84.size() + this->input_data_.high_zhandou_point_wgs84.size());
        for (const auto &pt : this->input_data_.leader_midway_point_wgs84) {
            this->output_data_.leader_show_points.emplace_back(WGS84Coord(pt.lon, pt.lat, pt.alt));
        }
        double leader_show_last_alt = 0.0;
        if (!this->output_data_.uav_leader_plane1.empty()) {
            leader_show_last_alt = this->output_data_.uav_leader_plane1.back().alt;
        } else if (!this->input_data_.leader_midway_point_wgs84.empty()) {
            leader_show_last_alt = this->input_data_.leader_midway_point_wgs84.back().alt;
        }
        for (const auto &pt : this->input_data_.high_zhandou_point_wgs84) {
            this->output_data_.leader_show_points.emplace_back(
                WGS84Coord(pt.lon, pt.lat, leader_show_last_alt + this->input_data_.leader_fly_high));
        }
    } else {
        // 非编队：战斗区点高度用 plane3 目标高度（若可用），否则用起点高度 + leader_fly_high
        double zhandou_alt = 0.0;
        if (!this->output_data_.uav_leader_plane3.empty()) {
            zhandou_alt = this->output_data_.uav_leader_plane3.front().alt;
        } else {
            double base_alt = this->input_data_.uav_leader_start_point_wgs84.alt;
            if ((!std::isfinite(base_alt) || std::abs(base_alt) < 1e-6) && !this->input_data_.leader_midway_point_wgs84.empty()) {
                base_alt = this->input_data_.leader_midway_point_wgs84.front().alt;
            }
            zhandou_alt = base_alt + this->input_data_.leader_fly_high;
        }
        this->output_data_.leader_show_points.reserve(this->input_data_.high_zhandou_point_wgs84.size());
        for (const auto &pt : this->input_data_.high_zhandou_point_wgs84) {
            this->output_data_.leader_show_points.emplace_back(WGS84Coord(pt.lon, pt.lat, zhandou_alt));
        }
    }

    {
        std::cout << "Constructed using_uav_list from input: [";
        for (size_t i = 0; i < this->output_data_.using_uav_list.size(); ++i) {
            std::cout << this->output_data_.using_uav_list[i];
            if (i + 1 < this->output_data_.using_uav_list.size()) std::cout << ",";
        }
        std::cout << "]" << std::endl;
    }

    // check_change：使用结构化数据接口（InputData/OutputData）
    if (!this->check_change(this->input_data_, this->output_data_)) {
        std::cerr << "check_change failed." << std::endl;
    } else {
        std::cerr << "check_change succeeded." << std::endl;
    }
/*   -------------------规划结果的辅助输出，主要给前端/可视化/标注用，不影响路径求解本身-----------------------   */
    //输出结果增加 midway_point_num（编队结果里才有）
    // 含义：把你输入的 leader_midway_point_wgs84（长机中途点）映射到最终生成的航迹点序列上，记录“每个中途点在航迹里的下标”。
    // 结构：list[int]，长度≈leader_midway_point_wgs84 的长度；找不到就填 -1。
    // 例如"midway_point_num": [
    //         6,
    //         44,
    //         57,
    //         92,
    //         103
    //     ]
    if (formation_enabled) {
        json work_json;
        this->outputDataToJson(this->output_data_, work_json);
        json midway = buildMidwayPointNum(this->input_data_, work_json);
        this->output_data_.midway_point_num.clear();
        if (midway.is_array()) {
            for (const auto &v : midway) {
                if (v.is_number_integer()) {
                    this->output_data_.midway_point_num.push_back(v.get<int>());
                }
            }
        }
    } else {
        this->output_data_.midway_point_num.clear();
    }

    // 最终统一转换为 json 输出
    this->outputDataToJson(this->output_data_, output_json);
    return true;
}

//生成僚机编队轨迹
json UavPathPlanner::generateFollowerTrajectories(const InputData &input_data,
                                  const std::vector<ENUPoint> &Trajectory_ENU,
                                  const std::vector<WGS84Point> &Trajectory_WGS84) {
    (void)Trajectory_WGS84;
    json plane_array = json::array();

    if (input_data.formation_using != 1) {
        return plane_array;
    }

    if (input_data.uavs_id.empty() || input_data.uav_start_point_wgs84.empty()) {
        return plane_array; // empty
    }

    json uavs_ids = json::array();
    for (int id : input_data.uavs_id) {
        uavs_ids.push_back(id);
    }

    json uav_starts = json::array();
    for (const auto &p : input_data.uav_start_point_wgs84) {
        json pt = json::array();
        pt.push_back(p.lon);
        pt.push_back(p.lat);
        pt.push_back(p.alt);
        uav_starts.push_back(pt);
    }

    // 计算 leader 起始 ENU（使用相同的 origin）
    WGS84Point leader_start_wgs{input_data.uav_leader_start_point_wgs84.lon,
                                input_data.uav_leader_start_point_wgs84.lat,
                                input_data.uav_leader_start_point_wgs84.alt};
    ENUPoint leader_start_enu = this->wgs84ToENU(leader_start_wgs, this->origin_);

    // 计算 leader 初始航向（使用采样轨迹的前两点）
    double leader_initial_heading = 0.0;
    if (Trajectory_ENU.size() >= 2) {
        double dx = Trajectory_ENU[1].east - Trajectory_ENU[0].east;
        double dy = Trajectory_ENU[1].north - Trajectory_ENU[0].north;
        leader_initial_heading = atan2(dy, dx);
    }
    double cos0 = cos(leader_initial_heading);
    double sin0 = sin(leader_initial_heading);
    Eigen::Matrix2d R0;
    R0 << cos0, -sin0,
          sin0,  cos0;

    size_t N = Trajectory_ENU.size();
    std::vector<Eigen::Vector2d> leader_xy; leader_xy.reserve(N);
    for (const auto &p : Trajectory_ENU) leader_xy.emplace_back(p.east, p.north);

    // 预计算并平滑航向角，避免转弯时编队变形或冲突
    std::vector<double> leader_headings(N);
    for (size_t t = 0; t < N; ++t) {
        double dx = 0.0, dy = 0.0;
        if (t == 0) {
            if (N > 1) {
                dx = leader_xy[1].x() - leader_xy[0].x();
                dy = leader_xy[1].y() - leader_xy[0].y();
            } else {
                dx = cos(leader_initial_heading);
                dy = sin(leader_initial_heading);
            }
        } else if (t == N - 1) {
            dx = leader_xy[N - 1].x() - leader_xy[N - 2].x();
            dy = leader_xy[N - 1].y() - leader_xy[N - 2].y();
        } else {
            // 使用中心差分计算切线方向，比前向差分更平滑
            dx = leader_xy[t + 1].x() - leader_xy[t - 1].x();
            dy = leader_xy[t + 1].y() - leader_xy[t - 1].y();
        }
        leader_headings[t] = atan2(dy, dx);
    }

    // 对航向角进行滑动窗口平滑
    if (N > 5) {
        std::vector<double> smoothed = leader_headings;
        int window = 10; // 窗口大小可调整
        for (size_t t = 0; t < N; ++t) {
            double sum_sin = 0.0, sum_cos = 0.0;
            int count = 0;
            for (int k = -window; k <= window; ++k) {
                int idx = (int)t + k;
                if (idx >= 0 && idx < (int)N) {
                    sum_sin += sin(leader_headings[idx]);
                    sum_cos += cos(leader_headings[idx]);
                    count++;
                }
            }
            if (count > 0) {
                smoothed[t] = atan2(sum_sin, sum_cos);
            }
        }
        leader_headings = smoothed;
    }

    // 编队距离：改用 formation_distance
    double formation_distance = 50.0; // 默认值（若 config.yaml/输入未提供，则沿用旧默认）

    // 竖一字形编队：每一列最多支持的行数（无人机数量）
    int uav_formation_max_row = 8; // 默认值，最小为 1

    // 编队距离下限约束：如果 formation_distance < (2*position_misalignment + uav_R)*1.41421 则自动扩大
    double position_misalignment = 0.0; // 无人机定位误差（米）
    double uav_R = 2.0;                 // 机体半径/安全包络（米）

    // 使用构造时加载的 config.yaml 参数（不在这里重复读取）
    formation_distance = this->config_.path_planning.formation_distance;
    position_misalignment = this->config_.path_planning.position_misalignment;
    uav_formation_max_row = this->config_.path_planning.uav_formation_max_row;
    uav_R = this->config_.altitude_optimization.uav_R;

    // 允许从输入 JSON 覆盖
    if (input_data.formation_distance > 0.0) formation_distance = input_data.formation_distance;
    if (input_data.position_misalignment >= 0.0) position_misalignment = input_data.position_misalignment;
    if (input_data.uav_R > 0.0) uav_R = input_data.uav_R;
    if (input_data.uav_formation_max_row > 0) uav_formation_max_row = input_data.uav_formation_max_row;

    if (uav_formation_max_row < 1) uav_formation_max_row = 1;

    const double min_formation_distance = (2.0 * position_misalignment + uav_R) * 1.41421;
    if (formation_distance < min_formation_distance) {
        std::cout << "formation_distance is too small (" << formation_distance << "), clamped to min "
                  << min_formation_distance << " (position_misalignment=" << position_misalignment
                  << ", uav_R=" << uav_R << ")" << std::endl;
        formation_distance = min_formation_distance;
    }

    // std::cout << "Formation Model: " << input_data.formation_model << std::endl;

    if (input_data.formation_model == 1) {
        std::cout << "Formation Model: 人字形编队" << std::endl;
        plane_array = generateVShapeTrajectories(uavs_ids, uav_starts, leader_start_enu, R0, leader_xy, leader_headings, Trajectory_ENU, formation_distance);
    } else if (input_data.formation_model == 2) {
        std::cout << "Formation Model: 一字形编队(横)" << std::endl;
        plane_array = generateLineShapeTrajectories(uavs_ids, uav_starts, leader_start_enu, R0, leader_xy, leader_headings, Trajectory_ENU, formation_distance);
    } else if (input_data.formation_model == 3) {
        std::cout << "Formation Model: 一字形编队(竖)" << std::endl;
        plane_array = generateVerticalLineShapeTrajectories(uavs_ids, uav_starts, leader_start_enu, R0, leader_xy, leader_headings, Trajectory_ENU, uav_formation_max_row, formation_distance);
    } else if (input_data.formation_model == 4) {
        std::cout << "Formation Model: 三角形编队" << std::endl;
        plane_array = generateTriangleShapeTrajectories(uavs_ids, uav_starts, leader_start_enu, R0, leader_xy, leader_headings, Trajectory_ENU, formation_distance);
    } else {
        // 默认使用人字形 (V-shape)
        std::cout << "Formation Model: 人字形编队" << std::endl;
        plane_array = generateVShapeTrajectories(uavs_ids, uav_starts, leader_start_enu, R0, leader_xy, leader_headings, Trajectory_ENU, formation_distance);
    }

    return plane_array;
}
//生成人字形编队轨迹
json UavPathPlanner::generateVShapeTrajectories(const json &uavs_ids, const json &uav_starts, 
                                                const ENUPoint &leader_start_enu, const Eigen::Matrix2d &R0,
                                                const std::vector<Eigen::Vector2d> &leader_xy, 
                                                const std::vector<double> &leader_headings,
                                                const std::vector<ENUPoint> &Trajectory_ENU,
                                                double safety_distance) {
    json plane_array = json::array();
    size_t N = Trajectory_ENU.size();

    for (size_t idx = 0; idx < uavs_ids.size(); ++idx) {
        int uid = uavs_ids[idx].get<int>();
        if (idx >= uav_starts.size()) break;

        // Get start point
        auto s = uav_starts[idx];
        WGS84Point start_wp;
        start_wp.lon = s[0].get<double>();
        start_wp.lat = s[1].get<double>();
        start_wp.alt = s.size() >= 3 ? s[2].get<double>() : 0.0;
        
        // Hardcoded V-Shape Geometry
        // Index 0: Left, Index 1: Right, Index 2: Left, ...
        int row = (idx / 2) + 1;
        int side = (idx % 2 == 0) ? 1 : -1; // 1 for Left, -1 for Right
        
        // 45 degree V-shape: dx = dy
        double dx = -row * safety_distance; // Behind
        double dy = side * row * safety_distance; // Lateral
        
        Eigen::Vector2d rel_body(dx, dy);
        double rel_up = 0.0; // Same altitude

        // 构造 follower 的 entry
        json follower_entry = json::array();
        follower_entry.push_back(uid);

        for (size_t t = 0; t < N; ++t) {
            // if (t == 0) {
            //     json pa = json::array();
            //     pa.push_back(start_wp.lon);
            //     pa.push_back(start_wp.lat);
            //     pa.push_back(start_wp.alt);
            //     follower_entry.push_back(pa);
            //     continue;
            // }

            double heading = leader_headings[t];
            double c = cos(heading);
            double s_ = sin(heading);
            Eigen::Matrix2d Rt; Rt << c, -s_,
                                      s_,  c;
            Eigen::Vector2d offset_global = Rt * rel_body;

            ENUPoint fp;
            fp.east = leader_xy[t].x() + offset_global.x();
            fp.north = leader_xy[t].y() + offset_global.y();
            fp.up = Trajectory_ENU[t].up + rel_up;

            WGS84Point wp = this->enuToWGS84(fp, this->origin_);
            json pa = json::array();
            pa.push_back(wp.lon);
            pa.push_back(wp.lat);
            pa.push_back(wp.alt);
            follower_entry.push_back(pa);
        }

        plane_array.push_back(follower_entry);
    }
    return plane_array;
}
//生成一字形编队轨迹
json UavPathPlanner::generateLineShapeTrajectories(const json &uavs_ids, const json &uav_starts, 
                                                   const ENUPoint &leader_start_enu, const Eigen::Matrix2d &R0,
                                                   const std::vector<Eigen::Vector2d> &leader_xy, 
                                                   const std::vector<double> &leader_headings,
                                                   const std::vector<ENUPoint> &Trajectory_ENU,
                                                   double safety_distance) {
    json plane_array = json::array();
    size_t N = Trajectory_ENU.size();
    double leader_start_alt = (N > 0) ? Trajectory_ENU.front().up : 0.0;

    for (size_t idx = 0; idx < uavs_ids.size(); ++idx) {
        int uid = uavs_ids[idx].get<int>();
        if (idx >= uav_starts.size()) break;

        // Get start point
        auto s = uav_starts[idx];
        WGS84Point start_wp;
        start_wp.lon = s[0].get<double>();
        start_wp.lat = s[1].get<double>();
        start_wp.alt = s.size() >= 3 ? s[2].get<double>() : 0.0;
        
        // Hardcoded Line-Shape (Abreast) Geometry
        // Left and Right of leader.
        int row = (idx / 2) + 1;
        int side = (idx % 2 == 0) ? 1 : -1; // 1 for Left, -1 for Right
        
        double dx = 0.0; // Abreast
        double dy = side * row * safety_distance; // Lateral
        
        Eigen::Vector2d rel_body(dx, dy);
        double rel_up = 0.0;

        // 构造 follower 的 entry
        json follower_entry = json::array();
        follower_entry.push_back(uid);

        for (size_t t = 0; t < N; ++t) {
            if (t == 0) {
                json pa = json::array();
                pa.push_back(start_wp.lon);
                pa.push_back(start_wp.lat);
                pa.push_back(leader_start_alt);
                follower_entry.push_back(pa);
                continue;
            }

            double heading = leader_headings[t];
            double c = cos(heading);
            double s_ = sin(heading);
            Eigen::Matrix2d Rt; Rt << c, -s_,
                                      s_,  c;
            Eigen::Vector2d offset_global = Rt * rel_body;

            ENUPoint fp;
            fp.east = leader_xy[t].x() + offset_global.x();
            fp.north = leader_xy[t].y() + offset_global.y();
            fp.up = Trajectory_ENU[t].up + rel_up;

            WGS84Point wp = this->enuToWGS84(fp, this->origin_);
            json pa = json::array();
            pa.push_back(wp.lon);
            pa.push_back(wp.lat);
            pa.push_back(wp.alt);
            follower_entry.push_back(pa);
        }

        plane_array.push_back(follower_entry);
    }
    return plane_array;
}

//生成一字形编队轨迹(竖)
json UavPathPlanner::generateVerticalLineShapeTrajectories(const json &uavs_ids, const json &uav_starts,
                                                           const ENUPoint &leader_start_enu, const Eigen::Matrix2d &R0,
                                                           const std::vector<Eigen::Vector2d> &leader_xy,
                                                           const std::vector<double> &leader_headings,
                                                           const std::vector<ENUPoint> &Trajectory_ENU,
                                                           int uav_formation_max_row,
                                                           double safety_distance)
{
    (void)leader_start_enu;
    (void)R0;

    if (uav_formation_max_row < 1) uav_formation_max_row = 1;

    json plane_array = json::array();
    size_t N = Trajectory_ENU.size();
    double leader_start_alt = (N > 0) ? Trajectory_ENU.front().up : 0.0;

    for (size_t idx = 0; idx < uavs_ids.size(); ++idx) {
        int uid = uavs_ids[idx].get<int>();
        if (idx >= uav_starts.size()) break;

        // Get start point
        auto s = uav_starts[idx];
        WGS84Point start_wp;
        start_wp.lon = s[0].get<double>();
        start_wp.lat = s[1].get<double>();
        start_wp.alt = s.size() >= 3 ? s[2].get<double>() : 0.0;

        // Vertical line (trail) geometry: behind leader along body -x
        // 每一列最多 uav_formation_max_row 架；超出则换到下一列。
        int col = static_cast<int>(idx) / uav_formation_max_row;      // 0,1,2...
        int row_in_col = static_cast<int>(idx) % uav_formation_max_row; // 0..max_row-1

        // rows start from 1 step behind leader
        double dx = -(row_in_col + 1) * safety_distance;

        // columns expand laterally (body +y is left): 0, +1, -1, +2, -2, ...
        double dy = 0.0;
        if (col > 0) {
            int side = (col % 2 == 1) ? 1 : -1;          // odd -> left, even -> right
            int level = (col + 1) / 2;                   // 1,1,2,2,...
            dy = side * level * safety_distance;
        }

        Eigen::Vector2d rel_body(dx, dy);
        double rel_up = 0.0;

        json follower_entry = json::array();
        follower_entry.push_back(uid);

        for (size_t t = 0; t < N; ++t) {
            if (t == 0) {
                json pa = json::array();
                pa.push_back(start_wp.lon);
                pa.push_back(start_wp.lat);
                pa.push_back(leader_start_alt);
                follower_entry.push_back(pa);
                continue;
            }

            double heading = leader_headings[t];
            double c = cos(heading);
            double s_ = sin(heading);
            Eigen::Matrix2d Rt;
            Rt << c, -s_,
                  s_,  c;
            Eigen::Vector2d offset_global = Rt * rel_body;

            ENUPoint fp;
            fp.east = leader_xy[t].x() + offset_global.x();
            fp.north = leader_xy[t].y() + offset_global.y();
            fp.up = Trajectory_ENU[t].up + rel_up;

            WGS84Point wp = this->enuToWGS84(fp, this->origin_);
            json pa = json::array();
            pa.push_back(wp.lon);
            pa.push_back(wp.lat);
            pa.push_back(wp.alt);
            follower_entry.push_back(pa);
        }

        plane_array.push_back(follower_entry);
    }

    return plane_array;
}

//生成三角形编队轨迹
json UavPathPlanner::generateTriangleShapeTrajectories(const json &uavs_ids, const json &uav_starts,
                                                       const ENUPoint &leader_start_enu, const Eigen::Matrix2d &R0,
                                                       const std::vector<Eigen::Vector2d> &leader_xy,
                                                       const std::vector<double> &leader_headings,
                                                       const std::vector<ENUPoint> &Trajectory_ENU,
                                                       double safety_distance)
{
    (void)leader_start_enu;
    (void)R0;

    json plane_array = json::array();
    size_t N = Trajectory_ENU.size();
    double leader_start_alt = (N > 0) ? Trajectory_ENU.front().up : 0.0;

    for (size_t idx = 0; idx < uavs_ids.size(); ++idx) {
        int uid = uavs_ids[idx].get<int>();
        if (idx >= uav_starts.size()) break;

        // Get start point
        auto s = uav_starts[idx];
        WGS84Point start_wp;
        start_wp.lon = s[0].get<double>();
        start_wp.lat = s[1].get<double>();
        start_wp.alt = s.size() >= 3 ? s[2].get<double>() : 0.0;

        // Triangle (Delta) layers for followers:
        // row 1: 2 UAVs at (-d, ±d)
        // row 2: 3 UAVs at (-2d, -2d/0/+2d)
        // ... row r: (r+1) UAVs
        int k = static_cast<int>(idx) + 1; // 1-based follower index
        int row = 1;
        int prev_count = 0;
        while (prev_count + (row + 1) < k) {
            prev_count += (row + 1);
            row++;
        }
        int pos_in_row = k - prev_count - 1; // 0..row

        // Use body frame: +x forward, +y left
        double dx = -row * safety_distance;
        double center = (double)row / 2.0;
        // Make the first UAV on the left side (consistent with V-shape idx%2==0 => left)
        double dy = (center - pos_in_row) * 2.0 * safety_distance;

        Eigen::Vector2d rel_body(dx, dy);
        double rel_up = 0.0;

        json follower_entry = json::array();
        follower_entry.push_back(uid);

        for (size_t t = 0; t < N; ++t) {
            if (t == 0) {
                json pa = json::array();
                pa.push_back(start_wp.lon);
                pa.push_back(start_wp.lat);
                pa.push_back(leader_start_alt);
                follower_entry.push_back(pa);
                continue;
            }

            double heading = leader_headings[t];
            double c = cos(heading);
            double s_ = sin(heading);
            Eigen::Matrix2d Rt;
            Rt << c, -s_,
                  s_,  c;
            Eigen::Vector2d offset_global = Rt * rel_body;

            ENUPoint fp;
            fp.east = leader_xy[t].x() + offset_global.x();
            fp.north = leader_xy[t].y() + offset_global.y();
            fp.up = Trajectory_ENU[t].up + rel_up;

            WGS84Point wp = this->enuToWGS84(fp, this->origin_);
            json pa = json::array();
            pa.push_back(wp.lon);
            pa.push_back(wp.lat);
            pa.push_back(wp.alt);
            follower_entry.push_back(pa);
        }

        plane_array.push_back(follower_entry);
    }

    return plane_array;
}

// Minisnap_EN: 仅在平面上(东/北)进行轨迹优化，忽略高度(将高度固定为起点高度)
std::vector<ENUPoint> UavPathPlanner::Minisnap_EN (std::vector<ENUPoint> Enu_waypoint_, double distance_, double v_avg_override)
{
    std::vector<ENUPoint> result{};
    int dot_num = Enu_waypoint_.size();  //路径点个数
    if (dot_num < 2) return result;

    // 构造仅包含 XY 的路线（Z 设为 0），让最小 snap 规划器只优化平面轨迹
    Eigen::MatrixXd route(dot_num, 3);
    for (int i = 0; i < dot_num; i++) {
        route(i, 0) = Enu_waypoint_[i].east;   // x -> east
        route(i, 1) = Enu_waypoint_[i].north;  // y -> north
        route(i, 2) = 0.0;                     // z -> 0 (忽略高度)
    }

    // use global generator (loaded from path_planning.minimum_snap_config_file)
    if (distance_ > 0) {
        std::cerr << "Using input JSON distance_points as sampling distance: " << distance_ << " m" << std::endl;
    }
    if (v_avg_override > 0) {
        std::cerr << "Using leader_speed as average speed: " << v_avg_override << " m/s" << std::endl;
    }

    Eigen::MatrixXd sampled = this->generator_.GenerateTrajectoryMatrix(route, this->config_.minimum_snap, distance_, v_avg_override);

    // 将采样的平面点转换回 ENUPoint，并把高度固定为起点高度（高度可保持不变或任意）
    double up_value = Enu_waypoint_[0].up; // 固定高度为起点高度
    for (int i = 0; i < sampled.rows(); ++i) {
        ENUPoint enu_point;
        enu_point.east = sampled(i, 0);
        enu_point.north = sampled(i, 1);
        enu_point.up = up_value;
        result.push_back(enu_point);
    }

    std::cout << "Generated 2D trajectory point number:" << result.size() << std::endl;
    return result;
}

// Minisnap_3D: 三维轨迹优化（保留高度），如果未找到专门配置同样使用默认 YAML 路径
std::vector<ENUPoint> UavPathPlanner::Minisnap_3D(std::vector<ENUPoint> Enu_waypoint_, double distance_, double v_avg_override)
{
    std::vector<ENUPoint> result{};
    int dot_num = Enu_waypoint_.size();  //路径点个数
    if (dot_num < 2) return result;

    Eigen::MatrixXd route(dot_num, 3);
    for (int i = 0; i < dot_num; i++) {
        route(i, 0) = Enu_waypoint_[i].east;   // x -> east
        route(i, 1) = Enu_waypoint_[i].north;  // y -> north  
        route(i, 2) = Enu_waypoint_[i].up;     // z -> up
    }

    // use global generator (loaded from path_planning.minimum_snap_config_file)
    if (distance_ > 0) {
        std::cerr << "Using input JSON distance_points as sampling distance: " << distance_ << " m" << std::endl;
    }
    if (v_avg_override > 0) {
        std::cerr << "Using leader_speed as average speed: " << v_avg_override << " m/s" << std::endl;
    }

    Eigen::MatrixXd sampled = this->generator_.GenerateTrajectoryMatrix(route, this->config_.minimum_snap, distance_, v_avg_override);

    // 将采样点转换为 ENUPoint 向量返回
    for (int i = 0; i < sampled.rows(); ++i) {
        ENUPoint enu_point;
        enu_point.east = sampled(i, 0);
        enu_point.north = sampled(i, 1);
        enu_point.up = sampled(i, 2);
        result.push_back(enu_point);
    }

    std::cout << "Generated 3D trajectory point number:" << result.size() << std::endl;
    return result;
}

// Bezier_3D: 贝塞尔曲线轨迹生成,未使用
std::vector<ENUPoint> UavPathPlanner::Bezier_3D(std::vector<ENUPoint> Enu_waypoint_, double distance_, double V_avg_override, double min_radius) //转弯半径无效
{
    std::vector<ENUPoint> result{};
    int dot_num = Enu_waypoint_.size();
    if (dot_num < 2) return result;

    Eigen::MatrixXd route(dot_num, 3);
    for (int i = 0; i < dot_num; i++) {
        route(i, 0) = Enu_waypoint_[i].east;
        route(i, 1) = Enu_waypoint_[i].north;
        route(i, 2) = Enu_waypoint_[i].up;
    }

    math_util::Bezier bezier;
    math_util::BezierConfig config;
    if (min_radius > 0) {
        config.min_radius = 300;
    }
    bezier.SetConfig(config);

    // yaml_path is not used by Bezier currently, pass empty string or dummy
    Eigen::MatrixXd sampled = bezier.GenerateTrajectoryMatrix(route, "", distance_, V_avg_override);

    for (int i = 0; i < sampled.rows(); ++i) {
        ENUPoint enu_point;
        enu_point.east = sampled(i, 0);
        enu_point.north = sampled(i, 1);
        enu_point.up = sampled(i, 2);
        result.push_back(enu_point);
    }

    std::cout << "Generated Bezier 3D trajectory point number:" << result.size() << std::endl;
    return result;
}
//将input json数据加载到input_data结构体中
bool UavPathPlanner::loadData(InputData &input_data, json &input_json)
{
    input_data = InputData{};

    // 路径采样间距（单位：米）
    bool has_json_distance = false;
    if (input_json.contains("distance_points") && !input_json["distance_points"].empty()) {
        try {
            if (input_json["distance_points"].is_number()) {
                input_data.distance_points = input_json["distance_points"].get<double>();
                has_json_distance = true;
            } else if (input_json["distance_points"].is_array() && !input_json["distance_points"].empty() && input_json["distance_points"][0].is_number()) {
                input_data.distance_points = input_json["distance_points"][0].get<double>();
                has_json_distance = true;
            }
        } catch (...) {
            has_json_distance = false;
        }
    }
    if (!has_json_distance) {
        input_data.distance_points = this->config_.path_planning.distance_points;
    }

    if (input_json.contains("leader_speed") && input_json["leader_speed"].is_number()) {
        input_data.leader_speed = input_json["leader_speed"].get<double>();
    }
    if (input_json.contains("leader_fly_high") && input_json["leader_fly_high"].is_number()) {
        input_data.leader_fly_high = input_json["leader_fly_high"].get<double>();
    }
    if (input_json.contains("formation_model") && input_json["formation_model"].is_number_integer()) {
        input_data.formation_model = input_json["formation_model"].get<int>();
    }
    if (input_json.contains("formation_using") && input_json["formation_using"].is_number_integer()) {
        input_data.formation_using = input_json["formation_using"].get<int>();
    }

    auto parse_wgs84_list = [&](const char *key, std::vector<WGS84Coord> &out, double default_alt = 0.0) {
        if (!input_json.contains(key) || !input_json[key].is_array()) return;
        for (const auto &item : input_json[key]) {
            WGS84Point p;
            if (!parseWGS84PointValue(item, p)) continue;
            if (!item.is_array() || item.size() < 3) p.alt = default_alt;
            out.emplace_back(WGS84Coord(p.lon, p.lat, p.alt));
        }
    };

    auto parse_int_list = [&](const char *key, std::vector<int> &out) {
        if (!input_json.contains(key)) return;
        const auto &v = input_json[key];
        if (v.is_array()) {
            for (const auto &x : v) {
                if (x.is_number_integer()) out.push_back(x.get<int>());
            }
        } else if (v.is_number_integer()) {
            out.push_back(v.get<int>());
        }
    };

    parse_wgs84_list("leader_midway_point_wgs84", input_data.leader_midway_point_wgs84, 0.0);
    parse_wgs84_list("high_zhandou_point_wgs84", input_data.high_zhandou_point_wgs84, 0.0);
    input_data.ready_zone.zone_id = 0;
    input_data.ready_zone.zone_type = "ready_zone";
    parse_wgs84_list("ready_zone", input_data.ready_zone.polygon, 0.0);
    if(input_json.contains("ready_high_list") && input_json["ready_high_list"].is_array() && input_json["ready_high_list"].size() >= 2) {
        input_data.ready_zone.height_range = {
            input_json["ready_high_list"][0].get<double>(),
            input_json["ready_high_list"][1].get<double>()
        };
    }

    parse_wgs84_list("uav_start_point_wgs84", input_data.uav_start_point_wgs84, 0.0);

    parse_int_list("uavs_id", input_data.uavs_id);
    parse_int_list("ready_id", input_data.ready_id);
    parse_int_list("uav_leader_id", input_data.uav_leader_ids);
    parse_int_list("using_uav_list", input_data.using_uav_list);

    // battle_zones 字段解析
    input_data.battle_zone_list.clear();
    if (input_json.contains("battle_zone_list") && input_json["battle_zone_list"].is_array()) {
        for (const auto &zone_id : input_json["battle_zone_list"]) {
            if (zone_id.is_number_integer()) {
                input_data.battle_zone_list.push_back(zone_id.get<int>());
            }
        }
    }
    
    std::vector<double> temp_battle_high_list;
    if (input_json.contains("battle_high_list") && input_json["battle_high_list"].is_array()) {
        for (const auto &h : input_json["battle_high_list"]) {
            if (h.is_number()) temp_battle_high_list.push_back(h.get<double>());
        }
    }
    
    std::vector<int> temp_battle_zone_link_flag;
    parse_int_list("battle_zone_link_flag", temp_battle_zone_link_flag);
    
    if (input_json.contains("battle_zone_wgs84") && input_json["battle_zone_wgs84"].is_array()) {
        int idx = 0;
        for (const auto &poly : input_json["battle_zone_wgs84"]) {
            if (!poly.is_array()) continue;
            FlightZone bz;
            bz.zone_id = idx;
            bz.zone_type = "battle_zone";
            for (const auto &pt : poly) {
                WGS84Point p;
                if (!parseWGS84PointValue(pt, p)) continue;
                bz.polygon.emplace_back(WGS84Coord(p.lon, p.lat, p.alt));
            }
            if (bz.polygon.size() >= 3) {
                if(idx < temp_battle_high_list.size()) {
                    bz.height_range = {temp_battle_high_list[idx], temp_battle_high_list[idx]}; // 或者根据实际使用情况赋值
                }
                if(idx < temp_battle_zone_link_flag.size()) {
                    bz.link_flag = temp_battle_zone_link_flag[idx];
                }
                input_data.battle_zones.push_back(std::move(bz));
            }
            idx++;
        }
    }
    if (!input_data.uav_leader_ids.empty()) {
        input_data.uav_leader_id = input_data.uav_leader_ids.front();
    }

    if (input_json.contains("uav_leader_start_point_wgs84") && input_json["uav_leader_start_point_wgs84"].is_array() &&
        !input_json["uav_leader_start_point_wgs84"].empty()) {
        WGS84Point p;
        if (parseWGS84PointValue(input_json["uav_leader_start_point_wgs84"][0], p)) {
            input_data.uav_leader_start_point_wgs84 = WGS84Coord(p.lon, p.lat, p.alt);
        }
    }

    if (input_json.contains("uavs_plane_data") && input_json["uavs_plane_data"].is_array()) {
        for (const auto &iter : input_json["uavs_plane_data"]) {
            if (!iter.is_array() || iter.size() < 3) continue;
            if (!iter[0].is_number_integer() || !iter[1].is_number_integer() || !iter[2].is_number_integer()) continue;
            int uid = iter[0].get<int>();
            int seg = iter[1].get<int>();
            int idx = iter[2].get<int>();
            input_data.uavs_plane_data = std::make_pair(uid, std::make_pair(seg, idx));
            input_data.uavs_plane_data_list.push_back({uid, seg, idx});
        }
    }

    auto parse_zone_key = [&](const char *key, std::vector<ProhibitedZone> &out, bool &has_zone_flag) {
        has_zone_flag = false;
        if (!input_json.contains(key) || !input_json[key].is_array() || input_json[key].empty()) return;

        for (const auto &zone : input_json[key]) {
            ProhibitedZone pz;
            pz.height_range = {-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};

            if (zone.is_array()) {
                if (zone.size() < 3) continue;

                bool has_height_range = false;
                double hmin = pz.height_range.first;
                double hmax = pz.height_range.second;
                if (zone.size() >= 4 && parseHeightRangeArray(zone.back(), hmin, hmax)) {
                    pz.height_range = {hmin, hmax};
                    has_height_range = true;
                }

                size_t pt_end = has_height_range ? (zone.size() - 1) : zone.size();
                for (size_t i = 0; i < pt_end; ++i) {
                    WGS84Point p;
                    if (parseWGS84PointValue(zone[i], p)) {
                        pz.polygon.emplace_back(WGS84Coord(p.lon, p.lat, p.alt));
                    }
                }
            } else if (zone.is_object()) {
                if (zone.contains("height_range")) {
                    double hmin = pz.height_range.first;
                    double hmax = pz.height_range.second;
                    if (parseHeightRangeArray(zone["height_range"], hmin, hmax)) {
                        pz.height_range = {hmin, hmax};
                    }
                }

                const char *poly_keys[] = {"polygon", "points", "zone"};
                for (const char *poly_key : poly_keys) {
                    if (!zone.contains(poly_key) || !zone[poly_key].is_array()) continue;
                    for (const auto &pt : zone[poly_key]) {
                        WGS84Point p;
                        if (parseWGS84PointValue(pt, p)) {
                            pz.polygon.emplace_back(WGS84Coord(p.lon, p.lat, p.alt));
                        }
                    }
                    if (!pz.polygon.empty()) break;
                }
            }

            if (pz.polygon.size() >= 3) {
                out.push_back(std::move(pz));
            }
        }

        has_zone_flag = !out.empty();
    };

    parse_zone_key("prohibited_zone_wgs84", input_data.prohibited_zones, input_data.has_prohibited_zone);
    parse_zone_key("check_prohibited_zone_wgs84", input_data.check_prohibited_zones, input_data.has_check_prohibited_zone);

    if (input_json.contains("ready_high_list") && input_json["ready_high_list"].is_array() && input_json["ready_high_list"].size() >= 2 &&
        input_json["ready_high_list"][0].is_number() && input_json["ready_high_list"][1].is_number()) {
        input_data.ready_zone.height_range = std::make_pair(input_json["ready_high_list"][0].get<double>(), input_json["ready_high_list"][1].get<double>());
    }
    if (input_json.contains("high_list") && input_json["high_list"].is_array() && input_json["high_list"].size() >= 2 &&
        input_json["high_list"][0].is_number() && input_json["high_list"][1].is_number()) {
        input_data.height_list = std::make_pair(input_json["high_list"][0].get<double>(), input_json["high_list"][1].get<double>());
    }

    if (input_json.contains("min_turning_radius") && input_json["min_turning_radius"].is_number()) {
        input_data.min_turning_radius = input_json["min_turning_radius"].get<double>();
    }
    if (input_data.min_turning_radius <= 0.0 && this->config_.path_planning.min_turning_radius > 0.0) {
        input_data.min_turning_radius = this->config_.path_planning.min_turning_radius;
    }

    input_data.using_midway_lines.clear();
    input_data.existing_midway_lines.clear();
    if (input_json.contains("using_midway_lines") && input_json["using_midway_lines"].is_array()) {
        for (const auto &line : input_json["using_midway_lines"]) {
            if (!line.is_array() || line.size() <= 2) {
                continue;
            }
            if (!line[0].is_number_integer() || !line[1].is_number_integer()) {
                continue;
            }

            UavTrajectoryLine typed_line;
            typed_line.uav_id = line[0].get<int>();
            typed_line.segment_id = line[1].get<int>();

            for (size_t i = 2; i < line.size(); ++i) {
                WGS84Point p;
                if (parseWGS84PointValue(line[i], p)) {
                    typed_line.points.emplace_back(WGS84Coord(p.lon, p.lat, p.alt));
                    input_data.existing_midway_lines.emplace_back(WGS84Coord(p.lon, p.lat, p.alt));
                }
            }

            input_data.using_midway_lines.push_back(std::move(typed_line));
        }
    }

    // 编队参数覆盖
    if (input_json.contains("formation_distance") && input_json["formation_distance"].is_number()) {
        input_data.formation_distance = input_json["formation_distance"].get<double>();
    } else if (input_json.contains("safety_distance") && input_json["safety_distance"].is_number()) {
        input_data.formation_distance = input_json["safety_distance"].get<double>();
    }
    if (input_json.contains("position_misalignment") && input_json["position_misalignment"].is_number()) {
        input_data.position_misalignment = input_json["position_misalignment"].get<double>();
    }
    if (input_json.contains("uav_R") && input_json["uav_R"].is_number()) {
        input_data.uav_R = input_json["uav_R"].get<double>();
    }
    if (input_json.contains("uav_formation_max_row") && input_json["uav_formation_max_row"].is_number_integer()) {
        input_data.uav_formation_max_row = input_json["uav_formation_max_row"].get<int>();
    }

    // 高度优化参数覆盖
    if (input_json.contains("uav_R") && input_json["uav_R"].is_number()) input_data.ao_uav_R = input_json["uav_R"].get<double>();
    if (input_json.contains("safe_distance") && input_json["safe_distance"].is_number()) input_data.ao_safe_distance = input_json["safe_distance"].get<double>();
    if (input_json.contains("lambda_follow") && input_json["lambda_follow"].is_number()) input_data.ao_lambda_follow = input_json["lambda_follow"].get<double>();
    if (input_json.contains("lambda_smooth") && input_json["lambda_smooth"].is_number()) input_data.ao_lambda_smooth = input_json["lambda_smooth"].get<double>();
    if (input_json.contains("max_climb_rate") && input_json["max_climb_rate"].is_number()) input_data.ao_max_climb_rate = input_json["max_climb_rate"].get<double>();

    return true;
}

//创建局部代价地图（ENU坐标系）
void UavPathPlanner::buildLocalENUCostMap(double margin, double resolution)//margin:边距参数,扩大局部地图覆盖范围
{
    if (!elev_cost_map_ || !elev_cost_map_->isElevationValid() || Trajectory_ENU.empty()) return;

  // 1. Determine bounding box of the trajectory in ENU
  double min_e = std::numeric_limits<double>::infinity();
  double max_e = -std::numeric_limits<double>::infinity();
  double min_n = std::numeric_limits<double>::infinity();
  double max_n = -std::numeric_limits<double>::infinity();

  for (const auto& p : Trajectory_ENU) {
      if (p.east < min_e) min_e = p.east;
      if (p.east > max_e) max_e = p.east;
      if (p.north < min_n) min_n = p.north;
      if (p.north > max_n) max_n = p.north;
  }

  // Add margin
  min_e -= margin;
  max_e += margin;
  min_n -= margin;
  max_n += margin;

  // 2. Create CostMap grid in ENU
  int w = static_cast<int>(std::ceil((max_e - min_e) / resolution));
  int h = static_cast<int>(std::ceil((max_n - min_n) / resolution));
  
  // Ensure at least 1x1
  if (w <= 0) w = 1;
  if (h <= 0) h = 1;

    elev_cost_map_->createCostMap(w, h, resolution, min_e, max_n); // Origin is top-left (min_e, max_n)

  std::cerr << std::fixed << std::setprecision(15)
            // << "Building local ENU CostMap: " << w << "x" << h 
            // << " res=" << resolution 
            // << " origin=(" << min_e << "," << max_n << ")" 
            << " Ref(Lon,Lat)=(" << this->origin_.lon << "," << this->origin_.lat << ")"
            << std::endl;
  std::cerr.unsetf(std::ios::fixed);

  // 3. Fill CostMap by sampling ElevationMap
  // Iterate over CostMap pixels (ENU), convert to WGS84, query ElevationMap
  for (int r = 0; r < h; ++r) {
      for (int c = 0; c < w; ++c) {
          // Center of pixel in ENU
          double enu_e = min_e + (c + 0.5) * resolution;
          double enu_n = max_n - (r + 0.5) * resolution; // Y decreases down
          
          ENUPoint enu_pt;
          enu_pt.east = enu_e;
          enu_pt.north = enu_n;
          enu_pt.up = 0.0; // Height doesn't matter for conversion to LatLon

          WGS84Point wgs_pt = this->enuToWGS84(enu_pt, this->origin_);
          
          double elev_val = -std::numeric_limits<double>::infinity();
          // Query ElevationMap (which is in WGS84)
          // Note: getElevationAt expects x,y in the ElevationMap's CRS (WGS84 Lon/Lat)
          if (elev_cost_map_->getElevationAt(wgs_pt.lon, wgs_pt.lat, elev_val)) {
              elev_cost_map_->setCost(c, r, static_cast<float>(elev_val));
          } else {
              // If out of bounds of elevation map, set to -inf or some default
              elev_cost_map_->setCost(c, r, -std::numeric_limits<float>::infinity());
          }
      }
  }
  std::cerr << "Local ENU CostMap built." << std::endl;
}

// 计算轨迹最小转弯半径
double UavPathPlanner::calculateMinTurningRadius(const std::vector<ENUPoint>& path) {
    if (path.size() < 3) {
        return -1.0; // Not enough points
    }

    double min_radius = std::numeric_limits<double>::max();

    for (size_t i = 0; i < path.size() - 2; ++i) {
        const auto& p1 = path[i];
        const auto& p2 = path[i+1];
        const auto& p3 = path[i+2];

        double a = std::sqrt(std::pow(p2.east - p3.east, 2) + std::pow(p2.north - p3.north, 2) + std::pow(p2.up - p3.up, 2));
        double b = std::sqrt(std::pow(p1.east - p3.east, 2) + std::pow(p1.north - p3.north, 2) + std::pow(p1.up - p3.up, 2));
        double c = std::sqrt(std::pow(p1.east - p2.east, 2) + std::pow(p1.north - p2.north, 2) + std::pow(p1.up - p2.up, 2));

        if (a < 1e-3 || b < 1e-3 || c < 1e-3) continue; // Skip if points are too close

        double s = (a + b + c) / 2.0;
        double area_sq = s * (s - a) * (s - b) * (s - c);
        
        if (area_sq < 1e-6) continue; // Collinear or nearly collinear

        double radius = (a * b * c) / (4.0 * std::sqrt(area_sq));
        if (radius < min_radius) {
            min_radius = radius;
        }
    }

    if (min_radius == std::numeric_limits<double>::max()) {
        return -1.0; // Could not calculate (e.g. straight line)
    }

    return min_radius;
}

// 计算第二段过渡轨迹（切圆切入优化）并更新第三段巡逻轨迹
void UavPathPlanner::computeTransitionAndRotatePatrol(const ENUPoint& p0, double heading0, double minR, double resolution, 
                                                      const std::vector<ENUPoint>& Patrol_Path, OutputData &output_data)
{
    std::cout<<"开始计算第二段任务区域过渡轨迹 (切圆切入优化)"<<std::endl;
    std::vector<ENUPoint> Transition_Path;
    std::vector<ENUPoint> Rotated_Patrol;
    this->buildTransitionAndRotatePatrol(p0, heading0, minR, resolution, Patrol_Path,
                                         Transition_Path, Rotated_Patrol);

    if (!Rotated_Patrol.empty()) {
        this->enforceTransitionClimbRateAndBorrowPatrolPrefix(Transition_Path, Rotated_Patrol,
                                                              "leader plane2(formation)");
    }

    if (!Rotated_Patrol.empty()) {
        std::vector<WGS84Point> Rotated_Patrol_WGS84 = this->enuToWGS84_Batch(Rotated_Patrol, this->origin_);
        this->writeLeaderSegment(output_data.uav_leader_plane3, output_data, 3, Rotated_Patrol_WGS84);
    }

    std::vector<WGS84Point> Transition_Path_WGS84 = this->enuToWGS84_Batch(Transition_Path,this->origin_);

    this->writeLeaderSegment(output_data.uav_leader_plane2, output_data, 2, Transition_Path_WGS84);
}
//避开禁飞区域,高度避障或者绕障
std::vector<ENUPoint> UavPathPlanner::avoidProhibitedZones(const std::vector<ENUPoint>& path) {
    if (input_data_.prohibited_zones.empty() || path.size() < 2) {
        return path;
    }

    // Pre-process zones to ENU polygons
    struct ENUZone {
        Polygon2d poly;
        double min_h, max_h;
    };
    std::vector<ENUZone> enu_zones;
    for (const auto& zone : input_data_.prohibited_zones) {
        std::vector<Vec2d> points;
        for (const auto& p : zone.polygon) {
            WGS84Point pt_wgs84;
            pt_wgs84.lon = p.lon;
            pt_wgs84.lat = p.lat;
            pt_wgs84.alt = p.alt;
            ENUPoint enu = wgs84ToENU(pt_wgs84, origin_);
            points.emplace_back(enu.east, enu.north);
        }
        if (points.size() < 3) continue;
        enu_zones.push_back({Polygon2d(points), zone.height_range.first, zone.height_range.second});
    }

    std::vector<ENUPoint> current_path = path;
    bool collision_found = true;
    int max_iterations = 5;
    int iter = 0;

    while (collision_found && iter < max_iterations) {
        collision_found = false;
        std::vector<ENUPoint> next_path;
        next_path.push_back(current_path[0]);
        iter++;
        
        // std::cout << "Avoidance iteration " << iter << ", path size: " << current_path.size() << std::endl;

        for (size_t i = 0; i < current_path.size() - 1; ++i) {
            ENUPoint p1 = next_path.back(); 
            ENUPoint p2 = current_path[i+1];
            
            LineSegment2d segment(Vec2d(p1.east, p1.north), Vec2d(p2.east, p2.north));
            double seg_min_h = std::min(p1.up, p2.up);
            double seg_max_h = std::max(p1.up, p2.up);

            int intersected_zone_idx = -1;
            
            for (size_t z = 0; z < enu_zones.size(); ++z) {
                const auto& zone = enu_zones[z];
                if (seg_max_h < zone.min_h || seg_min_h > zone.max_h) continue;

                if (zone.poly.DistanceTo(segment) < this->config_.path_planning.prohibited_zone_conflict_distance) { // 判定为冲突的距离阈值（米）
                    intersected_zone_idx = z;
                    break; 
                }
            }

            if (intersected_zone_idx != -1) {
                collision_found = true;
                std::cout << "Avoidance: Segment intersects prohibited zone " << intersected_zone_idx << " (iter " << iter << ")" << std::endl;
                const auto& zone = enu_zones[intersected_zone_idx];
                
                // --- 1. Calculate Horizontal Detour Cost ---
                std::vector<Vec2d> nodes;
                nodes.push_back(Vec2d(p1.east, p1.north)); // Node 0
                nodes.push_back(Vec2d(p2.east, p2.north)); // Node 1
                
                Vec2d center(0,0);
                for(const auto& pt : zone.poly.points()) center += pt;
                center /= zone.poly.points().size();

                for (const auto& pt : zone.poly.points()) {
                    Vec2d dir = pt - center;
                    dir.Normalize();
                    nodes.push_back(pt + dir * 100.0); // Expand 100m for safety
                }

                int n = nodes.size();
                std::vector<double> dist(n, std::numeric_limits<double>::max());
                std::vector<int> parent(n, -1);
                dist[0] = 0;
                
                std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<>> pq;
                pq.push({0, 0});

                while(!pq.empty()) {
                    double d = pq.top().first;
                    int u = pq.top().second;
                    pq.pop();

                    if (d > dist[u]) continue;
                    if (u == 1) break; 

                    for (int v = 0; v < n; ++v) {
                        if (u == v) continue;
                        
                        Vec2d mid = (nodes[u] + nodes[v]) / 2.0;
                        // Check if strictly inside
                        if (zone.poly.IsPointIn(mid) && zone.poly.DistanceToBoundary(mid) > 0.1) continue; 

                        double weight = nodes[u].DistanceTo(nodes[v]);
                        if (dist[u] + weight < dist[v]) {
                            dist[v] = dist[u] + weight;
                            parent[v] = u;
                            pq.push({dist[v], v});
                        }
                    }
                }

                double horizontal_cost = dist[1];
                if (horizontal_cost != std::numeric_limits<double>::max()) {
                    horizontal_cost += std::abs(p2.up - p1.up);
                }

                // --- 2. Calculate Vertical Avoidance Cost ---
                double safe_alt = zone.max_h + 50.0; // 20m buffer
                double target_h = std::max(safe_alt, std::max(p1.up, p2.up));
                
                Vec2d overlap_start, overlap_end;
                bool has_overlap = zone.poly.GetOverlap(segment, &overlap_start, &overlap_end);
                
                double vertical_cost = std::numeric_limits<double>::max();
                
                if (has_overlap) {
                    // Ensure start is closer to p1
                    if (overlap_start.DistanceSquareTo(Vec2d(p1.east, p1.north)) > overlap_end.DistanceSquareTo(Vec2d(p1.east, p1.north))) {
                        std::swap(overlap_start, overlap_end);
                    }
                    
                    // Cost = p1 -> start_high -> end_high -> p2
                    // We approximate the flight path distance
                    double d1 = std::hypot(p1.east - overlap_start.x(), p1.north - overlap_start.y()); 
                    double h1 = std::abs(target_h - p1.up);
                    double leg1 = std::hypot(d1, h1); 
                    
                    double d2 = overlap_start.DistanceTo(overlap_end);
                    double leg2 = d2; // at constant height target_h
                    
                    double d3 = std::hypot(p2.east - overlap_end.x(), p2.north - overlap_end.y());
                    double h3 = std::abs(target_h - p2.up);
                    double leg3 = std::hypot(d3, h3);
                    
                    vertical_cost = leg1 + leg2 + leg3;
                } else {
                    // Fallback if no strict overlap but close
                    double dist_2d = std::hypot(p1.east - p2.east, p1.north - p2.north);
                    vertical_cost = std::abs(target_h - p1.up) + dist_2d + std::abs(target_h - p2.up);
                }

                // --- 3. Choose Best Strategy ---
                if (horizontal_cost != std::numeric_limits<double>::max() && horizontal_cost <= vertical_cost) {
                    // Apply Horizontal Detour
                    std::cout << "  Strategy: Horizontal Detour (Cost: " << horizontal_cost << " vs Vertical: " << vertical_cost << ")" << std::endl;
                    std::vector<ENUPoint> detour;
                    int curr = 1;
                    while (curr != 0) {
                        detour.push_back({nodes[curr].x(), nodes[curr].y(), 0.0}); 
                        curr = parent[curr];
                    }
                    std::reverse(detour.begin(), detour.end());
                    
                    for (size_t k = 0; k < detour.size(); ++k) {
                        if (k == detour.size() - 1) {
                            detour[k].up = p2.up;
                        } else {
                            detour[k].up = p1.up;
                        }
                    }
                    next_path.insert(next_path.end(), detour.begin(), detour.end());

                } else {
                    // Apply Vertical Avoidance
                    std::cout << "  Strategy: Vertical Avoidance (Cost: " << vertical_cost << " vs Horizontal: " << horizontal_cost << ")" << std::endl;
                    
                    if (has_overlap) {
                        // Ensure start is closer to p1 (re-check as swap was local to if block scope if I didn't use references, but here I used local vars)
                        // Actually I swapped the local vars overlap_start/end, so they are correct now.
                        
                        next_path.push_back({overlap_start.x(), overlap_start.y(), target_h});
                        next_path.push_back({overlap_end.x(), overlap_end.y(), target_h});
                        next_path.push_back(p2);
                    } else {
                        next_path.push_back({p1.east, p1.north, target_h});
                        next_path.push_back({p2.east, p2.north, target_h});
                        next_path.push_back(p2);
                    }
                }

            } else {
                next_path.push_back(p2);
            }
        }
        current_path = next_path;
    }
    
    if (collision_found) {
        std::cerr << "Avoidance: Max iterations reached, path might still intersect zones." << std::endl;
    }

    return current_path;
}
bool UavPathPlanner::outputDataToJson(const OutputData &output_data, json &output_json) {
    output_json = json::object();

    output_json["abnormal_uav_plane"] = output_data.abnormal_uav_plane;
    output_json["using_uav_list"] = output_data.using_uav_list;
    output_json["ready_id"] = output_data.ready_id;
    output_json["midway_point_num"] = output_data.midway_point_num;

    output_json["leader_show_points"] = json::array();
    for (const auto &p : output_data.leader_show_points) {
        json pa = json::array();
        pa.push_back(p.lon);
        pa.push_back(p.lat);
        pa.push_back(p.alt);
        output_json["leader_show_points"].push_back(pa);
    }

    output_json["uav_leader_plane1"] = json::array();
    for (const auto &p : output_data.uav_leader_plane1) {
        json pa = json::array();
        pa.push_back(p.lon);
        pa.push_back(p.lat);
        pa.push_back(p.alt);
        output_json["uav_leader_plane1"].push_back(pa);
    }

    output_json["uav_leader_plane2"] = json::array();
    for (const auto &p : output_data.uav_leader_plane2) {
        json pa = json::array();
        pa.push_back(p.lon);
        pa.push_back(p.lat);
        pa.push_back(p.alt);
        output_json["uav_leader_plane2"].push_back(pa);
    }

    output_json["uav_leader_plane3"] = json::array();
    for (const auto &p : output_data.uav_leader_plane3) {
        json pa = json::array();
        pa.push_back(p.lon);
        pa.push_back(p.lat);
        pa.push_back(p.alt);
        output_json["uav_leader_plane3"].push_back(pa);
    }

    output_json["uav_plane1"] = json::array();
    for (const auto &line : output_data.uav_plane1) {
        json entry = json::array();
        entry.push_back(line.uav_id);
        for (const auto &p : line.points) {
            json pa = json::array();
            pa.push_back(p.lon);
            pa.push_back(p.lat);
            pa.push_back(p.alt);
            entry.push_back(pa);
        }
        output_json["uav_plane1"].push_back(entry);
    }

    output_json["uav_plane2"] = json::array();
    for (const auto &line : output_data.uav_plane2) {
        json entry = json::array();
        entry.push_back(line.uav_id);
        for (const auto &p : line.points) {
            json pa = json::array();
            pa.push_back(p.lon);
            pa.push_back(p.lat);
            pa.push_back(p.alt);
            entry.push_back(pa);
        }
        output_json["uav_plane2"].push_back(entry);
    }

    output_json["uav_plane3"] = json::array();
    for (const auto &line : output_data.uav_plane3) {
        json entry = json::array();
        entry.push_back(line.uav_id);
        for (const auto &p : line.points) {
            json pa = json::array();
            pa.push_back(p.lon);
            pa.push_back(p.lat);
            pa.push_back(p.alt);
            entry.push_back(pa);
        }
        output_json["uav_plane3"].push_back(entry);
    }

    output_json["using_midway_lines"] = json::array();
    for (const auto &line : output_data.using_midway_lines) {
        json entry = json::array();
        entry.push_back(line.uav_id);
        entry.push_back(line.segment_id);
        for (const auto &p : line.points) {
            json pa = json::array();
            pa.push_back(p.lon);
            pa.push_back(p.lat);
            pa.push_back(p.alt);
            entry.push_back(pa);
        }
        output_json["using_midway_lines"].push_back(entry);
    }

    return true;
}
