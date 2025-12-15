#include "uavPathPlanning.hpp"
#include "math_util/minimum_snap.hpp"
#include "math_util/bezier.hpp"
#include <yaml-cpp/yaml.h>
#include <cmath>
#include <cpl_string.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <sys/stat.h>
#include <iomanip> // Added for std::setprecision

#ifdef HAVE_GDAL
#include <gdal_priv.h>
#include <cpl_conv.h>
#endif

// (moved into UavPathPlanner::generator_)

// 配置命名空间：将时间分配相关默认值和 YAML 加载器内联到此文件中，避免使用单独的 _config.h
namespace TimeAllocConfig {
    inline double V_avg_default = 10.0;

    inline double angle_threshold = 0.2;        // 弧度
    inline double min_time_s = 40.0;            // 最小分配时间 s
    inline double prev_influence = 0.5;         // 前一段影响系数

    inline double base_pow = 0.65;
    inline double base_divider = 50.0;

    inline double long_threshold = 10.0;
    inline double med_threshold = 3.0;
    inline double small_threshold = 1.0;
    inline double tiny_threshold = 0.1;

    inline double long_adj = 0.7;
    inline double med_adj = 1.0;
    inline double small_adj = 1.1;
    inline double tiny_adj = 1.3;

    inline double min_alloc = 1.0;
    inline double max_alloc = 300.0;

    // 从 YAML 文件加载配置，若字段缺失则保留默认值
    inline bool loadFromYAML(const std::string &yaml_path) {
        try {
            YAML::Node config = YAML::LoadFile(yaml_path);
            if (config["V_avg_default"]) V_avg_default = config["V_avg_default"].as<double>();

            if (config["angle_threshold"]) angle_threshold = config["angle_threshold"].as<double>();
            if (config["min_time_s"]) min_time_s = config["min_time_s"].as<double>();
            if (config["prev_influence"]) prev_influence = config["prev_influence"].as<double>();

            if (config["base_pow"]) base_pow = config["base_pow"].as<double>();
            if (config["base_divider"]) base_divider = config["base_divider"].as<double>();

            if (config["long_threshold"]) long_threshold = config["long_threshold"].as<double>();
            if (config["med_threshold"]) med_threshold = config["med_threshold"].as<double>();
            if (config["small_threshold"]) small_threshold = config["small_threshold"].as<double>();
            if (config["tiny_threshold"]) tiny_threshold = config["tiny_threshold"].as<double>();

            if (config["long_adj"]) long_adj = config["long_adj"].as<double>();
            if (config["med_adj"]) med_adj = config["med_adj"].as<double>();
            if (config["small_adj"]) small_adj = config["small_adj"].as<double>();
            if (config["tiny_adj"]) tiny_adj = config["tiny_adj"].as<double>();

            if (config["min_alloc"]) min_alloc = config["min_alloc"].as<double>();
            if (config["max_alloc"]) max_alloc = config["max_alloc"].as<double>();

            return true;
        } catch (const std::exception &e) {
            std::cerr << "Failed to load time allocation YAML: " << e.what() << std::endl;
            return false;
        }
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

// 生成 arc-line-arc（相切圆 - 直线 - 相切圆）路径实现
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

    // Try combinations of left/right circle sides to find a feasible tangent
    bool found = false;
    std::pair<double,double> C1, C2; // centers
    std::pair<double,double> T1, T2; // tangent points
    int best_s0 = 0;
    int best_s1 = 0;

    auto tangent_at = [&](double theta, int sign)->std::pair<double,double>{
        if (sign > 0) return std::make_pair(-sin(theta), cos(theta));
        else return std::make_pair(sin(theta), -cos(theta));
    };

    // Determine which side p1 is relative to p0's heading
    double dx_chk = p1.east - p0.east;
    double dy_chk = p1.north - p0.north;
    double cross_chk = std::cos(h0) * dy_chk - std::sin(h0) * dx_chk;
    int forced_s0 = (cross_chk >= 0) ? 1 : -1;

    for (int s0 : {forced_s0}) {
        auto n0 = rotate90(std::cos(h0), std::sin(h0), s0);
        std::pair<double,double> c1 = {p0.east + radius * n0.first, p0.north + radius * n0.second};
        int s1 = s0;
        {
            auto n1 = rotate90(std::cos(h1), std::sin(h1), s1);
            std::pair<double,double> c2 = {p1.east + radius * n1.first, p1.north + radius * n1.second};

            double vx = c2.first - c1.first; double vy = c2.second - c1.second;
            double d = std::hypot(vx, vy);
            
            if (d >= 1e-6) {
                struct TangentLine { std::pair<double,double> t1, t2; };
                std::vector<TangentLine> candidates;

                // External tangents only (since s1 == s0)
                for (int sign : {1, -1}) {
                    auto vperp = rotate90(vx/d, vy/d, sign);
                    candidates.push_back({ 
                        {c1.first + radius * vperp.first, c1.second + radius * vperp.second},
                        {c2.first + radius * vperp.first, c2.second + radius * vperp.second}
                    });
                }

                for (const auto& line : candidates) {
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

                     // Found valid
                     found = true;
                     C1 = c1; C2 = c2; T1 = line.t1; T2 = line.t2;
                     best_s0 = s0; best_s1 = s1;
                     break;
                }
            }
        }
        if (found) break;
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
bool UavPathPlanner::runAltitudeOptimization(const std::string &elev_file, json &output_json, const json &input_json)
{
    try {
        if (elev_file.empty()) return false;
        if (this->Trajectory_ENU.empty()) {
            std::cerr << "runAltitudeOptimization: Trajectory_ENU is empty, nothing to optimize\n";
            return false;
        }

        // Filter Trajectory_ENU logic removed as per request (filtering input points instead)

        if (!this->loadElevationData(elev_file)) {
            std::cerr << "AltitudeOptimizer: failed to load elevation file: " << elev_file << "\n";
            return false;
        }

        // Initialize CostMap from ElevationMap
        // this->initCostMapFromElevation();
        // Instead of copying the whole map (which is in WGS84), we build a local ENU costmap
        // around the trajectory. This solves the coordinate system mismatch.
        // Use a reasonable resolution (e.g. 10m) and margin (e.g. 200m)
        this->buildLocalENUCostMap(500.0, 10.0);

        // build waypoint vector for optimizer from member Trajectory_ENU
        std::vector<Eigen::Vector3d> wpts;
        wpts.reserve(this->Trajectory_ENU.size());
        for (size_t i = 0; i < this->Trajectory_ENU.size(); ++i) {
            const auto &p = this->Trajectory_ENU[i];
            wpts.emplace_back(p.east, p.north, p.up);
        }

        AltitudeParams params;
        // Try to load from config.yaml
        std::vector<std::string> config_paths = {"config.yaml", "../config.yaml", "../../config.yaml"};
        bool config_loaded = false;
        for (const auto& config_path : config_paths) {
            struct stat buffer;
            if (stat(config_path.c_str(), &buffer) == 0) {
                 try {
                    YAML::Node config = YAML::LoadFile(config_path);
                    if (config["altitude_optimization"]) {
                        auto alt_conf = config["altitude_optimization"];
                        if (alt_conf["uav_R"]) params.uav_R = alt_conf["uav_R"].as<double>();
                        if (alt_conf["safe_distance"]) params.safe_distance = alt_conf["safe_distance"].as<double>();
                        if (alt_conf["lambda_follow"]) params.lambda_follow = alt_conf["lambda_follow"].as<double>();
                        if (alt_conf["lambda_smooth"]) params.lambda_smooth = alt_conf["lambda_smooth"].as<double>();
                        if (alt_conf["max_climb_rate"]) params.max_climb_rate = alt_conf["max_climb_rate"].as<double>();
                        std::cerr << "AltitudeOptimizer: loaded params from " << config_path << "\n";
                        config_loaded = true;
                        break;
                    }
                 } catch (const std::exception &e) {
                     std::cerr << "AltitudeOptimizer: failed to parse " << config_path << ": " << e.what() << std::endl;
                 }
            }
        }
        
        if (!config_loaded) {
             std::cerr << "AltitudeOptimizer: config.yaml not found or invalid, using defaults.\n";
        }

        std::vector<double> out_z;
        if (this->optimizeHeights(wpts, params, out_z)) {
            // First pass applied to Trajectory_ENU
            for (size_t i = 0; i < out_z.size() && i < this->Trajectory_ENU.size(); ++i) {
                this->Trajectory_ENU[i].up = out_z[i];
            }
            std::cerr << "AltitudeOptimizer: first pass optimized heights applied (" << out_z.size() << " points)\n";

            // Second pass: Global Smoothing
            // std::cerr << "AltitudeOptimizer: running second pass (global smoothing)...\n";
            std::vector<double> out_z_smooth;
            
            // Adjust params for second pass
            AltitudeParams params_smooth = params;
            params_smooth.lambda_smooth *= 10.0; // Increase smoothing weight
            params_smooth.max_climb_rate *= 0.5; // Stricter climb rate
            
            if (this->optimizeHeightsGlobalSmooth(out_z, wpts, params_smooth, out_z_smooth)) {
                 out_z = out_z_smooth;
                 for (size_t i = 0; i < out_z.size() && i < this->Trajectory_ENU.size(); ++i) {
                    this->Trajectory_ENU[i].up = out_z[i];
                 }
                //  std::cerr << "AltitudeOptimizer: second pass completed and applied.\n";
            } else {
                 std::cerr << "AltitudeOptimizer: second pass failed, keeping first pass results.\n";
            }

            // std::cout << "--- Optimized Trajectory Heights ---" << std::endl;
            // for (size_t i = 0; i < out_z.size(); ++i) {
            //     std::cout << "Point " << i << ": " << out_z[i] << " m" << std::endl;
            // }
            // std::cout << "------------------------------------" << std::endl;

            // Update output_json
            // 1. Convert updated ENU trajectory to WGS84
            std::vector<WGS84Point> Trajectory_WGS84 = this->enuToWGS84_Batch(this->Trajectory_ENU, this->origin_);
            
            // 2. Update leader trajectory in output_json
            this->putWGS84ToJson(output_json, "uav_leader_plane1", Trajectory_WGS84);
            
            // 3. Re-generate follower trajectories
            InputData input_data;
            json input_json_copy = input_json; // Make a copy because loadData takes non-const ref
            if (this->loadData(input_data, input_json_copy)) {
                 try {
                    json plane_array = this->generateFollowerTrajectories(input_json_copy, input_data, this->Trajectory_ENU, Trajectory_WGS84);
                    output_json["uav_plane1"] = plane_array;
                    std::cerr << "AltitudeOptimizer: updated output_json with new trajectories.\n";
                } catch (const std::exception &e) {
                    std::cerr << "AltitudeOptimizer: failed to update follower trajectories: " << e.what() << std::endl;
                }
            } else {
                 std::cerr << "AltitudeOptimizer: failed to parse input_json for follower update.\n";
            }

            return true;
        } else {
            std::cerr << "AltitudeOptimizer: optimization failed\n";
            return false;
        }
    } catch (const std::exception &e) {
        std::cerr << "runAltitudeOptimization exception: " << e.what() << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

bool UavPathPlanner::loadElevationData(const std::string &path)
{
  std::string path_to_load = path;
  
  struct stat stat_buf;
  if (stat(path.c_str(), &stat_buf) == 0) {
      long long size = stat_buf.st_size;
      long long limit = 200LL * 1024 * 1024; // 200MB
      
      if (size > limit) {
          std::string ovr_path = path + ".ovr";
          struct stat stat_ovr;
          if (stat(ovr_path.c_str(), &stat_ovr) == 0) {
              std::cout << "Elevation file is large (" << size << " bytes). Found .ovr file, using it instead: " << ovr_path << std::endl;
              path_to_load = ovr_path;
          }
      }
  }

#ifdef HAVE_GDAL
  GDALAllRegister();
  std::cerr << "ElevationMap: opening TIFF: " << path_to_load << std::endl;
  // open readonly first to probe size and geotransform
  GDALDataset *poDS = (GDALDataset*)GDALOpen(path_to_load.c_str(), GA_ReadOnly);
  if (!poDS) {
    cerr << "ElevationMap: GDALOpen failed for " << path_to_load << "\n";
    return false;
  }
  int full_w = poDS->GetRasterXSize();
  int full_h = poDS->GetRasterYSize();
  std::cerr << "ElevationMap: dataset size: " << full_w << " x " << full_h << std::endl;
  double gt[6];
  if (poDS->GetGeoTransform(gt) == CE_None) {
    elev_geo_.origin_x = gt[0];
    elev_geo_.pixel_w = gt[1];
    elev_geo_.rot_x = gt[2];
    elev_geo_.origin_y = gt[3];
    elev_geo_.rot_y = gt[4];
    elev_geo_.pixel_h = gt[5];
    std::cerr << "ElevationMap: geotransform: origin(" << elev_geo_.origin_x << ", " << elev_geo_.origin_y << ") pixel_w=" << elev_geo_.pixel_w << " pixel_h=" << elev_geo_.pixel_h << std::endl;
  } else {
    // default geotransform (pixel-as-meter)
    elev_geo_.origin_x = 0.0; elev_geo_.pixel_w = 1.0; elev_geo_.rot_x = 0.0;
    elev_geo_.origin_y = 0.0; elev_geo_.rot_y = 0.0; elev_geo_.pixel_h = -1.0;
  }

  // estimate memory needed (float32 per pixel)
  const uint64_t bytes_needed = static_cast<uint64_t>(full_w) * static_cast<uint64_t>(full_h) * sizeof(float);
   // target maximum bytes after downsampling: 1 GB (user request)
   const uint64_t TARGET_BYTES = 200ULL * 1024ULL * 1024ULL;
//    std::cerr << "ElevationMap: estimated memory for full raster (bytes): " << bytes_needed << " target_bytes=" << TARGET_BYTES << std::endl;

  // Check for existing overviews first
  GDALRasterBand *band = poDS->GetRasterBand(1);
  int existing_ov = band ? band->GetOverviewCount() : 0;
  bool found_suitable_ov = false;
  if (existing_ov > 0) {
    for (int i = 0; i < existing_ov; ++i) {
      GDALRasterBand *ob = band->GetOverview(i);
      if (ob) {
        uint64_t ow = static_cast<uint64_t>(ob->GetXSize());
        uint64_t oh = static_cast<uint64_t>(ob->GetYSize());
        if (ow * oh * sizeof(float) <= TARGET_BYTES) {
          found_suitable_ov = true;
          break;
        }
      }
    }
  }

   // If file is large and no suitable overview exists, perform in-code downsampling using MAX pooling
   if (bytes_needed > TARGET_BYTES && !found_suitable_ov) {
     bool ok = performDownsampling(poDS, full_w, full_h, bytes_needed, TARGET_BYTES, path_to_load);
     GDALClose(poDS);
     if (ok) {
       printElevationInfo();
       return true;
     } else {
       return false;
     }
   }  // If file is large and no overviews, try to create external overviews (.ovr) by reopening in update mode
  if (bytes_needed > TARGET_BYTES && existing_ov == 0) {
    std::cerr << "ElevationMap: large raster and no overviews; attempting to build external overviews..." << std::endl;
  GDALClose(poDS);
  // try to open in update mode to build overviews
  // Force creation of external overviews (.ovr) by setting GDAL config option
  CPLSetConfigOption("GDAL_TIFF_OVR", "YES");
  std::cerr << "ElevationMap: CPLSetConfigOption GDAL_TIFF_OVR=YES to force external overviews" << std::endl;
  GDALDataset *poDSup = (GDALDataset*)GDALOpen(path_to_load.c_str(), GA_Update);
    if (poDSup) {
      // build overview levels 2,4,8,... until smallest dimension < 256
      std::vector<int> ov_levels;
      int w = full_w, h = full_h;
      int level = 2;
      while ((w / level) >= 25600 || (h / level) >= 25600) {
        ov_levels.push_back(level);
        level *= 2;
      }
      std::cerr << "ElevationMap: requested overview levels: ";
      for (size_t i=0;i<ov_levels.size();++i) std::cerr << ov_levels[i] << (i+1<ov_levels.size()?",":"\n");
  if (!ov_levels.empty()) {
        // BuildOverviews signature expects: resampling, nListCount, panOverviewList, nBands, panBands, pfnProgress, pProgressData
        int nBands = 1;
        int panBands[1] = {1};
        // Try to build overviews using MAX sampling (take maximum of covered pixels) to preserve peaks.
        const char *resampling_try = "MAX";
        CPLErr berr = poDSup->BuildOverviews(resampling_try, static_cast<int>(ov_levels.size()), ov_levels.data(), nBands, panBands, nullptr, nullptr);
        if (berr != CE_None) {
          std::cerr << "ElevationMap: BuildOverviews with resampling=MAX failed, falling back to AVERAGE for " << path_to_load << "\n";
          resampling_try = "AVERAGE";
          berr = poDSup->BuildOverviews(resampling_try, static_cast<int>(ov_levels.size()), ov_levels.data(), nBands, panBands, nullptr, nullptr);
        }
        if (berr != CE_None) {
          cerr << "ElevationMap: BuildOverviews failed for " << path_to_load << "\n";
        } else {
          std::cerr << "ElevationMap: external overviews built for " << path_to_load << " using resampling=" << resampling_try << "\n";
        }
      } else {
        std::cerr << "ElevationMap: no overview levels chosen (small raster)" << std::endl;
      }
      GDALClose(poDSup);
    } else {
      // cannot open update; warn and fall back to read-only
      std::cerr << "ElevationMap: cannot open " << path_to_load << " in update mode to build overviews; continuing" << std::endl;
    }
    // clear the config option we set earlier
    CPLSetConfigOption("GDAL_TIFF_OVR", nullptr);
    // reopen readonly to continue
    poDS = (GDALDataset*)GDALOpen(path_to_load.c_str(), GA_ReadOnly);
    if (!poDS) {
      cerr << "ElevationMap: GDALOpen failed for " << path_to_load << " after overview attempt\n";
      return false;
    }
    band = poDS->GetRasterBand(1);
    existing_ov = band ? band->GetOverviewCount() : 0;
  }

  // If overviews exist, try to pick a lower-resolution overview that fits memory threshold
  GDALRasterBand *readBand = band;
  int use_ov_index = -1;
  if (band && existing_ov > 0) {
    // iterate from largest overview (highest index) to smallest
    for (int i = existing_ov - 1; i >= 0; --i) {
      GDALRasterBand *ob = band->GetOverview(i);
      if (!ob) continue;
      uint64_t ow = static_cast<uint64_t>(ob->GetXSize());
      uint64_t oh = static_cast<uint64_t>(ob->GetYSize());
  uint64_t need = ow * oh * sizeof(float);
  if (need <= TARGET_BYTES) {
        readBand = ob;
        use_ov_index = i;
        elev_width_ = static_cast<int>(ow);
        elev_height_ = static_cast<int>(oh);
        // adjust geotransform: pixel size scales by factor = full_w / ow (assume integer)
        double scale_x = static_cast<double>(full_w) / static_cast<double>(ow);
        double scale_y = static_cast<double>(full_h) / static_cast<double>(oh);
        elev_geo_.pixel_w = elev_geo_.pixel_w * scale_x;
        elev_geo_.pixel_h = elev_geo_.pixel_h * scale_y;
        break;
      }
    }
    // std::cerr << "ElevationMap: existing overviews count=" << existing_ov << " chosen overview index=" << use_ov_index << std::endl;
  }

  // If we didn't pick an overview, read full resolution
  if (use_ov_index == -1) {
    elev_width_ = full_w;
    elev_height_ = full_h;
  }

  // allocate and read from the chosen band (overview or base)
  elev_data_.assign(static_cast<size_t>(elev_width_) * static_cast<size_t>(elev_height_), std::numeric_limits<float>::quiet_NaN());
  CPLErr err;
  if (readBand) {
    // readBand may be an overview; use RasterIO on the band
    err = readBand->RasterIO(GF_Read, 0, 0, elev_width_, elev_height_, elev_data_.data(), elev_width_, elev_height_, GDT_Float32, 0, 0);
  } else {
    err = CE_Failure;
  }
  if (err != CE_None) {
    cerr << "ElevationMap: GDAL RasterIO failed\n";
    GDALClose(poDS);
    return false;
  }
  GDALClose(poDS);
  elev_valid_ = true;
  std::cerr << "ElevationMap: finished reading raster (width=" << elev_width_ << " height=" << elev_height_ << ")" << std::endl;
  // print a short summary of the map
//   printElevationInfo();
  return true;
#else
  cerr << "ElevationMap: GDAL not available at build time. Cannot load TIFF.\n";
  return false;
#endif
}

bool UavPathPlanner::performDownsampling(void* poDS_ptr, int full_w, int full_h, uint64_t bytes_needed, uint64_t target_bytes, const std::string& path) {
#ifdef HAVE_GDAL
    GDALDataset *poDS = (GDALDataset*)poDS_ptr;
    std::cerr << "ElevationMap: raster exceeds target size and no suitable overview found, performing in-code max downsample\n";
    GDALRasterBand *band = poDS->GetRasterBand(1);
    if (!band) {
      std::cerr << "ElevationMap: no raster band found\n";
      return false;
    }
    int hasNoData = 0;
    double noDataVal = band->GetNoDataValue(&hasNoData);

    // compute initial reduction factor so that output fits within TARGET_BYTES
    double scale = std::sqrt((double)bytes_needed / (double)target_bytes);
    int factor = static_cast<int>(std::ceil(scale));
    if (factor < 1) factor = 1;

    // We'll try decreasing the factor if the resulting downsample has too little valid data
    const double MIN_VALID_FRACTION = 0.01; // require at least 1% valid pixels
    const int MAX_ITER = 8;
    int iter = 0;
    bool done = false;
    std::vector<float> best_data;
    int best_w = 0, best_h = 0;
    while (!done && iter < MAX_ITER) {
      ++iter;
      int out_w = static_cast<int>((full_w + factor - 1) / factor);
      int out_h = static_cast<int>((full_h + factor - 1) / factor);
      std::cerr << "ElevationMap: attempt " << iter << " downsample factor=" << factor << " out_size=" << out_w << "x" << out_h << std::endl;

      std::vector<float> outdata(static_cast<size_t>(out_w) * static_cast<size_t>(out_h), std::numeric_limits<float>::quiet_NaN());
      std::vector<char> hasValTmp(static_cast<size_t>(out_w) * static_cast<size_t>(out_h), 0);
      std::vector<float> rowbuf(full_w);
      for (int y = 0; y < full_h; ++y) {
        CPLErr rerr = band->RasterIO(GF_Read, 0, y, full_w, 1, rowbuf.data(), full_w, 1, GDT_Float32, 0, 0);
        if (rerr != CE_None) {
          std::cerr << "ElevationMap: RasterIO row read failed at y=" << y << "\n";
          return false;
        }
        int oy = y / factor;
        for (int x = 0; x < full_w; ++x) {
          float v = rowbuf[x];
          if (!std::isfinite(v)) continue;
          bool isNoData = false;
          if (hasNoData) {
            if (v == static_cast<float>(noDataVal)) isNoData = true;
          } else {
            if (v == -32767.0f || v == -32768.0f || v == -9999.0f || v == -99999.0f) isNoData = true;
          }
          if (isNoData) continue;
          int ox = x / factor;
          size_t idx = static_cast<size_t>(oy) * static_cast<size_t>(out_w) + static_cast<size_t>(ox);
          if (!hasValTmp[idx]) {
            outdata[idx] = v;
            hasValTmp[idx] = 1;
          } else {
            if (v > outdata[idx]) outdata[idx] = v;
          }
        }
      }
      size_t total_out = static_cast<size_t>(out_w) * static_cast<size_t>(out_h);
      size_t valid_count = 0;
      for (size_t i = 0; i < total_out; ++i) if (hasValTmp[i]) ++valid_count;
      double valid_frac = total_out ? (double)valid_count / (double)total_out : 0.0;
      std::cerr << "ElevationMap: valid fraction=" << valid_frac << " (" << valid_count << "/" << total_out << ")" << std::endl;

      if (valid_frac >= MIN_VALID_FRACTION || factor == 1) {
        // accept this result
        best_data.swap(outdata);
        best_w = out_w; best_h = out_h;
        done = true;
        break;
      } else {
        // keep as fallback, but try to reduce factor to get more valid pixels
        best_data.swap(outdata);
        best_w = out_w; best_h = out_h;
        int new_factor = std::max(1, factor / 2);
        if (new_factor == factor) { done = true; break; }
        factor = new_factor;
        std::cerr << "ElevationMap: too few valid pixels, reducing factor to " << factor << " and retrying" << std::endl;
      }
    }

    if (best_data.empty()) {
      std::cerr << "ElevationMap: downsample produced no data" << std::endl;
      return false;
    }

    // move best_data into data_
    elev_width_ = best_w; elev_height_ = best_h;
    elev_data_.swap(best_data);
    // update geo transform based on final factor estimate: compute effective factor = full_w / width_
    double eff_scale_x = static_cast<double>(full_w) / static_cast<double>(elev_width_);
    double eff_scale_y = static_cast<double>(full_h) / static_cast<double>(elev_height_);
    elev_geo_.pixel_w = elev_geo_.pixel_w * eff_scale_x;
    elev_geo_.pixel_h = elev_geo_.pixel_h * eff_scale_y;

    // Save the downsampled data to .ovr file for future use
    std::string ovrPath = path + ".ovr";
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (poDriver) {
        char **papszOptions = nullptr;
        papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");
        papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
        GDALDataset *poOvrDS = poDriver->Create(ovrPath.c_str(), elev_width_, elev_height_, 1, GDT_Float32, papszOptions);
        if (poOvrDS) {
             double new_gt[6];
             new_gt[0] = elev_geo_.origin_x;
             new_gt[1] = elev_geo_.pixel_w;
             new_gt[2] = elev_geo_.rot_x;
             new_gt[3] = elev_geo_.origin_y;
             new_gt[4] = elev_geo_.rot_y;
             new_gt[5] = elev_geo_.pixel_h;
             poOvrDS->SetGeoTransform(new_gt);
             poOvrDS->SetProjection(poDS->GetProjectionRef());
             GDALRasterBand *ovBand = poOvrDS->GetRasterBand(1);
             if (hasNoData) ovBand->SetNoDataValue(noDataVal);
             if (ovBand->RasterIO(GF_Write, 0, 0, elev_width_, elev_height_, elev_data_.data(), elev_width_, elev_height_, GDT_Float32, 0, 0) != CE_None) {
                 std::cerr << "ElevationMap: failed to write data to overview file " << ovrPath << std::endl;
             }
             GDALClose(poOvrDS);
             std::cerr << "ElevationMap: saved downsampled data to " << ovrPath << std::endl;
        } else {
             std::cerr << "ElevationMap: failed to create overview file " << ovrPath << std::endl;
        }
        CSLDestroy(papszOptions);
    }

    elev_valid_ = true;
    std::cerr << "ElevationMap: in-code max downsample finished (width=" << elev_width_ << " height=" << elev_height_ << ")" << std::endl;
    return true;
#else
    return false;
#endif
}

void UavPathPlanner::printElevationInfo() const
{
  if (!elev_valid_) {
    std::cerr << "ElevationMap: no valid data loaded" << std::endl;
    return;
  }
  std::cerr << "ElevationMap: info -> width=" << elev_width_ << " height=" << elev_height_ << " resolution_x=" << elev_geo_.pixel_w << " resolution_y=" << elev_geo_.pixel_h << std::endl;
  // compute basic stats (skip NaNs)
  double minv = std::numeric_limits<double>::infinity();
  double maxv = -std::numeric_limits<double>::infinity();
  uint64_t nan_count = 0;
  size_t total = static_cast<size_t>(elev_width_) * static_cast<size_t>(elev_height_);
  for (size_t i = 0; i < total; ++i) {
    float v = elev_data_[i];
    if (std::isnan(v)) { ++nan_count; continue; }
    if (v < minv) minv = v;
    if (v > maxv) maxv = v;
  }
  if (nan_count == total) {
    std::cerr << "ElevationMap: all values are NaN" << std::endl;
  } else {
    std::cerr << "ElevationMap: data min=" << minv << " max=" << maxv << " NaN_count=" << nan_count << std::endl;
  }
  // print a few sample values
  std::cerr << "ElevationMap: sample values (first up to 8): ";
  for (size_t i = 0; i < std::min<size_t>(8, total); ++i) {
    if (i) std::cerr << ", ";
    std::cerr << elev_data_[i];
  }
  std::cerr << std::endl;
}

bool UavPathPlanner::getElevationAt(double x, double y, double &elev) const
{
  if (!elev_valid_) return false;
  // Convert world x,y to pixel indices (assuming geo_.origin is top-left)
  // pixel center at origin + (i+0.5)*pixel_w
  double px = (x - elev_geo_.origin_x) / elev_geo_.pixel_w - 0.5;
  double py = (y - elev_geo_.origin_y) / elev_geo_.pixel_h - 0.5; // pixel_h may be negative
  // Note: if pixel_h negative, dividing by negative flips sign; above handles both
  // Bilinear interpolation
  int ix = static_cast<int>(floor(px));
  int iy = static_cast<int>(floor(py));
  double fx = px - ix;
  double fy = py - iy;
  if (ix < 0 || iy < 0 || ix+1 >= elev_width_ || iy+1 >= elev_height_) return false;
  int i00 = iy * elev_width_ + ix;
  int i10 = iy * elev_width_ + (ix+1);
  int i01 = (iy+1) * elev_width_ + ix;
  int i11 = (iy+1) * elev_width_ + (ix+1);
  double v00 = elev_data_[i00];
  double v10 = elev_data_[i10];
  double v01 = elev_data_[i01];
  double v11 = elev_data_[i11];
  elev = v00*(1-fx)*(1-fy) + v10*(fx)*(1-fy) + v01*(1-fx)*(fy) + v11*(fx)*(fy);
  return true;
}

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
    if (this->getCostAt(wx, wy, celev)) {
        elev = static_cast<double>(celev);
        has_elev = true;
    } else {
        // Fallback to ElevationMap (need to convert ENU wx,wy to WGS84)
        ENUPoint enu_pt;
        enu_pt.east = wx;
        enu_pt.north = wy;
        enu_pt.up = 0.0;
        WGS84Point wgs_pt = this->enuToWGS84(enu_pt, this->origin_);
        
        if (this->getElevationAt(wgs_pt.lon, wgs_pt.lat, elev)) {
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
      
      if (this->getCostAt(wx, wy, celev)) {
          min_h = static_cast<double>(celev) + p.safe_distance;
      } else {
          double elev;
          // Fallback to ElevationMap (need to convert ENU wx,wy to WGS84)
          ENUPoint enu_pt;
          enu_pt.east = wx;
          enu_pt.north = wy;
          enu_pt.up = 0.0;
          WGS84Point wgs_pt = this->enuToWGS84(enu_pt, this->origin_);

          if (this->getElevationAt(wgs_pt.lon, wgs_pt.lat, elev)) {
              min_h = elev + p.safe_distance;
          }
      }
      
      if (min_h > -std::numeric_limits<double>::infinity() && out_z[i] < min_h) {
          out_z[i] = min_h;
      }
  }

  return true;
}

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


//get the input data of uav waypoints
// 读取JSON文件中的路径点并直接返回东北天坐标
// 从JSON对象中读取路径点并直接返回东北天坐标

std::vector<ENUPoint> UavPathPlanner::getENUFromJSON(const json& j, const std::string& key1, const std::string& key2) 
{
    std::vector<ENUPoint> enu_path;
    
    try {
        // 检查键是否存在
        if (!j.contains(key1)) {
            throw std::runtime_error("JSON中未找到键: " + key1);
        }
        if (!j.contains(key2)) {
            throw std::runtime_error("JSON中未找到键: " + key2);
        }

        auto points_array = j[key1];
        auto add_points_array = j[key2];

        double last_alt = 0.0;
        std::vector<WGS84Point> wgs84_points;
        if (!points_array.empty())         
        {// 读取所有WGS84点
            
            for (const auto& point_array : points_array) {
                if (point_array.size() >= 3) {
                    WGS84Point point;
                    point.lon = point_array[0].get<double>();
                    point.lat = point_array[1].get<double>();
                    point.alt = point_array[2].get<double>();
                    last_alt = point.alt;
                    wgs84_points.push_back(point);
                }
            }
        }
        
        // //准备区的边界点需要排序,第一个点距离上次保存的最后一个点最近,然后由方向判断顺时针还是逆时针排序
        std::vector<WGS84Point> add_wgs84_points;
        if (!add_points_array.empty()) {
            // 读取所有WGS84点
            for (const auto& add_point_array : add_points_array) {
                if (add_point_array.size() >= 2) {
                    WGS84Point point;
                    point.lon = add_point_array[0].get<double>();
                    point.lat = add_point_array[1].get<double>();
                    point.alt = last_alt;
                    add_wgs84_points.push_back(point);
                }
            }
        }
        // wgs84_points.push_back(add_wgs84_points[0]); // Close the loop
        if (!wgs84_points.empty() && !add_wgs84_points.empty()) {
            WGS84Point last_pt = wgs84_points.back();
            
            // Find closest point
            int min_idx = -1;
            double min_dist = 1.79769e+308;
            
            for (size_t i = 0; i < add_wgs84_points.size(); ++i) {
                ENUPoint enu = this->wgs84ToENU(add_wgs84_points[i], last_pt);
                double dist = enu.east * enu.east + enu.north * enu.north + enu.up * enu.up;
                if (dist < min_dist) {
                    min_dist = dist;
                    min_idx = i;
                }
            }
            
            if (min_idx != -1) {
                // Rotate to make min_idx the first element
                std::rotate(add_wgs84_points.begin(), add_wgs84_points.begin() + min_idx, add_wgs84_points.end());
                
                // Determine direction (CW or CCW)
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
                // add_wgs84_points.push_back(add_wgs84_points[0]); // Close the loop
                // add_wgs84_points.push_back(add_wgs84_points[1]); //规划保持结束点的朝向,适应循环飞行
            }
        }
        
        wgs84_points.insert(wgs84_points.end(), add_wgs84_points.begin(), add_wgs84_points.end());
        

        if (wgs84_points.empty()) {
            return enu_path;
        }
        
        // 设置ENU坐标系原点（使用第一个点作为原点）
        this->origin_ = wgs84_points[0];
        this->origin_.alt = 0;
        enu_path = this->wgs84ToENU_Batch(wgs84_points, this->origin_);
    } catch (const std::exception& e) {
        std::cerr << "从JSON获取ENU坐标错误: " << e.what() << std::endl;
    }
    return enu_path;
}

double UavPathPlanner::getDistanceFromJSON(const json& j, const std::string& key) 
{
    double distance_;
    
    try {
        // 检查键是否存在
        if (!j.contains(key)) {
            throw std::runtime_error("JSON中未找到键: " + key);
        }
        
        auto distance_json = j[key];
        if (distance_json.empty()) {
            return distance_;
        }
        distance_ = distance_json;               
    } catch (const std::exception& e) {
        std::cerr << "从JSON获取DISTANCE错误: " << e.what() << std::endl;
    }
    return distance_;
}

//

bool UavPathPlanner::getPlan(json &input_json, json &output_json, bool use3D, std::string algorithm)
{
    std::vector<ENUPoint> Enu_waypoint;
    int midwaypoint_num=0;
    int zhandoupoint_num=0;
    double distance;
    InputData input_data;
    if (!loadData(input_data, input_json))
    {
        std::cerr << "Failed to load intput json data." << std::endl;
    }
    else
    {
        std::cerr << "Successfully load intput json data." << std::endl;
    }
    //记录路径点和区域边界点数量, 用不同的路径优化方法规划
    if (input_json.contains("leader_midway_point_wgs84")) {
        midwaypoint_num = input_json["leader_midway_point_wgs84"].size();
        std::cout << "leader_midway_point_wgs84 count: " << input_json["leader_midway_point_wgs84"].size() << std::endl;
    }
    if (input_json.contains("high_zhandou_point_wgs84")) {
        zhandoupoint_num = input_json["high_zhandou_point_wgs84"].size();
        std::cout << "high_zhandou_point_wgs84 count: " << input_json["high_zhandou_point_wgs84"].size() << std::endl;
    }

    //Enu_waypoint包含了路径点和区域边界点
    Enu_waypoint = this->getENUFromJSON(input_json,"leader_midway_point_wgs84","high_zhandou_point_wgs84"); //
    
    //增加打印leader_midway_point_wgs84点数和high_zhandou_point_wgs84点数
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
                std::cerr << "getPlan: merging waypoint " << i << " to next (dist=" << dist2d << "m)\n";
            }
               }
        // filtered_wpts.push_back(Enu_waypoint.back()); // Always keep the last point
        // 将剩余的点（包括最后一个 midway 点和所有 zhandou 点）加入
        int start_idx = (midwaypoint_num > 0) ? (midwaypoint_num - 1) : 0;
        for (size_t i = start_idx; i < Enu_waypoint.size(); ++i) {
            filtered_wpts.push_back(Enu_waypoint[i]);
        }
        
        if (filtered_wpts.size() < midwaypoint_num) {
            std::cerr << "getPlan: filtered waypoints from " << midwaypoint_num 
                      << " to " << filtered_wpts.size() << " points.\n";
            Enu_waypoint = filtered_wpts;
        }
    }

    distance  = this->getDistanceFromJSON(input_json,"distance_points");
    // 支持从输入 JSON 中指定平均速度（key: "leader_speed"），兼容旧键 "V_avg" / "v_avg"
    double leader_speed = -1.0;
    try {
        if (input_json.contains("leader_speed") && !input_json["leader_speed"].empty()) {
            leader_speed = input_json["leader_speed"].get<double>();
        }
    } catch(...) { leader_speed = -1.0; }
    // 将 leader_speed 写回 input_data 以便其他函数使用
    if (leader_speed > 0) input_data.leader_speed = leader_speed;

    
    // 提取用于轨迹生成的点（排除末尾的 zhandou points）
    std::vector<ENUPoint> planning_waypoints;
    if (Enu_waypoint.size() >= zhandoupoint_num) {
        planning_waypoints.assign(Enu_waypoint.begin(), Enu_waypoint.end() - zhandoupoint_num );
    } else {
        planning_waypoints = Enu_waypoint;
    }

    if (algorithm == "minimum_snap") {
        if (use3D) {
            Trajectory_ENU = this->Minisnap_3D(planning_waypoints, distance, input_data.leader_speed);
        } else {
            Trajectory_ENU = this->Minisnap_EN(planning_waypoints, distance, input_data.leader_speed);
        }
    } else if (algorithm == "bspline") {
        std::cerr << "bspline algorithm not implemented yet." << std::endl;
        return false;
    } else if (algorithm == "bezier") {
        Trajectory_ENU = this->Bezier_3D(planning_waypoints, distance, input_data.leader_speed, input_data.min_turning_radius);
    } else {
        std::cerr << "Unknown algorithm: " << algorithm << ". Please use one of: minimum_snap, bspline, bezier" << std::endl;
        return false;
    }
    // Altitude optimization is no longer performed automatically here.
    // Use UavPathPlanner::runAltitudeOptimization(elev_file) externally to run
    // a separate altitude-only optimization step when desired.
    // std::cout << "trajectory point number:" << Trajectory_ENU.size() << std::endl;
    //东北天坐标转换为经纬高
    std::vector<WGS84Point> Trajectory_WGS84 = this->enuToWGS84_Batch(Trajectory_ENU,this->origin_);
    // 将经纬高路径写入 output_json —— leader
    this->putWGS84ToJson(output_json, "uav_leader_plane1", Trajectory_WGS84);
    // uav_plane1 用于保存编队（followers）列表，每个子数组以 follower id 开头后跟点列表
    json plane_array = json::array();

    // 计算并输出最小转弯半径
    double min_radius = this->calculateMinTurningRadius(Trajectory_ENU);
    if (min_radius > 0) {
        std::cout << "Minimum turning radius: " << min_radius << " m" << std::endl;
    } else {
        std::cout << "Minimum turning radius: N/A (Straight line or insufficient points)" << std::endl;
    }

    // 计算 Trajectory_ENU 最终点的朝向
    double final_heading = 0.0;
    if (Trajectory_ENU.size() >= 2) {
        const auto& p_last = Trajectory_ENU.back();
        const auto& p_prev = Trajectory_ENU[Trajectory_ENU.size() - 2];
        final_heading = std::atan2(p_last.north - p_prev.north, p_last.east - p_prev.east);
    } else if (!planning_waypoints.empty()) {
        // 如果轨迹点不足，尝试使用规划点的最后两个点
        if (planning_waypoints.size() >= 2) {
             const auto& p_last = planning_waypoints.back();
             const auto& p_prev = planning_waypoints[planning_waypoints.size() - 2];
             final_heading = std::atan2(p_last.north - p_prev.north, p_last.east - p_prev.east);
        }
    }
    // std::cout << "Final heading: " << final_heading << " rad (" << rad2deg(final_heading) << " deg)" << std::endl;

    // 生成并写入跟随者轨迹
    try {
        json plane_array = this->generateFollowerTrajectories(input_json, input_data, Trajectory_ENU, Trajectory_WGS84);
        output_json["uav_plane1"] = plane_array;
    } catch (const std::exception &e) {
        std::cerr << "调用 generateFollowerTrajectories 出错: " << e.what() << std::endl;
    }

    std::vector<ENUPoint> Patrol_Path; // 提升作用域以供第二段轨迹计算使用

/*                   第三段段任务区域巡逻轨迹计算               */
if(input_json["high_zhandou_point_wgs84"].size()!=0)   //存在高战斗区域点且有前序轨迹
{
    std::cout<<"开始计算第三段任务区域巡逻轨迹"<<std::endl;
    
    // 提取用于区域边界点的点（仅包含 zhandou points）
    std::vector<ENUPoint> Patrol_waypoints={};
    if (Enu_waypoint.size() >= zhandoupoint_num) {
        Patrol_waypoints.assign(Enu_waypoint.end() - zhandoupoint_num ,Enu_waypoint.end());
    } 
    Patrol_waypoints.push_back(Patrol_waypoints[0]); //闭合巡逻路径

    // 为了保证闭合处的切向连续（即第一个点和最后一个点的朝向相同），
    // 我们在末尾再多加一个点（即第二个点），使路径变为 P0 -> P1 -> ... -> P0 -> P1
    // 然后截取 P0 -> ... -> P0 部分
    if (Patrol_waypoints.size() > 2) {
        Patrol_waypoints.push_back(Patrol_waypoints[1]);
    }

    std::vector<ENUPoint> Patrol_Path_Full = Minisnap_3D(Patrol_waypoints, distance, input_data.leader_speed);
    
    if (Patrol_waypoints.size() > 2 && !Patrol_Path_Full.empty()) {
        // 目标截断点是倒数第二个航点 (P0)
        ENUPoint target_p = Patrol_waypoints[Patrol_waypoints.size() - 2];
        
        // 从后往前搜索最接近 target_p 的点
        // 限制搜索范围为后半段，避免匹配到起点的 P0
        size_t best_idx = Patrol_Path_Full.size() - 1;
        double min_dist = std::numeric_limits<double>::max();
        size_t search_start = Patrol_Path_Full.size() / 2;
        
        for (size_t i = Patrol_Path_Full.size() - 1; i >= search_start; --i) {
            double dx = Patrol_Path_Full[i].east - target_p.east;
            double dy = Patrol_Path_Full[i].north - target_p.north;
            double dz = Patrol_Path_Full[i].up - target_p.up;
            double d = dx*dx + dy*dy + dz*dz;
            
            if (d < min_dist) {
                min_dist = d;
                best_idx = i;
            }
        }
        
        // 截取到 best_idx
        if (best_idx < Patrol_Path_Full.size()) {
            Patrol_Path.assign(Patrol_Path_Full.begin(), Patrol_Path_Full.begin() + best_idx + 1);
        } else {
            Patrol_Path = Patrol_Path_Full;
        }
    } else {
        Patrol_Path = Patrol_Path_Full;
    }
    Patrol_Path.push_back(Patrol_Path[0]); //闭合巡逻路径
    std::vector<WGS84Point> Patrol_Path_WGS84 = this->enuToWGS84_Batch(Patrol_Path,this->origin_);
    this->putWGS84ToJson(output_json, "uav_leader_plane3", Patrol_Path_WGS84);

}


/*                   第二段任务区域过渡轨迹计算      input_json["high_zhandou_point_wgs84"].size()           */
if(input_json["high_zhandou_point_wgs84"].size()!=0 && !Trajectory_ENU.empty() && !Patrol_Path.empty())   //存在高战斗区域点且有前序轨迹
{
    std::cout<<"开始计算第二段任务区域过渡轨迹 (切圆切入优化)"<<std::endl;
    ENUPoint p0 = Trajectory_ENU.back(); // 起点 ENU
    double heading0 = final_heading; // 起点朝向

    // 计算并更新跟随者轨迹
    json plane_array = this->generateFollowerTrajectories(input_json, input_data, Trajectory_ENU, Trajectory_WGS84);
    output_json["uav_plane1"] = plane_array;

    // 计算第二段过渡轨迹（切圆切入优化）并更新第三段巡逻轨迹
    computeTransitionAndRotatePatrol(p0, heading0, input_data.min_turning_radius, distance, Patrol_Path, output_json);
}


    return true;
}


bool UavPathPlanner::putWGS84ToJson(json &j, const std::string &key, const std::vector<WGS84Point> &traj) {
    try {
        // 创建二维数组
        json trajectory_array = json::array();
        
        for (const auto& point : traj) {
            // 每个点创建一维数组 [经度, 纬度, 高度]
            json point_array = json::array();
            point_array.push_back(point.lon);
            point_array.push_back(point.lat);
            point_array.push_back(point.alt);
            
            trajectory_array.push_back(point_array);
        }        
        // 将轨迹数组添加到 JSON 对象的指定键
        j[key] = trajectory_array;        
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error saving WGS84 to JSON: " << e.what() << std::endl;
        return false;
    }
}

json UavPathPlanner::generateFollowerTrajectories(const json &input_json, const InputData &input_data,
                                  const std::vector<ENUPoint> &Trajectory_ENU,
                                  const std::vector<WGS84Point> &Trajectory_WGS84) {
    json plane_array = json::array();

    if (input_data.formation_using != 1) {
        return plane_array;
    }

    if (!(input_json.contains("uavs_id") && input_json.contains("uav_start_point_wgs84"))) {
        return plane_array; // empty
    }

    auto uavs_ids = input_json["uavs_id"];
    auto uav_starts = input_json["uav_start_point_wgs84"];

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

    // 获取安全距离
    double safety_distance = 50.0; // 默认值
    if (input_json.contains("safety_distance")) {
        safety_distance = input_json["safety_distance"].get<double>();
    }

    // std::cout << "Formation Model: " << input_data.formation_model << std::endl;

    if (input_data.formation_model == 1) {
        std::cout << "Formation Model: 人字形编队" << std::endl;
        plane_array = generateVShapeTrajectories(uavs_ids, uav_starts, leader_start_enu, R0, leader_xy, leader_headings, Trajectory_ENU, safety_distance);
    } else if (input_data.formation_model == 2) {
        std::cout << "Formation Model: 一字形编队" << std::endl;
        plane_array = generateLineShapeTrajectories(uavs_ids, uav_starts, leader_start_enu, R0, leader_xy, leader_headings, Trajectory_ENU, safety_distance);
    } else {
        // 默认使用人字形 (V-shape)
        std::cout << "Formation Model: 人字形编队" << std::endl;
        plane_array = generateVShapeTrajectories(uavs_ids, uav_starts, leader_start_enu, R0, leader_xy, leader_headings, Trajectory_ENU, safety_distance);
    }

    return plane_array;
}

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
            if (t == 0) {
                json pa = json::array();
                pa.push_back(start_wp.lon);
                pa.push_back(start_wp.lat);
                pa.push_back(start_wp.alt);
                follower_entry.push_back(pa);
                continue;
            }

            double heading = leader_headings[t];
            double c = cos(heading);
            double s_ = sin(heading);
            Eigen::Matrix2d Rt; Rt << c, -s_, s_, c;
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

json UavPathPlanner::generateLineShapeTrajectories(const json &uavs_ids, const json &uav_starts, 
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
                pa.push_back(start_wp.alt);
                follower_entry.push_back(pa);
                continue;
            }

            double heading = leader_headings[t];
            double c = cos(heading);
            double s_ = sin(heading);
            Eigen::Matrix2d Rt; Rt << c, -s_, s_, c;
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

    // use global generator
    // 与 Minisnap_3D 相同的配置查找策略
    std::vector<std::string> try_paths = {
        "math_util/minimum_snap_config.ymal",
        "../math_util/minimum_snap_config.ymal",
        "../../math_util/minimum_snap_config.ymal"
    };

    std::string yaml_cfg;
    bool found = false;
    for (const auto &p : try_paths) {
        std::ifstream ifs(p);
        if (ifs.good()) {
            yaml_cfg = p;
            found = true;
            break;
        }
    }
    if (!found) {
        yaml_cfg = "math_util/minimum_snap_config.ymal";
        std::cerr << "Warning: cannot find config in usual locations; will try '" << yaml_cfg << "' and fall back to defaults if unreadable." << std::endl;
    } else {
        std::cerr << "Using YAML config: " << yaml_cfg << std::endl;
    }

    if (distance_ > 0) {
        std::cerr << "Using input JSON distance_points as sampling distance: " << distance_ << " m" << std::endl;
    }
    if (v_avg_override > 0) {
        std::cerr << "Using leader_speed as average speed: " << v_avg_override << " m/s" << std::endl;
    }

    Eigen::MatrixXd sampled = this->generator_.GenerateTrajectoryMatrix(route, yaml_cfg, distance_, v_avg_override);

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

    // use global generator
    std::vector<std::string> try_paths = {
        "math_util/minimum_snap_config.ymal",
        "../math_util/minimum_snap_config.ymal",
        "../../math_util/minimum_snap_config.ymal"
    };

    std::string yaml_cfg;
    bool found = false;
    for (const auto &p : try_paths) {
        std::ifstream ifs(p);
        if (ifs.good()) {
            yaml_cfg = p;
            found = true;
            break;
        }
    }
    if (!found) {
        yaml_cfg = "math_util/minimum_snap_config.ymal";
        std::cerr << "Warning: cannot find config in usual locations; will try '" << yaml_cfg << "' and fall back to defaults if unreadable." << std::endl;
    } else {
        std::cerr << "Using YAML config: " << yaml_cfg << std::endl;
    }

    if (distance_ > 0) {
        std::cerr << "Using input JSON distance_points as sampling distance: " << distance_ << " m" << std::endl;
    }
    if (v_avg_override > 0) {
        std::cerr << "Using leader_speed as average speed: " << v_avg_override << " m/s" << std::endl;
    }

    Eigen::MatrixXd sampled = this->generator_.GenerateTrajectoryMatrix(route, yaml_cfg, distance_, v_avg_override);

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

// Bezier_3D: 贝塞尔曲线轨迹生成
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

bool UavPathPlanner::loadData(InputData &input_data, json &input_json)
{
    if (input_json["battle_high_list"].size() != 0)
    {
        auto value = input_json["battle_high_list"];
    }
    else
    {
        std::cout << "battle_high_list is empty." << std::endl;
    }
    if (input_json["battle_zone_list"].size() != 0)
    {
        auto value = input_json["battle_zone_list"];
    }
    else
    {
        std::cout << "battle_zone_list is empty." << std::endl;
    }
    if (input_json["distance_points"].size() != 0)
    {
        auto value = input_json["distance_points"];
        input_data.distance_points = value;
    }
    else
    {
        std::cout << "distance_points is empty." << std::endl;
    }
    // 允许从输入 JSON 指定平均速度 V_avg（m/s）
    if (input_json.contains("leader_speed") && !input_json["leader_speed"].empty()) {
        try {
            input_data.leader_speed = input_json["leader_speed"].get<double>();
        } catch(...) {
            // ignore parse errors, keep default
        }
    }
    if (input_json["leader_fly_high"].size() != 0)
    {
        auto value = input_json["leader_fly_high"];
        input_data.leader_fly_high = value;
    }
    else
    {
        std::cout << "leader_fly_high is empty." << std::endl;
    }
    if (input_json["formation_model"].size() != 0)
    {
        auto value = input_json["formation_model"];
        input_data.formation_model = value;
    }
    else
    {
        std::cout << "formation_model is empty." << std::endl;
    }
    if (input_json["formation_using"].size() != 0)
    {
        auto value = input_json["formation_using"];
        input_data.formation_using = value;
    }
    else
    {
        std::cout << "formation_using is empty." << std::endl;
    }
    if (input_json["uav_leader_id"].size() != 0)
    {
        auto value = input_json["uav_leader_id"];
        input_data.uav_leader_id = value[0];
    }
    else
    {
        std::cout << "uav_leader_id is empty." << std::endl;
    }
    if (input_json["leader_midway_point_wgs84"].size() != 0)
    {
        auto value = input_json["leader_midway_point_wgs84"];
        for (const auto &iter : value)
        {
            input_data.leader_midway_point_wgs84.emplace_back(WGS84Coord(iter[0], iter[1], iter[2]));
        }
    }
    else
    {
        std::cout << "leader_midway_point_wgs84 is empty." << std::endl;
    }
    if (input_json["ready_zone"].size() != 0)
    {
               auto value = input_json["ready_zone"];
        for (const auto &iter : value)
        {
            input_data.ready_zone.emplace_back(WGS84Coord(iter[0], iter[1], 0.0));
        }
    }
    else
    {
        std::cout << "ready_zone is empty." << std::endl;
    }
    if (input_json["uav_leader_start_point_wgs84"].size() != 0)
    {
        auto value = input_json["uav_leader_start_point_wgs84"];
        input_data.uav_leader_start_point_wgs84 = WGS84Coord(value[0][0], value[0][1], 0.0);
    }
    else
    {
        std::cout << "uav_leader_start_point_wgs84 is empty." << std::endl;
    }
    if (input_json["uavs_plane_data"].size() != 0)
    {
        auto value = input_json["uavs_plane_data"];
        for (const auto &iter : value)
        {
            input_data.uavs_plane_data = std::make_pair(iter[0], std::make_pair(iter[1], iter[2]));
        }
    }
    else
    {
        std::cout << "uavs_plane_data is empty." << std::endl;
    }
    if (input_json["uav_leader_id"].size() != 0)
    {
        auto value = input_json["uav_leader_id"];
    }
    else
    {
        std::cout << "uav_leader_id is empty." << std::endl;
    }
    if (input_json["ready_high_list"].size() != 0)
    {
        auto value = input_json["ready_high_list"];
        input_data.ready_high_list = std::make_pair(value[0], value[1]);
    }
    else
    {
        std::cout << "ready_high_list is empty." << std::endl;
    }
    if (input_json["high_list"].size() != 0)
    {
        auto value = input_json["high_list"];
        input_data.height_list = std::make_pair(value[0], value[1]);
    }
    else
    {
        std::cout << "high_list is empty." << std::endl;
    }
    if (input_json.contains("min_turning_radius"))
    {
        input_data.min_turning_radius = input_json["min_turning_radius"];
    }
    else
    {
        // std::cout << "min_turning_radius is empty." << std::endl;
    }
    // Load config from yaml if available to set min_turning_radius if not provided in input
    std::vector<std::string> config_paths = {"config.yaml", "../config.yaml", "../../config.yaml"};
    for (const auto& config_path : config_paths) {
        struct stat buffer;
        if (stat(config_path.c_str(), &buffer) == 0) {
             try {
                YAML::Node config = YAML::LoadFile(config_path);
                if (config["path_planning"]) {
                    auto pp_conf = config["path_planning"];
                    if (pp_conf["min_turning_radius"]) {
                        double r = pp_conf["min_turning_radius"].as<double>();
                        if (input_data.min_turning_radius <= 0.0) {
                             input_data.min_turning_radius = r;
                             std::cerr << "Loaded min_turning_radius from " << config_path << ": " << r << std::endl;
                        }
                    }
                }
             } catch (const std::exception &e) {
                 std::cerr << "Failed to parse " << config_path << ": " << e.what() << std::endl;
             }
             break; 
        }
    }

    return true;
}

// CostMap implementation merged into UavPathPlanner
void UavPathPlanner::createCostMap(int width, int height, double resolution, double origin_x, double origin_y)
{
  cost_width_ = width; 
  cost_height_ = height; 
  cost_resolution_ = resolution; 
  cost_origin_x_ = origin_x; 
  cost_origin_y_ = origin_y;
  cost_data_.assign(cost_width_ * cost_height_, -std::numeric_limits<float>::infinity());
}

void UavPathPlanner::setCost(int x, int y, float c)
{
  if (x < 0 || y < 0 || x >= cost_width_ || y >= cost_height_) return;
  cost_data_[y * cost_width_ + x] = c;
}

float UavPathPlanner::getCost(int x, int y) const
{
  if (x < 0 || y < 0 || x >= cost_width_ || y >= cost_height_) return -std::numeric_limits<float>::infinity();
  return cost_data_[y * cost_width_ + x];
}

bool UavPathPlanner::getCostAt(double x, double y, float &val) const
{
  // Assume top-left origin, x increases right, y decreases down (standard map)
  int c = static_cast<int>(floor((x - cost_origin_x_) / cost_resolution_));
  int r = static_cast<int>(floor((cost_origin_y_ - y) / cost_resolution_));
  
  if (c < 0 || c >= cost_width_ || r < 0 || r >= cost_height_) return false;
  val = cost_data_[r * cost_width_ + c];
  return true;
}

void UavPathPlanner::buildLocalENUCostMap(double margin, double resolution)//margin:边距参数,扩大局部地图覆盖范围
{
  if (!elev_valid_ || Trajectory_ENU.empty()) return;

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

  createCostMap(w, h, resolution, min_e, max_n); // Origin is top-left (min_e, max_n)

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
          if (this->getElevationAt(wgs_pt.lon, wgs_pt.lat, elev_val)) {
              setCost(c, r, static_cast<float>(elev_val));
          } else {
              // If out of bounds of elevation map, set to -inf or some default
              setCost(c, r, -std::numeric_limits<float>::infinity());
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
                                                      const std::vector<ENUPoint>& Patrol_Path, json& output_json)
{
    std::cout<<"开始计算第二段任务区域过渡轨迹 (切圆切入优化)"<<std::endl;

    // 3. 寻找最佳切入点 (遍历左右转弯及所有巡逻点)
    double best_score = std::numeric_limits<double>::max();
    size_t best_idx = 0;
    double best_arc_len = 0;
    double best_line_len = 0;
    double best_theta_end = 0;
    int best_s = 0;
    double best_cx = 0;
    double best_cy = 0;
    double best_theta_start = 0;
    bool found_any = false;

    // 尝试左右两种转弯方向，寻找最优解
    for (int s : {1, -1}) {
        // 计算切圆圆心
        // 左转(s=1): (-sin, cos); 右转(s=-1): (sin, -cos)
        double cx = p0.east - s * minR * std::sin(heading0);
        double cy = p0.north + s * minR * std::cos(heading0);
        double theta_start = std::atan2(p0.north - cy, p0.east - cx);

        for(size_t i = 0; i < Patrol_Path.size(); ++i) {
            const auto& pt = Patrol_Path[i];
            // 计算该点处的巡逻方向 (指向下一点)
            const auto& next_pt = Patrol_Path[(i + 1) % Patrol_Path.size()];
            double patrol_dx = next_pt.east - pt.east;
            double patrol_dy = next_pt.north - pt.north;
            double patrol_len = std::hypot(patrol_dx, patrol_dy);
            if(patrol_len < 1e-3) continue;
            patrol_dx /= patrol_len;
            patrol_dy /= patrol_len;

            // 计算从圆心到目标点的向量及距离
            double v_cx = pt.east - cx;
            double v_cy = pt.north - cy;
            double dist_cp = std::hypot(v_cx, v_cy);
            
            if(dist_cp <= minR) continue; // 目标点在圆内，无法做切线

            // 计算切点
            // 我们需要的是“沿圆弧运动方向”飞出的切线
            double alpha = std::atan2(v_cy, v_cx);
            double beta = std::acos(minR / dist_cp);
            
            // 两个可能的切点角度
            double candidates[2] = {alpha + beta, alpha - beta};
            for(double theta : candidates) {
                // 切点坐标
                double tx = cx + minR * std::cos(theta);
                double ty = cy + minR * std::sin(theta);
                
                // 切线向量 (T -> P)
                double lx = pt.east - tx;
                double ly = pt.north - ty;
                double l_len = std::hypot(lx, ly);
                if(l_len < 1e-3) continue;
                double l_dx = lx / l_len;
                double l_dy = ly / l_len;

                // 圆在该切点处的切线方向 (运动方向)
                double tan_x = -s * std::sin(theta);
                double tan_y = s * std::cos(theta);

                // 检查1: 直线方向必须与圆的切线方向一致 (cos theta > 0)
                if(tan_x * l_dx + tan_y * l_dy < 0.99) continue; 

                // 检查2: 切入角度与巡逻方向的对齐程度
                double alignment = l_dx * patrol_dx + l_dy * patrol_dy;
                // 严格限制切入角度，例如必须小于 36 度 (cos > 0.8)
                if(alignment < 0.8) continue;

                // 计算弧长
                double d_theta = theta - theta_start;
                if (s > 0) { // Left, CCW
                    while (d_theta <= 0) d_theta += 2 * M_PI;
                    while (d_theta > 2 * M_PI) d_theta -= 2 * M_PI;
                } else { // Right, CW
                    while (d_theta >= 0) d_theta -= 2 * M_PI;
                    while (d_theta < -2 * M_PI) d_theta += 2 * M_PI;
                }
                double arc_len = std::abs(d_theta) * minR;
                
                // 评分: 路径长度 + 角度惩罚
                // 角度惩罚权重: 使得算法更倾向于切向切入 (alignment -> 1)
                // 假设 1000m 的惩罚系数，意味着 10度偏差 (cos~0.985, 1-cos~0.015) => 15m 惩罚
                double penalty = 1000.0 * (1.0 - alignment);
                double total_cost = arc_len + l_len + penalty;
                
                if(total_cost < best_score) {
                    best_score = total_cost;
                    best_idx = i;
                    best_arc_len = arc_len;
                    best_line_len = l_len;
                    best_theta_end = theta;
                    best_s = s;
                    best_cx = cx;
                    best_cy = cy;
                    best_theta_start = theta_start;
                    found_any = true;
                }
            }
        }
    }

    std::vector<ENUPoint> Transition_Path;
    if(found_any) {
        // 生成圆弧部分
        int steps_arc = std::max(1, (int)std::ceil(best_arc_len / resolution));
        double d_theta_total = (best_s > 0) ? (best_arc_len / minR) : -(best_arc_len / minR);
        
        for(int i=0; i<=steps_arc; ++i) {
            double t = (double)i / steps_arc;
            double ang = best_theta_start + d_theta_total * t;
            ENUPoint pt;
            pt.east = best_cx + minR * std::cos(ang);
            pt.north = best_cy + minR * std::sin(ang);
            // 高度线性插值
            pt.up = p0.up + (Patrol_Path[best_idx].up - p0.up) * (t * best_arc_len / (best_arc_len + best_line_len));
            Transition_Path.push_back(pt);
        }

        // 生成直线部分
        ENUPoint T_end = Transition_Path.back();
        ENUPoint P_target = Patrol_Path[best_idx];
        int steps_line = std::max(1, (int)std::ceil(best_line_len / resolution));
        for(int i=1; i<=steps_line; ++i) {
            double t = (double)i / steps_line;
            ENUPoint pt;
            pt.east = T_end.east + t * (P_target.east - T_end.east);
            pt.north = T_end.north + t * (P_target.north - T_end.north);
            pt.up = T_end.up + t * (P_target.up - T_end.up);
            Transition_Path.push_back(pt);
        }
        
        // 关键步骤: 重新调整巡逻路径 (Plane 3) 的起点，使其从 best_idx 开始
        // 这样无人机飞完 Plane 2 后能无缝衔接 Plane 3
        std::vector<ENUPoint> Rotated_Patrol;
        Rotated_Patrol.reserve(Patrol_Path.size());
        for(size_t i=0; i<Patrol_Path.size(); ++i) {
            Rotated_Patrol.push_back(Patrol_Path[(best_idx + i) % Patrol_Path.size()]);
        }
        // 闭合
        Rotated_Patrol.push_back(Rotated_Patrol[0]);
        
        std::vector<WGS84Point> Rotated_Patrol_WGS84 = this->enuToWGS84_Batch(Rotated_Patrol, this->origin_);
        this->putWGS84ToJson(output_json, "uav_leader_plane3", Rotated_Patrol_WGS84);
        
    } else {
        std::cerr << "Warning: Failed to find valid tangent transition, falling back to straight line." << std::endl;
        // Fallback: 直接连线到最近点
        ENUPoint p1 = Patrol_Path[0];
        double dx = p1.east - p0.east;
        double dy = p1.north - p0.north;
        double dist = std::hypot(dx, dy);
        int steps = std::max(1, (int)std::ceil(dist / resolution));
        for (int i = 0; i <= steps; ++i) {
            double t = (double)i / steps;
            ENUPoint q{p0.east + t * dx, p0.north + t * dy, p0.up + t * (p1.up - p0.up)};
            Transition_Path.push_back(q);
        }
    }

    std::vector<WGS84Point> Transition_Path_WGS84 = this->enuToWGS84_Batch(Transition_Path,this->origin_);
    this->putWGS84ToJson(output_json, "uav_leader_plane2", Transition_Path_WGS84);
}
