#include "uavPathPlanning.hpp"
#include "math_util/minimum_snap.hpp"
#include "math_util/altitude_optimizer.hpp"
#include <yaml-cpp/yaml.h>

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

// UavPathPlanner::runAltitudeOptimization: 现在只接受一个高程文件路径参数 (.tif 或 PGM)
bool UavPathPlanner::runAltitudeOptimization(const std::string &elev_file)
{
    try {
        if (elev_file.empty()) return false;
        if (this->Trajectory_ENU.empty()) {
            std::cerr << "runAltitudeOptimization: Trajectory_ENU is empty, nothing to optimize\n";
            return false;
        }

        math_util::AltitudeOptimizer altopt;
        if (!altopt.load_elevation_data(elev_file)) {
            std::cerr << "AltitudeOptimizer: failed to load elevation file: " << elev_file << "\n";
            return false;
        }

        const auto &emap = altopt.getElevationMap();

        // compute elev min/max (naive scan)
        double elev_min = 1e300, elev_max = -1e300;
        for (int r = 0; r < emap.getHeight(); ++r) {
            for (int c = 0; c < emap.getWidth(); ++c) {
                double wx = emap.getGeoTransform().origin_x + (c + 0.5) * emap.getGeoTransform().pixel_w;
                double wy = emap.getGeoTransform().origin_y + (r + 0.5) * emap.getGeoTransform().pixel_h;
                double ev;
                if (emap.getElevationAt(wx, wy, ev)) {
                    elev_min = std::min(elev_min, ev);
                    elev_max = std::max(elev_max, ev);
                }
            }
        }
        if (elev_min > elev_max) return false;

        math_util::CostMap cmap;
        cmap.fromElevationMap(emap, elev_min, elev_max);

        // build waypoint vector for optimizer from member Trajectory_ENU
        std::vector<Eigen::Vector3d> wpts;
        wpts.reserve(this->Trajectory_ENU.size());
        for (size_t i = 0; i < this->Trajectory_ENU.size(); ++i) {
            const auto &p = this->Trajectory_ENU[i];
            wpts.emplace_back(p.east, p.north, p.up);
        }

        math_util::AltitudeOptimizer::Params params; // use defaults

        std::vector<double> out_z;
        if (altopt.optimizeHeights(wpts, emap, cmap, params, out_z)) {
            for (size_t i = 0; i < out_z.size() && i < this->Trajectory_ENU.size(); ++i) {
                this->Trajectory_ENU[i].up = out_z[i];
            }
            std::cerr << "AltitudeOptimizer: optimized heights applied (" << out_z.size() << " points)\n";
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

//get the input data of uav waypoints
// 读取JSON文件中的路径点并直接返回东北天坐标
// 从JSON对象中读取路径点并直接返回东北天坐标

std::vector<ENUPoint> UavPathPlanner::getENUFromJSON(const json& j, const std::string& key) 
{
    std::vector<ENUPoint> enu_path;
    
    try {
        // 检查键是否存在
        if (!j.contains(key)) {
            throw std::runtime_error("JSON中未找到键: " + key);
        }
        
        auto points_array = j[key];
        if (points_array.empty()) {
            return enu_path;
        }
        
        // 读取所有WGS84点
        std::vector<WGS84Point> wgs84_points;
        for (const auto& point_array : points_array) {
            if (point_array.size() >= 3) {
                WGS84Point point;
                point.lon = point_array[0].get<double>();
                point.lat = point_array[1].get<double>();
                point.alt = point_array[2].get<double>();
                wgs84_points.push_back(point);
            }
        }
        
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

bool UavPathPlanner::getPlan(json &input_json, json &output_json, bool use3D)
{
    std::vector<ENUPoint> Enu_waypoint;
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
    Enu_waypoint = this->getENUFromJSON(input_json,"leader_midway_point_wgs84"); //
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

    
    if (use3D) {
        Trajectory_ENU = this->Minisnap_3D(Enu_waypoint, distance, input_data.leader_speed);
    } else {
        Trajectory_ENU = this->Minisnap_EN(Enu_waypoint, distance, input_data.leader_speed);
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

    // 生成并写入跟随者轨迹
    try {
        json plane_array = this->generateFollowerTrajectories(input_json, input_data, Trajectory_ENU, Trajectory_WGS84);
        output_json["uav_plane1"] = plane_array;
    } catch (const std::exception &e) {
        std::cerr << "调用 generateFollowerTrajectories 出错: " << e.what() << std::endl;
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

    for (size_t idx = 0; idx < uavs_ids.size(); ++idx) {
        int uid = uavs_ids[idx].get<int>();
        if (idx >= uav_starts.size()) break;
        auto s = uav_starts[idx];
        WGS84Point follower_start_wgs{s[0].get<double>(), s[1].get<double>(), s.size()>=3 ? s[2].get<double>() : 0.0};
    ENUPoint follower_start_enu = this->wgs84ToENU(follower_start_wgs, this->origin_);

        Eigen::Vector2d rel_global(follower_start_enu.east - leader_start_enu.east,
                                   follower_start_enu.north - leader_start_enu.north);
        double rel_up = follower_start_enu.up - leader_start_enu.up;
        Eigen::Vector2d rel_body = R0.transpose() * rel_global;

        // 构造 follower 的 entry
        json follower_entry = json::array();
        follower_entry.push_back(uid);

        for (size_t t = 0; t < N; ++t) {
            double heading = 0.0;
            if (t + 1 < N) {
                Eigen::Vector2d diff = leader_xy[t+1] - leader_xy[t];
                heading = atan2(diff.y(), diff.x());
            } else if (t > 0) {
                Eigen::Vector2d diff = leader_xy[t] - leader_xy[t-1];
                heading = atan2(diff.y(), diff.x());
            } else {
                heading = leader_initial_heading;
            }
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
    return true;
}
