#ifndef _MINIMUM_SNAP_HPP_
#define _MINIMUM_SNAP_HPP_

#include <Eigen/Dense>
#include <vector>
#include <string>

// Trajectory generator for minimum-snap / minimum-jerk closed-form solution
class TrajectoryGeneratorTool {
private:
    int Factorial(int x);

public:
    TrajectoryGeneratorTool() = default;
    ~TrajectoryGeneratorTool() = default;

    Eigen::MatrixXd SolveQPClosedForm(
            int order,
            const Eigen::MatrixXd &Path,
            const Eigen::MatrixXd &Vel,
            const Eigen::MatrixXd &Acc,
            const Eigen::VectorXd &Time);

    // 生成完整轨迹并返回采样点矩阵 (N x 3)。
    // Path: (M x 3) 路径点矩阵
    // yaml_path: 配置文件路径，包含 order, V_avg, min_time_s, sample_distance, start/end vel/acc
    // 如果 sample_distance_override > 0 则覆盖 YAML 中的 sample_distance 值（单位：米）
    Eigen::MatrixXd GenerateTrajectoryMatrix(const Eigen::MatrixXd &Path, const std::string &yaml_path, double sample_distance_override = -1.0);
};

// --- Planner public API merged here for convenience ---

// 纯C++的Pose和Path结构体定义（用于兼容旧接口）
struct Pose
{
    double position_x;
    double position_y;
    double position_z;
    Pose(double x = 0.0, double y = 0.0, double z = 0.0)
        : position_x(x), position_y(y), position_z(z) {}
};

struct Path
{
    std::vector<Pose> poses;
    void clear() { poses.clear(); }
    size_t size() const { return poses.size(); }
    void push_back(const Pose& pose) { poses.push_back(pose); }
};

class planner {
private:
    int dot_num;
    int poly_coeff_num;
    int mode = 4; // 或者从配置中获取
    Eigen::MatrixXd route;
    Eigen::VectorXd time;
    Eigen::MatrixXd poly_coeff;

public:
    void getparam(void);
    Eigen::MatrixXd getcoeff(void);
    Eigen::MatrixXd getcoeff(Eigen::MatrixXd route_,Eigen::VectorXd time_);
    Eigen::Vector3d getPosPoly(Eigen::MatrixXd polyCoeff, int k, double t);
    Path trajectory_path(void);
    Path trajectory_path(Eigen::MatrixXd route_,Eigen::VectorXd time_); //1秒采样一个路径点
    Path trajectory_path(Eigen::MatrixXd route_,Eigen::VectorXd time_,double dis); //dis距离采样点
    Path trajectory_{};
};

#endif