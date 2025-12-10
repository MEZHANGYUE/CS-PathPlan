#ifndef _MINIMUM_SNAP_HPP_
#define _MINIMUM_SNAP_HPP_

#include <Eigen/Dense>
#include <vector>
#include <string>

// Trajectory generator for minimum-snap / minimum-jerk closed-form solution
class TrajectoryGeneratorTool {
private:
    int Factorial(int x);
    // 读取是否启用路径偏差惩罚（使每段收敛到直线）的权重
    double path_weight = 0.0;
public:
    TrajectoryGeneratorTool() = default;
    ~TrajectoryGeneratorTool() = default;

    Eigen::MatrixXd SolveQPClosedForm(
            int order,
            const Eigen::MatrixXd &Path,
            const Eigen::MatrixXd &Vel,
            const Eigen::MatrixXd &Acc,
        const Eigen::VectorXd &Time,
        double path_weight = 0.0,
        double vel_zero_weight = 0.0,
        double *max_deviation = nullptr);

    // 生成完整轨迹并返回采样点矩阵 (N x 3)。
    // Path: (M x 3) 路径点矩阵
    // yaml_path: 配置文件路径，包含 order, V_avg, min_time_s, sample_distance, start/end vel/acc
    // 如果 sample_distance_override > 0 则覆盖 YAML 中的 sample_distance 值（单位：米）
    // sample_distance_override: if >0 overrides YAML sample_distance
    // v_avg_override: if >0 overrides YAML V_avg (average speed in m/s)
    Eigen::MatrixXd GenerateTrajectoryMatrix(const Eigen::MatrixXd &Path, const std::string &yaml_path, double sample_distance_override = -1.0, double v_avg_override = -1.0);
};



#endif