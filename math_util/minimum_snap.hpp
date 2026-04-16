#ifndef _MINIMUM_SNAP_HPP_
#define _MINIMUM_SNAP_HPP_

#include <Eigen/Dense>
#include <vector>
#include <string>

// minimum-snap 配置（合并自 minimum_snap_config.ymal / config.yaml.minimum_snap）
struct MinimumSnapConfig {
    // 多项式导数阶数（d_order），min-snap 常用 order=3
    int order = 3;

    // 偏离直线段惩罚权重
    double path_weight = 0.0;

    // 每个路径点的 0 速度约束权重（经过路径点减速）
    double vel_zero_weight = 0.0;

    // 平均巡航速度 (m/s)
    double V_avg = 5.0;

    // 单段最小时长 (s)
    double min_time_s = 0.1;

    // 轨迹点采样间隔距离 (m)
    double sample_distance = 1.0;

    // 起始/终止速度与加速度（可选）
    Eigen::Vector3d start_vel = Eigen::Vector3d::Zero();
    Eigen::Vector3d end_vel = Eigen::Vector3d::Zero();
    Eigen::Vector3d start_acc = Eigen::Vector3d::Zero();
    Eigen::Vector3d end_acc = Eigen::Vector3d::Zero();
};

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
    // cfg: 最小 snap 配置（已从 config.yaml 读入）
    // 如果 sample_distance_override > 0 则覆盖 cfg.sample_distance 值（单位：米）
    // v_avg_override: if >0 overrides cfg.V_avg (average speed in m/s)
    Eigen::MatrixXd GenerateTrajectoryMatrix(const Eigen::MatrixXd &Path, const MinimumSnapConfig &cfg,
                                            double sample_distance_override = -1.0, double v_avg_override = -1.0);

};



#endif