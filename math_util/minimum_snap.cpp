
#include "minimum_snap.hpp"
#include <iostream>
#include <yaml-cpp/yaml.h>

using namespace std;
using namespace Eigen;

/*!
 * 计算x的阶乘
 * @param x
 * @return x!
 */
int TrajectoryGeneratorTool::Factorial(int x) {
    int fac = 1;
    for (int i = x; i > 0; i--)
        fac = fac * i;
    return fac;
}

Eigen::MatrixXd TrajectoryGeneratorTool::GenerateTrajectoryMatrix(const Eigen::MatrixXd &Path, const std::string &yaml_path, double sample_distance_override, double v_avg_override) {
    // 配置默认值
    int order = 3; // 默认最小化 snap 对应 d_order=3 (可在yaml中覆盖)
    double V_avg = 5.0; // m/s
    double min_time_s = 0.1;
    double sample_distance = 1.0; // 轨迹采样距离 m

    // 起始/终止速度与加速度（可在yaml中配置为数组 [x,y,z]）
    Eigen::MatrixXd Vel = Eigen::MatrixXd::Zero(2, 3);
    Eigen::MatrixXd Acc = Eigen::MatrixXd::Zero(2, 3);

    // 读取 YAML 配置（若文件存在且字段可用则覆盖默认值）
    try {
        YAML::Node cfg = YAML::LoadFile(yaml_path);
        if (cfg["order"]) order = cfg["order"].as<int>();
        if (cfg["V_avg"]) V_avg = cfg["V_avg"].as<double>();
        if (cfg["min_time_s"]) min_time_s = cfg["min_time_s"].as<double>();
        if (cfg["sample_distance"]) sample_distance = cfg["sample_distance"].as<double>();

        if (cfg["start_vel"] && cfg["start_vel"].IsSequence() && cfg["start_vel"].size()>=3) {
            Vel(0,0) = cfg["start_vel"][0].as<double>();
            Vel(0,1) = cfg["start_vel"][1].as<double>();
            Vel(0,2) = cfg["start_vel"][2].as<double>();
        }
        if (cfg["end_vel"] && cfg["end_vel"].IsSequence() && cfg["end_vel"].size()>=3) {
            Vel(1,0) = cfg["end_vel"][0].as<double>();
            Vel(1,1) = cfg["end_vel"][1].as<double>();
            Vel(1,2) = cfg["end_vel"][2].as<double>();
        }
        if (cfg["start_acc"] && cfg["start_acc"].IsSequence() && cfg["start_acc"].size()>=3) {
            Acc(0,0) = cfg["start_acc"][0].as<double>();
            Acc(0,1) = cfg["start_acc"][1].as<double>();
            Acc(0,2) = cfg["start_acc"][2].as<double>();
        }
        if (cfg["end_acc"] && cfg["end_acc"].IsSequence() && cfg["end_acc"].size()>=3) {
            Acc(1,0) = cfg["end_acc"][0].as<double>();
            Acc(1,1) = cfg["end_acc"][1].as<double>();
            Acc(1,2) = cfg["end_acc"][2].as<double>();
        }
    } catch (const std::exception &e) {
        std::cerr << "TrajectoryGeneratorTool: failed to load YAML config '"<< yaml_path << "': " << e.what() << ". Using defaults." << std::endl;
    }

    // 如果调用者传入了覆盖采样距离的参数，优先使用它
    if (sample_distance_override > 0.0) {
        sample_distance = sample_distance_override;
    }
    // 如果调用者传入了覆盖平均速度参数，优先使用它
    if (v_avg_override > 0.0) {
        V_avg = v_avg_override;
    }

    // 打印实际使用的采样距离，便于调试
    // std::cerr << "TrajectoryGeneratorTool: sampling distance used = " << sample_distance << " m" << std::endl;

    // 输入路径检查
    if (Path.rows() < 2 || Path.cols() < 3) {
        std::cerr << "TrajectoryGeneratorTool::GenerateTrajectoryMatrix: Path must be (N>=2 x 3)" << std::endl;
        return Eigen::MatrixXd();
    }

    const int num_points = Path.rows();
    const int num_segments = num_points - 1;

    // 计算每段长度并分配时间（简单：长度 / V_avg，最小为 min_time_s）
    Eigen::VectorXd Time(num_segments);
    for (int i = 0; i < num_segments; ++i) {
        double dx = Path(i+1,0) - Path(i,0);
        double dy = Path(i+1,1) - Path(i,1);
        double dz = Path(i+1,2) - Path(i,2);
        double len = std::sqrt(dx*dx + dy*dy + dz*dz);
        double t = (V_avg > 1e-6) ? (len / V_avg) : min_time_s;
        if (t < min_time_s) t = min_time_s;
        Time(i) = t;
    }

    // 调用闭式解求解多项式系数
    Eigen::MatrixXd polyCoeff = SolveQPClosedForm(order, Path, Vel, Acc, Time);

    // 多项式次数与每段系数个数
    const int p_order = 2 * order - 1;
    const int p_num1d = p_order + 1;

    // 采样：精确的固定距离采样
    // 使用一个小的时间步长步进，但在每个步进内通过线性插值多次插入满足 sample_distance 的点，保证记录点间距≈sample_distance
    // double dt_default = (V_avg > 1e-6) ? std::max(0.01, sample_distance / V_avg) : 0.01;
    double dt_default = 0.1;
    std::vector<Eigen::Vector3d> samples;
    samples.reserve(1000); //预分配内存,减少动态扩展次数

    auto eval_poly_at = [&](int seg, double t)->Eigen::Vector3d {
        Eigen::Vector3d pt(0,0,0);
        for (int dim = 0; dim < 3; ++dim) {
            Eigen::VectorXd coeff = polyCoeff.row(seg).segment(dim * p_num1d, p_num1d);
            double val = 0.0;
            for (int k = 0; k < p_num1d; ++k) {
                double c = coeff(k);
                int exp = p_num1d - 1 - k;
                val += c * std::pow(t, exp);
            }
            pt(dim) = val;
        }
        return pt;
    };

    Eigen::Vector3d last_recorded(0,0,0);
    bool has_last = false;
    Eigen::Vector3d prev_pt(0,0,0);

    for (int seg = 0; seg < num_segments; ++seg) {
        double T = Time(seg);
        double dt = dt_default;
        if (dt > T/10.0) dt = T/10.0;    //每段至少取10个点

        // 记录本段开始前已有的采样点数，用于统计本段新增的采样点
        size_t samples_before = samples.size();

        // 在段起始点处进行一次评估以获得连续性（但不在非首段重复写入）
        Eigen::Vector3d t0_pt = eval_poly_at(seg, 0.0);
        if (!has_last) {
            samples.push_back(t0_pt);
            last_recorded = t0_pt;
            has_last = true;
        }
        prev_pt = t0_pt;

        for (double t = dt; t <= T + 1e-12; t += dt) {
            double tt = std::min(t, T);
            Eigen::Vector3d cur_pt = eval_poly_at(seg, tt);
            Eigen::Vector3d seg_vec = cur_pt - prev_pt;
            double seg_len = seg_vec.norm();
            if(seg_len>=sample_distance)

            // 将 prev_pt 移动到 cur_pt 以供下一步使用
            {
                prev_pt = cur_pt;
                samples.push_back(cur_pt);
            }
        }
    // 输出本段统计信息：段索引与本段新增采样点数量
    size_t samples_added = samples.size() - samples_before;
    std::cout << "segment: " << seg << " samples_added: " << samples_added << std::endl;
        // 在最后一段时，确保终点被加入（避免重复）
        if (seg == num_segments - 1) {
            Eigen::Vector3d endpt = eval_poly_at(seg, T);
            if (samples.empty() || (samples.back() - endpt).norm() > 1e-6) samples.push_back(endpt);
        }
    }

    // 转换为矩阵输出
    Eigen::MatrixXd out(samples.size(), 3);
    for (size_t i = 0; i < samples.size(); ++i) {
        out(i,0) = samples[i](0);
        out(i,1) = samples[i](1);
        out(i,2) = samples[i](2);
    }

    return out;
}

/*!
 * 通过闭式求解QP，得到每段拟合轨迹的多项式系数
 * @param order 导数阶数。例如最小化jerk，则需要求解三次导数，则 d_order=3
 * @param Path 航迹点的空间坐标(3D)
 * @param Vel 航迹点对应的速度(中间点速度是待求的未知量)
 * @param Acc 航迹点对应的加速度(中间点加速度是待求的未知量)
 * @param Time 每段轨迹对应的时间周期
 * @return 轨迹x,y,z三个方向上的多项式系数
 *
 * 返回矩阵(PolyCoeff)的数据格式：每一行是一段轨迹，第一列是x方向上的多项式次数，越左次数越高
 * 第一段轨迹三个方向上的系数 | px_i px_(i-1) px_(i-2) ... px_1 px_0 | y ... | z ... |
 * 第二段轨迹三个方向上的系数 | px_i px_(i-1) px_(i-2) ... px_1 px_0 | y ... | z ... |
 *                          ........
 *
 * 注意：给定起始点和终点的速度加速度，更高阶的导数设置为0
 */
Eigen::MatrixXd TrajectoryGeneratorTool::SolveQPClosedForm(
        int order,
        const Eigen::MatrixXd &Path,
        const Eigen::MatrixXd &Vel,
        const Eigen::MatrixXd &Acc,
        const Eigen::VectorXd &Time) {

    const int p_order = 2 * order - 1;//多项式的最高次数 p^(p_order)t^(p_order) + ...
    const int p_num1d = p_order + 1;//每一段轨迹的变量个数，对于五阶多项式为：p5, p4, ... p0
    std::cout << "input points number : " << Path.rows() << std::endl;
    const int number_segments = Time.size();
    //每一段都有x,y,z三个方向，每一段多项式的系数的个数有3*p_num1d
    MatrixXd PolyCoeff = MatrixXd::Zero(number_segments, 3 * p_num1d);
    //整条轨迹在ｘ,y,z方向上共多少个未知系数
    const int number_coefficients = p_num1d * number_segments;
    VectorXd Px(number_coefficients), Py(number_coefficients), Pz(number_coefficients);

    const int M_block_rows = order * 2;
    const int M_block_cols = p_num1d;
    //M：转换矩阵，将系数向量转换为方程的微分量
    MatrixXd M = MatrixXd::Zero(number_segments * M_block_rows, number_segments * M_block_cols);
    for (int i = 0; i < number_segments; ++i) {
        int row = i * M_block_rows, col = i * M_block_cols;
        MatrixXd sub_M = MatrixXd::Zero(M_block_rows, M_block_cols);

        for (int j = 0; j < order; ++j) {
            for (int k = 0; k < p_num1d; ++k) {
                if (k < j)
                    continue;

                sub_M(j, p_num1d - 1 - k) = Factorial(k) / Factorial(k - j) * pow(0, k - j);
                sub_M(j + order, p_num1d - 1 - k) = Factorial(k) / Factorial(k - j) * pow(Time(i), k - j);
            }
        }

        M.block(row, col, M_block_rows, M_block_cols) = sub_M;
    }

    //构造选择矩阵C的过程非常复杂，但是只要多花点时间探索一些规律，举几个例子，应该是能写出来的!!
    const int number_valid_variables = (number_segments + 1) * order;
    const int number_fixed_variables = 2 * order + (number_segments - 1);
    //C_T：选择矩阵，用于分离未知量和已知量
    MatrixXd C_T = MatrixXd::Zero(number_coefficients, number_valid_variables);
    for (int i = 0; i < number_coefficients; ++i) {
        if (i < order) {
            C_T(i, i) = 1;
            continue;
        }

        if (i >= number_coefficients - order) {
            const int delta_index = i - (number_coefficients - order);
            C_T(i, number_fixed_variables - order + delta_index) = 1;
            continue;
        }

        if ((i % order == 0) && (i / order % 2 == 1)) {
            const int index = i / (2 * order) + order;
            C_T(i, index) = 1;
            continue;
        }

        if ((i % order == 0) && (i / order % 2 == 0)) {
            const int index = i / (2 * order) + order - 1;
            C_T(i, index) = 1;
            continue;
        }

        if ((i % order != 0) && (i / order % 2 == 1)) {
            const int temp_index_0 = i / (2 * order) * (2 * order) + order;
            const int temp_index_1 = i / (2 * order) * (order - 1) + i - temp_index_0 - 1;
            C_T(i, number_fixed_variables + temp_index_1) = 1;
            continue;
        }

        if ((i % order != 0) && (i / order % 2 == 0)) {
            const int temp_index_0 = (i - order) / (2 * order) * (2 * order) + order;
            const int temp_index_1 = (i - order) / (2 * order) * (order - 1) + (i - order) - temp_index_0 - 1;
            C_T(i, number_fixed_variables + temp_index_1) = 1;
            continue;
        }
    }

    // Q：二项式的系数矩阵
    MatrixXd Q = MatrixXd::Zero(number_coefficients, number_coefficients);
    for (int k = 0; k < number_segments; ++k) {
        MatrixXd sub_Q = MatrixXd::Zero(p_num1d, p_num1d);
        for (int i = 0; i <= p_order; ++i) {
            for (int l = 0; l <= p_order; ++l) {
                if (p_num1d - i <= order || p_num1d - l <= order)
                    continue;

                sub_Q(i, l) = (Factorial(p_order - i) / Factorial(p_order - order - i)) *
                              (Factorial(p_order - l) / Factorial(p_order - order - l)) /
                              (p_order - i + p_order - l - (2 * order - 1)) *
                              pow(Time(k), p_order - i + p_order - l - (2 * order - 1));
            }
        }

        const int row = k * p_num1d;
        Q.block(row, row, p_num1d, p_num1d) = sub_Q;
    }

    MatrixXd R = C_T.transpose() * M.transpose().inverse() * Q * M.inverse() * C_T;

    for (int axis = 0; axis < 3; ++axis) {
        VectorXd d_selected = VectorXd::Zero(number_valid_variables);
        for (int i = 0; i < number_coefficients; ++i) {
            if (i == 0) {
                d_selected(i) = Path(0, axis);
                continue;
            }

            if (i == 1 && order >= 2) {
                d_selected(i) = Vel(0, axis);
                continue;
            }

            if (i == 2 && order >= 3) {
                d_selected(i) = Acc(0, axis);
                continue;
            }

            if (i == number_coefficients - order + 2 && order >= 3) {
                d_selected(number_fixed_variables - order + 2) = Acc(1, axis);
                continue;
            }

            if (i == number_coefficients - order + 1 && order >= 2) {
                d_selected(number_fixed_variables - order + 1) = Vel(1, axis);
                continue;
            }

            if (i == number_coefficients - order) {
                d_selected(number_fixed_variables - order) = Path(number_segments, axis);
                continue;
            }

            if ((i % order == 0) && (i / order % 2 == 0)) {
                const int index = i / (2 * order) + order - 1;
                d_selected(index) = Path(i / (2 * order), axis);
                continue;
            }
        }

        MatrixXd R_PP = R.block(number_fixed_variables, number_fixed_variables,
                                number_valid_variables - number_fixed_variables,
                                number_valid_variables - number_fixed_variables);
        VectorXd d_F = d_selected.head(number_fixed_variables);
        MatrixXd R_FP = R.block(0, number_fixed_variables, number_fixed_variables,
                                number_valid_variables - number_fixed_variables);

        MatrixXd d_optimal = -R_PP.inverse() * R_FP.transpose() * d_F;

        d_selected.tail(number_valid_variables - number_fixed_variables) = d_optimal;
        VectorXd d = C_T * d_selected;

        if (axis == 0)
            Px = M.inverse() * d;

        if (axis == 1)
            Py = M.inverse() * d;

        if (axis == 2)
            Pz = M.inverse() * d;
    }

    for (int i = 0; i < number_segments; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (j == 0) {
                PolyCoeff.block(i, j * p_num1d, 1, p_num1d) =
                        Px.block(p_num1d * i, 0, p_num1d, 1).transpose();
                continue;
            }

            if (j == 1) {
                PolyCoeff.block(i, j * p_num1d, 1, p_num1d) =
                        Py.block(p_num1d * i, 0, p_num1d, 1).transpose();
                continue;
            }

            if (j == 2) {
                PolyCoeff.block(i, j * p_num1d, 1, p_num1d) =
                        Pz.block(p_num1d * i, 0, p_num1d, 1).transpose();
                continue;
            }
        }
    }

    return PolyCoeff;
}

