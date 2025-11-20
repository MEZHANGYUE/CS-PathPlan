#include <vector>
#include <Eigen/Dense>
#include "planner.h"



/*
*获取配置文件参数 - 需要重写为从文件读取或其他配置源
*/
void planner::getparam(void)
{
    // 这里需要根据你的实际配置源重写
    // 示例硬编码数据，实际应用中应该从文件读取
    dot_num = 4; // 示例：4个路径点
    
    route.resize(dot_num, 3);
    // 示例路径点数据
    route << 0.0, 0.0, 0.0,
             1.0, 0.0, 0.0,
             1.0, 1.0, 0.0,
             0.0, 1.0, 0.0;
    
    time.resize(3); // 示例：3个时间段
    time << 2.0, 2.0, 2.0; // 每个段2秒
}

/*
*获取路径参数矩阵
*/
//从文件中获取参数,见getparam();
Eigen::MatrixXd planner::getcoeff(void)
{
    Eigen::MatrixXd polycoeff;
    Eigen::MatrixXd vel = Eigen::MatrixXd::Zero(2, 3);
    Eigen::MatrixXd acc = Eigen::MatrixXd::Zero(2, 3);
    TrajectoryGeneratorTool TrajectoryGeneratorTool;
    getparam();
    polycoeff = TrajectoryGeneratorTool.SolveQPClosedForm(mode, route, vel, acc, time);
    return polycoeff; 
}
//传入路径和时间分配参数
Eigen::MatrixXd planner::getcoeff(Eigen::MatrixXd route_,Eigen::VectorXd time_)
{
    Eigen::MatrixXd polycoeff;
    Eigen::MatrixXd vel = Eigen::MatrixXd::Zero(2, 3);
    Eigen::MatrixXd acc = Eigen::MatrixXd::Zero(2, 3);
    TrajectoryGeneratorTool TrajectoryGeneratorTool;
    //路径和时间采用传入的参数
    polycoeff = TrajectoryGeneratorTool.SolveQPClosedForm(mode, route_, vel, acc, time_);
    return polycoeff; 
}
/*!
 * 求解第k个轨迹段t时刻对应的位置
 * @param polyCoeff 多项式系数矩阵
 * @param k 轨迹段序号
 * @param t 时刻
 * @return [x,y,z]^T
 */
Eigen::Vector3d planner::getPosPoly(Eigen::MatrixXd polyCoeff, int k, double t) 
{
    Eigen::Vector3d pt;
    poly_coeff_num = 2 * mode;

    for (int dim = 0; dim < 3; dim++) 
    {
        Eigen::VectorXd coeff;
        coeff.resize(poly_coeff_num);
        
        coeff = (polyCoeff.row(k)).segment(dim * poly_coeff_num, poly_coeff_num);
        Eigen::VectorXd times = Eigen::VectorXd::Zero(poly_coeff_num);
        
        for (int j = 0; j < poly_coeff_num; j++)
            if (j == 0)
                times(j) = 1.0;
            else
                times(j) = pow(t, j);
        
        double temp_pose = 0.0;
        for (int i = 0; i < times.rows(); ++i) 
        {
            temp_pose = temp_pose + coeff(i) * times(times.rows() - i - 1);
        }
        pt(dim) = temp_pose;
    }

    return pt;
}

/*
* 生成轨迹路径
*/
Path planner::trajectory_path(void)  //参数文件读取时间分配 (手动)
{
    Path trajectory;
    Eigen::Vector3d pos;

    for (int i = 0; i < time.size(); i++) 
    {
        for (double t = 0.0; t < time(i); t += 0.01) 
        {
            pos = getPosPoly(poly_coeff, i, t);
            trajectory.push_back(Pose(pos(0), pos(1), pos(2)));
        }
    }
    
    return trajectory;
}
Path planner::trajectory_path(Eigen::MatrixXd route_,Eigen::VectorXd time_)  //输入时间分配,1s一个点
{
    Path trajectory;
    Eigen::Vector3d pos;
    poly_coeff = getcoeff(route_,time_);
    for (int i = 0; i < time_.size(); i++) 
    {
        for (double t = 0.0; t < time_(i); t += 1) 
        {
            pos = getPosPoly(poly_coeff, i, t);
            trajectory.push_back(Pose(pos(0), pos(1), pos(2)));
        }
    }
    // std::cout << "trajectory point number:" << trajectory.size() << std::endl;
    return trajectory;
}
Path planner::trajectory_path(Eigen::MatrixXd route_, Eigen::VectorXd time_, double dis)  //输入时间分配,dis距离取一个点
{
    Path trajectory;
    Eigen::Vector3d pos, prev_pos;
    bool first_point = true;
    poly_coeff = getcoeff(route_, time_);
    
    for (int i = 0; i < time_.size(); i++) 
    {
        double accumulated_distance = 0.0;
        prev_pos = getPosPoly(poly_coeff, i, 0.0); // 获取段起始点
        
        for (double t = 0.0; t <= time_(i); t += 0.1) // 使用较小时间步长进行精细采样
        {
            pos = getPosPoly(poly_coeff, i, t);
            
            if (first_point) {
                trajectory.push_back(Pose(pos(0), pos(1), pos(2)));
                prev_pos = pos;
                first_point = false;
                continue;
            }
            
            // 计算当前点与前一个点的距离
            double segment_distance = (pos - prev_pos).norm();
            accumulated_distance += segment_distance;
            
            // 如果累积距离达到或超过目标距离，添加路径点
            if (accumulated_distance >= dis) {
                trajectory.push_back(Pose(pos(0), pos(1), pos(2)));
                accumulated_distance = 0.0; // 重置累积距离
            }
            
            prev_pos = pos;
        }
        
        // 确保添加段的终点
        pos = getPosPoly(poly_coeff, i, time_(i));
        trajectory.push_back(Pose(pos(0), pos(1), pos(2)));
    }
    
    // std::cout << "trajectory point number:" << trajectory.size() << std::endl;
    return trajectory;
}

