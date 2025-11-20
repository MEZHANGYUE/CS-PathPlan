#ifndef PLANNER_H
#define PLANNER_H

#include <vector>
#include <Eigen/Dense>
#include"trajectory_generator.h"
// 前向声明
// 纯C++的Pose和Path结构体定义
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

class planner 
{
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