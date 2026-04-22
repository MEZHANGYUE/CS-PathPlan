#pragma once

#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <limits>
#include <cmath>
#include <numeric> 

using namespace std;

namespace yingji {
    // 使用 OpenCV 图像矩计算多边形重心（轮廓闭合图形适用）
    cv::Point2f computeMomentCentroid(const std::vector<cv::Point>& polygon);

    // 匈牙利算法（最小化匹配成本）
    // 适用于船数 >= 区域数（将补齐方阵）
    std::vector<int> hungarianAlgorithm(const std::vector<std::vector<double>>& cost);

    // 主函数：计算分配关系
    std::vector<int> assignShipsToRegions(const std::vector<cv::Point>& start_pixels,
                                        const std::vector<std::vector<cv::Point>>& patrol_regions);

    std::vector<int> assignShipsToRegions_SortedDirectly(
    const std::vector<cv::Point>& start_pixels,
    const std::vector<std::vector<cv::Point>>& patrol_regions);

}
