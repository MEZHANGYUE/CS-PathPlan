#pragma once

#include <opencv2/opencv.hpp>
#include <random>
#include <iostream>
#include "clipper.hpp"

using namespace ClipperLib;
using namespace std;

namespace yingji {
	// 缩放因子（Clipper 使用整数）
	const double SCALE = 1000.0;

    // 转换函数
    Path toClipperPath(const std::vector<cv::Point2f>& pts);

	Path toClipperPath(const std::vector<cv::Point>& pts);

	std::vector<cv::Point> fromClipperPath(const Path& path);

	// 几何方法计算多边形质心（适用于凸/凹多边形）
    // 原理：基于多边形有向面积公式（Shoelace Formula）
    cv::Point2f computePolygonCentroid(const std::vector<cv::Point>& poly);

	// 初始化种子点
	std::vector<cv::Point2f> initSeedPoints(const std::vector<cv::Point>& polygon, int n);

	// 一轮 Lloyd 迭代
	std::vector<cv::Point2f> lloydIteration(
    	const std::vector<cv::Point2f>& seeds,
    	const std::vector<cv::Point>& boundary_polygon,
    	std::vector<std::vector<cv::Point>>& out_regions
	);

	std::vector<cv::Point> shrinkPolygon(const std::vector<cv::Point>& polygon, 
	double resolution, double shrink_meters);

}