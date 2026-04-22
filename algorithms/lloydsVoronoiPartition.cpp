#include "lloydsVoronoiPartition.h"

ClipperLib::Path yingji::toClipperPath(const std::vector<cv::Point2f>& pts) {
    Path path;
    for (const auto& pt : pts)
        path.emplace_back(IntPoint((cInt)(pt.x * SCALE), (cInt)(pt.y * SCALE)));
    return path;
}

ClipperLib::Path yingji::toClipperPath(const std::vector<cv::Point>& pts) {
    ClipperLib::Path path;
    for (const auto& pt : pts)
        path.emplace_back(IntPoint((cInt)(pt.x * SCALE), (cInt)(pt.y * SCALE)));
    return path;
}

std::vector<cv::Point> yingji::fromClipperPath(const ClipperLib::Path& path) {
    std::vector<cv::Point> result;
    for (const auto& pt : path)
        result.emplace_back((int)(pt.X / SCALE), (int)(pt.Y / SCALE));
    return result;
}

// 几何方法计算多边形质心（适用于凸/凹多边形）
// 原理：基于多边形有向面积公式（Shoelace Formula）
cv::Point2f yingji::computePolygonCentroid(const std::vector<cv::Point>& poly) {
    double cx = 0, cy = 0;
    double A = 0.0;
    int n = poly.size();
    for (int i = 0; i < n; ++i) {
        const cv::Point& p0 = poly[i];
        const cv::Point& p1 = poly[(i + 1) % n];
        double cross = static_cast<double>(p0.x) * p1.y - static_cast<double>(p1.x) * p0.y;
        A += cross;
        cx += (p0.x + p1.x) * cross;
        cy += (p0.y + p1.y) * cross;
    }
    A *= 0.5;
    if (std::abs(A) < 1e-5) return poly[0]; // 面积太小，退化为第一个点
    cx /= (6.0 * A);
    cy /= (6.0 * A);
    return cv::Point2f(static_cast<float>(cx), static_cast<float>(cy));
}

// 初始化种子点
std::vector<cv::Point2f> yingji::initSeedPoints(const std::vector<cv::Point>& polygon, int n) {
    std::vector<cv::Point2f> seeds;
    cv::Rect bbox = cv::boundingRect(polygon);
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<> dist_x(bbox.x, bbox.x + bbox.width);
    std::uniform_real_distribution<> dist_y(bbox.y, bbox.y + bbox.height);

    while ((int)seeds.size() < n) {
        cv::Point2f pt(dist_x(rng), dist_y(rng));
        if (cv::pointPolygonTest(polygon, pt, false) >= 0) {
            seeds.push_back(pt);
        }
    }
    return seeds;
}

// 一轮 Lloyd 迭代
std::vector<cv::Point2f> yingji::lloydIteration(
    const std::vector<cv::Point2f>& seeds,
    const std::vector<cv::Point>& boundary_polygon,
    std::vector<std::vector<cv::Point>>& out_regions
) {
    cv::Rect rect = cv::boundingRect(boundary_polygon);
    rect.x -= 50; rect.y -= 50; rect.width += 100; rect.height += 100;

    cv::Subdiv2D subdiv(rect);
    for (const auto& pt : seeds) subdiv.insert(pt);

    std::vector<std::vector<cv::Point2f>> facets;
    std::vector<cv::Point2f> centers;
    subdiv.getVoronoiFacetList({}, facets, centers);

    Path boundary = toClipperPath(boundary_polygon);
    std::vector<cv::Point2f> new_seeds;
    out_regions.clear();

    for (const auto& facet : facets) {
        Clipper clipper;
        Path subj = toClipperPath(facet);
        clipper.AddPath(subj, ptSubject, true);
        clipper.AddPath(boundary, ptClip, true);
        Paths solution;
        clipper.Execute(ctIntersection, solution, pftNonZero, pftNonZero);
        if (!solution.empty()) {
            std::vector<cv::Point> region = fromClipperPath(solution[0]);
            out_regions.push_back(region);
            new_seeds.push_back(computePolygonCentroid(region));
        }
    }
    return new_seeds;
}

std::vector<cv::Point> yingji::shrinkPolygon(const std::vector<cv::Point>& polygon, 
	double resolution, double shrink_meters) {
    // 像素收缩距离 → 缩放后的整数距离
    double shrink_pixels = shrink_meters / resolution;
    double offset = -shrink_pixels * SCALE;  // 向内收缩

    // 转换为 Clipper 路径（整数坐标）
    Path subj;
    for (const auto& pt : polygon) {
        subj.push_back(IntPoint(static_cast<cInt>(pt.x * SCALE), static_cast<cInt>(pt.y * SCALE)));
    }

    ClipperOffset co;
    co.AddPath(subj, jtMiter, etClosedPolygon);

    Paths solution;
    co.Execute(solution, offset);

    // 返回收缩后的多边形（只取第一个，如果有多个可按需合并）
    std::vector<cv::Point> result;
    if (!solution.empty()) {
        for (const auto& pt : solution[0]) {
            result.emplace_back(static_cast<int>(pt.X / SCALE), static_cast<int>(pt.Y / SCALE));
        }
    }

    return result;
}