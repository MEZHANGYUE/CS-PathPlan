#include "hungarianAlgorithm.h"
#include "csLog/csLog.h"

// 使用 OpenCV 图像矩计算多边形重心（轮廓闭合图形适用）
cv::Point2f yingji::computeMomentCentroid(const std::vector<cv::Point>& polygon) {
    cv::Moments m = cv::moments(polygon);
    if (m.m00 == 0) return polygon[0]; // 或者 {0,0}，根据需求调整
    return cv::Point2f(
        static_cast<float>(m.m10 / m.m00),
        static_cast<float>(m.m01 / m.m00)
    );
}

// 匈牙利算法（最小化匹配成本）
// 适用于船数 >= 区域数（将补齐方阵）
std::vector<int> yingji::hungarianAlgorithm(const std::vector<std::vector<double>>& cost) {
    int n = cost.size();       // 行数（船数）
    int m = cost[0].size();    // 列数（区域数）
    int size = std::max(n, m); // 转为方阵

    // 构造方阵，初始化为 INF
    std::vector<std::vector<double>> cost_matrix(size, std::vector<double>(size, 1e9));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            cost_matrix[i][j] = cost[i][j];

    std::vector<double> u(size), v(size);
    std::vector<int> p(size + 1), way(size + 1);

    for (int i = 0; i < size; ++i) {
        p[0] = i;
        int j0 = 0;
        std::vector<double> minv(size, 1e9);
        std::vector<bool> used(size, false);
        do {
            used[j0] = true;
            int i0 = p[j0];
            int j1 = -1;
            double delta = 1e9;

            for (int j = 0; j < size; ++j) {
                if (!used[j]) {
                    double cur = cost_matrix[i0][j] - u[i0] - v[j];
                    if (cur < minv[j]) {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if (minv[j] < delta) {
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }

            if (j1 == -1) {
                VRLOG_DEBUG_F << "Hungarian warning: no candidate found during this iteration (i=" << i << "). Possibly due to degenerate 1x1 case."
                << std::endl;
                break; 
            }

            for (int j = 0; j < size; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
        } while (p[j0] != 0);


        do {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }

    // p[j] -> i，表示第 i 行分配给第 j 列
    std::vector<int> result(m, -1); // 区域 -> 船
    for (int j = 1; j <= size; ++j) {
        if (p[j] < n) {
            result[j - 1] = p[j]; // 第 j-1 个区域分配给第 p[j] 船
        }
    }
    return result;
}

// 主函数：计算分配关系
std::vector<int> yingji::assignShipsToRegions(const std::vector<cv::Point>& start_pixels,
                                      const std::vector<std::vector<cv::Point>>& patrol_regions) {
    int num_ships = start_pixels.size();
    int num_regions = patrol_regions.size();

    std::vector<std::vector<double>> cost(num_ships, std::vector<double>(num_regions));

    for (int i = 0; i < num_ships; ++i) {
        for (int j = 0; j < num_regions; ++j) {
            cv::Point2f centroid = computeMomentCentroid(patrol_regions[j]);
            cv::Point2f start_pixel = 
            {static_cast<float>(start_pixels[i].x), static_cast<float>(start_pixels[i].y)};
            cost[i][j] = cv::norm(start_pixel - centroid);
        }
    }

    return hungarianAlgorithm(cost); // 返回区域 -> 船编号
}

std::vector<int> yingji::assignShipsToRegions_SortedDirectly(
    const std::vector<cv::Point>& start_pixels,
    const std::vector<std::vector<cv::Point>>& patrol_regions) {

    int num_ships = start_pixels.size();
    int num_regions = patrol_regions.size();

    // 1. 计算区域重心
    std::vector<cv::Point2f> centroids(num_regions);
    for (int j = 0; j < num_regions; ++j)
        centroids[j] = computeMomentCentroid(patrol_regions[j]);

    // 2. 判断使用横向（x）还是纵向（y）排序
    cv::Rect ship_bbox = cv::boundingRect(start_pixels);
    bool horizontal_split = (ship_bbox.width > ship_bbox.height);

    // 3. 排序索引
    std::vector<int> ship_indices(num_ships);
    std::vector<int> region_indices(num_regions);
    std::iota(ship_indices.begin(), ship_indices.end(), 0);
    std::iota(region_indices.begin(), region_indices.end(), 0);

    if (horizontal_split) {
        std::sort(ship_indices.begin(), ship_indices.end(),
                  [&](int a, int b) { return start_pixels[a].x < start_pixels[b].x; });
        std::sort(region_indices.begin(), region_indices.end(),
                  [&](int a, int b) { return centroids[a].x < centroids[b].x; });
    } else {
        std::sort(ship_indices.begin(), ship_indices.end(),
                  [&](int a, int b) { return start_pixels[a].y < start_pixels[b].y; });
        std::sort(region_indices.begin(), region_indices.end(),
                  [&](int a, int b) { return centroids[a].y < centroids[b].y; });
    }

    // 4. 一一匹配：排序第 i 船 → 排序第 i 区
    std::vector<int> final_result(num_regions, -1);
    int n = std::min(num_ships, num_regions);
    for (int i = 0; i < n; ++i) {
        final_result[region_indices[i]] = ship_indices[i];
    }

    return final_result;  // 区域索引 -> 对应船编号
}

