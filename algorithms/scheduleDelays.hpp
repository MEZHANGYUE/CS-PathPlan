#pragma once

#include <vector>
#include <opencv2/core.hpp>
#include <cmath>
#include "../pathPlanningCommon.h"

inline std::vector<double> scheduleDelays(
    const std::vector<std::vector<cv::Point2f>>& ships,
    double dt,             // 时间步长
    double safetyR,        // 安全半径（单位：像素）
    double tauSafe = 0.0   // 最小时间间隔（单位：秒）
) {
    int N = ships.size();
    std::vector<double> d(N, 0.0);  // 输出每艘船的延时（单位：秒）

    for (int i = 0; i < N; ++i) {
        double delay = 0.0;
        bool ok = false;

        while (!ok) {
            ok = true;

            for (int j = 0; j < i && ok; ++j) {
                const auto& path_i = ships[i];
                const auto& path_j = ships[j];

                for (size_t k = 0; k < path_i.size() && ok; ++k) {
                    double t_i = k * dt + delay;

                    for (size_t l = 0; l < path_j.size(); ++l) {
                        double t_j = l * dt + d[j];
                        double tij = std::abs(t_i - t_j);

                        if (tij < tauSafe &&
                            cv::norm(path_i[k] - path_j[l]) < safetyR) {
                            ok = false;
                            break;
                        }
                    }
                }
            }

            if (!ok) delay += dt;  // 若冲突则延迟一格
        }

        d[i] = delay;
    }

    return d;
}


// 按固定时间间隔重新采样路径
inline std::vector<cv::Point2f> resamplePathByTime(
    const std::vector<Point>& patrol_path,
    const std::vector<double>& speeds,   // 每个点的速度，单位 m/s
    double resolution,                   // 像素分辨率，单位 m/pixel
    double dt)                           // 采样时间步长，单位 s
{
    std::vector<cv::Point2f> float_path;
    int N = patrol_path.size();
    for (const auto& pt : patrol_path)
        float_path.emplace_back(static_cast<float>(pt.x), static_cast<float>(pt.y));

    std::vector<double> segment_times;  // 每段耗时
    std::vector<double> cum_times = {0.0};

    for (int i = 1; i < N; ++i) {
        double dx = float_path[i].x - float_path[i - 1].x;
        double dy = float_path[i].y - float_path[i - 1].y;
        double dist_m = std::sqrt(dx * dx + dy * dy) * resolution;
        double avg_speed = 0.5 * (speeds[i - 1] + speeds[i]);
        double time = dist_m / avg_speed;
        segment_times.push_back(time);
        cum_times.push_back(cum_times.back() + time);
    }

    double total_time = cum_times.back();
    std::vector<cv::Point2f> result;

    for (double t = 0.0; t <= total_time; t += dt) {
        int i = 0;
        while (i + 1 < cum_times.size() && cum_times[i + 1] < t)
            ++i;
        if (i + 1 >= cum_times.size()) break;

        double t0 = cum_times[i];
        double t1 = cum_times[i + 1];
        double alpha = (t - t0) / (t1 - t0);

        const auto& p0 = float_path[i];
        const auto& p1 = float_path[i + 1];
        cv::Point2f interp = p0 + static_cast<float>(alpha) * (p1 - p0);
        result.push_back(interp);
    }

    return result;
}