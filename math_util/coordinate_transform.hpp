// #pragma once
#ifndef P_PLANNING_H
#define P_PLANNING_H
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES // 用于启用M_PI常量
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
namespace math_util
{

struct WGS84Coord {
    double lon; // 经度（度）
    double lat; // 纬度（度）
    double alt; // 高度（米）

    // 构造函数
    WGS84Coord(double lon, double lat, double alt) : lon(lon), lat(lat), alt(alt) {}

    // 重载小于运算符，用于排序
    bool operator < (const WGS84Coord& p) const {
        return lon < p.lon || (lon == p.lon && lat < p.lat);
    }
    
    // 重载等于运算符
    bool operator == (const WGS84Coord& p) const {
        return lon == p.lon && lat == p.lat;
    }

};

struct ECEFCoord {
    double x;
    double y;
    double z;
    ECEFCoord(double x, double y, double z) : x(x), y(y), z(z) {}
};

struct ENUCoord {
    double e;
    double n;
    double u;
    ENUCoord(double e, double n, double u) : e(e), n(n), u(u) {}
};

// WGS84坐标系转ECEF（地心地固坐标系）
inline ECEFCoord wgs84_to_ecef(const WGS84Coord &coord) {
    const double a = 6378137.0;          // WGS84椭球长半轴
    const double e_sq = 0.00669437999013; // 第一偏心率的平方

    double lon_rad = coord.lon * M_PI / 180.0;
    double lat_rad = coord.lat * M_PI / 180.0;
    double sin_lat = sin(lat_rad);
    double cos_lat = cos(lat_rad);
    double sin_lon = sin(lon_rad);
    double cos_lon = cos(lon_rad);

    double N = a / sqrt(1 - e_sq * sin_lat * sin_lat);

    double x = (N + coord.alt) * cos_lat * cos_lon;
    double y = (N + coord.alt) * cos_lat * sin_lon;
    double z = (N * (1 - e_sq) + coord.alt) * sin_lat;

    return {x, y, z};
}

// ECEF坐标系转WGS84
inline WGS84Coord ecef_to_wgs84(const ECEFCoord &ecef) {
    const double a = 6378137.0;
    const double e_sq = 0.00669437999013;
    const double epsilon = 1e-12; // 收敛阈值
    const int max_iter = 100;

    double p = sqrt(ecef.x * ecef.x + ecef.y * ecef.y);
    double lon_rad = atan2(ecef.y, ecef.x);

    double lat_rad = atan2(ecef.z, p * (1 - e_sq));
    double h = 0.0;
    double diff = 1.0;
    int iter = 0;

    while (diff > epsilon && iter < max_iter) {
        double sin_lat = sin(lat_rad);
        double N = a / sqrt(1 - e_sq * sin_lat * sin_lat);
        h = p / cos(lat_rad) - N;

        double lat_new = atan2(ecef.z, p * (1 - e_sq * N / (N + h)));
        diff = fabs(lat_new - lat_rad);
        lat_rad = lat_new;
        iter++;
    }

    return {lon_rad * 180.0 / M_PI, lat_rad * 180.0 / M_PI, h};
}

// ECEF坐标系转ENU（东北天坐标系）
inline ENUCoord ecef_to_enu(const ECEFCoord &target, 
                    const ECEFCoord &ref_ecef,
                    const WGS84Coord &ref_wgs84) {
    double dx = target.x - ref_ecef.x;
    double dy = target.y - ref_ecef.y;
    double dz = target.z - ref_ecef.z;

    double lon_rad = ref_wgs84.lon * M_PI / 180.0;
    double lat_rad = ref_wgs84.lat * M_PI / 180.0;

    double sin_lon = sin(lon_rad);
    double cos_lon = cos(lon_rad);
    double sin_lat = sin(lat_rad);
    double cos_lat = cos(lat_rad);

    // 计算ENU坐标
    double e = -sin_lon * dx + cos_lon * dy;
    double n = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz;
    double u = cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz;

    return {e, n, u};
}

// ENU坐标系转ECEF
inline ECEFCoord enu_to_ecef(const ENUCoord &enu,
                     const ECEFCoord &ref_ecef,
                     const WGS84Coord &ref_wgs84) {
    double lon_rad = ref_wgs84.lon * M_PI / 180.0;
    double lat_rad = ref_wgs84.lat * M_PI / 180.0;

    double sin_lon = sin(lon_rad);
    double cos_lon = cos(lon_rad);
    double sin_lat = sin(lat_rad);
    double cos_lat = cos(lat_rad);

    // 计算ECEF坐标增量
    double dx = -sin_lon * enu.e - sin_lat * cos_lon * enu.n + cos_lat * cos_lon * enu.u;
    double dy =  cos_lon * enu.e - sin_lat * sin_lon * enu.n + cos_lat * sin_lon * enu.u;
    double dz =  cos_lat * enu.n + sin_lat * enu.u;

    return {ref_ecef.x + dx, ref_ecef.y + dy, ref_ecef.z + dz};
}

}

#endif