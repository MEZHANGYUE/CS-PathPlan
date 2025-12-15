#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <algorithm>
#include <chrono>
#include <yaml-cpp/yaml.h>
#include "uavPathPlanning.hpp"
#include "math_util/minimum_snap.hpp"
int main(int argc, char *argv[])
{
    /****** test code ******/
    // 支持通过命令行参数指定输入标识符，例如: ./uavPathPlanningTest 32_0
    std::string token;
    if (argc > 1) token = argv[1];
    // 搜索 data 目录下匹配的 input 文件（文件名包含 token 并包含 "input"）
    std::filesystem::path data_dir = std::filesystem::path("../data");
    std::string input_path;
    if (token.empty()) {
        // 默认文件（向后兼容）
        input_path = (data_dir / "uav32_0_planning_input.json").string();
    } else {
        for (auto &entry : std::filesystem::directory_iterator(data_dir)) {
            if (!entry.is_regular_file()) continue;
            auto p = entry.path();
            if (p.extension() == ".json") {
                std::string fname = p.filename().string();
                // 区分大小写查找 token 与 "input"
                if (fname.find(token) != std::string::npos && fname.find("input") != std::string::npos) {
                    input_path = p.string();
                    break;
                }
            }
        }
        if (input_path.empty()) {
            std::cerr << "未在 ../data 中找到匹配 '" << token << "' 且包含 'input' 的 json 文件，使用默认文件。" << std::endl;
            input_path = (data_dir / "uav32_0_planning_input.json").string();
        }
    }

    std::cout << "Using input file: " << input_path << std::endl;

    // 输出文件名：将文件名中的第一个 "input" 替换为 "output"
    std::string output_path = input_path;
    auto pos = output_path.rfind("input");
    if (pos != std::string::npos) {
        output_path.replace(pos, 5, "output");
    } else {
        // fallback: append _output.json
        std::filesystem::path op(output_path);
        output_path = (op.parent_path() / (op.stem().string() + "_output" )).string() + op.extension().string();
    }

    std::cout << "Will write output to: " << output_path << std::endl;

    // 读取并解析JSON文件内容到json对象
    json obj_total, result_json;
    try {
        std::ifstream file(input_path);
        if (!file.is_open()) {
            std::cerr << "无法打开文件: " << input_path << std::endl;
            return 1;
        }
        file >> obj_total;
        file.close();
    } catch (const std::exception &e) {
        std::cerr << "JSON 解析错误: " << e.what() << std::endl;
        return 1;
    }
    std::string json_string = obj_total.dump();
    std::cout << "input:" << std::endl;
    std::cout << json_string << std::endl;

    UavPathPlanner planner;
    if (!planner.getPlan(obj_total, result_json,true,"minimum_snap"))
    {
        std::cerr << "Failed to plan!" << std::endl;
    }

    // Load configuration
    std::string config_path = "../config.yaml";
    bool altitude_opt_enabled = false;
    std::string elevation_file;

    try {
        if (std::filesystem::exists(config_path)) {
            YAML::Node config = YAML::LoadFile(config_path);
            if (config["altitude_optimization"]) {
                if (config["altitude_optimization"]["enabled"]) {
                    altitude_opt_enabled = config["altitude_optimization"]["enabled"].as<bool>();
                    std::cout << "altitude_opt_enabled: " << altitude_opt_enabled << std::endl;
                }
                if (config["altitude_optimization"]["elevation_file"]) {
                    elevation_file = config["altitude_optimization"]["elevation_file"].as<std::string>();
                } else {
                    std::cerr << "Warning: 'elevation_file' parameter missing in " << config_path << std::endl;
                }
            }
            std::cout << "Loaded config from " << config_path << std::endl;
        } else {
             std::cout << "Config file " << config_path << " not found." << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Warning: Failed to load config.yaml (" << e.what() << ")" << std::endl;
    }
    //需要高度优化
    if (altitude_opt_enabled) {
        if (!elevation_file.empty()) {
            std::cout << "Running Altitude Optimization with file: " << elevation_file << std::endl;
            auto start_time = std::chrono::high_resolution_clock::now();
            if (!planner.runAltitudeOptimization(elevation_file, result_json, obj_total))
            {
                std::cerr << "Failed to Altitude Optimization!" << std::endl;
            }
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            std::cout << "Altitude optimization time: " << elapsed.count() << "s" << std::endl;
        } else {
             std::cerr << "Altitude optimization enabled but no elevation file specified." << std::endl;
        }
    } else {
        std::cout << "Altitude optimization skipped." << std::endl;
    }

    // 从输入 JSON 的 uav_leader_id 和 uavs_id 构建 using_uav_list（包含长机与僚机，去重）
    result_json["using_uav_list"] = json::array();
    if (obj_total.contains("uav_leader_id") && obj_total["uav_leader_id"].is_array()) {
        for (const auto &id : obj_total["uav_leader_id"]) {
            if (id.is_number()) result_json["using_uav_list"].push_back(id);
        }
    }
    if (obj_total.contains("uavs_id") && obj_total["uavs_id"].is_array()) {
        for (const auto &id : obj_total["uavs_id"]) {
            if (!id.is_number()) continue;
            bool found = false;
            for (const auto &existing : result_json["using_uav_list"]) {
                if (existing == id) { found = true; break; }
            }
            if (!found) result_json["using_uav_list"].push_back(id);
        }
    }
    result_json["ready_id"] = json::array();
    if (obj_total.contains("ready_id") && obj_total["ready_id"].is_array())
    result_json["ready_id"] = obj_total["ready_id"];

    result_json["leader_show_points"]= json::array();
    if (obj_total.contains("leader_midway_point_wgs84") && obj_total["leader_midway_point_wgs84"].is_array())
        result_json["leader_show_points"] = obj_total["leader_midway_point_wgs84"];
    // 如果输入中包含 high_zhandou_point_wgs84，则将这些点追加到 leader_show_points，
    // 并把高度设置为 leader_show_points 中最后一个点的高度（尝试常见字段名：z、altitude、height，或数组第3个元素）
    if (obj_total.contains("high_zhandou_point_wgs84") && obj_total["high_zhandou_point_wgs84"].is_array()) {
        double last_alt = 0.0;
        if (result_json["leader_show_points"].is_array() && !result_json["leader_show_points"].empty()) {
            auto last_pt = result_json["leader_show_points"].back();
            if (last_pt.is_object()) {
                if (last_pt.contains("z") && last_pt["z"].is_number()) last_alt = last_pt["z"].get<double>();
                else if (last_pt.contains("altitude") && last_pt["altitude"].is_number()) last_alt = last_pt["altitude"].get<double>();
                else if (last_pt.contains("height") && last_pt["height"].is_number()) last_alt = last_pt["height"].get<double>();
            } else if (last_pt.is_array() && last_pt.size() >= 3 && last_pt[2].is_number()) {
                last_alt = last_pt[2].get<double>();
            }
        }
        for (const auto &pt : obj_total["high_zhandou_point_wgs84"]) {
            if (pt.is_object()) {
                json newpt = pt;
                if (newpt.contains("z")) newpt["z"] = last_alt;
                else if (newpt.contains("altitude")) newpt["altitude"] = last_alt;
                else if (newpt.contains("height")) newpt["height"] = last_alt;
                else newpt["z"] = last_alt;
                result_json["leader_show_points"].push_back(newpt);
            } else if (pt.is_array() && pt.size() >= 2) {
                json newpt = pt;
                if (pt.size() >= 3) {
                    newpt[2] = last_alt;
                } else {
                    newpt.push_back(last_alt);
                }
                result_json["leader_show_points"].push_back(newpt);
            }
        }
    }
    
    std::cout << "Constructed using_uav_list from input: " << result_json["using_uav_list"].dump() << std::endl;

    // 保存到自动推断的输出路径
    planner.saveJsonToFile(result_json, output_path);
    std::cout << "Output JSON updated in memory." << std::endl;
    
}
