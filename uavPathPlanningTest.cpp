#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <algorithm>
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
    if (!planner.getPlan(obj_total, result_json))
    {
        std::cerr << "Failed to plan!" << std::endl;
    }

    // if (!planner.runAltitudeOptimization("../data/neimeng.tif.ovr"))
    // {
    //     std::cerr << "Failed to Altitude Optimization!" << std::endl;
    // }
    // else
    // {
    //     std::cout << "output:" << std::endl;
    //     std::string out_string = result_json.dump();
    //     std::cout << out_string << std::endl;
    // }
    // 保存到自动推断的输出路径
    planner.saveJsonToFile(result_json, output_path);
    
}
