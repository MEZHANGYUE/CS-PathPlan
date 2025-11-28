#include <vector>
#include <iostream>
#include <fstream>
#include "uavPathPlanning.hpp"
#include "math_util/minimum_snap.hpp"
int main(int argc, char *argv[])
{
    /****** test code ******/
    // 打开JSON文件
    std::ifstream file("../data/uav34_0_planning_input.json");
    if (!file.is_open())
    {
        std::cerr << "无法打开文件" << std::endl;
        return 1;
    }
    // 读取并解析JSON文件内容到json对象
    json obj_total, result_json;
    file >> obj_total; // 或者使用 json::parse(file) 来代替file >> j;
    file.close();      // 关闭文件流，但不是必须的，因为file会在作用域结束时自动关闭。    
    std::string json_string = obj_total.dump();
    std::cout << "input:" << std::endl;
    std::cout << json_string << std::endl;

    if (!getPlan(obj_total, result_json))
    {
        std::cerr << "Failed to plan!" << std::endl;
    }
    // else
    // {
    //     std::cout << "output:" << std::endl;
    //     std::string out_string = result_json.dump();
    //     std::cout << out_string << std::endl;
    // }
    saveJsonToFile(result_json, "../data/uav34_0_planning_output.json");
    
}
