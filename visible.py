import json
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import glob
from matplotlib import font_manager

def setup_chinese_font():
    """设置中文字体支持"""
    try:
        # 尝试使用系统中已有的中文字体
        font_dirs = ['/usr/share/fonts/', '/usr/local/share/fonts/', os.path.expanduser('~/.fonts/')]
        font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
        
        chinese_fonts = []
        for font_file in font_files:
            try:
                prop = font_manager.FontProperties(fname=font_file)
                font_name = prop.get_name()
                # 检查是否是中文字体
                if any(char in font_name for char in ['宋体', '黑体', '微软', 'Microsoft', 'Sim', 'Kai', 'Hei', 'Song']):
                    chinese_fonts.append(font_file)
            except:
                continue
        
        if chinese_fonts:
            # 使用找到的第一个中文字体
            font_prop = font_manager.FontProperties(fname=chinese_fonts[0])
            plt.rcParams['font.family'] = font_prop.get_name()
        else:
            # 如果没有找到中文字体，使用默认字体并禁用中文显示
            print("警告：未找到中文字体，将使用英文显示")
            plt.rcParams['font.family'] = ['DejaVu Sans', 'Arial', 'sans-serif']
            
    except Exception as e:
        print(f"字体设置警告: {e}")
        plt.rcParams['font.family'] = ['DejaVu Sans', 'Arial', 'sans-serif']

def read_json_file(filename):
    """读取JSON文件"""
    try:
        with open(filename, 'r', encoding='utf-8') as file:
            data = json.load(file)
        return data
    except FileNotFoundError:
        print(f"错误：文件 {filename} 未找到")
        return None
    except json.JSONDecodeError:
        print(f"错误：文件 {filename} 不是有效的JSON格式")
        return None

def extract_coordinates(data, path_key):
    """从数据中提取坐标点"""
    coordinates = []
    
    if isinstance(data, dict) and path_key in data:
        points = data[path_key]
        if isinstance(points, list):
            for point in points:
                if isinstance(point, (list, tuple)) and len(point) >= 2:
                    # 假设坐标格式为 [经度, 纬度] 或 [x, y]
                    coordinates.append((point[0], point[1]))
                elif isinstance(point, dict):
                    # 如果坐标是对象形式，如 {"lat": xx, "lng": xx}
                    if 'lat' in point and 'lng' in point:
                        coordinates.append((point['lng'], point['lat']))
                    elif 'latitude' in point and 'longitude' in point:
                        coordinates.append((point['longitude'], point['latitude']))
                    elif 'x' in point and 'y' in point:
                        coordinates.append((point['x'], point['y']))
    
    return coordinates

def plot_path_and_trajectory(waypoints, trajectory_points, title="Path and Trajectory Visualization"):
    """绘制路径点和轨迹点"""
    plt.figure(figsize=(12, 8))
    
    # 提取坐标
    if waypoints:
        wp_x, wp_y = zip(*waypoints)
        plt.scatter(wp_x, wp_y, c='red', s=100, marker='o', label='Waypoints', zorder=5)
        plt.plot(wp_x, wp_y, 'r--', alpha=0.7, linewidth=2, label='Planned Path')
        
        # 在路径点上添加编号
        for i, (x, y) in enumerate(waypoints):
            plt.annotate(f'{i}', (x, y), xytext=(5, 5), textcoords='offset points', 
                        fontsize=8, color='red')
    
    if trajectory_points:
        traj_x, traj_y = zip(*trajectory_points)
        plt.scatter(traj_x, traj_y, c='blue', s=30, marker='.', alpha=0.6, label='Trajectory Points')
        plt.plot(traj_x, traj_y, 'b-', alpha=0.8, linewidth=1.5, label='Actual Trajectory')
    
    plt.xlabel('Longitude / X Coordinate')
    plt.ylabel('Latitude / Y Coordinate')
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 设置等比例显示
    plt.axis('equal')
    plt.tight_layout()
    plt.show()

def find_matching_files(file_pattern):
    """智能查找匹配的文件"""
    # 可能的文件后缀模式
    patterns = [
        f"{file_pattern}*input*.json",
        f"{file_pattern}*output*.json",
        f"{file_pattern}*planning*input*.json",
        f"{file_pattern}*planning*output*.json",
        f"{file_pattern}_*input*.json",
        f"{file_pattern}_*output*.json"
    ]
    
    input_files = []
    output_files = []
    
    for pattern in patterns:
        matches = glob.glob(pattern)
        for match in matches:
            filename_lower = match.lower()
            if 'input' in filename_lower and 'output' not in filename_lower:
                input_files.append(match)
            elif 'output' in filename_lower and 'input' not in filename_lower:
                output_files.append(match)
    
    return input_files, output_files

def auto_detect_files(base_path):
    """自动检测输入输出文件"""
    # 如果直接提供了完整的文件路径
    if base_path.endswith('.json'):
        if 'input' in base_path.lower():
            input_file = base_path
            # 尝试找到对应的输出文件
            output_file = base_path.replace('input', 'output')
            if not os.path.exists(output_file):
                output_file = base_path.replace('Input', 'Output')
        elif 'output' in base_path.lower():
            output_file = base_path
            # 尝试找到对应的输入文件
            input_file = base_path.replace('output', 'input')
            if not os.path.exists(input_file):
                input_file = base_path.replace('Output', 'Input')
        else:
            print("Error: Cannot identify file type, please ensure filename contains 'input' or 'output'")
            return None, None
    else:
        # 使用模式匹配查找文件
        input_files, output_files = find_matching_files(base_path + '*')
        
        if not input_files and not output_files:
            # 如果没找到，尝试在data目录下查找
            input_files, output_files = find_matching_files(f"./data/{os.path.basename(base_path)}*")
        
        # 选择最匹配的文件
        input_file = input_files[0] if input_files else None
        output_file = output_files[0] if output_files else None
    
    return input_file, output_file

def main(file_path):
    """主函数，智能匹配输入输出文件"""
    
    # 设置中文字体（在绘图前调用）
    setup_chinese_font()
    
    # 自动检测文件
    input_file, output_file = auto_detect_files(file_path)
    
    if not input_file:
        print("Error: Input file not found")
        print("Please check file path or use one of the following formats:")
        print("  python3 visible.py ./data/uav34_0_planning")
        print("  python3 visible.py ./data/uav34_0")
        print("  python3 visible.py ./data/uav34")
        print("  python3 visible.py ./data/uav34_0_planning_input.json")
        return
    
    if not output_file:
        print("Error: Output file not found")
        print(f"Input file: {input_file}")
        return
    
    print(f"Found input file: {input_file}")
    print(f"Found output file: {output_file}")
    
    # 读取输入文件（路径点）
    print("Reading input file...")
    input_data = read_json_file(input_file)
    if input_data is None:
        return
    
    # 读取输出文件（轨迹点）
    print("Reading output file...")
    output_data = read_json_file(output_file)
    if output_data is None:
        return
    
    # 提取路径点
    waypoints = extract_coordinates(input_data, "leader_midway_point_wgs84")
    print(f"Found {len(waypoints)} waypoints")
    
    # 提取轨迹点
    trajectory_points = extract_coordinates(output_data, "uav_leader_plane1")
    print(f"Found {len(trajectory_points)} trajectory points")
    
    if not waypoints and not trajectory_points:
        print("No valid waypoint or trajectory data found")
        return
    
    # 从文件名中提取 UAV 编号用于标题
    uav_id = os.path.basename(file_path).split('_')[0]
    
    # 可视化
    print("Generating visualization...")
    plot_path_and_trajectory(waypoints, trajectory_points, 
                           f"{uav_id} Path Planning and Execution Trajectory")

def analyze_data_structure(filename):
    """分析JSON文件数据结构（用于调试）"""
    data = read_json_file(filename)
    if data:
        print(f"\n{filename} data structure analysis:")
        print_json_structure(data)

def print_json_structure(obj, indent=0):
    """递归打印JSON结构"""
    prefix = "  " * indent
    
    if isinstance(obj, dict):
        for key, value in obj.items():
            if isinstance(value, (dict, list)):
                print(f"{prefix}{key}: {type(value).__name__}")
                print_json_structure(value, indent + 1)
            else:
                print(f"{prefix}{key}: {type(value).__name__}")
    elif isinstance(obj, list):
        if obj:
            print(f"{prefix}[0]: {type(obj[0]).__name__}")
            if len(obj) > 1:
                print(f"{prefix}... ({len(obj)} elements)")
        else:
            print(f"{prefix}[] (empty list)")

def print_usage():
    """打印使用说明"""
    print("Usage: python3 visible.py <file path or prefix>")
    print("Examples:")
    print("  python3 visible.py ./data/uav34_0_planning")
    print("  python3 visible.py ./data/uav34_0")
    print("  python3 visible.py ./data/uav34")
    print("  python3 visible.py ./data/uav34_0_planning_input.json")
    print("Auto-matching rules:")
    print("  - Find corresponding files containing 'input' and 'output'")
    print("  - Case-insensitive matching supported")

if __name__ == "__main__":
    # 检查命令行参数
    if len(sys.argv) != 2:
        print("Error: Please provide file path as argument")
        print_usage()
        sys.exit(1)
    
    file_path = sys.argv[1]
    
    # 如果需要分析数据结构，取消下面的注释
    # analyze_data_structure(f"{file_path}_input.json")
    # analyze_data_structure(f"{file_path}_output.json")
    
    main(file_path)