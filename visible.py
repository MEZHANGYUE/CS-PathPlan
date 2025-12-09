import json
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import glob
import re
from matplotlib import font_manager
import math
from PIL import Image
Image.MAX_IMAGE_PIXELS = None

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


def load_tiff_elevation(path, max_dim=2000):
    """Load TIFF elevation data, optionally downsample if too large, and extract extent."""
    try:
        img = Image.open(path)
        
        # Extract GeoTIFF tags
        tags = img.tag_v2
        extent = None
        if 33550 in tags and 33922 in tags:
            scale = tags[33550]
            tiepoint = tags[33922]
            
            origin_x = tiepoint[3]
            origin_y = tiepoint[4]
            pixel_w = scale[0]
            pixel_h = -scale[1] 
            
            w, h = img.size
            xmin = origin_x
            xmax = origin_x + pixel_w * w
            ymin = origin_y + pixel_h * h
            ymax = origin_y 
            
            if ymin > ymax:
                ymin, ymax = ymax, ymin
                
            extent = [xmin, xmax, ymin, ymax]
            print(f"Loaded TIFF elevation with extent: {extent}, original size: {w}x{h}")
            
        # Resize if too large to save memory/time
        w, h = img.size
        if w > max_dim or h > max_dim:
            ratio = min(max_dim / w, max_dim / h)
            new_w = int(w * ratio)
            new_h = int(h * ratio)
            img = img.resize((new_w, new_h), Image.BILINEAR)
            print(f"Resized elevation data to {new_w}x{new_h}")
            
        data = np.array(img)
        return data, extent
    except Exception as e:
        print(f"Error loading TIFF elevation: {e}")
        return None, None

def load_pgm_image(path):
    """简单的 PGM (P5/P2) 读取器，返回 numpy 数组（H x W）或 None。"""
    try:
        with open(path, 'rb') as f:
            # 读取魔数（如 P5 或 P2）并跳过注释
            magic = f.readline().strip()
            if not magic:
                return None
            magic = magic.decode('ascii') if isinstance(magic, bytes) else magic
            # 跳过注释和空白行，读取尺寸
            def read_token():
                tok = b''
                while True:
                    ch = f.read(1)
                    if not ch:
                        break
                    if ch.isspace():
                        if tok:
                            break
                        else:
                            continue
                    if ch == b'#':
                        # 跳过整行注释
                        f.readline()
                        continue
                    tok += ch
                return tok.decode('ascii') if tok else None

            w = int(read_token())
            h = int(read_token())
            maxv = int(read_token())
            if magic == 'P5':
                # 二进制数据
                # consume single whitespace/newline after maxv if any
                # 已由 read_token 停留在下一个字节
                data = f.read(w * h)
                arr = np.frombuffer(data, dtype=np.uint8)
                if arr.size != w * h:
                    # 有时像素为 2 字节（maxv > 255）
                    f.seek(0)
                    # 回退并使用更稳健的读取方式
                    f = open(path, 'r')
                    content = f.read()
                    # fallback: use matplotlib.imread
                    try:
                        img = plt.imread(path)
                        if img.ndim == 3:
                            img = img[:, :, 0]
                        return img.astype(np.float32)
                    except Exception:
                        return None
                arr = arr.reshape((h, w))
                return arr.astype(np.float32)
            elif magic == 'P2':
                # ASCII 格式
                toks = []
                while len(toks) < w * h:
                    t = read_token()
                    if t is None:
                        break
                    toks.append(int(t))
                arr = np.array(toks, dtype=np.float32).reshape((h, w))
                return arr
            else:
                return None
    except Exception:
        return None


def parse_aux_xml_for_geotransform(aux_path):
    """尝试从 .aux.xml 文件中解析 GeoTransform（返回 6 元素列表或 None）。"""
    if not os.path.exists(aux_path):
        return None
    try:
        text = open(aux_path, 'r', encoding='utf-8', errors='ignore').read()
        # 搜索 <GeoTransform> .. numbers .. </GeoTransform>
        m = re.search(r"<GeoTransform>([\s\S]*?)</GeoTransform>", text)
        if not m:
            return None
        nums = re.findall(r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?", m.group(1))
        if len(nums) >= 6:
            gt = [float(x) for x in nums[:6]]
            return gt
    except Exception:
        return None
    return None

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
                # 列表或元组形式: [lon, lat] 或 [lon, lat, alt]
                if isinstance(point, (list, tuple)) and len(point) >= 2:
                    if len(point) >= 3:
                        coordinates.append((point[0], point[1], point[2]))
                    else:
                        coordinates.append((point[0], point[1]))
                # 字典形式，支持多种键名
                elif isinstance(point, dict):
                    if 'lat' in point and 'lng' in point:
                        if 'alt' in point or 'altitude' in point:
                            coordinates.append((point.get('lng'), point.get('lat'), point.get('alt') or point.get('altitude')))
                        else:
                            coordinates.append((point.get('lng'), point.get('lat')))
                    elif 'latitude' in point and 'longitude' in point:
                        if 'altitude' in point:
                            coordinates.append((point.get('longitude'), point.get('latitude'), point.get('altitude')))
                        else:
                            coordinates.append((point.get('longitude'), point.get('latitude')))
                    elif 'x' in point and 'y' in point:
                        if 'z' in point:
                            coordinates.append((point.get('x'), point.get('y'), point.get('z')))
                        else:
                            coordinates.append((point.get('x'), point.get('y')))
    
    return coordinates


def extract_all_plane_trajectories(data, prefix_regex=r'^uav_plane'):
    """解析 JSON 中所有以 uav_plane 开头的键，返回列表 [(id, [(lon,lat), ...]), ...]

    处理的数据格式示例：
    "uav_plane1": [
        [33, [lon, lat, alt], [lon, lat, alt], ...],
        [34, [lon, lat, alt], ...]
    ]
    """
    plane_trajs = []
    if not isinstance(data, dict):
        return plane_trajs

    for key, val in data.items():
        if re.match(prefix_regex, key):
            if isinstance(val, list):
                for entry in val:
                    # entry 应为列表，第一项是 id，后面是点
                    if isinstance(entry, list) and len(entry) >= 2:
                        try:
                            pid = entry[0]
                            pts = []
                            for p in entry[1:]:
                                if isinstance(p, (list, tuple)) and len(p) >= 2:
                                    # 支持 [lon, lat] 或 [lon, lat, alt]
                                    if len(p) >= 3:
                                        pts.append((p[0], p[1], p[2]))
                                    else:
                                        pts.append((p[0], p[1]))
                                elif isinstance(p, dict):
                                    if 'lat' in p and 'lng' in p:
                                        if 'alt' in p or 'altitude' in p:
                                            pts.append((p.get('lng'), p.get('lat'), p.get('alt') or p.get('altitude')))
                                        else:
                                            pts.append((p['lng'], p['lat']))
                                    elif 'latitude' in p and 'longitude' in p:
                                        if 'altitude' in p:
                                            pts.append((p.get('longitude'), p.get('latitude'), p.get('altitude')))
                                        else:
                                            pts.append((p['longitude'], p['latitude']))
                            if pts:
                                plane_trajs.append((pid, pts))
                        except Exception:
                            # 忽略单条格式错误的 entry
                            continue

    return plane_trajs

def plot_path_and_trajectory(waypoints=None, leader_traj=None, plane_trajs=None, title="Path and Trajectory Visualization", save_path=None, show_plot=True, plot_3d=None, bg_img=None, bg_extent=None, bg_cmap='gray', elevation_data=None, elevation_extent=None):
    """绘制路径点、leader 轨迹（可选）以及多条 plane 轨迹（plane_trajs 为 [(id, [(lon,lat),...]), ...]）。

    - waypoints: 列表 [(lon,lat), ...]
    - leader_traj: 列表 [(lon,lat), ...]
    - plane_trajs: 列表 [(id, [(lon,lat), ...]), ...]
    """
    # 自动检测是否需要 3D 绘图：如果任何坐标包含第三维度，则切换为 3D
    need_3d = False
    # def has_z(points):
    #     if not points:
    #         return False
    #     for p in points:
    #         if isinstance(p, (list, tuple)) and len(p) >= 3:
    #             return True
    #     return False

    # plot_3d 参数优先：None = 自动检测；True/False = 强制
    if plot_3d is None:
        # if has_z(waypoints) or has_z(leader_traj) or (plane_trajs and any(has_z(pts) for (_id, pts) in plane_trajs)):
            need_3d = False
    else:
        need_3d = bool(plot_3d)

    fig = plt.figure(figsize=(12, 8))
    if need_3d:
        ax = fig.add_subplot(111, projection='3d')
        
        # 绘制3D高程地形
        if elevation_data is not None and elevation_extent is not None:
            try:
                xmin, xmax, ymin, ymax = elevation_extent
                h, w = elevation_data.shape
                
                # 计算轨迹包围盒
                all_coords = []
                if waypoints: all_coords.extend(waypoints)
                if leader_traj: all_coords.extend(leader_traj)
                if plane_trajs:
                    for _, pts in plane_trajs:
                        all_coords.extend(pts)
                
                crop_xmin, crop_xmax, crop_ymin, crop_ymax = xmin, xmax, ymin, ymax
                
                if all_coords:
                    xs = [p[0] for p in all_coords]
                    ys = [p[1] for p in all_coords]
                    if xs and ys:
                        t_xmin, t_xmax = min(xs), max(xs)
                        t_ymin, t_ymax = min(ys), max(ys)
                        # 增加 10% 边距
                        # pad_x = max((t_xmax - t_xmin) * 0.1, 0.002)
                        # pad_y = max((t_ymax - t_ymin) * 0.1, 0.002)
                        pad_x = (t_xmax - t_xmin)*0.8
                        pad_y = (t_ymax - t_ymin)*0.8
                        # 与地图范围求交集
                        crop_xmin = max(xmin, t_xmin - pad_x)
                        crop_xmax = min(xmax, t_xmax + pad_x)
                        crop_ymin = max(ymin, t_ymin - pad_y)
                        crop_ymax = min(ymax, t_ymax + pad_y)

                if crop_xmax > crop_xmin and crop_ymax > crop_ymin:
                    # 计算对应的索引范围
                    # X: xmin -> xmax
                    ix_start = int((crop_xmin - xmin) / (xmax - xmin) * w)
                    ix_end = int((crop_xmax - xmin) / (xmax - xmin) * w)
                    
                    # Y: ymax -> ymin (图像行索引 0 对应 ymax)
                    iy_start = int((ymax - crop_ymax) / (ymax - ymin) * h)
                    iy_end = int((ymax - crop_ymin) / (ymax - ymin) * h)
                    
                    # 限制索引范围并确保至少 2x2
                    ix_start = max(0, min(w - 2, ix_start))
                    ix_end = max(ix_start + 2, min(w, ix_end))
                    
                    iy_start = max(0, min(h - 2, iy_start))
                    iy_end = max(iy_start + 2, min(h, iy_end))
                    
                    # 重新计算实际的 crop 范围，以匹配调整后的索引
                    real_crop_xmin = xmin + (ix_start / w) * (xmax - xmin)
                    real_crop_xmax = xmin + (ix_end / w) * (xmax - xmin)
                    real_crop_ymax = ymax - (iy_start / h) * (ymax - ymin)
                    real_crop_ymin = ymax - (iy_end / h) * (ymax - ymin)

                    if ix_end > ix_start and iy_end > iy_start:
                        sub_data = elevation_data[iy_start:iy_end, ix_start:ix_end]
                        
                        # 生成网格
                        sub_x = np.linspace(real_crop_xmin, real_crop_xmax, ix_end - ix_start)
                        sub_y = np.linspace(real_crop_ymax, real_crop_ymin, iy_end - iy_start)
                        X, Y = np.meshgrid(sub_x, sub_y)
                        
                        # 使用 terrain colormap
                        ax.plot_surface(X, Y, sub_data, cmap='terrain', alpha=0.5, linewidth=0, antialiased=False, zorder=0)
                        print(f"Plotted 3D elevation surface (grid: {sub_data.shape}, range: x[{real_crop_xmin:.4f}, {real_crop_xmax:.4f}], y[{real_crop_ymin:.4f}, {real_crop_ymax:.4f}])")
                    else:
                        print("Warning: Cropped elevation area is empty.")
                else:
                    print("Warning: Trajectory is outside elevation map extent.")

            except Exception as e:
                print(f"Warning: failed to plot 3D elevation: {e}")
                
    else:
        ax = fig.add_subplot(111)

    # 绘制背景高程图（如果提供）
    if bg_img is not None and bg_extent is not None and not need_3d:
        try:
            # 确保 bg_img 为二维数组
            img = bg_img
            if img.ndim == 3:
                img = img[:, :, 0]
            extent = bg_extent  # [xmin, xmax, ymin, ymax]
            ax.imshow(img, cmap=bg_cmap, extent=extent, origin='upper', alpha=0.7, zorder=0)
        except Exception as e:
            print(f"Warning: failed to draw background elevation image: {e}")

    # 绘制路径点
    if waypoints:
        # 支持 (x,y) 或 (x,y,z)
        if need_3d:
            wp_x = [p[0] for p in waypoints]
            wp_y = [p[1] for p in waypoints]
            wp_z = [p[2] if len(p) >= 3 else 0 for p in waypoints]
            ax.scatter(wp_x, wp_y, wp_z, c='red', s=80, marker='o', label='Waypoints', zorder=5)
            ax.plot(wp_x, wp_y, wp_z, 'r--', alpha=0.7, linewidth=2, label='Planned Path')
            for i, (x, y, *rest) in enumerate(waypoints):
                z = rest[0] if rest else 0
                try:
                    ax.text(x, y, z, f'{i}', fontsize=8, color='red')
                except Exception:
                    pass
        else:
            wp_x, wp_y = zip(*[(p[0], p[1]) for p in waypoints])
            ax.scatter(wp_x, wp_y, c='red', s=80, marker='o', label='Waypoints', zorder=5)
            ax.plot(wp_x, wp_y, 'r--', alpha=0.7, linewidth=2, label='Planned Path')
            for i, (x, y, z) in enumerate(waypoints):  #高度不用,但需要读取
                try:
                    ax.annotate(f'{i}', (x, y), xytext=(5, 5), textcoords='offset points', fontsize=8, color='red')
                except Exception:
                    pass

    # 绘制 leader 轨迹（优先用蓝色）
    if leader_traj:
        try:
            if need_3d:
                traj_x = [p[0] for p in leader_traj]
                traj_y = [p[1] for p in leader_traj]
                traj_z = [p[2] if len(p) >= 3 else 0 for p in leader_traj]
                ax.scatter(traj_x, traj_y, traj_z, c='blue', s=30, marker='.', alpha=0.6, label='Leader Trajectory')
                ax.plot(traj_x, traj_y, traj_z, 'b-', alpha=0.8, linewidth=1.5)
            else:
                traj_x, traj_y = zip(*[(p[0], p[1]) for p in leader_traj])
                ax.scatter(traj_x, traj_y, c='blue', s=30, marker='.', alpha=0.6, label='Leader Trajectory')
                ax.plot(traj_x, traj_y, 'b-', alpha=0.8, linewidth=1.5)
        except Exception:
            pass

    # 绘制 plane_trajs（多条），使用 colormap 区分
    if plane_trajs:
        cmap = plt.get_cmap('tab10')
        for idx, (pid, pts) in enumerate(plane_trajs):
            if not pts:
                continue
            color = cmap(idx % 10)
            if need_3d:
                xs = [p[0] for p in pts]
                ys = [p[1] for p in pts]
                zs = [p[2] if len(p) >= 3 else 0 for p in pts]
                ax.plot(xs, ys, zs, '-', color=color, linewidth=1.5, alpha=0.9, label=f'Plane {pid}')
                ax.scatter(xs, ys, zs, c=[color], s=20, marker='.', alpha=0.9)
                try:
                    ax.text(xs[0], ys[0], zs[0], f'id:{pid}', fontsize=8, color=color)
                except Exception:
                    pass
            else:
                xs, ys = zip(*[(p[0], p[1]) for p in pts])
                ax.plot(xs, ys, '-', color=color, linewidth=1.5, alpha=0.9, label=f'Plane {pid}')
                ax.scatter(xs, ys, c=[color], s=20, marker='.', alpha=0.9)
                # 在首点处标注 id
                try:
                    ax.annotate(f'id:{pid}', (xs[0], ys[0]), xytext=(4, 4), textcoords='offset points', fontsize=8, color=color)
                except Exception:
                    pass

    ax.set_xlabel('Longitude / X Coordinate')
    ax.set_ylabel('Latitude / Y Coordinate')
    if need_3d:
        try:
            ax.set_zlabel('Altitude / Z')
        except Exception:
            pass
    plt.title(title)
    plt.legend()
    try:
        ax.grid(True, alpha=0.3)
    except Exception:
        pass

    # 3D 不支持 axis('equal') 的通用实现，这里仅在 2D 时使用
    if not need_3d:
        try:
            ax.set_aspect('equal', 'box')
        except Exception:
            pass

    plt.tight_layout()
    # 如果指定了保存路径，先保存图像
    if save_path:
        try:
            plt.savefig(save_path, dpi=200)
            print(f"Saved plot to: {save_path}")
        except Exception as e:
            print(f"Warning: failed to save plot to {save_path}: {e}")

    if show_plot:
        try:
            plt.show()
        except Exception:
            # 在无 GUI 环境中，plt.show() 可能失败，已保存图像则可忽略
            print("Note: unable to show plot (no display). If save_path provided, image was saved.")
            plt.close()
    else:
        plt.close()

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

    # 提取路径点（leader 中的中间航点）
    waypoints = extract_coordinates(input_data, "leader_midway_point_wgs84")
    print(f"Found {len(waypoints)} waypoints")

    # 提取 leader 轨迹（如果存在）
    leader_traj = extract_coordinates(output_data, "uav_leader_plane1")
    print(f"Found {len(leader_traj)} leader trajectory points")

    # 提取所有 uav_plane* 的轨迹（支持多条，每条带 id）
    plane_trajs = extract_all_plane_trajectories(output_data)
    print(f"Found {len(plane_trajs)} plane trajectories (uav_plane*)")

    if not waypoints and not leader_traj and not plane_trajs:
        print("No valid waypoint or trajectory data found")
        return

    # 从文件名中提取 UAV 编号用于标题
    uav_id = os.path.basename(file_path).split('_')[0]

    # 可视化
    print("Generating visualization...")
    # 默认把图像保存到与输出 JSON 同名的 PNG 文件
    try:
        save_path = os.path.splitext(output_file)[0] + '.png'
    except Exception:
        save_path = None
    if save_path:
        print(f"Will save plot to: {save_path}")
    # 尝试加载 data 目录下的 neimeng 高程 PGM，并解析 aux.xml 获取 GeoTransform
    bg_img = None
    bg_extent = None
    elevation_data = None
    elevation_extent = None
    
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    # fallback to ./data if __file__ resolution fails
    if not os.path.isdir(data_dir):
        data_dir = './data'
    
    # Load .tif.ovr for 3D elevation
    ovr_path = os.path.join(data_dir, 'neimeng.tif.ovr')
    if os.path.exists(ovr_path):
        print(f"Loading 3D elevation data from: {ovr_path}")
        elevation_data, elevation_extent = load_tiff_elevation(ovr_path)
    
    pgm_path = os.path.join(data_dir, 'neimeng.tif.band1.pgm')
    aux_path = os.path.join(data_dir, 'neimeng.tif.aux.xml')
    if os.path.exists(pgm_path):
        print(f"Loading elevation image from: {pgm_path}")
        img = load_pgm_image(pgm_path)
        if img is not None:
            bg_img = img
            gt = parse_aux_xml_for_geotransform(aux_path) if os.path.exists(aux_path) else None
            h, w = img.shape
            if gt is not None:
                origin_x = gt[0]
                pixel_w = gt[1]
                origin_y = gt[3]
                pixel_h = gt[5]
                xmin = origin_x
                xmax = origin_x + pixel_w * w
                # 根据 pixel_h 正负计算 ymin/ymax
                if pixel_h < 0:
                    ymax = origin_y
                    ymin = origin_y + pixel_h * h
                else:
                    ymin = origin_y
                    ymax = origin_y + pixel_h * h
                bg_extent = [xmin, xmax, ymin, ymax]
                print(f"Parsed GeoTransform from aux: extent={bg_extent}")
            else:
                # 根据数据中的点范围来推断 extent（略放大一点）
                all_coords = []
                if waypoints:
                    all_coords += [(p[0], p[1]) for p in waypoints]
                if leader_traj:
                    all_coords += [(p[0], p[1]) for p in leader_traj]
                for (_id, pts) in plane_trajs:
                    all_coords += [(p[0], p[1]) for p in pts]
                if all_coords:
                    xs = [c[0] for c in all_coords]
                    ys = [c[1] for c in all_coords]
                    xmin, xmax = min(xs), max(xs)
                    ymin, ymax = min(ys), max(ys)
                    # 扩展一些边距
                    padx = max(1e-6, (xmax - xmin) * 0.1)
                    pady = max(1e-6, (ymax - ymin) * 0.1)
                    bg_extent = [xmin - padx, xmax + padx, ymin - pady, ymax + pady]
                    print(f"Using coordinates-derived extent for background: {bg_extent}")
                else:
                    # 无法确定坐标系，使用像素坐标范围
                    bg_extent = [0, w, 0, h]
                    print("No trajectory coordinates found; using pixel coordinates for elevation extent")
        else:
            print(f"Warning: failed to load PGM image: {pgm_path}")
    else:
        print(f"Elevation PGM not found at {pgm_path}, skipping background image")

    plot_path_and_trajectory(waypoints=waypoints, leader_traj=leader_traj if leader_traj else None, plane_trajs=plane_trajs, title=f"{uav_id} Path Planning and Execution Trajectory", save_path=save_path, bg_img=bg_img, bg_extent=bg_extent, elevation_data=elevation_data, elevation_extent=elevation_extent)

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
    # 支持可选的第二个参数来强制 3D/2D：
    #   python3 visible.py <file>            -> 自动检测 2D/3D
    #   python3 visible.py <file> 3d|--3d    -> 强制 3D
    #   python3 visible.py <file> 2d|--2d    -> 强制 2D
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Error: Please provide file path as argument")
        print_usage()
        sys.exit(1)

    file_path = sys.argv[1]
    force_3d = None
    if len(sys.argv) == 3:
        flag = sys.argv[2].lower()
        if flag in ('3d', '--3d'):
            force_3d = True
        elif flag in ('2d', '--2d'):
            force_3d = False
        else:
            print(f"Warning: unknown flag '{sys.argv[2]}', ignoring and using auto-detect")
            force_3d = None

    # 如果不强制指定，则直接使用 main 的自动检测流程
    if force_3d is None:
        # 如果需要分析数据结构，取消下面的注释
        # analyze_data_structure(f"{file_path}_input.json")
        # analyze_data_structure(f"{file_path}_output.json")
        main(file_path)
    else:
        # 强制 3D/2D：读取同 main 的文件并传入 plot_3d 标志
        setup_chinese_font()
        input_file, output_file = auto_detect_files(file_path)
        if not input_file:
            print("Error: Input file not found")
            print_usage()
            sys.exit(1)
        if not output_file:
            print("Error: Output file not found")
            sys.exit(1)

        input_data = read_json_file(input_file)
        output_data = read_json_file(output_file)
        waypoints = extract_coordinates(input_data, "leader_midway_point_wgs84")
        leader_traj = extract_coordinates(output_data, "uav_leader_plane1")
        plane_trajs = extract_all_plane_trajectories(output_data)

        uav_id = os.path.basename(file_path).split('_')[0]
        try:
            save_path = os.path.splitext(output_file)[0] + '.png'
        except Exception:
            save_path = None

        # Load elevation data if needed (especially for 3D)
        elevation_data = None
        elevation_extent = None
        bg_img = None
        bg_extent = None
        
        data_dir = os.path.join(os.path.dirname(__file__), 'data')
        if not os.path.isdir(data_dir):
            data_dir = './data'
            
        if force_3d:
            ovr_path = os.path.join(data_dir, 'neimeng.tif.ovr')
            if os.path.exists(ovr_path):
                print(f"Loading 3D elevation data from: {ovr_path}")
                elevation_data, elevation_extent = load_tiff_elevation(ovr_path)
        
        # Also load PGM for 2D fallback or if needed
        pgm_path = os.path.join(data_dir, 'neimeng.tif.band1.pgm')
        if os.path.exists(pgm_path):
             bg_img = load_pgm_image(pgm_path)

        print(f"Generating {'3D' if force_3d else '2D'} visualization (forced)...")
        plot_path_and_trajectory(waypoints=waypoints, leader_traj=leader_traj if leader_traj else None, plane_trajs=plane_trajs, title=f"{uav_id} Path Planning and Execution Trajectory", save_path=save_path, plot_3d=force_3d, elevation_data=elevation_data, elevation_extent=elevation_extent)