# -*- coding: utf-8 -*-
# ##################################################################
#
#   可控制纤维体积分数和多空隙缺陷的Abaqus三维RVE建模脚本
#
# 功能说明:
# 1. 根据设定的纤维体积分数(Vf)和空隙率(Vp)自动生成三维RVE几何模型
# 2. 算法流程:
#    a. [纤维优先] 采用RSA播种 + 锚定松弛算法排布纤维, 保证纤维间(F-F)间距 (2D)
#    b. [空隙后置] 采用3D RSA算法排布空隙, 播种时同时检查:
#       i.  空隙-空隙 (V-V) 间距 (3D 中心-中心 27盒子检测)
#       ii. 空隙-纤维 (F-V) 间距 (2D 中心-中心 9盒子(XY)检测)
# 3. 自动创建材料、赋予截面、划分网格
# 4. 施加三维周期性边界条件(PBCs) - [通过ABAQUS Micromechanics Plugin]
# 5. 在纤维-基体界面插入三维cohesive单元(COH3D8/COH3D6)
# 6. 导出纤维中心和空隙信息为CSV文件
#
# 参考文献:
# [1] Li, X., et al. (2025). A novel method for constructing 3D void RVE
#     elements and rapid homogenization of composite materials.
#     Composite Structures, 360, 119040.
#
# 作者: 刘正鹏 (Liu Zhengpeng)
# 版本: v1.1 (MP-PBC + 优化的排布逻辑)
# 创建日期: 2025-09-30
# 适用软件: ABAQUS 2023 (需安装 Micromechanics Plugin)
# Python版本: 2.7 (ABAQUS内置)
# 技术交流: GitHub/CSDN/知乎 @小盆i
# 联系方式: 1370872708@qq.com / Zhengpeng0105@gmail.com
#
# ##################################################################

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import math
import random as rd
import time
import os
import sys

from part import *
from material import *
from section import *
from sketch import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from visualization import *
from connectorBehavior import *

executeOnCaeStartup()


# =================================================================
#                 3D RSA 算法 (用于空隙排布)
# =================================================================

def _get_3d_periodic_distance_sq(p1, p2, rveSize):
    """计算三维周期性空间中两个点之间的最短距离的平方"""
    dx = abs(p1[0] - p2[0])
    dy = abs(p1[1] - p2[1])
    dz = abs(p1[2] - p2[2])

    if dx > rveSize[0] / 2.0:
        dx = rveSize[0] - dx
    if dy > rveSize[1] / 2.0:
        dy = rveSize[1] - dy
    if dz > rveSize[2] / 2.0:
        dz = rveSize[2] - dz

    return dx * dx + dy * dy + dz * dz


def _generate_void_coordinates_rsa(num_voids, rveSize, min_void_dist_3D_VV,
                                   fiber_centers, min_fiber_void_dist_2D_FV,
                                   max_attempts=5000):
    """
    使用3D RSA(随机顺序吸附)算法生成空隙
    V-V 检查: 3D 中心-中心 27盒子干涉检测
    F-V 检查: 2D 中心-中心 9盒子(XY)干涉检测

    参数:
        num_voids: 目标空隙数量
        rveSize: RVE尺寸 [L, W, H]
        min_void_dist_3D_VV: 空隙中心点的最小允许距离 (3D, V-V)
        fiber_centers: 已固定的纤维中心 [(x, y), ...]
        min_fiber_void_dist_2D_FV: 纤维-空隙最小中心距离 (2D, F-V)
        max_attempts: 每次放置尝试的最大次数

    返回:
        void_list: [(cx, cy, cz, 0.0, 0.0), ...] (球体不需要旋转)
    """
    print("--- Initializing 3D RSA for void placement (V-V 3D-Center & F-V 2D-Center) ---")
    print("    Target voids: %d" % num_voids)
    print("    Min V-V center distance (3D): %.6f" % min_void_dist_3D_VV)
    print("    Min F-V center distance (2D): %.6f" % min_fiber_void_dist_2D_FV)

    void_list = []
    L, W, H = rveSize[0], rveSize[1], rveSize[2]

    # V-V 检查
    min_dist_sq_VV = min_void_dist_3D_VV ** 2
    existing_void_centers = []  # 用于 V-V 检查

    # F-V 检查
    min_dist_sq_FV = min_fiber_void_dist_2D_FV ** 2
    offsets_2d_xy = []
    for dx in [-L, 0, L]:
        for dy in [-W, 0, W]:
            offsets_2d_xy.append((dx, dy))

    for i in range(num_voids):
        placed = False
        for attempt in range(max_attempts):
            # 1. 生成新空隙候选
            cx_new = rd.uniform(0, L)
            cy_new = rd.uniform(0, W)
            cz_new = rd.uniform(0, H)
            center_new = (cx_new, cy_new, cz_new)

            is_too_close = False

            # 2. V-V 检查 (3D Center-to-Center)
            for center_existing in existing_void_centers:
                dist_sq = _get_3d_periodic_distance_sq(center_new, center_existing, rveSize)
                if dist_sq < min_dist_sq_VV:
                    is_too_close = True
                    break
            if is_too_close:
                continue

            # 3. F-V 检查 (2D Center-to-Center, 9-box XY)
            for fx, fy in fiber_centers:
                # 检查新空隙(cx_new, cy_new)与纤维的9个XY平面镜像(fx_mir, fy_mir)的距离
                for dx, dy in offsets_2d_xy:
                    fx_mir = fx + dx
                    fy_mir = fy + dy

                    # 2D distance check
                    dist_sq_2d = (cx_new - fx_mir) ** 2 + (cy_new - fy_mir) ** 2

                    if dist_sq_2d < min_dist_sq_FV:
                        is_too_close = True
                        break  # (offsets_2d_xy loop)
                if is_too_close:
                    break  # (fiber_centers loop)
            if is_too_close:
                continue  # 尝试新的 attempt

            # 4. 放置
            if not is_too_close:
                existing_void_centers.append(center_new)
                # 球体不需要旋转
                void_list.append((cx_new, cy_new, cz_new, 0.0, 0.0))
                placed = True
                break

        if not placed:
            print("    WARNING: RSA void placement congested. Placed %d of %d voids." % (len(void_list), num_voids))
            break

        if (i + 1) % 50 == 0 or i == num_voids - 1:
            print("    ... Placed %d voids." % (i + 1))

    print("--- Void placement complete. Total voids: %d ---" % len(void_list))
    return void_list


# =================================================================
#                 3D 空间间距检测辅助
# =================================================================

def determineMirrorOffsets_3D(center, rveSize, threshold):
    """(智能策略) 判断空隙需要在哪些方向上创建周期性镜像"""
    centx, centy, centz = center
    L, W, H = rveSize

    # 独立判断每个轴是否靠近边界
    near_x = centx < threshold or centx > L - threshold
    near_y = centy < threshold or centy > W - threshold
    near_z = centz < threshold or centz > H - threshold

    if not (near_x or near_y or near_z):
        return []

    x_offsets = [-L, 0, L] if near_x else [0]
    y_offsets = [-W, 0, W] if near_y else [0]
    z_offsets = [-H, 0, H] if near_z else [0]

    mirror_offsets = []
    for dx in x_offsets:
        for dy in y_offsets:
            for dz in z_offsets:
                if dx == 0 and dy == 0 and dz == 0:
                    continue
                mirror_offsets.append((dx, dy, dz))

    return mirror_offsets


# =================================================================
#           纤维排布算法 (F-F 约束)
# =================================================================

def _relax_fibers_ff_only(initial_coords, seeding_count, fiberCount,
                          rveSize, minFiberDistance):
    """
    锚定松弛算法
    仅处理 纤维-纤维 (F-F) 间距约束
    """
    print("--- Initializing Anchored Relaxation (F-F only) ---")

    coords = [list(c) for c in initial_coords]
    max_iterations = 20000
    movement_factor = 0.5
    anchor_damping_factor = 0.05
    min_dist_sq_ff = minFiberDistance ** 2

    rveW, rveH = rveSize[0], rveSize[1]

    for iter_num in range(max_iterations):
        max_movement_sq = 0.0
        net_forces = [[0.0, 0.0] for _ in range(fiberCount)]
        any_overlap_this_iteration = False

        # 1. 纤维间的相互作用力 (2D 周期性)
        for i in range(fiberCount):
            for j in range(i + 1, fiberCount):
                dx = coords[j][0] - coords[i][0]
                dy = coords[j][1] - coords[i][1]
                if dx > rveW / 2: dx -= rveW
                if dx < -rveW / 2: dx += rveW
                if dy > rveH / 2: dy -= rveH
                if dy < -rveH / 2: dy += rveH
                dist_sq = dx * dx + dy * dy
                if dist_sq < min_dist_sq_ff:
                    any_overlap_this_iteration = True
                    dist = math.sqrt(dist_sq) if dist_sq > 0 else 1e-9
                    overlap = minFiberDistance - dist
                    force_magnitude = overlap
                    force_x = force_magnitude * (dx / dist)
                    force_y = force_magnitude * (dy / dist)
                    net_forces[i][0] -= force_x
                    net_forces[i][1] -= force_y
                    net_forces[j][0] += force_x
                    net_forces[j][1] += force_y

        # 2. 应用力和检查收敛
        max_movement_sq = 0.0  # 在循环内重置
        for i in range(fiberCount):
            move_x = net_forces[i][0] * movement_factor
            move_y = net_forces[i][1] * movement_factor
            if i < seeding_count:
                move_x *= anchor_damping_factor
                move_y *= anchor_damping_factor
            coords[i][0] = (coords[i][0] + move_x) % rveW
            coords[i][1] = (coords[i][1] + move_y) % rveH
            current_movement_sq = move_x ** 2 + move_y ** 2
            if current_movement_sq > max_movement_sq:
                max_movement_sq = current_movement_sq

        if (iter_num + 1) % 100 == 0:
            print("... F-F Relaxation Iteration %d, Max movement: %.2e" %
                  (iter_num + 1, math.sqrt(max_movement_sq)))

        # 检查收敛 (必须同时满足: 1. 移动量足够小 AND 2. 没有任何重叠)
        if (max_movement_sq < (1e-6 * minFiberDistance) ** 2) and (not any_overlap_this_iteration):
            print("--- F-F Relaxation Converged after %d iterations (Movement settled AND no overlaps). ---" % (
                    iter_num + 1))
            break

    if iter_num == max_iterations - 1:
        print("--- Max F-F relaxation iterations reached. Proceeding to validation. ---")

    # ==================== 最终验证 (带微小公差) ====================
    print("--- Initializing F-F Final Validation Check ---")
    final_violations = []

    # 检查纤维 (2D)
    for i in range(fiberCount):
        for j in range(i + 1, fiberCount):
            dx = coords[j][0] - coords[i][0];
            dy = coords[j][1] - coords[i][1]
            if dx > rveW / 2: dx -= rveW
            if dx < -rveW / 2: dx += rveW
            if dy > rveH / 2: dy -= rveH
            if dy < -rveH / 2: dy += rveH
            dist_sq = dx * dx + dy * dy

            # 引入一个微小的数值公差 (epsilon) 来处理浮点数精度问题
            if dist_sq < min_dist_sq_ff * (1.0 - 1e-12):
                final_violations.append("Fibers %d-%d too close: %.6f (Required: %.6f)" %
                                        (i, j, math.sqrt(dist_sq), minFiberDistance))

    if final_violations:
        error_msg = ("FATAL: Cannot satisfy F-F constraints after %d iterations!\n"
                     "Violations:\n  %s\n"
                     "Please try:\n"
                     "  1. Reduce target fiber volume fraction\n"
                     "  2. Increase RVE size" %
                     (max_iterations, "\n  ".join(final_violations[:10])))
        raise Exception(error_msg)

    print("--- F-F Final Validation PASSED. ---")
    return coords


def _generate_fiber_layout(fiberCount, rveSize, minFiberDistance, rsa_seeding_ratio):
    """
    生成纤维布局的主函数 (仅 F-F 约束)
    1. RSA 播种
    2. 随机放置
    3. 松弛算法
    """
    rsa_max_attempts = 500
    rveW, rveH = rveSize[0], rveSize[1]
    min_dist_sq_ff = minFiberDistance ** 2

    if fiberCount == 0:
        return []

    print("\n--- Stage 1: RSA Seeding (F-F Only Constraints) ---")
    seeding_count = int(fiberCount * rsa_seeding_ratio)
    seeded_coords = []
    for i in range(seeding_count):
        placed = False
        for _ in range(rsa_max_attempts):
            xt = rd.uniform(0, rveW)
            yt = rd.uniform(0, rveH)
            is_too_close = False

            # 检查1: 纤维-纤维 间距 (2D 周期性)
            for xc, yc in seeded_coords:
                dx = abs(xt - xc);
                dy = abs(yt - yc)
                if dx > rveW / 2: dx = rveW - dx
                if dy > rveH / 2: dy = rveH - dy
                if dx * dx + dy * dy < min_dist_sq_ff:
                    is_too_close = True
                    break
            if is_too_close:
                continue  # 尝试新的 (xt, yt)

            # 安全, 放置
            seeded_coords.append((xt, yt))
            placed = True
            break

        if not placed:
            print("RSA seeding congested at %d fibers (Target: %d)." % (len(seeded_coords), seeding_count))
            break

    seeding_count_actual = len(seeded_coords)
    print("RSA placed %d anchor fibers." % seeding_count_actual)

    print("\n--- Stage 2: Random Placement ---")
    remaining = fiberCount - seeding_count_actual
    initial_coords = list(seeded_coords)
    for _ in range(remaining):
        initial_coords.append((rd.uniform(0, rveW), rd.uniform(0, rveH)))

    print("\n--- Stage 3: F-F Relaxation & Enforcement ---")
    try:
        fiber_centers = _relax_fibers_ff_only(
            initial_coords, seeding_count_actual, fiberCount,
            rveSize, minFiberDistance
        )
        return fiber_centers
    except Exception as e:
        raise e  # 将异常向上传递


# =================================================================
#                 CSV坐标输出模块
# =================================================================

def exportVoidGeometryToCSV(filename, void_list, void_radius, rveVolume, target_porosity):
    """将多个空隙的几何信息导出为CSV文件"""
    try:
        work_dir = os.getcwd()
        filepath = os.path.join(work_dir, filename)

        num_voids = len(void_list)
        # 体积是 4/3 * pi * R^3
        vol_per_void = (4.0 / 3.0) * math.pi * (void_radius ** 3)
        actual_volume = num_voids * vol_per_void
        actual_porosity = actual_volume / rveVolume

        with open(filepath, 'w') as f:
            f.write("# 3D RVE Multi-Void Geometry Information\n")
            f.write("# Generated: %s\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
            f.write("# Void Type: Sphere\n")
            f.write("# Target Porosity: %.4f%%\n" % (target_porosity * 100))
            f.write("# RVE Volume: %.8f\n" % rveVolume)
            f.write("# --- Per Void --- \n")
            f.write("# Radius: %.6f\n" % (void_radius))
            f.write("# Volume per Void: %.8f\n" % vol_per_void)
            f.write("# --- Total --- \n")
            f.write("# Void Count (Placed): %d\n" % num_voids)
            f.write("# Total Void Volume (Actual): %.8f\n" % actual_volume)
            f.write("# Actual Porosity: %.4f%%\n" % (actual_porosity * 100))
            f.write("#\n")
            f.write("Void_ID,Center_X,Center_Y,Center_Z,Radius\n")
            for i, (cx, cy, cz, _, _) in enumerate(void_list, start=1):
                f.write("%d,%.8f,%.8f,%.8f,%.8f\n" % (i, cx, cy, cz, void_radius))

        print("\n" + "=" * 60)
        print("SUCCESS: Void geometry exported to CSV")
        print("File location: %s" % filepath)
        print("=" * 60 + "\n")
        return {'actual_porosity': actual_porosity, 'count': num_voids, 'actual_volume': actual_volume}
    except Exception as e:
        print("\nWARNING: Failed to export void geometry: %s\n" % str(e))
        return None


def exportFiberCentersToCSV(fiber_centers, filename, rveSize, fiberRadius,
                            depth, target_Vf, void_info=None):
    """将纤维中心坐标导出为CSV文件"""
    try:
        work_dir = os.getcwd()
        filepath = os.path.join(work_dir, filename)
        with open(filepath, 'w') as f:
            f.write("# 3D RVE Fiber Center Coordinates (Unidirectional)\n")
            f.write("# Generated: %s\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
            f.write("# RVE Size (Width x Height x Depth): %.6f x %.6f x %.6f\n" %
                    (rveSize[0], rveSize[1], depth))
            f.write("# Fiber Radius: %.6f\n" % fiberRadius)
            f.write("# Fiber Length: %.6f\n" % depth)
            f.write("# Target Fiber Volume Fraction: %.4f\n" % target_Vf)
            f.write("# Total Fiber Count: %d\n" % len(fiber_centers))
            if void_info:
                f.write("# Void Present: Yes\n")
                f.write("# Void Porosity (Actual): %.4f%%\n" % (void_info.get('actual_porosity', 0.0) * 100))
                f.write("# Void Count: %d\n" % void_info.get('count', 0))
            else:
                f.write("# Void Present: No\n")
            f.write("#\n")
            f.write("Fiber_ID,X_Coordinate,Y_Coordinate,Z_Start,Z_End\n")
            for i, (x, y) in enumerate(fiber_centers, start=1):
                f.write("%d,%.8f,%.8f,%.8f,%.8f\n" % (i, x, y, 0.0, depth))
        print("\n" + "=" * 60)
        print("SUCCESS: Fiber coordinates exported to CSV")
        print("File location: %s" % filepath)
        print("Total fibers: %d" % len(fiber_centers))
        if void_info:
            print("Actual porosity: %.4f%%" % (void_info.get('actual_porosity', 0.0) * 100))
        print("=" * 60 + "\n")
        return filepath
    except Exception as e:
        print("\nWARNING: Failed to export coordinates: %s\n" % str(e))
        return None


# =================================================================
#              距离验证模块
# =================================================================

def verifyMinimumFiberDistance3D(fiber_centers, rveSize, fiberRadius, minFiberDistance):
    """验证纤维间的最小距离(XY平面)"""
    fiberCount = len(fiber_centers)

    print("\n" + "=" * 70)
    print("MINIMUM FIBER DISTANCE VERIFICATION (XY Plane)")
    print("=" * 70)
    print("RVE Size: %.6f x %.6f x %.6f" % (rveSize[0], rveSize[1], rveSize[2]))
    print("Fiber Radius: %.6f mm" % fiberRadius)
    print("Required Minimum Center Distance: %.6f mm" % minFiberDistance)
    print("Total Fibers: %d" % fiberCount)

    min_distance_found = float('inf')
    violating_pairs = 0

    # 验证最小距离的平方
    minFiberDistance_sq = minFiberDistance ** 2

    if fiberCount < 2:
        print("Not enough fibers to verify.")
        return True, {}

    for i in range(fiberCount):
        for j in range(i + 1, fiberCount):
            dx = fiber_centers[j][0] - fiber_centers[i][0]
            dy = fiber_centers[j][1] - fiber_centers[i][1]
            if dx > rveSize[0] / 2.0: dx -= rveSize[0]
            if dx < -rveSize[0] / 2.0: dx += rveSize[0]
            if dy > rveSize[1] / 2.0: dy -= rveSize[1]
            if dy < -rveSize[1] / 2.0: dy += rveSize[1]

            distance_sq = dx * dx + dy * dy
            if distance_sq < min_distance_found:
                min_distance_found = distance_sq

            # 引入一个微小的数值公差 (epsilon) 来处理浮点数精度问题
            if distance_sq < minFiberDistance_sq * (1.0 - 1e-12):
                violating_pairs += 1

    min_distance_found = math.sqrt(min_distance_found)
    print("\n" + "-" * 70)
    print("DISTANCE STATISTICS")
    print("-" * 70)
    print("  Minimum Distance Found: %.6f mm" % min_distance_found)

    if violating_pairs == 0:
        print("VERIFICATION RESULT: PASSED")
        verification_passed = True
    else:
        print("VERIFICATION RESULT: FAILED (%d pairs violation)" % violating_pairs)
        verification_passed = False
    print("=" * 70 + "\n")

    stats = {
        'fiber_count': fiberCount,
        'min_distance': min_distance_found,
        'required_distance': minFiberDistance,
        'violations_count': violating_pairs
    }
    return verification_passed, stats


# =================================================================
#       三维纤维-基体分类算法
# =================================================================

def buildAllPeriodicCenters3D(centers, rveSize, radius):
    """构建包含周期性镜像的完整组分中心列表 (XY平面)"""
    all_centers = []
    rveW, rveH = rveSize[0], rveSize[1]
    for xt, yt in centers:
        all_centers.append((xt, yt))
        if xt < radius: all_centers.append((xt + rveW, yt))
        if xt > rveW - radius: all_centers.append((xt - rveW, yt))
        if yt < radius: all_centers.append((xt, yt + rveH))
        if yt > rveH - radius: all_centers.append((xt, yt - rveH))
        if xt < radius and yt < radius and math.sqrt(xt ** 2 + yt ** 2) < radius:
            all_centers.append((xt + rveW, yt + rveH))
        if xt < radius and yt > rveH - radius and math.sqrt(xt ** 2 + (rveH - yt) ** 2) < radius:
            all_centers.append((xt + rveW, yt - rveH))
        if xt > rveW - radius and yt > rveH - radius and math.sqrt((rveW - xt) ** 2 + (rveH - yt) ** 2) < radius:
            all_centers.append((xt - rveW, yt - rveH))
        if xt > rveW - radius and yt < radius and math.sqrt((rveW - xt) ** 2 + yt ** 2) < radius:
            all_centers.append((xt - rveW, yt + rveH))
    return all_centers


def getCellCenterFromVertices(cell):
    """通过体单元的顶点坐标计算几何中心"""
    try:
        vertices = cell.getVertices()
        if not vertices:
            return None
        x_coords, y_coords, z_coords = [], [], []
        for vertex in vertices:
            try:
                if hasattr(vertex, 'pointOn'):
                    coord = vertex.pointOn[0]
                    x_coords.append(coord[0])
                    y_coords.append(coord[1])
                    z_coords.append(coord[2])
            except:
                continue
        if x_coords:
            return (sum(x_coords) / len(x_coords),
                    sum(y_coords) / len(y_coords),
                    sum(z_coords) / len(z_coords))
        return None
    except Exception as e:
        return None


def classifyCellsImproved(all_cells, rveSize, rveVolume,
                          fiber_centers, fiberRadius):
    """
    改进的体单元分类算法: 识别纤维、基体
    (注: 孔隙是通过布尔Cut实现的, 它们不是独立的cell, 故不在此处寻找)
    """
    print("\n" + "=" * 70)
    print("CELL CLASSIFICATION")
    print("=" * 70)
    print("Total cells to classify: %d" % len(all_cells))

    # 步骤1: 按体积排序 (最大的通常是基体)
    try:
        sorted_cells = sorted(all_cells, key=lambda c: c.getSize(), reverse=True)
    except Exception as e:
        print("FATAL ERROR: Failed to sort cells by size. %s" % str(e))
        return [], []

    matrix_cells_list = [sorted_cells[0]]
    potential_cells = sorted_cells[1:]
    fiber_cells_list = []
    print("    Main matrix cell identified (Volume: %.6e)" % sorted_cells[0].getSize())
    print("    Cells to classify: %d" % len(potential_cells))

    # 步骤2: 构建完整的周期性纤维中心列表
    all_fiber_centers = buildAllPeriodicCenters3D(fiber_centers, rveSize, fiberRadius)
    print("    Total fiber centers (with periodicity): %d" % len(all_fiber_centers))

    # 步骤3: 详细分类
    print("\n  Step 3: Detailed classification using geometry")
    matrix_fragment_count = 0
    unclassified_cells = []

    for idx, cell in enumerate(potential_cells):
        cell_center = None
        try:
            if hasattr(cell, 'pointOn') and cell.pointOn:
                point = cell.pointOn[0]
                cell_center = (point[0], point[1], point[2])
        except:
            pass
        if cell_center is None:
            try:
                cell_centroid = cell.getCentroid()
                if cell_centroid and len(cell_centroid) >= 3:
                    cell_center = (cell_centroid[0], cell_centroid[1], cell_centroid[2])
            except:
                pass
        if cell_center is None:
            try:
                cell_center = getCellCenterFromVertices(cell)
            except:
                pass
        if cell_center is None:
            print("  WARNING: Unable to determine cell center for cell %d (idx). Skipping." % idx)
            unclassified_cells.append(cell)
            continue

        cell_x, cell_y = cell_center[0], cell_center[1]

        # 计算到最近纤维中心的XY距离
        min_dist_fib = float('inf')
        for fc_x, fc_y in all_fiber_centers:
            dist = math.sqrt((cell_x - fc_x) ** 2 + (cell_y - fc_y) ** 2)
            if dist < min_dist_fib:
                min_dist_fib = dist

        # 判断归属 (孔隙是Cut掉的, 不会在这里被找到)
        if min_dist_fib < fiberRadius:
            fiber_cells_list.append(cell)
        else:
            matrix_cells_list.append(cell)
            matrix_fragment_count += 1

    # 步骤4: 输出分类结果
    print("\n" + "-" * 70)
    print("CLASSIFICATION RESULTS")
    print("-" * 70)
    print("  Fiber cells identified: %d" % len(fiber_cells_list))
    print("  Matrix cells total: %d (1 main + %d fragments)" %
          (len(matrix_cells_list), matrix_fragment_count))
    if unclassified_cells:
        print("  Unclassified cells: %d (will be treated as matrix)" % len(unclassified_cells))
        matrix_cells_list.extend(unclassified_cells)

    # 步骤5: 验证体积分数
    fiber_total_volume = sum([c.getSize() for c in fiber_cells_list])
    actual_Vf_fib = fiber_total_volume / rveVolume
    target_Vf_fib_count = (len(fiber_centers) * math.pi * fiberRadius ** 2 * rveSize[2]) / rveVolume

    print("\n  Fiber Vf Validation:")
    print("    Target Vf (from input): %.4f" % target_Vf_fib_count)
    print("    Actual Vf (from cells): %.4f" % actual_Vf_fib)

    if len(fiber_cells_list) != len(all_fiber_centers):
        print("\n  WARNING: Classified fiber cell count (%d) does not match"
              " periodic fiber count (%d)." % (len(fiber_cells_list), len(all_fiber_centers)))
        print("  This may be due to boolean operation failures or classification errors.")

    print("=" * 70 + "\n")
    return fiber_cells_list, matrix_cells_list


# =================================================================
#         主建模函数(集成多空隙功能)
# =================================================================

def create3DRVEModelWithVoids(modelName='RVE_3D_with_Voids',
                              rveSize=[1.0, 1.0, 1.0],
                              fiberRadius=0.05,
                              target_Vf=0.50,
                              fiber_min_dist_factor=2.05,
                              globalSeedSize=0.02,
                              deviationFactor=0.1,
                              minSizeFactor=0.1,
                              rsa_seeding_ratio=0.3,
                              export_coordinates=True,
                              csv_filename=None,
                              # 空隙相关参数
                              enable_void=True,
                              void_porosity=0.02,
                              void_sphere_radius=0.01,
                              min_void_dist_factor=2.0,
                              fiber_void_min_dist_factor=1.05,
                              void_boundary_threshold_factor=1.0,
                              export_void_geometry=True,
                              # 纤维材料参数 (GPa)
                              fiber_E1=230.0,
                              fiber_E2=15.0,
                              fiber_G12=15.0,
                              fiber_G23=7.0,
                              fiber_nu12=0.2,
                              # 基体材料参数 (MPa)
                              matrix_E=3170.0,
                              matrix_nu=0.35,
                              matrix_friction_angle=16.0,
                              matrix_flow_stress_ratio=1.0,
                              matrix_dilation_angle=16.0,
                              matrix_hardening_yield=106.4,
                              matrix_hardening_plastic_strain=0.0,
                              matrix_damage_strain=0.01,
                              matrix_damage_stress_triax=0.0,
                              matrix_damage_strain_rate=0.0,
                              matrix_damage_displacement=5e-05,
                              # 界面材料参数
                              cohesive_K_nn=1e8,
                              cohesive_K_ss=1e8,
                              cohesive_K_tt=1e8,
                              cohesive_t_n=44.0,
                              cohesive_t_s=82.0,
                              cohesive_t_t=82.0,
                              cohesive_GIC=0.001,
                              cohesive_GIIC=0.002,
                              cohesive_GIIIC=0.002,
                              cohesive_eta=1.5,
                              cohesive_stab_coeff=0.0001,
                              # 作业和插件参数
                              job_name='Job-1',
                              num_cpus=8,
                              abaqus_plugin_path='d:/Abaqus2023/Plugins'):
    """创建带多空隙缺陷的三维RVE模型主函数"""

    print("\n" + "=" * 70)
    print("Starting 3D RVE Model Generation (WITH MULTI-VOID DEFECTS)")
    print("Algorithm: Fibers First (F-F), Voids Second (V-V 3D-Center & F-V 2D-Center).")
    print("PBC: ABAQUS Micromechanics Plugin.")
    print("=" * 70)

    # ==================== 步骤 1: 计算几何参数 (纤维) ====================
    print("\nStep 1: Calculating geometric parameters (Fibers)...")

    rveVolume = rveSize[0] * rveSize[1] * rveSize[2]
    depth = rveSize[2]

    # --- 纤维参数 ---
    rveArea = rveSize[0] * rveSize[1]
    fiberArea = math.pi * fiberRadius ** 2
    fiberCount = int(round((target_Vf * rveArea) / fiberArea))
    minFiberDistance = fiber_min_dist_factor * fiberRadius

    print("  Target Vf: %.4f" % target_Vf)
    print("  Calculated Fiber Count: %d" % fiberCount)
    print("  Min Fiber-Fiber Center Distance (2D): %.6f" % minFiberDistance)

    # ==================== 步骤 2: 生成纤维布局 (F-F) ====================
    print("\nStep 2: Generating fiber layout (F-F constraints)...")

    fiber_centers = []
    try:
        fiber_centers = _generate_fiber_layout(
            fiberCount, rveSize, minFiberDistance, rsa_seeding_ratio
        )
    except Exception as e:
        print("\n" + "#" * 70);
        print(str(e));
        print("#" * 70 + "\n");
        return  # 退出主函数

    if not fiber_centers and fiberCount > 0:
        print("FATAL: Fiber generation failed.")
        return

    # (可选) 双重检查
    if fiber_centers:
        verification_passed, _ = verifyMinimumFiberDistance3D(
            fiber_centers, rveSize, fiberRadius, minFiberDistance
        )
        if not verification_passed:
            raise Exception("Minimum fiber distance requirement not satisfied!")

    print("Step 2 Complete. Fiber layout generated.")

    # ==================== 步骤 3: 计算几何参数 (空隙) ====================
    print("\nStep 3: Calculating geometric parameters (Voids)...")

    void_list = []
    void_info = None

    # --- 空隙参数 ---
    # 根据 球体 参数 [Radius]
    void_radius = void_sphere_radius

    # --- F-V 间距 (2D 中心-中心) ---
    minFiberVoidDistance_2D_FV = fiber_void_min_dist_factor * (fiberRadius + void_radius)

    # --- V-V 间距 (3D 中心-中心) ---
    minVoidDistance_3D_VV = min_void_dist_factor * (void_radius)

    # 边界检查 (使用半径)
    void_r_max_overall = void_radius

    num_voids_target = 0
    if enable_void and void_porosity > 0:
        target_total_void_volume = void_porosity * rveVolume
        # 体积是 4/3 * pi * R^3
        vol_per_void = (4.0 / 3.0) * math.pi * (void_radius ** 3)
        if vol_per_void == 0:
            raise ValueError("Void radius cannot be zero.")
        num_voids_target = int(round(target_total_void_volume / vol_per_void))

        print("  Target Void Porosity: %.4f%%" % (void_porosity * 100))
        print("  Void Radius: %.6f" % (void_radius))
        print("  Calculated Number of Voids: %d" % num_voids_target)
        print(
            "  Min Fiber-Void Center-to-Center Distance (2D): %.6f" % minFiberVoidDistance_2D_FV)
        print("  Min Void-Void Center-to-Center Distance (3D): %.6f" % minVoidDistance_3D_VV)
    else:
        print("  Void generation SKIPPED (enable_void=False or porosity=0)")
        enable_void = False

    # ==================== 步骤 4: 生成空隙坐标 (V-V & F-V) ====================
    if enable_void:
        print("\nStep 4: Generating void coordinates (V-V & F-V constraints)...")

        void_list = _generate_void_coordinates_rsa(
            num_voids_target, rveSize,
            minVoidDistance_3D_VV,
            fiber_centers,
            minFiberVoidDistance_2D_FV
        )

        if not void_list:
            print("  WARNING: No voids were placed (congestion or 0 target).")
            enable_void = False  # 如果没生成空隙, 则禁用后续步骤
        else:
            if export_void_geometry:
                void_csv_filename = "VoidGeometry_3D_Por%d_%s.csv" % (
                    int(void_porosity * 10000), time.strftime("%Y%m%d_%H%M%S"))
                void_info = exportVoidGeometryToCSV(void_csv_filename, void_list,
                                                    void_radius, rveVolume, void_porosity)
        print("Step 4 Complete.")
    else:
        print("\nStep 4: Void generation skipped.")

    # ==================== 步骤 5: 导出坐标 (可选) ====================
    print("\nStep 5: Exporting coordinates...")
    if export_coordinates and fiber_centers:
        if csv_filename is None:
            if enable_void and void_info:
                csv_filename = "FiberCenters_3D_Vf%d_Void%d_%s.csv" % (
                    int(target_Vf * 100), int(void_info['actual_porosity'] * 10000),
                    time.strftime("%Y%m%d_%H%M%S"))
            else:
                csv_filename = "FiberCenters_3D_Vf%d_%s.csv" % (
                    int(target_Vf * 100), time.strftime("%Y%m%d_%H%M%S"))
        exportFiberCentersToCSV(fiber_centers, csv_filename,
                                rveSize, fiberRadius, depth, target_Vf, void_info)

    # ==================== 步骤 6: 创建周期性镜像 (纤维) ====================
    print("\nStep 6: Creating periodic mirror coordinates (Fibers)...")
    # 调用函数构建包含周期性镜像的完整纤维中心列表
    all_fiber_centers_for_geom = buildAllPeriodicCenters3D(fiber_centers, rveSize, fiberRadius)
    xCoords = [c[0] for c in all_fiber_centers_for_geom]
    yCoords = [c[1] for c in all_fiber_centers_for_geom]
    print("Total fiber coordinates (including mirrors): %d" % len(xCoords))

    # ==================== 步骤 7: 创建几何 ====================
    print("\nStep 7: Creating Abaqus 3D geometry...")

    if modelName in mdb.models: del mdb.models[modelName]
    mdb.Model(name=modelName, modelType=STANDARD_EXPLICIT)
    model = mdb.models[modelName]
    if 'Model-1' in mdb.models and modelName != 'Model-1' and len(mdb.models['Model-1'].parts) == 0:
        del mdb.models['Model-1']

    # 创建基体部件
    print("  Creating matrix part...")
    s_matrix = model.ConstrainedSketch(name='MatrixSketch', sheetSize=max(rveSize) * 2)
    s_matrix.rectangle(point1=(0.0, 0.0), point2=(rveSize[0], rveSize[1]))
    p_matrix = model.Part(name='Matrix', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p_matrix.BaseSolidExtrude(sketch=s_matrix, depth=depth)

    # 创建纤维部件
    if fiberCount > 0:
        print("  Creating fiber part...")
        s_fiber = model.ConstrainedSketch(name='FiberSketch', sheetSize=max(rveSize) * 3)
        for i in range(len(xCoords)):
            s_fiber.CircleByCenterPerimeter(center=(xCoords[i], yCoords[i]),
                                            point1=(xCoords[i] + fiberRadius, yCoords[i]))
        p_fiber = model.Part(name='Fiber', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        p_fiber.BaseSolidExtrude(sketch=s_fiber, depth=depth)

    # ==================== 步骤 8: 布尔运算 (Merge, Cut) ====================
    print("\nStep 8: Creating assembly using Boolean operations...")
    assembly = model.rootAssembly
    assembly.DatumCsysByDefault(CARTESIAN)

    # 创建基体和纤维的实例
    inst_matrix = assembly.Instance(name='Matrix-1', part=p_matrix, dependent=OFF)
    instances_to_merge = [inst_matrix]

    if fiberCount > 0:
        inst_fiber = assembly.Instance(name='Fiber-1', part=p_fiber, dependent=OFF)
        instances_to_merge.append(inst_fiber)

    # 第一步: 合并基体和纤维
    if len(instances_to_merge) > 1:
        print("  Step 8a: Merging Matrix and Fibers...")
        assembly.InstanceFromBooleanMerge(
            name='RVE-Temp-WithFibers',
            instances=tuple(instances_to_merge),
            keepIntersections=ON,
            originalInstances=DELETE
        )
        print("    Merge completed.")
        part_to_be_cut_instance = assembly.instances['RVE-Temp-WithFibers-1']
        part_to_be_cut_name = 'RVE-Temp-WithFibers'
    else:
        part_to_be_cut_instance = inst_matrix
        part_to_be_cut_name = 'Matrix'
        assembly.features.changeKey(fromName='Matrix-1', toName='RVE-Temp-WithFibers-1')

    # 第二步: 创建空隙实例(含3D周期性)并执行Cut
    if enable_void:  # 此时 enable_void 必须为 True 且 void_list 必须有内容
        print("  Step 8b: Creating void instances (with 3D periodicity) for cutting...")

        # 验证空隙尺寸是否合理
        min_void_size = void_radius
        if min_void_size < 0.0005:  # 小于0.5微米
            print("    WARNING: Void dimensions are very small (min: %.6f mm)" % min_void_size)
            print("    This may cause geometry creation issues in Abaqus.")
            print("    Recommendation: Use void dimensions >= 0.001 mm")

        cutting_instances = []
        temp_void_part_names = []  # 存储名称以供以后清理
        inst_counter = 0

        # 空隙的最大尺寸，用于边界判断
        void_threshold = void_r_max_overall * void_boundary_threshold_factor

        for i, void_data in enumerate(void_list):
            cx, cy, cz, _, _ = void_data
            center = (cx, cy, cz)

            # (智能策略) 判断该空隙需要在哪些镜像位置创建副本
            mirror_offsets = determineMirrorOffsets_3D(center, rveSize, void_threshold)

            # 创建原始空隙实例 (中心RVE内) 和所有必要的镜像实例
            positions_to_create = [(0, 0, 0)]
            positions_to_create.extend(mirror_offsets)

            for offset in positions_to_create:
                inst_counter += 1
                inst_name = 'Void-Inst-%d' % inst_counter

                # 1. 为每个空隙创建独立的球体Part
                temp_part_name = 'VoidPart-%d' % inst_counter
                temp_void_part_names.append(temp_part_name)
                p_void_sphere = model.Part(name=temp_part_name, dimensionality=THREE_D, type=DEFORMABLE_BODY)

                # 2. 创建球体几何 (通过旋转半圆)
                sketch_size = max(4.0 * void_radius, 0.02)
                s_sphere = model.ConstrainedSketch(name='__sphere_profile_%d__' % inst_counter,
                                                   sheetSize=sketch_size)

                # 绘制旋转轴 (Y轴)
                s_sphere.ConstructionLine(point1=(0.0, -2.0 * void_radius), point2=(0.0, 2.0 * void_radius))
                s_sphere.assignCenterline(line=s_sphere.geometry.findAt((0.0, 0.0)))

                # 绘制半圆弧 (右侧)
                s_sphere.ArcByCenterEnds(center=(0.0, 0.0),
                                         point1=(0.0, void_radius),
                                         point2=(0.0, -void_radius),
                                         direction=CLOCKWISE)
                # 封闭直线
                s_sphere.Line(point1=(0.0, void_radius), point2=(0.0, -void_radius))

                try:
                    p_void_sphere.BaseSolidRevolve(sketch=s_sphere, angle=360.0, flipRevolveDirection=OFF)
                    del model.sketches['__sphere_profile_%d__' % inst_counter]
                except Exception as e:
                    print("      ERROR: Failed to create void sphere #%d. Error: %s" % (inst_counter, str(e)))
                    try:
                        del model.sketches['__sphere_profile_%d__' % inst_counter]
                    except:
                        pass
                    try:
                        del model.parts[temp_part_name]
                    except:
                        pass
                    continue

                # 3. 创建实例
                inst_void = assembly.Instance(name=inst_name, part=p_void_sphere, dependent=OFF)

                # 4. 平移到最终位置（原始中心 + 周期性偏移）
                # (球体是旋转不变的, 不需要旋转步骤)
                final_pos = (cx + offset[0], cy + offset[1], cz + offset[2])
                assembly.translate(instanceList=(inst_name,), vector=final_pos)

                cutting_instances.append(inst_void)

        print("    Created %d total void instances (including 3D mirrors)." % len(cutting_instances))

        if cutting_instances:
            assembly.InstanceFromBooleanCut(
                name='RVE-3D',
                instanceToBeCut=part_to_be_cut_instance,
                cuttingInstances=tuple(cutting_instances),
                originalInstances=DELETE
            )
            print("    Cut completed.")
        else:
            print("    No void instances to cut. Renaming part.")
            assembly.features.changeKey(fromName='RVE-Temp-WithFibers-1', toName='RVE-3D-1')
            model.parts.changeKey(fromName=part_to_be_cut_name, toName='RVE-3D')

    else:
        print("  Step 8b: No voids to cut. Renaming final part.")
        assembly.features.changeKey(fromName='RVE-Temp-WithFibers-1', toName='RVE-3D-1')
        model.parts.changeKey(fromName=part_to_be_cut_name, toName='RVE-3D')

    print("Step 8 Complete.")

    # ==================== 步骤 9: 裁剪 ====================
    print("\nStep 9: Trimming to RVE boundaries...")
    p_rve = model.parts['RVE-3D']
    faces = p_rve.faces
    edges = p_rve.edges

    top_face = None
    try:
        zMax = rveSize[2]
        tol = 1e-6
        top_face_candidates = p_rve.faces.getByBoundingBox(
            xMin=-tol, xMax=rveSize[0] + tol, yMin=-tol, yMax=rveSize[1] + tol,
            zMin=zMax - tol, zMax=zMax + tol)
        top_face = max(top_face_candidates, key=lambda f: f.getSize())
    except:
        for face in faces:
            try:
                if abs(face.pointOn[0][2] - depth) < 1e-6:
                    top_face = face
                    break
            except:
                continue
    if top_face is None: top_face = faces[0]

    top_edge = None
    tol = 1e-6
    try:
        edge_candidates = p_rve.edges.getByBoundingBox(
            xMin=-tol, xMax=tol, yMin=-tol, yMax=rveSize[1] + tol, zMin=-tol, zMax=tol)
        if edge_candidates:
            top_edge = max(edge_candidates, key=lambda e: e.getSize())
    except:
        pass
    if top_edge is None: top_edge = edges[0]

    t = p_rve.MakeSketchTransform(sketchPlane=top_face,
                                  sketchUpEdge=top_edge,
                                  sketchPlaneSide=SIDE1, sketchOrientation=RIGHT,
                                  origin=(rveSize[0] / 2, rveSize[1] / 2, depth))
    s_cut = model.ConstrainedSketch(name='CutSketch', sheetSize=max(rveSize) * 3, transform=t)
    p_rve.projectReferencesOntoSketch(sketch=s_cut, filter=COPLANAR_EDGES)
    margin = 2 * (fiberRadius + (void_r_max_overall if enable_void else 0))
    s_cut.rectangle(point1=(-rveSize[0] / 2, -rveSize[1] / 2),
                    point2=(rveSize[0] / 2, rveSize[1] / 2))
    s_cut.rectangle(point1=(-rveSize[0] / 2 - margin, -rveSize[1] / 2 - margin),
                    point2=(rveSize[0] / 2 + margin, rveSize[1] / 2 + margin))
    p_rve.CutExtrude(sketchPlane=top_face, sketchUpEdge=top_edge,
                     sketchPlaneSide=SIDE1, sketchOrientation=RIGHT,
                     sketch=s_cut, flipExtrudeDirection=OFF)
    print("Trimming completed.")

    # ==================== 步骤 10: 创建几何集合 ====================
    print("\nStep 10: Creating geometry sets...")
    p_rve = model.parts['RVE-3D']
    all_cells = p_rve.cells
    p_rve.Set(cells=all_cells, name='set_AllCells')

    if len(all_cells) > 1:
        # 调用分类函数
        fiber_cells_list, matrix_cells_list = classifyCellsImproved(
            all_cells, rveSize, rveVolume,
            fiber_centers, fiberRadius
        )
        if matrix_cells_list:
            p_rve.Set(name='set_MatrixCell', cells=CellArray(matrix_cells_list))
        else:
            p_rve.Set(name='set_MatrixCell', cells=p_rve.cells[0:0])

        if fiber_cells_list:
            p_rve.Set(name='set_FiberCells', cells=CellArray(fiber_cells_list))
        else:
            p_rve.Set(name='set_FiberCells', cells=p_rve.cells[0:0])

    else:
        p_rve.Set(name='set_MatrixCell', cells=all_cells)
        p_rve.Set(name='set_FiberCells', cells=p_rve.cells[0:0])

    print("  Fiber cells set: %d" % len(p_rve.sets['set_FiberCells'].cells))
    print("  Matrix cells set: %d" % len(p_rve.sets['set_MatrixCell'].cells))

    # 创建纤维-基体界面集合 (为Cohesive做准备)
    print("  Finding Fiber-Matrix interface...")
    matrix_faces_list = []
    if 'set_MatrixCell' in p_rve.sets and len(p_rve.sets['set_MatrixCell'].cells) > 0:
        all_part_faces = p_rve.faces
        for cell in p_rve.sets['set_MatrixCell'].cells:
            face_indices = cell.getFaces()
            for index in face_indices:
                matrix_faces_list.append(all_part_faces[index])
    if matrix_faces_list:
        p_rve.Set(name='set_AllMatrixFaces', faces=FaceArray(matrix_faces_list))
    else:
        p_rve.Set(name='set_AllMatrixFaces', faces=p_rve.faces[0:0])

    fiber_faces_list = []
    if 'set_FiberCells' in p_rve.sets and len(p_rve.sets['set_FiberCells'].cells) > 0:
        all_part_faces = p_rve.faces
        for cell in p_rve.sets['set_FiberCells'].cells:
            face_indices = cell.getFaces()
            for index in face_indices:
                fiber_faces_list.append(all_part_faces[index])
    if fiber_faces_list:
        p_rve.Set(name='set_AllFiberFaces', faces=FaceArray(fiber_faces_list))
    else:
        p_rve.Set(name='set_AllFiberFaces', faces=p_rve.faces[0:0])

    if ('set_AllFiberFaces' in p_rve.sets and 'set_AllMatrixFaces' in p_rve.sets and
            len(p_rve.sets['set_AllFiberFaces'].faces) > 0 and
            len(p_rve.sets['set_AllMatrixFaces'].faces) > 0):
        try:
            p_rve.SetByBoolean(name='set_FiberMatrixInterface',
                               sets=(p_rve.sets['set_AllFiberFaces'], p_rve.sets['set_AllMatrixFaces']),
                               operation=INTERSECTION)
        except Exception as e:
            print("  ERROR: Interface intersection failed: %s" % str(e))
            p_rve.Set(name='set_FiberMatrixInterface', faces=p_rve.faces[0:0])
    else:
        print("  Matrix or Fiber face sets are empty, no interface to find.")
        p_rve.Set(name='set_FiberMatrixInterface', faces=p_rve.faces[0:0])

    print("  Fiber-Matrix interface faces (for cohesive): %d" % len(p_rve.sets['set_FiberMatrixInterface'].faces))
    print("Step 10 Complete.")

    # ==================== 步骤 11: 材料和截面 ====================
    print("\nStep 11: Defining materials and sections...")

    # 基体材料 (Drucker-Prager)
    mat_matrix = model.Material(name='Matrix-Material')
    mat_matrix.Elastic(table=((matrix_E, matrix_nu),))
    mat_matrix.DruckerPrager(table=((matrix_friction_angle, matrix_flow_stress_ratio, matrix_dilation_angle),))
    mat_matrix.druckerPrager.DruckerPragerHardening(table=((matrix_hardening_yield, matrix_hardening_plastic_strain),))
    mat_matrix.DuctileDamageInitiation(
        table=((matrix_damage_strain, matrix_damage_stress_triax, matrix_damage_strain_rate),))
    mat_matrix.ductileDamageInitiation.DamageEvolution(type=DISPLACEMENT, table=((matrix_damage_displacement,),))

    # 纤维材料 (正交各向异性, MPa)
    mat_fiber = model.Material(name='Fiber-Material')
    mat_fiber.Elastic(type=ENGINEERING_CONSTANTS,
                      table=((fiber_E1 * 1000, fiber_E2 * 1000, fiber_E2 * 1000,
                              fiber_nu12, fiber_nu12, 0.25,  # 假设 nu_23 = 0.25
                              fiber_G12 * 1000, fiber_G12 * 1000, fiber_G23 * 1000),))

    # 界面材料 (Cohesive)
    mat_cohesive = model.Material(name='Cohesive-Material')
    mat_cohesive.Elastic(type=TRACTION, table=((cohesive_K_nn, cohesive_K_ss, cohesive_K_tt),))
    mat_cohesive.QuadsDamageInitiation(table=((cohesive_t_n, cohesive_t_s, cohesive_t_t),))
    mat_cohesive.quadsDamageInitiation.DamageEvolution(type=ENERGY, mixedModeBehavior=BK, power=cohesive_eta,
                                                       table=((cohesive_GIC, cohesive_GIIC, cohesive_GIIIC),))
    # 粘性稳定系数
    mat_cohesive.quadsDamageInitiation.DamageStabilizationCohesive(cohesiveCoeff=cohesive_stab_coeff)

    # 截面定义
    model.HomogeneousSolidSection(name='FiberSection', material='Fiber-Material')
    model.HomogeneousSolidSection(name='MatrixSection', material='Matrix-Material')
    model.CohesiveSection(name='CohesiveSection', material='Cohesive-Material',
                          response=TRACTION_SEPARATION)

    # 截面指派
    if 'set_FiberCells' in p_rve.sets and p_rve.sets['set_FiberCells'].cells:
        p_rve.SectionAssignment(region=p_rve.sets['set_FiberCells'], sectionName='FiberSection')
    if 'set_MatrixCell' in p_rve.sets and p_rve.sets['set_MatrixCell'].cells:
        p_rve.SectionAssignment(region=p_rve.sets['set_MatrixCell'], sectionName='MatrixSection')

    # 纤维材料方向指派 (Direction-1沿Z轴)
    if 'set_FiberCells' in p_rve.sets and len(p_rve.sets['set_FiberCells'].cells) > 0:
        print("  Assigning fiber material orientation (Dir-1 along Z-axis)...")
        try:
            p_rve.MaterialOrientation(
                region=p_rve.sets['set_FiberCells'],
                orientationType=SYSTEM,
                axis=AXIS_3,
                localCsys=None,
                additionalRotationType=ROTATION_NONE,
                stackDirection=STACK_1
            )
        except Exception as e:
            print("  ERROR: Failed to assign fiber orientation: %s" % str(e))

    print("Step 11 Complete.")

    # ==================== 步骤 12: 网格 ====================
    print("\nStep 12: Meshing RVE...")
    p_rve = model.parts['RVE-3D']
    p_rve.seedPart(size=globalSeedSize, deviationFactor=deviationFactor,
                   minSizeFactor=minSizeFactor, constraint=FREE)

    # 定义体单元类型
    elemType_bulk_tet = ElemType(elemCode=C3D10, elemLibrary=STANDARD,
                                 elemDeletion=ON, maxDegradation=0.99)
    elemType_bulk_wedge = ElemType(elemCode=C3D15, elemLibrary=STANDARD,
                                   elemDeletion=ON, maxDegradation=0.99)
    elemType_bulk_hex = ElemType(elemCode=C3D8R, elemLibrary=STANDARD,
                                 elemDeletion=ON, maxDegradation=0.99)

    # 将体单元类型仅赋予基体和纤维
    bulk_cells_list = []
    if 'set_MatrixCell' in p_rve.sets:
        bulk_cells_list.extend(p_rve.sets['set_MatrixCell'].cells)
    if 'set_FiberCells' in p_rve.sets:
        bulk_cells_list.extend(p_rve.sets['set_FiberCells'].cells)

    if bulk_cells_list:
        if enable_void:
            print("  Using TET-dominated mesh (suitable for void-containing geometry)...")
            p_rve.setMeshControls(regions=CellArray(bulk_cells_list),
                                  elemShape=TET, technique=FREE)
            # 必须同时允许TET和WEDGE，因为TET网格划分器会生成WEDGE作为过渡
            p_rve.setElementType(regions=(CellArray(bulk_cells_list),),
                                 elemTypes=(elemType_bulk_tet, elemType_bulk_wedge))
        else:
            print("  Using HEX-dominated mesh with SWEEP technique...")
            p_rve.setMeshControls(regions=CellArray(bulk_cells_list),
                                  elemShape=HEX_DOMINATED, technique=SWEEP)
            # 允许六面体和楔形
            p_rve.setElementType(regions=(CellArray(bulk_cells_list),),
                                 elemTypes=(elemType_bulk_hex, elemType_bulk_wedge))
    else:
        print("  WARNING: No bulk cells (Matrix/Fiber) found to assign element type.")

    print("  Generating mesh for all regions...")
    p_rve.generateMesh()
    print("  Total elements (pre-insertion): %d" % len(p_rve.elements))
    print("Step 12 Complete.")

    # ==================== 步骤 13: 粘接单元处理 ====================
    print("\nStep 13: Inserting cohesive elements...")

    if 'set_FiberMatrixInterface' in p_rve.sets and len(p_rve.sets['set_FiberMatrixInterface'].faces) > 0:
        try:
            p_rve.insertElements(faces=p_rve.sets['set_FiberMatrixInterface'])
            print("  Cohesive elements inserted at Fiber-Matrix interface.")
        except Exception as e:
            print("  ERROR: Cohesive element insertion failed: %s" % str(e))
            print("  This often happens with complex geometry or poor mesh quality.")
    else:
        print("  No Fiber-Matrix interface found. Skipping cohesive insertion.")

    # 创建单元集合
    all_elements = p_rve.elements
    p_rve.Set(elements=all_elements, name='set_AllElements')

    # 重新获取几何集合 (插入Cohesive后, 单元归属可能变化)
    if 'set_FiberCells' in p_rve.sets and len(p_rve.sets['set_FiberCells'].elements) > 0:
        p_rve.Set(elements=p_rve.sets['set_FiberCells'].elements, name='set_FiberElements')
    else:
        p_rve.Set(elements=p_rve.elements[0:0], name='set_FiberElements')

    if 'set_MatrixCell' in p_rve.sets and len(p_rve.sets['set_MatrixCell'].elements) > 0:
        p_rve.Set(elements=p_rve.sets['set_MatrixCell'].elements, name='set_MatrixElements')
    else:
        p_rve.Set(elements=p_rve.elements[0:0], name='set_MatrixElements')

    # 识别Cohesive单元
    p_rve.SetByBoolean(name='set_CohesiveElements',
                       sets=(p_rve.sets['set_AllElements'],
                             p_rve.sets['set_FiberElements'],
                             p_rve.sets['set_MatrixElements']),
                       operation=DIFFERENCE)

    print("\n  Element Classification (Post-Insertion):")
    print("    Total elements: %d" % len(all_elements))
    print("    Fiber elements: %d" % len(p_rve.sets['set_FiberElements'].elements))
    print("    Matrix elements: %d" % len(p_rve.sets['set_MatrixElements'].elements))
    print("    Cohesive elements: %d" % len(p_rve.sets['set_CohesiveElements'].elements))

    # 指派Cohesive单元类型和截面
    if len(p_rve.sets['set_CohesiveElements'].elements) > 0:
        # 定义 Cohesive 单元类型 (8节点 和 6节点)
        elemType_coh_hex = ElemType(elemCode=COH3D8, elemLibrary=STANDARD,
                                    elemDeletion=ON, maxDegradation=0.99)
        elemType_coh_wedge = ElemType(elemCode=COH3D6, elemLibrary=STANDARD,
                                      elemDeletion=ON, maxDegradation=0.99)

        # 同时分配 COH3D8 和 COH3D6 ***
        # 因为TET网格(C3D4)的面是三角形, 生成COH3D6
        # SWEEP网格(C3D8/C3D6)的面是四边形或三角形, 生成COH3D8或COH3D6
        p_rve.setElementType(regions=(p_rve.sets['set_CohesiveElements'].elements,),
                             elemTypes=(elemType_coh_hex, elemType_coh_wedge))

        p_rve.SectionAssignment(region=p_rve.sets['set_CohesiveElements'],
                                sectionName='CohesiveSection')
        print("  Cohesive elements configured (COH3D8 and COH3D6).")

    print("Step 13 Complete.")

    # ==================== 步骤 14: 创建分析作业并施加周期性边界条件 (MP) ====================
    print("\nStep 14: Creating Analysis Job and Applying PBCs (Micromechanics Plugin)...")

    # 关键步骤：删除旧的几何实例，并从已划分网格的Part创建新实例
    rootAssembly = model.rootAssembly
    if 'RVE-3D-1' in rootAssembly.instances:
        del rootAssembly.instances['RVE-3D-1']
        print("  Deleted old geometric instance 'RVE-3D-1'.")

    # 确保我们获取的是最新的、已网格化的part
    p_rve = model.parts['RVE-3D']
    rveInstance = rootAssembly.Instance(name='RVE-3D-1', part=p_rve, dependent=ON)
    print("  Created new mesh instance 'RVE-3D-1' in assembly.")

    # 导入Micromechanics Plugin模块
    try:
        if abaqus_plugin_path not in sys.path:
            sys.path.insert(0, abaqus_plugin_path)
        import microMechanics
        from microMechanics.mmpBackend import Interface
        from microMechanics.mmpBackend.mmpInterface.mmpRVEConstants import *
    except ImportError as e:
        print("\n  FATAL ERROR: Failed to import Micromechanics Plugin.")
        print("  Please ensure the plugin is installed and the path is correct:")
        print("  PATH: %s" % abaqus_plugin_path)
        print("  Error details: %s" % str(e))
        raise e

    try:
        # 1. 插件要求Job必须预先存在
        print("  Creating analysis job '%s' with %d CPUs..." % (job_name, num_cpus))
        if job_name in mdb.jobs:
            del mdb.jobs[job_name]

        mdb.Job(name=job_name, model=modelName, description='RVE Homogenization Analysis',
                type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
                memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
                explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
                modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
                scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=num_cpus,
                numDomains=num_cpus, numGPUs=0)
        print("  Job '%s' created successfully." % job_name)

        # 2. 调用插件函数施加PBC
        print("  Calling Micromechanics Plugin to apply PBC...")
        # --- 定义载荷历史 ---
        # 这是一个12个元素的元组: (E11, E22, E33, E12, E13, E23, E11_Amp, E22_Amp, E33_Amp, E12_Amp, E13_Amp, E23_Amp)
        field_history_tuple = (
            '',  # E11_val [0]
            '',  # E22_val [1]
            '',  # E33_val [2]
            '0.02',  # E12_val [3]
            '',  # E13_val [4]
            '',  # E23_val [5]
            '',  # E11_amp [6]
            '',  # E22_amp [7]
            '',  # E33_amp [8]
            'Default',  # E12_amp [9]
            '',  # E13_amp [10]
            ''  # E23_amp [11]
            # 只有12个元素, 'fieldHistory1'。如果需要温度, 将使用 'fieldHistory2'
        )
        # --- 载荷历史定义结束 ---

        Interface.Loading.MechanicalModelMaker(
            modelName=modelName,
            jobName=job_name,
            # constraintType='PERIODIC',
            constraintType='TAYLOR',
            drivenField='STRAIN',
            mechanicalHistoryType=LOADUSERDEFINED,  # 对应 "User-Defined"
            fieldHistory1=field_history_tuple,  # 对应 E12=0.02, Amp=Default
            totalHistoryTime=1.0,  # 对应 Total Time: 1
            homogenizationIntervals=0,  # 对应 Frequency: Initial
            homogenizeProperties=(True, False, False),  # (Elastic, CTE, Density)
            doNotSubmit=True  # 对应 "Do not Submit Job"
        )

        print("  Micromechanics Plugin executed successfully.")
        print("  PBC, steps, and user-defined load history configured in job '%s'." % job_name)

    except Exception as e:
        print("\n  FATAL ERROR: Failed to apply PBC using Micromechanics Plugin.")
        print("  Error details: %s" % str(e))
        raise e

    print("\nStep 14 Complete.")

    # ==================== 步骤 15: 清理多余Part ====================
    print("\nStep 15: Cleaning up temporary parts...")
    parts_to_delete = [p for p in ['Fiber', 'Matrix', 'RVE-Temp-WithFibers'] if p in model.parts]

    # 将临时的空部件添加到删除列表中
    if enable_void and 'temp_void_part_names' in locals():
        for part_name in temp_void_part_names:
            if part_name in model.parts:
                parts_to_delete.append(part_name)

    for part_name in parts_to_delete:
        try:
            del model.parts[part_name]
            print("  Deleted part: %s" % part_name)
        except:
            print("  Could not delete part: %s" % part_name)
    print("Step 15 Complete.")

    # ==================== 完成 ====================
    print("\n" + "=" * 70)
    print("3D RVE MODEL GENERATION COMPLETED SUCCESSFULLY")
    print("=" * 70)
    print("\nModel Summary:")
    print("  Model Name: %s" % modelName)
    print("  RVE Dimensions: %.6f x %.6f x %.6f" % (rveSize[0], rveSize[1], depth))
    print("  Fiber Count: %d (Target Vf: %.4f)" % (fiberCount, target_Vf))
    if enable_void and void_info:
        print("  Void Count: %d" % void_info['count'])
        print("  Void Porosity: %.4f%% (Actual)" % (void_info['actual_porosity'] * 100))
    print("  Cohesive Elements: %d" % len(p_rve.sets['set_CohesiveElements'].elements))
    print("  Periodic BCs: Applied using the ABAQUS Micromechanics Plugin")
    print("  Analysis Job: '%s' created." % job_name)


# =================================================================
#                 主执行入口
# =================================================================

if __name__ == '__main__':

    # ========== 全局配置参数 ==========
    # RVE尺寸[宽度X, 高度Y, 深度Z](单位:mm)
    RVE_SIZE = [0.057, 0.057, 0.01]
    # 纤维半径(单位:mm)
    FIBER_RADIUS = 0.0035
    # 目标纤维体积分数(0-1)
    TARGET_VF = 0.5
    # RSA播种比例 (0.0 ~ 1.0)，高值(如0.9)排布更均匀，低值(如0.1)更接近物理堆积, 速度较慢
    RSA_SEEDING_RATIO = 0.9

    # ========== 间距参数 ==========
    # 纤维-纤维 最小中心距因子 (min_dist = 因子 * fiberRadius)
    FIBER_MIN_DIST_FACTOR = 2.05

    # 纤维-空隙 最小间距因子 (min_dist = 因子 * (fiberRadius + void_radius))
    # 这是一个2D投影的中心-中心距离
    FIBER_VOID_MIN_DIST_FACTOR = 1.05

    # 空隙-空隙 最小中心距因子 (min_dist = 因子 * void_radius)
    # 这是一个3D中心-中心距离
    MIN_VOID_DIST_FACTOR = 2.1

    # ========== 空隙缺陷参数 ==========
    ENABLE_VOID = True  # 是否启用空隙缺陷
    # 目标空隙率(0-1)
    VOID_POROSITY = 0.01

    # 单个球体空隙的半径 (单位:mm)
    # 脚本将自动计算所需数量以达到 VOID_POROSITY
    VOID_SPHERE_RADIUS = 0.0015

    # 空隙3D周期性边界阈值因子 (threshold = 因子 * Radius)
    VOID_BOUNDARY_THRESHOLD_FACTOR = 1.0

    EXPORT_VOID_GEOMETRY = True  # 是否导出空隙几何信息到CSV

    # ========== 网格参数 ==========
    GLOBAL_SEED_SIZE = 0.0008  # 全局网格尺寸(单位:mm)
    DEVIATION_FACTOR = 0.5  # 网格偏离因子: 允许网格偏离真实几何的程度
    # 最小网格尺寸因子: (最小单元 / 全局尺寸)。 必须足够小以捕捉到最薄的基体/空隙间隙。
    MIN_SIZE_FACTOR = 0.5

    # ========== 纤维材料参数 (GPa) ==========
    FIBER_E1 = 230.0  # 纤维轴向模量(单位:GPa)
    FIBER_E2 = 15.0  # 纤维横向模量(单位:GPa)
    FIBER_G12 = 15.0  # 纤维纵向剪切模量(单位:GPa)
    FIBER_G23 = 7.0  # 纤维横向剪切模量(单位:GPa)
    FIBER_NU12 = 0.2  # 纤维主泊松比

    # ========== 基体材料参数 (MPa) ==========
    MATRIX_E = 3170.0  # 基体弹性模量(单位:MPa)
    MATRIX_NU = 0.35  # 基体泊松比
    MATRIX_FRICTION_ANGLE = 16.0  # Drucker-Prager摩擦角(单位:度)
    MATRIX_FLOW_STRESS_RATIO = 1.0  # Drucker-Prager流动应力比
    MATRIX_DILATION_ANGLE = 16.0  # Drucker-Prager膨胀角(单位:度)
    MATRIX_HARDENING_YIELD = 106.4  # 硬化屈服应力(单位:MPa)
    MATRIX_HARDENING_PLASTIC_STRAIN = 0.0  # 硬化对应塑性应变
    MATRIX_DAMAGE_STRAIN = 0.01  # 韧性损伤起始应变
    MATRIX_DAMAGE_STRESS_TRIAX = 0.0  # 应力三轴度
    MATRIX_DAMAGE_STRAIN_RATE = 0.0  # 应变率
    MATRIX_DAMAGE_DISPLACEMENT = 5e-05  # 损伤演化位移(单位:mm)

    # ========== 界面材料参数 ==========
    COHESIVE_K_NN = 1e8  # 界面法向刚度(单位:N/mm³)
    COHESIVE_K_SS = 1e8  # 界面第一切向刚度(单位:N/mm³)
    COHESIVE_K_TT = 1e8  # 界面第二切向刚度(单位:N/mm³)
    COHESIVE_T_N = 44.0  # 界面法向强度(单位:MPa)
    COHESIVE_T_S = 82.0  # 界面第一切向强度(单位:MPa)
    COHESIVE_T_T = 82.0  # 界面第一切向强度(单位:MPa)
    COHESIVE_GIC = 0.001  # I型断裂能(单位:N/mm)
    COHESIVE_GIIC = 0.002  # II型断裂能(单位:N/mm)
    COHESIVE_GIIIC = 0.002  # III型断裂能(单位:N/mm)
    COHESIVE_ETA = 1.5  # BK准则指数
    COHESIVE_STAB_COEFF = 0.0001  # 粘性稳定系数

    # ========== 输出控制 ==========
    EXPORT_COORDINATES = True  # 是否导出纤维中心坐标到CSV
    CSV_FILENAME = None  # CSV文件名(None则自动生成)

    # ========== 模型命名 ==========
    if ENABLE_VOID and VOID_POROSITY > 0:
        targetModelName = 'RVE_3D_Vf%d_Void%dpct' % (int(TARGET_VF * 100), int(VOID_POROSITY * 100))
    else:
        targetModelName = 'RVE_3D_Vf%d' % int(TARGET_VF * 100)
    tempModelName = targetModelName + '_TEMP_' + str(int(time.time()))

    # ========== 分析作业与插件参数 ==========
    JOB_NAME = 'Job-' + targetModelName  # 分析作业的名称 (与模型名关联)
    NUM_CPUS = 8  # 分析时使用的CPU核心数
    # ABAQUS插件路径 (请根据您的安装情况修改)
    ABAQUS_PLUGIN_PATH = 'd:/Abaqus2023/Plugins'

    print("\n" + "=" * 70)
    print("Starting 3D RVE Generation (Multi-Void / MP-PBC Version)...")
    print("=" * 70)
    print("Target Model Name: %s" % targetModelName)
    print("Target Job Name: %s" % JOB_NAME)
    print("Void enabled: %s" % ("Yes" if (ENABLE_VOID and VOID_POROSITY > 0) else "No"))
    if ENABLE_VOID and VOID_POROSITY > 0:
        print("Target Void porosity: %.2f%%" % (VOID_POROSITY * 100))
    print("Fiber Vf: %.2f%%" % (TARGET_VF * 100))
    print("Using temporary name: %s" % tempModelName)
    print("=" * 70)

    # ========== 执行建模 ==========
    create3DRVEModelWithVoids(
        modelName=tempModelName,
        rveSize=RVE_SIZE,
        fiberRadius=FIBER_RADIUS,
        target_Vf=TARGET_VF,
        fiber_min_dist_factor=FIBER_MIN_DIST_FACTOR,
        globalSeedSize=GLOBAL_SEED_SIZE,
        deviationFactor=DEVIATION_FACTOR,
        minSizeFactor=MIN_SIZE_FACTOR,
        rsa_seeding_ratio=RSA_SEEDING_RATIO,
        export_coordinates=EXPORT_COORDINATES,
        csv_filename=CSV_FILENAME,
        # 空隙参数
        enable_void=ENABLE_VOID,
        void_porosity=VOID_POROSITY,
        void_sphere_radius=VOID_SPHERE_RADIUS,
        min_void_dist_factor=MIN_VOID_DIST_FACTOR,
        fiber_void_min_dist_factor=FIBER_VOID_MIN_DIST_FACTOR,
        void_boundary_threshold_factor=VOID_BOUNDARY_THRESHOLD_FACTOR,
        export_void_geometry=EXPORT_VOID_GEOMETRY,
        # 纤维材料参数
        fiber_E1=FIBER_E1, fiber_E2=FIBER_E2, fiber_G12=FIBER_G12,
        fiber_G23=FIBER_G23, fiber_nu12=FIBER_NU12,
        # 基体材料参数
        matrix_E=MATRIX_E, matrix_nu=MATRIX_NU,
        matrix_friction_angle=MATRIX_FRICTION_ANGLE,
        matrix_flow_stress_ratio=MATRIX_FLOW_STRESS_RATIO,
        matrix_dilation_angle=MATRIX_DILATION_ANGLE,
        matrix_hardening_yield=MATRIX_HARDENING_YIELD,
        matrix_hardening_plastic_strain=MATRIX_HARDENING_PLASTIC_STRAIN,
        matrix_damage_strain=MATRIX_DAMAGE_STRAIN,
        matrix_damage_stress_triax=MATRIX_DAMAGE_STRESS_TRIAX,
        matrix_damage_strain_rate=MATRIX_DAMAGE_STRAIN_RATE,
        matrix_damage_displacement=MATRIX_DAMAGE_DISPLACEMENT,
        # 界面材料参数
        cohesive_K_nn=COHESIVE_K_NN, cohesive_K_ss=COHESIVE_K_SS, cohesive_K_tt=COHESIVE_K_TT,
        cohesive_t_n=COHESIVE_T_N, cohesive_t_s=COHESIVE_T_S, cohesive_t_t=COHESIVE_T_T,
        cohesive_GIC=COHESIVE_GIC, cohesive_GIIC=COHESIVE_GIIC, cohesive_GIIIC=COHESIVE_GIIIC,
        cohesive_eta=COHESIVE_ETA, cohesive_stab_coeff=COHESIVE_STAB_COEFF,
        # 作业和插件参数
        job_name=JOB_NAME,
        num_cpus=NUM_CPUS,
        abaqus_plugin_path=ABAQUS_PLUGIN_PATH
    )

    # ========== 后处理 ==========
    if tempModelName in mdb.models:
        print("\n" + "=" * 70)
        print("Starting Model Cleanup and Rename...")
        print("=" * 70)

        # 清理所有其他模型, 包括可能存在的同名旧模型
        models_to_delete = [m for m in mdb.models.keys()
                            if m != tempModelName]
        if targetModelName in mdb.models:
            models_to_delete.append(targetModelName)

        if models_to_delete:
            print("\nDeleting old/existing models:")
            for m in list(set(models_to_delete)):  # 去重
                if m in mdb.models:
                    print("  - %s" % m)
                    del mdb.models[m]
            print("Old models deleted.")
        else:
            print("\nNo old models to delete.")

        print("\nRenaming model:")
        print("  From: '%s'" % tempModelName)
        print("  To:   '%s'" % targetModelName)
        mdb.models.changeKey(fromName=tempModelName,
                             toName=targetModelName)

        # 更新Job使其指向重命名后的模型
        if JOB_NAME in mdb.jobs:
            try:
                mdb.jobs[JOB_NAME].setValues(model=targetModelName)
                print("\nUpdating Job-Model association:")
                print("  Job '%s' is now linked to model '%s'." % (JOB_NAME, targetModelName))
            except Exception as e:
                print("\nWARNING: Failed to update job association for '%s'." % JOB_NAME)
                print("  Error: %s" % str(e))

        print("\n" + "=" * 70)
        print("MODEL GENERATION COMPLETE")
        print("=" * 70)
    else:
        print("\n" + "!" * 70)
        print("ERROR: Model generation failed. Temporary model not found.")
        print("!" * 70 + "\n")

