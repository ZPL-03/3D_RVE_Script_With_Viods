# -*- coding: utf-8 -*-
# ##################################################################
#
#   可控制纤维体积分数和多空隙缺陷(椭球体)的Abaqus三维RVE建模脚本
#
# 功能说明:
# 1. 根据设定的纤维体积分数(Vf)和空隙率(Vp)自动生成三维RVE几何模型
# 2. 算法流程: [纤维优先，空隙后置]
#    a. [纤维优先] 采用RSA播种 + 锚定松弛算法排布纤维, 保证纤维间(F-F)间距 (2D)
#    b. [空隙后置] 采用3D RSA算法排布旋转椭球体空隙, 播种时同时检查:
#       i.  空隙-空隙 (V-V) 间距 (3D 轴线-轴线 27盒子检测, A轴)
#       ii. 空隙-纤维 (F-V) 间距 (2D 中心-中心 9盒子(XY)检测, 采用保守包围球半径)
# 3. 自动创建材料、赋予截面、划分网格
# 4. 施加三维周期性边界条件(PBCs) - [通过节点配对施加*Equation约束]
# 5. 在纤维-基体界面插入三维cohesive单元(COH3D8/COH3D6)
# 6. 导出纤维中心和空隙信息为CSV文件
#
# 参考文献:
# [1] Li, X., et al. (2025). A novel method for constructing 3D void RVE
#     elements and rapid homogenization of composite materials.
#     Composite Structures, 360, 119040.
#
# 作者: 刘正鹏 (Liu Zhengpeng)
# 版本: v1.0
# 创建日期: 2025-09-30
# 适用软件: ABAQUS 2023
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
#                 3D 空间直线/间距检测
# =================================================================

def shortest_distance_between_segments(p1, q1, p2, q2):
    """计算空间中两条线段 p1-q1 和 p2-q2 之间的最短距离 (解析法)"""
    # 向量表示
    u = (q1[0] - p1[0], q1[1] - p1[1], q1[2] - p1[2])
    v = (q2[0] - p2[0], q2[1] - p2[1], q2[2] - p2[2])
    w = (p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2])

    # 点积
    a = u[0] * u[0] + u[1] * u[1] + u[2] * u[2]  # dot(u, u)
    b = u[0] * v[0] + u[1] * v[1] + u[2] * v[2]  # dot(u, v)
    c = v[0] * v[0] + v[1] * v[1] + v[2] * v[2]  # dot(v, v)
    d = u[0] * w[0] + u[1] * w[1] + u[2] * w[2]  # dot(u, w)
    e = v[0] * w[0] + v[1] * w[1] + v[2] * w[2]  # dot(v, w)

    D = a * c - b * b
    sD, tD = D, D

    if D < 1e-7:  # 平行线
        sN = 0.0
        sD = 1.0
        tN = e
        tD = c
    else:
        sN = (b * e - c * d)
        tN = (a * e - b * d)
        if sN < 0.0:
            sN = 0.0
            tN = e
            tD = c
        elif sN > sD:
            sN = sD
            tN = e + b
            tD = c

    if tN < 0.0:
        tN = 0.0
        if -d < 0.0:
            sN = 0.0
        elif -d > a:
            sN = a
        else:
            sN = -d
        sD = a
    elif tN > tD:
        tN = tD
        if (-d + b) < 0.0:
            sN = 0.0
        elif (-d + b) > a:
            sN = a
        else:
            sN = -d + b
        sD = a

    sc = 0.0 if abs(sN) < 1e-7 else sN / sD
    tc = 0.0 if abs(tN) < 1e-7 else tN / tD

    # 计算距离向量
    dP_x = w[0] + (sc * u[0]) - (tc * v[0])
    dP_y = w[1] + (sc * u[1]) - (tc * v[1])
    dP_z = w[2] + (sc * u[2]) - (tc * v[2])

    return math.sqrt(dP_x ** 2 + dP_y ** 2 + dP_z ** 2)


def _get_ellipsoid_endpoints(cx, cy, cz, rot_z, rot_x, long_semi_axis_a):
    """根据椭球体中心和方向计算其长轴(A轴)端点"""
    L_half = long_semi_axis_a
    cos_rz = math.cos(rot_z)
    sin_rz = math.sin(rot_z)
    cos_rx = math.cos(rot_x)
    sin_rx = math.sin(rot_x)
    # 起点 (local y = -L_half)
    p1_x = L_half * sin_rz
    p1_y = -L_half * cos_rz
    p2_x = p1_x
    p2_y = p1_y * cos_rx
    p2_z = p1_y * sin_rx
    p_void = (cx + p2_x, cy + p2_y, cz + p2_z)
    # 终点 (local y = +L_half)
    q1_x = -L_half * sin_rz
    q1_y = L_half * cos_rz
    q2_x = q1_x
    q2_y = q1_y * cos_rx
    q2_z = q1_y * sin_rx
    q_void = (cx + q2_x, cy + q2_y, cz + q2_z)
    return p_void, q_void


# =================================================================
#                 3D RSA 算法 (用于空隙排布)
# =================================================================

def _generate_void_coordinates_rsa(num_voids, rveSize, void_long_semi_axis_a,
                                   min_void_dist_axis,
                                   fiber_centers, min_fiber_void_dist_2D,
                                   max_attempts=5000):
    """
    使用3D RSA(随机顺序吸附)算法生成空隙
    V-V 检查: 3D 轴线-轴线 27盒子干涉检测 (使用椭球体长轴A)
    F-V 检查: 2D 中心-中心 9盒子(XY)干涉检测 (使用保守包围球半径)

    参数:
        num_voids: 目标空隙数量
        rveSize: RVE尺寸 [L, W, H]
        void_long_semi_axis_a: 椭球体空隙的长半轴 (A)
        min_void_dist_axis: 空隙轴线间的最小允许距离 (3D, V-V)
        fiber_centers: 已固定的纤维中心 [(x, y), ...]
        min_fiber_void_dist_2D: 纤维-空隙最小中心距离 (2D保守, F-V)
        max_attempts: 每次放置尝试的最大次数

    返回:
        void_list: [(cx, cy, cz, rot_z, rot_x), ...]
    """
    print("--- Initializing 3D RSA for void placement (V-V Axis-to-Axis & F-V 2D-Conservative Center) ---")
    print("    Target voids: %d" % num_voids)
    print("    Min V-V Axis-to-Axis distance (3D): %.6f" % min_void_dist_axis)
    print("    Min F-V Center-to-Center distance (2D Conservative): %.6f" % min_fiber_void_dist_2D)

    void_list = []
    L, W, H = rveSize[0], rveSize[1], rveSize[2]

    # 定义 27 (V-V) 周期性偏移
    offsets_3d = []
    for dx in [-L, 0, L]:
        for dy in [-W, 0, W]:
            for dz in [-H, 0, H]:
                offsets_3d.append((dx, dy, dz))

    # 定义 9 (F-V) 周期性偏移
    offsets_2d_xy = []
    for dx in [-L, 0, L]:
        for dy in [-W, 0, W]:
            offsets_2d_xy.append((dx, dy))

    # F-V 检查使用距离的平方
    min_fv_dist_sq = min_fiber_void_dist_2D ** 2

    for i in range(num_voids):
        placed = False
        for attempt in range(max_attempts):
            # 1. 生成新空隙候选
            cx_new = rd.uniform(0, L)
            cy_new = rd.uniform(0, W)
            cz_new = rd.uniform(0, H)
            rot_z_new = rd.uniform(0, 2 * math.pi)
            rot_x_new = rd.uniform(0, 2 * math.pi)

            is_too_close = False

            # 2. V-V 检查 (Axis-to-Axis, 27-box)
            # 仅当存在已放置空隙时才需要计算新空隙的3D端点
            if void_list:
                p_new, q_new = _get_ellipsoid_endpoints(
                    cx_new, cy_new, cz_new, rot_z_new, rot_x_new, void_long_semi_axis_a
                )

                for void_data_existing in void_list:
                    (cx_ex, cy_ex, cz_ex, rot_z_ex, rot_x_ex) = void_data_existing

                    # 获取已存在空隙的轴线端点
                    p_existing, q_existing = _get_ellipsoid_endpoints(
                        cx_ex, cy_ex, cz_ex, rot_z_ex, rot_x_ex, void_long_semi_axis_a
                    )

                    # 检查新空隙(p_new, q_new)与已存在空隙的27个镜像(p_mir, q_mir)的距离
                    for dx, dy, dz in offsets_3d:
                        p_mir = (p_existing[0] + dx, p_existing[1] + dy, p_existing[2] + dz)
                        q_mir = (q_existing[0] + dx, q_existing[1] + dy, q_existing[2] + dz)

                        dist = shortest_distance_between_segments(p_new, q_new, p_mir, q_mir)

                        if dist < min_void_dist_axis:
                            is_too_close = True
                            break

                    if is_too_close:
                        break

            if is_too_close:
                continue  # 尝试新的 attempt

            # 3. F-V 检查 (Conservative 2D Center-to-Center, 9-box XY)
            # 此检查不关心3D朝向, 只关心(cx, cy)
            for fx, fy in fiber_centers:
                # 检查新空隙(cx_new, cy_new)与纤维的9个XY平面镜像(fx_mir, fy_mir)的距离
                for dx, dy in offsets_2d_xy:
                    fx_mir = fx + dx
                    fy_mir = fy + dy

                    dist_sq = (cx_new - fx_mir) ** 2 + (cy_new - fy_mir) ** 2

                    if dist_sq < min_fv_dist_sq:
                        is_too_close = True
                        break  # (offsets_2d_xy loop)

                if is_too_close:
                    break  # (fiber_centers loop)

            if is_too_close:
                continue  # 尝试新的 attempt

            # 4. 放置
            if not is_too_close:
                void_list.append((cx_new, cy_new, cz_new, rot_z_new, rot_x_new))
                placed = True
                break  # (attempt loop)

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

def exportVoidGeometryToCSV(filename, void_list, void_axes, rveVolume, target_porosity):
    """将多个空隙的几何信息导出为CSV文件"""
    try:
        work_dir = os.getcwd()
        filepath = os.path.join(work_dir, filename)

        num_voids = len(void_list)
        # void_axes 是 [A, C]
        # 体积是 4/3 * pi * A * C * C
        vol_per_void = (4.0 / 3.0) * math.pi * void_axes[0] * void_axes[1] * void_axes[1]
        actual_volume = num_voids * vol_per_void
        actual_porosity = actual_volume / rveVolume

        with open(filepath, 'w') as f:
            f.write("# 3D RVE Multi-Void Geometry Information\n")
            f.write("# Generated: %s\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
            f.write("# Void Type: Spheroid (Rotation-Ellipsoid)\n")
            f.write("# Target Porosity: %.4f%%\n" % (target_porosity * 100))
            f.write("# RVE Volume: %.8f\n" % rveVolume)
            f.write("# --- Per Void --- \n")
            f.write("# Semi-axes (A-long, C-short): %.6f, %.6f\n" % (
                void_axes[0], void_axes[1]))
            f.write("# Volume per Void: %.8f\n" % vol_per_void)
            f.write("# --- Total --- \n")
            f.write("# Void Count (Placed): %d\n" % num_voids)
            f.write("# Total Void Volume (Actual): %.8f\n" % actual_volume)
            f.write("# Actual Porosity: %.4f%%\n" % (actual_porosity * 100))
            f.write("#\n")
            f.write("Void_ID,Center_X,Center_Y,Center_Z,Rotation_Z(rad),Rotation_X(rad)\n")
            for i, (cx, cy, cz, rot_z, rot_x) in enumerate(void_list, start=1):
                f.write("%d,%.8f,%.8f,%.8f,%.8f,%.8f\n" % (i, cx, cy, cz, rot_z, rot_x))

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
#                 周期性边界条件 (PBCs) 辅助函数
# =================================================================

def createReferencePoints3D(model):
    """创建三个参考点用于三维周期性边界条件"""
    rootAssembly = model.rootAssembly
    p_X = model.Part(name='Ref-X-Part', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p_X.ReferencePoint(point=(0.0, 0.0, 0.0))
    p_Y = model.Part(name='Ref-Y-Part', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p_Y.ReferencePoint(point=(0.0, 0.0, 0.0))
    p_Z = model.Part(name='Ref-Z-Part', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p_Z.ReferencePoint(point=(0.0, 0.0, 0.0))
    rootAssembly.Instance(name='Ref-X-Instance', part=p_X, dependent=ON)
    rootAssembly.Instance(name='Ref-Y-Instance', part=p_Y, dependent=ON)
    rootAssembly.Instance(name='Ref-Z-Instance', part=p_Z, dependent=ON)
    rootAssembly.Set(name='set_RefPoint_X',
                     referencePoints=(rootAssembly.instances['Ref-X-Instance'].referencePoints[1],))
    rootAssembly.Set(name='set_RefPoint_Y',
                     referencePoints=(rootAssembly.instances['Ref-Y-Instance'].referencePoints[1],))
    rootAssembly.Set(name='set_RefPoint_Z',
                     referencePoints=(rootAssembly.instances['Ref-Z-Instance'].referencePoints[1],))


def getRVEDimensions3D(model, instanceName):
    """获取三维RVE的边界坐标"""
    nodes = model.rootAssembly.instances[instanceName].nodes
    xMin = min(n.coordinates[0] for n in nodes)
    xMax = max(n.coordinates[0] for n in nodes)
    yMin = min(n.coordinates[1] for n in nodes)
    yMax = max(n.coordinates[1] for n in nodes)
    zMin = min(n.coordinates[2] for n in nodes)
    zMax = max(n.coordinates[2] for n in nodes)
    return xMin, xMax, yMin, yMax, zMin, zMax


def getBoundaryNodes3D(model, instanceName, dimensions):
    """获取六个面上的所有节点"""
    nodes = model.rootAssembly.instances[instanceName].nodes
    xMin, xMax, yMin, yMax, zMin, zMax = dimensions
    tol = 1e-6
    nodes_left, nodes_right = [], []
    nodes_front, nodes_back = [], []
    nodes_bottom, nodes_top = [], []
    for n in nodes:
        x, y, z = n.coordinates[0], n.coordinates[1], n.coordinates[2]
        if abs(x - xMin) < tol: nodes_left.append(n)
        if abs(x - xMax) < tol: nodes_right.append(n)
        if abs(y - yMin) < tol: nodes_front.append(n)
        if abs(y - yMax) < tol: nodes_back.append(n)
        if abs(z - zMin) < tol: nodes_bottom.append(n)
        if abs(z - zMax) < tol: nodes_top.append(n)
    nodes_left.sort(key=lambda n: (n.coordinates[1], n.coordinates[2]))
    nodes_right.sort(key=lambda n: (n.coordinates[1], n.coordinates[2]))
    nodes_front.sort(key=lambda n: (n.coordinates[0], n.coordinates[2]))
    nodes_back.sort(key=lambda n: (n.coordinates[0], n.coordinates[2]))
    nodes_bottom.sort(key=lambda n: (n.coordinates[0], n.coordinates[1]))
    nodes_top.sort(key=lambda n: (n.coordinates[0], n.coordinates[1]))
    return nodes_left, nodes_right, nodes_front, nodes_back, nodes_bottom, nodes_top


def pairBoundaryNodes3D(slave_nodes, master_nodes, tolerance, coord_indices):
    """三维节点配对算法"""
    paired_nodes = []
    master_pool = list(master_nodes)
    slave_pool = list(slave_nodes)
    idx1, idx2 = coord_indices
    for s_node in slave_pool:
        s_coord = s_node.coordinates
        best_match_node = None
        min_dist = float('inf')
        for m_node in master_pool:
            m_coord = m_node.coordinates
            delta1 = abs(s_coord[idx1] - m_coord[idx1])
            delta2 = abs(s_coord[idx2] - m_coord[idx2])
            planar_dist = math.sqrt(delta1 ** 2 + delta2 ** 2)
            if planar_dist < tolerance and planar_dist < min_dist:
                min_dist = planar_dist
                best_match_node = m_node
        if best_match_node:
            paired_nodes.append((s_node, best_match_node))
            master_pool.remove(best_match_node)
    return paired_nodes


def applyPeriodicConstraints3D(model, instanceName, node_pairs, pair_type, constrained_nodes_set):
    """
    施加三维周期性约束 (避免重复约束节点自由度)

    Args:
        model: Abaqus 模型对象.
        instanceName: RVE 实例的名称.
        node_pairs: 包含 (从节点_node, 主节点_node) 元组的列表.
        pair_type: 字符串 ('Left-Right', 'Front-Back', 'Bottom-Top').
        constrained_nodes_set: 一个包含已被约束节点标签的集合 (set).
                               此函数将更新这个集合.
    """
    r_assy = model.rootAssembly
    inst = r_assy.instances[instanceName]

    # 根据配对类型确定参考点名称、节点集标签前缀和需要约束的自由度
    if pair_type == 'Left-Right':
        ref_point_name = 'set_RefPoint_X'
        tag1, tag2 = 'L', 'R'  # tag1 是从面(Slave), tag2 是主面(Master)
        dof_indices = [1, 2, 3]  # 基于X方向周期性，需要约束 U1, U2, U3
    elif pair_type == 'Front-Back':
        ref_point_name = 'set_RefPoint_Y'
        tag1, tag2 = 'F', 'B'  # tag1 是从面(Slave), tag2 是主面(Master)
        dof_indices = [1, 2, 3]  # 基于Y方向周期性，需要约束 U1, U2, U3
    elif pair_type == 'Bottom-Top':
        ref_point_name = 'set_RefPoint_Z'
        tag1, tag2 = 'Bo', 'T'  # tag1 是从面(Slave), tag2 是主面(Master)
        dof_indices = [1, 2, 3]  # 基于Z方向周期性，需要约束 U1, U2, U3
    else:
        print("    ERROR: Unknown pair_type in applyPeriodicConstraints3D: %s" % pair_type)
        return

    # 定义方程约束的系数: u_slave - u_master - u_ref = 0
    # Abaqus 方程项格式: (系数, 节点集名称, 自由度索引)
    # 注意: 参考点系数的符号取决于后续如何施加位移。
    # 假设参考点正位移对应正的平均应变:
    # U_right - U_left = Disp_X => U_left = U_right - Disp_X
    # 方程形式: 1.0 * U_left - 1.0 * U_right + 1.0 * Disp_X = 0
    # Y 和 Z 方向同理。
    coeffs = (1.0, -1.0, 1.0)  # 分别对应 (从节点, 主节点, 参考点) 的系数

    nodes_constrained_in_this_call = 0  # 记录本次调用中实际约束了多少个新的从节点
    for i, (node1, node2) in enumerate(node_pairs):  # node1 是从节点(slave), node2 是主节点(master)

        # --- 检查从节点是否已经被之前的约束处理过 ---
        if node1.label in constrained_nodes_set:
            # 如果从节点的标签已经在集合中，说明它的自由度已被约束，跳过此节点对，避免重复约束
            continue  # 跳到下一个节点对

        # --- 如果从节点未被约束，则创建临时节点集并施加约束 ---
        # 为当前节点对创建临时的节点集名称
        set1_name = 'set_Node-%s-%d' % (tag1, i + 1)  # 从节点集
        set2_name = 'set_Node-%s-%d' % (tag2, i + 1)  # 主节点集
        # 使用 sequenceFromLabels 从实例中获取节点并创建节点集
        # 注意：这里需要确保 node1.label 和 node2.label 能正确获取节点标签
        r_assy.Set(nodes=inst.nodes.sequenceFromLabels(labels=(node1.label,)), name=set1_name)
        r_assy.Set(nodes=inst.nodes.sequenceFromLabels(labels=(node2.label,)), name=set2_name)

        # 根据 dof_indices 列表，为需要的自由度创建 *Equation 约束
        if 1 in dof_indices:  # 约束 DOF 1 (X方向位移 U1)
            model.Equation(name='Eq-%s%s-X-%d' % (tag1, tag2, i + 1),
                           terms=((coeffs[0], set1_name, 1), (coeffs[1], set2_name, 1), (coeffs[2], ref_point_name, 1)))
        if 2 in dof_indices:  # 约束 DOF 2 (Y方向位移 U2)
            model.Equation(name='Eq-%s%s-Y-%d' % (tag1, tag2, i + 1),
                           terms=((coeffs[0], set1_name, 2), (coeffs[1], set2_name, 2), (coeffs[2], ref_point_name, 2)))
        if 3 in dof_indices:  # 约束 DOF 3 (Z方向位移 U3)
            model.Equation(name='Eq-%s%s-Z-%d' % (tag1, tag2, i + 1),
                           terms=((coeffs[0], set1_name, 3), (coeffs[1], set2_name, 3), (coeffs[2], ref_point_name, 3)))

        # --- 将处理过的从节点标签添加到集合中 ---
        constrained_nodes_set.add(node1.label)
        nodes_constrained_in_this_call += 1  # 增加计数

    # 打印本次调用实际约束了多少个新的从节点
    print("    Applied %s constraints to %d new slave nodes." % (pair_type, nodes_constrained_in_this_call))


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
    """构建包含周期性镜像的完整纤维中心列表 (XY平面)"""
    all_centers = []
    rveW, rveH = rveSize[0], rveSize[1]
    for xt, yt in centers:
        all_centers.append((xt, yt))
        if xt < radius: all_centers.append((xt + rveW, yt))
        if xt > rveW - radius: all_centers.append((xt - rveW, yt))
        if yt < radius: all_centers.append((xt, yt + rveH))
        if yt > rveH - radius: all_centers.append((xt, yt - rveH))
        # 角落检查
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
    all_fiber_centers_for_classify = buildAllPeriodicCenters3D(fiber_centers, rveSize, fiberRadius)
    print("    Total fiber centers (with periodicity): %d" % len(all_fiber_centers_for_classify))

    # 步骤3: 详细分类
    print("\n  Step 3: Detailed classification using geometry")
    matrix_fragment_count = 0
    unclassified_cells = []

    for idx, cell in enumerate(potential_cells):
        cell_center = None
        try:
            # 优先使用 pointOn, 它保证在几何体内
            if hasattr(cell, 'pointOn') and cell.pointOn:
                point = cell.pointOn[0]
                cell_center = (point[0], point[1], point[2])
        except:
            pass
        if cell_center is None:
            try:
                # 其次尝试 getCentroid
                cell_centroid = cell.getCentroid()
                if cell_centroid and len(cell_centroid) >= 3:
                    cell_center = (cell_centroid[0], cell_centroid[1], cell_centroid[2])
            except:
                pass
        if cell_center is None:
            try:
                # 最后尝试用顶点均值
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
        for fc_x, fc_y in all_fiber_centers_for_classify:
            dist = math.sqrt((cell_x - fc_x) ** 2 + (cell_y - fc_y) ** 2)
            if dist < min_dist_fib:
                min_dist_fib = dist

        # 判断归属 (孔隙是Cut掉的, 不会在这里被找到)
        # 引入1e-9的公差防止浮点数误差
        if min_dist_fib < (fiberRadius - 1e-9):
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

    if len(fiber_cells_list) != len(all_fiber_centers_for_classify):
        print("\n  WARNING: Classified fiber cell count (%d) does not match"
              " periodic fiber count (%d)." % (len(fiber_cells_list), len(all_fiber_centers_for_classify)))
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
                              pairingToleranceFactor=3.0,
                              rsa_seeding_ratio=0.3,
                              export_coordinates=True,
                              csv_filename=None,
                              # 空隙相关参数
                              enable_void=True,
                              void_porosity=0.02,
                              void_semi_axes_spheroid=[0.01, 0.005],
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
                              cohesive_stab_coeff=0.0001):
    """创建带多空隙(椭球体)缺陷的三维RVE模型主函数"""

    print("\n" + "=" * 70)
    print("Starting 3D RVE Model Generation (WITH MULTI-SPHEROID-VOID DEFECTS)")
    print("Algorithm: Fibers First (F-F), Voids Second (V-V 3D-Axis & F-V 2D-Conservative-Center).")
    print("PBC: Node-Pairing Equations.")
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
    # 根据 Spheroid (旋转椭球体) 参数 [A-长半轴, C-短半轴]
    void_axis_a = void_semi_axes_spheroid[0]
    void_axis_c = void_semi_axes_spheroid[1]
    void_axes = [void_axis_a, void_axis_c]  # 存储 [A, C] 格式用于导出

    # --- F-V 间距 (2D 保守) ---
    # 采用空隙的包围球半径作为2D检测的保守半径 (基于用户请求)
    void_bounding_radius = max(void_axis_a, void_axis_c)
    # 纤维-空隙 最小间距 (基于2D中心, 采用保守的包围球半径)
    minFiberVoidDistance_2D_conservative = fiber_void_min_dist_factor * (fiberRadius + void_bounding_radius)

    # --- V-V 间距 (3D 轴-轴) ---
    # 空隙-空隙 最小轴-轴间距 (基于短半轴直径)
    minVoidDistance_axis = min_void_dist_factor * (void_axis_c * 2.0)

    # 边界检查 (使用真实包围球半径)
    void_r_max_overall = void_bounding_radius

    num_voids_target = 0
    if enable_void and void_porosity > 0:
        target_total_void_volume = void_porosity * rveVolume
        # 体积是 4/3 * pi * A * C * C
        vol_per_void = (4.0 / 3.0) * math.pi * void_axis_a * void_axis_c * void_axis_c
        if vol_per_void == 0:
            raise ValueError("Void semi-axes cannot be zero.")
        num_voids_target = int(round(target_total_void_volume / vol_per_void))

        print("  Target Void Porosity: %.4f%%" % (void_porosity * 100))
        print("  Void Semi-Axes (A-long, C-short): %.6f, %.6f" % (void_axis_a, void_axis_c))
        print("  Calculated Number of Voids: %d" % num_voids_target)
        print(
            "  Min Fiber-Void Center-to-Center Distance (2D Conservative): %.6f" % minFiberVoidDistance_2D_conservative)
        print("  Min Void-Void Axis-to-Axis Distance (3D): %.6f" % minVoidDistance_axis)
    else:
        print("  Void generation SKIPPED (enable_void=False or porosity=0)")
        enable_void = False

    # ==================== 步骤 4: 生成空隙坐标 (V-V & F-V) ====================
    if enable_void:
        print("\nStep 4: Generating void coordinates (V-V & F-V constraints)...")

        void_list = _generate_void_coordinates_rsa(
            num_voids_target, rveSize, void_axis_a,
            minVoidDistance_axis,
            fiber_centers, minFiberVoidDistance_2D_conservative
        )

        if not void_list:
            print("  WARNING: No voids were placed (congestion or 0 target).")
            enable_void = False  # 如果没生成空隙, 则禁用后续步骤
        else:
            if export_void_geometry:
                void_csv_filename = "VoidGeometry_3D_Por%d_%s.csv" % (
                    int(void_porosity * 10000), time.strftime("%Y%m%d_%H%M%S"))
                void_info = exportVoidGeometryToCSV(void_csv_filename, void_list,
                                                    void_axes, rveVolume, void_porosity)
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
    if enable_void:
        print("  Step 8b: Creating void instances (with 3D periodicity) for cutting...")

        # 验证空隙尺寸是否合理（Abaqus对极小几何体处理不佳）
        min_void_size = min(void_axis_a, void_axis_c)
        if min_void_size < 0.0005:  # 小于0.5微米
            print("    WARNING: Void dimensions are very small (min: %.6f mm)" % min_void_size)
            print("    This may cause geometry creation issues in Abaqus.")
            print("    Recommendation: Use void semi-axes >= 0.001 mm")

        cutting_instances = []
        temp_void_part_names = []  # 存储名称以供以后清理
        inst_counter = 0

        # 空隙的最大尺寸，用于边界判断
        void_threshold = void_r_max_overall * void_boundary_threshold_factor

        for i, void_data in enumerate(void_list):
            cx, cy, cz, rot_z, rot_x = void_data
            center = (cx, cy, cz)

            # (智能策略) 判断该空隙需要在哪些镜像位置创建副本
            mirror_offsets = determineMirrorOffsets_3D(center, rveSize, void_threshold)

            # 创建原始空隙实例 (中心RVE内) 和所有必要的镜像实例
            positions_to_create = [(0, 0, 0)]
            positions_to_create.extend(mirror_offsets)

            for offset in positions_to_create:
                inst_counter += 1
                inst_name = 'Void-Inst-%d' % inst_counter

                # 1. 为每个空隙创建独立的椭球体Part
                temp_part_name = 'VoidPart-%d' % inst_counter
                temp_void_part_names.append(temp_part_name)  # Add to cleanup list
                p_void_ellipsoid = model.Part(name=temp_part_name, dimensionality=THREE_D, type=DEFORMABLE_BODY)

                # 2. 创建椭球体几何（通过旋转椭圆轮廓）
                # 注意: 此方法创建的是旋转椭球体(spheroid), A(长半轴)为旋转轴, C(短半轴)为径向轴
                sketch_size = max(4.0 * max(void_axis_a, void_axis_c), 0.02)
                s_ellipse = model.ConstrainedSketch(name='__ellipse_profile_%d__' % inst_counter,
                                                    sheetSize=sketch_size)

                # 绘制旋转轴（Y轴作为构造线）并指定为中心线
                # 这使得A轴(长半轴)初始时沿着Y轴
                s_ellipse.ConstructionLine(point1=(0.0, -2.0 * void_axis_a), point2=(0.0, 2.0 * void_axis_a))
                s_ellipse.assignCenterline(line=s_ellipse.geometry.findAt((0.0, 0.0)))

                # 绘制椭圆的右半部分轮廓 (x >= 0)
                # 径向半径使用 void_axis_c (短半轴)
                num_points = 50
                points = []
                for j in range(num_points + 1):  # 51 points total
                    theta = math.pi * (float(j) / float(num_points) - 0.5)
                    # 右边椭圆轮廓的参数方程（x>=0）
                    # Y半轴为a, X半轴为c
                    x = void_axis_c * math.cos(theta)
                    y = void_axis_a * math.sin(theta)
                    points.append((x, y))

                # 分段创建样条曲线，每段使用较少的点以避免精度问题
                segment_size = 5
                for j in range(0, len(points) - segment_size, segment_size):
                    segment_points = points[j:j + segment_size + 1]
                    if len(segment_points) >= 2:
                        try:
                            s_ellipse.Spline(points=segment_points)
                        except:
                            # 如果样条失败，使用直线段代替
                            for k in range(len(segment_points) - 1):
                                s_ellipse.Line(point1=segment_points[k], point2=segment_points[k + 1])

                # 添加封闭线段（沿Y轴）连接顶部和底部
                s_ellipse.Line(point1=(0.0, void_axis_a), point2=(0.0, -void_axis_a))

                # 3. 通过旋转创建椭球体（绕Y轴旋转360度）
                try:
                    p_void_ellipsoid.BaseSolidRevolve(sketch=s_ellipse, angle=360.0, flipRevolveDirection=OFF)
                    revolve_success = True
                except Exception as e:
                    print("      WARNING: Failed to create void ellipsoid #%d using revolve." % inst_counter)
                    print("      Error: %s" % str(e))
                    print("      Attempting alternative method (sphere approximation)...")
                    revolve_success = False

                # 如果旋转失败，使用球体近似
                if not revolve_success:
                    # 清理失败的草图
                    try:
                        del model.sketches['__ellipse_profile_%d__' % inst_counter]
                    except:
                        pass

                    # 使用平均半径创建球体作为备选方案
                    avg_radius = (void_axis_a + void_axis_c) / 2.0
                    s_sphere = model.ConstrainedSketch(name='__sphere_profile_%d__' % inst_counter,
                                                       sheetSize=sketch_size)
                    s_sphere.ConstructionLine(point1=(0.0, -2.0 * avg_radius), point2=(0.0, 2.0 * avg_radius))
                    s_sphere.assignCenterline(line=s_sphere.geometry.findAt((0.0, 0.0)))
                    # 绘制半圆弧
                    s_sphere.ArcByCenterEnds(center=(0.0, 0.0),
                                             point1=(0.0, avg_radius),
                                             point2=(0.0, -avg_radius),
                                             direction=CLOCKWISE)
                    s_sphere.Line(point1=(0.0, avg_radius), point2=(0.0, -avg_radius))
                    try:
                        p_void_ellipsoid.BaseSolidRevolve(sketch=s_sphere, angle=360.0, flipRevolveDirection=OFF)
                        print("      Successfully created void using sphere approximation.")
                    except Exception as e2:
                        print("      ERROR: Failed to create void even with sphere approximation: %s" % str(e2))
                        # 清理并跳过这个空隙
                        try:
                            del model.sketches['__sphere_profile_%d__' % inst_counter]
                        except:
                            pass
                        try:
                            del model.parts[temp_part_name]
                        except:
                            pass
                        continue

                # 清理临时草图（根据使用的方法删除对应的草图）
                try:
                    if revolve_success:
                        del model.sketches['__ellipse_profile_%d__' % inst_counter]
                    else:
                        del model.sketches['__sphere_profile_%d__' % inst_counter]
                except:
                    pass  # 草图可能已被删除

                # 4. 创建实例
                inst_void = assembly.Instance(name=inst_name, part=p_void_ellipsoid, dependent=OFF)

                # 5. 旋转（Z-X欧拉角）
                # 注意：旋转顺序很重要 - 先绕Z轴旋转，再绕X轴旋转
                if abs(rot_z) > 1e-6:
                    assembly.rotate(instanceList=(inst_name,),
                                    axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1.0),
                                    angle=math.degrees(rot_z))
                if abs(rot_x) > 1e-6:
                    assembly.rotate(instanceList=(inst_name,),
                                    axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0),
                                    angle=math.degrees(rot_x))

                # 6. 平移到最终位置（原始中心 + 周期性偏移）
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
    elemType_bulk_tet = ElemType(elemCode=C3D4, elemLibrary=STANDARD,
                                 elemDeletion=ON, maxDegradation=0.99)
    elemType_bulk_wedge = ElemType(elemCode=C3D6, elemLibrary=STANDARD,
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

    # ==================== 步骤 14: 周期性边界条件 (节点配对) ====================
    print("\nStep 14: Applying Periodic Boundary Conditions (Node-Pairing)...")

    rootAssembly = model.rootAssembly

    # 关键步骤：删除旧的几何实例，并从已划分网格的Part创建新实例
    if 'RVE-3D-1' in rootAssembly.instances:
        del rootAssembly.instances['RVE-3D-1']
        print("  Deleted old geometric instance 'RVE-3D-1'.")

    p_rve = model.parts['RVE-3D']
    rveInstance = rootAssembly.Instance(name='RVE-3D-1', part=p_rve, dependent=ON)
    print("  Created new mesh instance 'RVE-3D-1' in assembly.")

    createReferencePoints3D(model)

    dimensions = getRVEDimensions3D(model, 'RVE-3D-1')
    (nodes_left, nodes_right, nodes_front, nodes_back,
     nodes_bottom, nodes_top) = getBoundaryNodes3D(model, 'RVE-3D-1', dimensions)

    pairing_tolerance = globalSeedSize * pairingToleranceFactor
    print("  PBC pairing tolerance: %.6f" % pairing_tolerance)

    if len(nodes_left) <= len(nodes_right):
        lr_pairs = pairBoundaryNodes3D(nodes_left, nodes_right, pairing_tolerance, (1, 2))
    else:
        lr_pairs = [(p[1], p[0]) for p in pairBoundaryNodes3D(nodes_right, nodes_left, pairing_tolerance, (1, 2))]
    if len(nodes_front) <= len(nodes_back):
        fb_pairs = pairBoundaryNodes3D(nodes_front, nodes_back, pairing_tolerance, (0, 2))
    else:
        fb_pairs = [(p[1], p[0]) for p in pairBoundaryNodes3D(nodes_back, nodes_front, pairing_tolerance, (0, 2))]
    if len(nodes_bottom) <= len(nodes_top):
        bt_pairs = pairBoundaryNodes3D(nodes_bottom, nodes_top, pairing_tolerance, (0, 1))
    else:
        bt_pairs = [(p[1], p[0]) for p in pairBoundaryNodes3D(nodes_top, nodes_bottom, pairing_tolerance, (0, 1))]

    print("\n  Boundary Node Pairing Results:")
    print("    Left/Right (X): %d pairs from %d/%d nodes" % (len(lr_pairs), len(nodes_left), len(nodes_right)))
    print("    Front/Back (Y): %d pairs from %d/%d nodes" % (len(fb_pairs), len(nodes_front), len(nodes_back)))
    print("    Bottom/Top (Z): %d pairs from %d/%d nodes" % (len(bt_pairs), len(nodes_bottom), len(nodes_top)))

    unpaired_nodes_total = (len(nodes_left) - len(lr_pairs)) + (len(nodes_right) - len(lr_pairs)) + \
                           (len(nodes_front) - len(fb_pairs)) + (len(nodes_back) - len(fb_pairs)) + \
                           (len(nodes_bottom) - len(bt_pairs)) + (len(nodes_top) - len(bt_pairs))
    if unpaired_nodes_total > 0:
        print("    WARNING: %d nodes unpaired. Check mesh or tolerance." % unpaired_nodes_total)
    else:
        print("    SUCCESS: All boundary nodes paired successfully!")

    # 初始化一个空的集合，用于跟踪已经被约束过的从节点标签
    constrained_nodes = set()
    print("\n  Applying Left-Right Constraints...")
    applyPeriodicConstraints3D(model, 'RVE-3D-1', lr_pairs, 'Left-Right', constrained_nodes)
    print("\n  Applying Front-Back Constraints...")
    applyPeriodicConstraints3D(model, 'RVE-3D-1', fb_pairs, 'Front-Back', constrained_nodes)
    print("\n  Applying Bottom-Top Constraints...")
    applyPeriodicConstraints3D(model, 'RVE-3D-1', bt_pairs, 'Bottom-Top', constrained_nodes)

    # 提醒用户：此方法不自动创建分析步
    print("\n  Node-Pairing PBCs applied.")
    print("  NOTE: Analysis steps and loads are NOT created automatically with this method.")

    print("\nStep 14 Complete.")

    # ==================== 步骤 15: 清理多余Part ====================
    print("\nStep 15: Cleaning up temporary parts...")
    parts_to_delete = [p for p in ['Fiber', 'Matrix', 'RVE-Temp-WithFibers'] if p in model.parts]

    # Add the temporary void parts to the deletion list
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
    print("  Periodic BCs: Applied using Node-Pairing Equations.")
    if unpaired_nodes_total > 0:
        print("\n  NOTE: %d boundary nodes could not be paired." % unpaired_nodes_total)


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
    TARGET_VF = 0.3
    # RSA播种比例 (0.0 ~ 1.0)，高值(如0.9)排布更均匀，低值(如0.1)更接近物理堆积, 速度较慢
    RSA_SEEDING_RATIO = 0.9

    # ========== 间距参数 ==========
    # 纤维-纤维 最小中心距因子 (min_dist = 因子 * fiberRadius)
    FIBER_MIN_DIST_FACTOR = 2.05

    # 纤维-空隙 最小间距因子 (min_dist = 因子 * (fiberRadius + void_bounding_radius))
    # 这是一个2D保守检测 (void_bounding_radius = max(A, C))
    FIBER_VOID_MIN_DIST_FACTOR = 1.05

    # 空隙-空隙 最小轴-轴间距因子 (min_dist = 因子 * void_short_diameter (C*2))
    # 这是一个3D轴线-轴线距离 (基于短轴直径)
    MIN_VOID_DIST_FACTOR = 1.5

    # ========== 空隙缺陷参数 ==========
    ENABLE_VOID = True  # 是否启用空隙缺陷
    # 目标空隙率(0-1)
    VOID_POROSITY = 0.01

    # 单个旋转椭球体空隙的半轴尺寸 [A-长半轴, C-短半轴] (单位:mm)
    # 脚本将自动计算所需数量以达到 VOID_POROSITY
    VOID_SEMI_AXES_SPHEROID = [0.002, 0.001]

    # 空隙3D周期性边界阈值因子 (threshold = 因子 * max(A, C))
    VOID_BOUNDARY_THRESHOLD_FACTOR = 1.0

    EXPORT_VOID_GEOMETRY = True  # 是否导出空隙几何信息到CSV

    # ========== 网格参数 ==========
    GLOBAL_SEED_SIZE = 0.001  # 全局网格尺寸(单位:mm)
    DEVIATION_FACTOR = 0.1  # 网格偏离因子
    MIN_SIZE_FACTOR = 0.1  # 最小网格尺寸因子

    # ========== 周期性边界条件参数 ==========
    PAIRING_TOLERANCE_FACTOR = 1.0  # 节点配对容差因子(最终容差 = 因子 × 全局网格尺寸)

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
    COHESIVE_T_T = 82.0  # 界面第二切向强度(单位:MPa)
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

    print("\n" + "=" * 70)
    print("Starting 3D RVE Generation (Fiber First / V-V/F-V Axis Check / Node-Pairing-PBC Version)...")
    print("=" * 70)
    print("Target Model Name: %s" % targetModelName)
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
        pairingToleranceFactor=PAIRING_TOLERANCE_FACTOR,
        rsa_seeding_ratio=RSA_SEEDING_RATIO,
        export_coordinates=EXPORT_COORDINATES,
        csv_filename=CSV_FILENAME,
        # 空隙参数
        enable_void=ENABLE_VOID,
        void_porosity=VOID_POROSITY,
        void_semi_axes_spheroid=VOID_SEMI_AXES_SPHEROID,
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
        cohesive_eta=COHESIVE_ETA, cohesive_stab_coeff=COHESIVE_STAB_COEFF
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

        print("\n" + "=" * 70)
        print("MODEL GENERATION COMPLETE")
        print("=" * 70)
    else:
        print("\n" + "!" * 70)
        print("ERROR: Model generation failed. Temporary model not found.")
        print("!" * 70 + "\n")