# 3D RVE Model Generator with Voids | 含空隙缺陷的三维RVE模型生成器

## English

### Overview

This repository contains ABAQUS Python scripts for automatically generating 3D Representative Volume Element (RVE) models of fiber-reinforced composites with controllable void defects. The scripts implement advanced algorithms for fiber placement and void distribution, along with automatic material assignment, mesh generation, and periodic boundary conditions.

### Key Features

- **Controllable Fiber Volume Fraction (Vf)**: Automatically calculates the required number of fibers based on target Vf
- **Three Void Shapes**: Supports spherical, cylindrical, and ellipsoidal void geometries
- **Advanced Placement Algorithms**:
  - Hybrid algorithm for fiber placement (RSA seeding + anchored relaxation + forced correction)
  - 3D Random Sequential Adsorption (RSA) for void placement
- **Minimum Distance Control**: Ensures minimum spacing between:
  - Fiber-to-fiber (2D periodic)
  - Void-to-void (3D periodic)
  - Fiber-to-void (2D)
- **Automatic Modeling**:
  - Material property assignment
  - Section definition
  - Mesh generation
  - 3D periodic boundary conditions (PBCs)
  - Cohesive elements insertion at fiber-matrix interfaces (COH3D8/COH3D6)
- **Data Export**: Exports fiber centers and void geometry information to CSV files

### Available Scripts

| Script Name | Void Shape | Description |
|------------|------------|-------------|
| `3D_RVE_Model_with_Voids_Sphere.py` | Sphere | Spherical void defects with uniform radius |
| `3D_RVE_Model_with_Voids_Cylinder.py` | Cylinder | Cylindrical void defects with random 3D orientation |
| `3D_RVE_Model_with_Voids_Ellipsoid.py` | Ellipsoid | Ellipsoidal void defects with random 3D orientation |

### Void Shape Specifications

#### 1. Spherical Voids (`Sphere`)
- **Geometry**: Perfect spheres with uniform radius
- **Parameter**: `VOID_SPHERE_RADIUS` (mm)
- **Rotation**: No rotation needed
- **Use Case**: Ideal for modeling round pores or gas bubbles

#### 2. Cylindrical Voids (`Cylinder`)
- **Geometry**: Cylinders with length and radius
- **Parameters**: `VOID_CYLINDER_DIMS = [length, radius]` (mm)
- **Rotation**: Random 3D orientation (rot_z, rot_x)
- **Use Case**: Suitable for modeling elongated defects or microcracks

#### 3. Ellipsoidal Voids (`Ellipsoid`)
- **Geometry**: Prolate spheroid (ellipsoid of revolution)
- **Parameters**: `VOID_SEMI_AXES_SPHEROID = [A, C]` where A=long axis, C=short axis (mm)
- **Rotation**: Random 3D orientation (rot_z, rot_x)
- **Use Case**: Ideal for modeling naturally occurring voids with aspect ratio

### Requirements

- **Software**: ABAQUS 2023 (or compatible versions)
- **Python**: 2.7 (ABAQUS built-in)
- **Operating System**: Windows/Linux

### Installation

```bash
# Clone the repository
git clone https://github.com/ZPL-03/3D-RVE-Voids-Generator.git

# Navigate to the directory
cd 3D-RVE-Voids-Generator
```

### Usage

#### Method 1: Run from ABAQUS CAE

1. Open ABAQUS CAE
2. Go to `File` → `Run Script`
3. Select the desired script (Sphere/Cylinder/Ellipsoid version)
4. The model will be automatically generated

#### Method 2: Run from Command Line

```bash
# For Windows
abaqus cae noGUI=3D_RVE_Model_with_Voids_Sphere.py

# For Linux
abaqus cae nogui=3D_RVE_Model_with_Voids_Sphere.py
```

### Key Parameters

All scripts share the same parameter structure. Edit the configuration section at the bottom of each script:

#### RVE Geometry
```python
RVE_SIZE = [0.057, 0.057, 0.01]  # [Width_X, Height_Y, Depth_Z] in mm
FIBER_RADIUS = 0.0035            # Fiber radius in mm
TARGET_VF = 0.5                  # Target fiber volume fraction (0-1)
```

#### Void Parameters (Shape-Specific)

**For Sphere:**
```python
ENABLE_VOID = True
VOID_POROSITY = 0.01             # Target void porosity (0-1)
VOID_SPHERE_RADIUS = 0.0015      # Sphere radius in mm
```

**For Cylinder:**
```python
ENABLE_VOID = True
VOID_POROSITY = 0.01
VOID_CYLINDER_DIMS = [0.002, 0.001]  # [length, radius] in mm
```

**For Ellipsoid:**
```python
ENABLE_VOID = True
VOID_POROSITY = 0.01
VOID_SEMI_AXES_SPHEROID = [0.002, 0.001]  # [A-long axis, C-short axis] in mm
```

#### Spacing Control
```python
FIBER_MIN_DIST_FACTOR = 2.05          # Fiber-fiber spacing factor
FIBER_VOID_MIN_DIST_FACTOR = 1.05     # Fiber-void spacing factor
MIN_VOID_DIST_FACTOR = 2.05           # Void-void spacing factor
```

#### Mesh Settings
```python
GLOBAL_SEED_SIZE = 0.0005        # Global mesh size in mm
DEVIATION_FACTOR = 0.1           # Mesh deviation factor
MIN_SIZE_FACTOR = 0.1            # Minimum mesh size factor
```

### Material Properties

The scripts include predefined material properties for:
- **Fiber**: Transversely isotropic material (Carbon fiber)
- **Matrix**: Elastic-plastic material with Drucker-Prager plasticity and damage
- **Interface**: Cohesive zone model with traction-separation law

All material properties can be customized in the configuration section.

### Output Files

1. **ABAQUS Model File**: `.cae` file with complete 3D RVE model
2. **Fiber Coordinates**: `*_fiber_centers.csv` (optional)
3. **Void Geometry**: `*_void_geometry.csv` (optional)

### Algorithm Details

#### Fiber Placement
1. **RSA Seeding**: Initial fiber placement using Random Sequential Adsorption
2. **Anchored Relaxation**: Iterative adjustment to eliminate overlaps while maintaining initial distribution
3. **Forced Correction**: Final check and correction for any remaining violations

#### Void Placement
- **3D RSA**: Random Sequential Adsorption in 3D space with periodic boundary consideration
- **Orientation**: Random 3D rotation for cylindrical and ellipsoidal voids
- **Collision Detection**: Ensures minimum spacing with existing voids and fibers

### Reference

If you use this code in your research, please cite:

```bibtex
@article{Li2025,
  title={A novel method for constructing 3D void RVE elements and rapid homogenization of composite materials},
  author={Li, X. and others},
  journal={Composite Structures},
  volume={360},
  pages={119040},
  year={2025}
}
```

### Author

**Liu Zhengpeng (刘正鹏)**
- GitHub: [@ZPL-03](https://github.com/ZPL-03)
- CSDN/知乎: @小盆i
- Email: 
  - 1370872708@qq.com
  - Zhengpeng0105@gmail.com

### Version History

- **v3.0** (2025-10-23): Current version with three void shape options
- **v2.0**: Added void defect functionality
- **v1.0**: Initial release with fiber placement only

### License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Acknowledgments

- ABAQUS scripting interface documentation
- Composite materials research community

---

## 中文

### 概述

本仓库包含用于自动生成含可控空隙缺陷的纤维增强复合材料三维代表性体积单元（RVE）模型的ABAQUS Python脚本。脚本实现了先进的纤维排布和空隙分布算法，以及材料赋值、网格划分和周期性边界条件的自动化。

### 主要功能

- **可控纤维体积分数（Vf）**：根据目标Vf自动计算所需纤维数量
- **三种空隙形状**：支持球形、圆柱形和椭球形空隙几何
- **先进的排布算法**：
  - 纤维排布的混合算法（RSA播种 + 锚定松弛 + 强制校正）
  - 空隙排布的三维随机顺序吸附（RSA）算法
- **最小间距控制**：确保以下间距：
  - 纤维-纤维（2D周期性）
  - 空隙-空隙（3D周期性）
  - 纤维-空隙（2D）
- **自动化建模**：
  - 材料属性赋值
  - 截面定义
  - 网格生成
  - 三维周期性边界条件（PBCs）
  - 纤维-基体界面内聚力单元插入（COH3D8/COH3D6）
- **数据导出**：导出纤维中心和空隙几何信息至CSV文件

### 可用脚本

| 脚本名称 | 空隙形状 | 描述 |
|---------|---------|------|
| `3D_RVE_Model_with_Voids_Sphere.py` | 球形 | 均匀半径的球形空隙缺陷 |
| `3D_RVE_Model_with_Voids_Cylinder.py` | 圆柱形 | 随机三维方向的圆柱形空隙缺陷 |
| `3D_RVE_Model_with_Voids_Ellipsoid.py` | 椭球形 | 随机三维方向的椭球形空隙缺陷 |

### 空隙形状说明

#### 1. 球形空隙（`Sphere`）
- **几何形状**：均匀半径的完美球体
- **参数**：`VOID_SPHERE_RADIUS`（单位：mm）
- **旋转**：无需旋转
- **适用场景**：适合模拟圆形孔隙或气泡

#### 2. 圆柱形空隙（`Cylinder`）
- **几何形状**：具有长度和半径的圆柱体
- **参数**：`VOID_CYLINDER_DIMS = [长度, 半径]`（单位：mm）
- **旋转**：随机三维方向（rot_z, rot_x）
- **适用场景**：适合模拟细长缺陷或微裂纹

#### 3. 椭球形空隙（`Ellipsoid`）
- **几何形状**：旋转椭球体（长椭球）
- **参数**：`VOID_SEMI_AXES_SPHEROID = [A, C]`，其中A为长半轴，C为短半轴（单位：mm）
- **旋转**：随机三维方向（rot_z, rot_x）
- **适用场景**：适合模拟具有长径比的自然形成空隙

### 软件要求

- **软件**：ABAQUS 2023（或兼容版本）
- **Python**：2.7（ABAQUS内置）
- **操作系统**：Windows/Linux

### 安装

```bash
# 克隆仓库
git clone https://github.com/ZPL-03/3D-RVE-Voids-Generator.git

# 进入目录
cd 3D-RVE-Voids-Generator
```

### 使用方法

#### 方法1：在ABAQUS CAE中运行

1. 打开ABAQUS CAE
2. 选择 `文件` → `运行脚本`
3. 选择所需脚本（球形/圆柱形/椭球形版本）
4. 模型将自动生成

#### 方法2：命令行运行

```bash
# Windows系统
abaqus cae noGUI=3D_RVE_Model_with_Voids_Sphere.py

# Linux系统
abaqus cae nogui=3D_RVE_Model_with_Voids_Sphere.py
```

### 关键参数

所有脚本共享相同的参数结构。编辑每个脚本底部的配置部分：

#### RVE几何参数
```python
RVE_SIZE = [0.057, 0.057, 0.01]  # [宽度X, 高度Y, 深度Z]，单位：mm
FIBER_RADIUS = 0.0035            # 纤维半径，单位：mm
TARGET_VF = 0.5                  # 目标纤维体积分数（0-1）
```

#### 空隙参数（形状特定）

**球形：**
```python
ENABLE_VOID = True
VOID_POROSITY = 0.01             # 目标空隙率（0-1）
VOID_SPHERE_RADIUS = 0.0015      # 球体半径，单位：mm
```

**圆柱形：**
```python
ENABLE_VOID = True
VOID_POROSITY = 0.01
VOID_CYLINDER_DIMS = [0.002, 0.001]  # [长度, 半径]，单位：mm
```

**椭球形：**
```python
ENABLE_VOID = True
VOID_POROSITY = 0.01
VOID_SEMI_AXES_SPHEROID = [0.002, 0.001]  # [A-长半轴, C-短半轴]，单位：mm
```

#### 间距控制
```python
FIBER_MIN_DIST_FACTOR = 2.05          # 纤维-纤维间距因子
FIBER_VOID_MIN_DIST_FACTOR = 1.05     # 纤维-空隙间距因子
MIN_VOID_DIST_FACTOR = 2.05           # 空隙-空隙间距因子
```

#### 网格设置
```python
GLOBAL_SEED_SIZE = 0.0005        # 全局网格尺寸，单位：mm
DEVIATION_FACTOR = 0.1           # 网格偏离因子
MIN_SIZE_FACTOR = 0.1            # 最小网格尺寸因子
```

### 材料属性

脚本包含预定义的材料属性：
- **纤维**：横观各向同性材料（碳纤维）
- **基体**：弹塑性材料，含Drucker-Prager塑性和损伤
- **界面**：含牵引-分离定律的内聚力区模型

所有材料属性均可在配置部分自定义。

### 输出文件

1. **ABAQUS模型文件**：`.cae`文件，包含完整的三维RVE模型
2. **纤维坐标**：`*_fiber_centers.csv`（可选）
3. **空隙几何**：`*_void_geometry.csv`（可选）

### 算法详情

#### 纤维排布
1. **RSA播种**：使用随机顺序吸附算法初始排布纤维
2. **锚定松弛**：迭代调整以消除重叠，同时保持初始分布
3. **强制校正**：最终检查和修正任何剩余违规

#### 空隙排布
- **三维RSA**：考虑周期性边界的三维空间随机顺序吸附
- **方向**：圆柱形和椭球形空隙的随机三维旋转
- **碰撞检测**：确保与现有空隙和纤维的最小间距

### 参考文献

如果您在研究中使用此代码，请引用：

```bibtex
@article{Li2025,
  title={A novel method for constructing 3D void RVE elements and rapid homogenization of composite materials},
  author={Li, X. and others},
  journal={Composite Structures},
  volume={360},
  pages={119040},
  year={2025}
}
```

### 作者

**刘正鹏 (Liu Zhengpeng)**
- GitHub: [@ZPL-03](https://github.com/ZPL-03)
- CSDN/知乎: @小盆i
- 邮箱: 
  - 1370872708@qq.com
  - Zhengpeng0105@gmail.com

### 版本历史

- **v3.0** (2025-10-23)：当前版本，包含三种空隙形状选项
- **v2.0**：添加空隙缺陷功能
- **v1.0**：初始版本，仅含纤维排布

### 许可证

本项目采用MIT许可证 - 详见[LICENSE](LICENSE)文件。

### 致谢

- ABAQUS脚本接口文档
- 复合材料研究社区

---


**如有问题或建议，欢迎提Issue或Pull Request！**

**For questions or suggestions, feel free to open an Issue or Pull Request!**
