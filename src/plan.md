# src 代码整理优化计划

## 当前问题分析

### 代码分类现状

**求解器类**（3个）：

- `sph.cpp` - 球对称玻色星求解器（3场：psi, nu, phi）
- `rbs.cpp` - 轴对称旋转玻色星求解器（5场：nu, incA, incB, incbt, phi），两阶段求解
- `msol.cpp` - 参数扫描求解器（从已有解继续，改变 omega）

**数据工具类**（4个）：

- `reader.cpp` - 读取解并计算物理量（ADM质量、Komar质量、角动量等）
- `dataout.cpp` - 输出度量场到2D文本文件（r, theta, ap, A, B, bt, phi）
- `sphdataout.cpp` - 球对称数据输出（r, A, B, phi，仅 theta=0）
- `initdataout.cpp` - 初始数据输出（r, theta, rsint, phi）

**初始化工具类**（1个）：

- `in.cpp` - 生成初始 phi 场并保存

**测试工具类**（1个）：

- `test.cpp` - 测试代码

**备份文件**（2个）：

- `rbscopy.cpp` - rbs.cpp 的备份
- `msolcopy.cpp` - msol.cpp 的备份

### 主要问题

1. **I/O 格式分裂**：轴对称含 `lambda`，球对称未写 `lambda`（等价 `lambda=0`），导致 reader/dataout 崩溃或格式不统一。
2. **文件命名不规范**：`rbs`, `msol`, `sph` 等缩写不直观。
3. **代码重复**：多个文件有相似的 MPI 初始化、GPU 初始化、空间设置代码。
4. **功能混乱**：数据输出工具功能重叠（dataout, sphdataout, initdataout）。
5. **缺少文档**：代码注释不足，难以理解各文件用途。
6. **目录结构混乱**：所有文件都在 src 根目录，没有分类。
7. **备份文件冗余**：rbscopy.cpp 和 msolcopy.cpp 应该删除或归档。

## 优化方案

### 阶段一：目录重组

创建清晰的目录结构：

```
src/
├── solvers/              # 求解器
│   ├── spherical/        # 球对称求解器
│   │   └── boson_star_spherical.cpp
│   └── axisymmetric/     # 轴对称求解器
│       ├── boson_star_rotating.cpp
│       └── boson_star_scan.cpp
├── tools/                # 工具程序
│   ├── analysis/        # 分析工具
│   │   └── compute_quantities.cpp
│   ├── output/           # 输出工具
│   │   ├── export_metric_fields.cpp
│   │   └── export_initial_data.cpp
│   └── initialization/   # 初始化工具
│       └── generate_initial_phi.cpp
├── utils/                # 公共工具
│   ├── mpi_utils.hpp     # MPI 初始化封装
│   ├── gpu_utils.hpp     # GPU 初始化封装
│   └── space_utils.hpp   # 空间设置工具
└── archive/              # 归档（备份文件）
    ├── rbscopy.cpp
    └── msolcopy.cpp
```

### 阶段二：代码重构

#### 2.1 提取公共代码

**创建 `utils/mpi_utils.hpp`**：

```cpp
// MPI 初始化和清理的封装
struct MPIContext {
    int rank;
    int size;
    MPIContext(int argc, char** argv);
    ~MPIContext();
};
```

**创建 `utils/gpu_utils.hpp`**：

```cpp
// GPU 初始化的封装
struct GPUContext {
    GPUContext(int rank);
    ~GPUContext();
};
```

**创建 `utils/space_utils.hpp`**：

```cpp
// 空间设置的辅助函数
Space_polar create_polar_space(int resol, const Array<double>& bounds);
void setup_rsint(Space_polar& space, Scalar& rsint);
```

#### 2.2 统一求解器接口

所有求解器遵循统一模式：

1. 参数解析（命令行或配置文件）
2. 空间初始化
3. 场初始化
4. 系统设置
5. 求解
6. 结果保存

#### 2.3 统一 I/O（含 lambda，自动判型）

- 统一二进制格式：写入顺序为 `Space_*` →（轴对称含 `kk:int`，球对称可占位或省略）→ `omega:double` → `lambda:double` → 场数据。
- 提供公共接口 `save_solution / load_solution`（新建 `include/io_commons.hpp` 或 `src/utils/io_commons.hpp`），所有求解器与工具复用。
- 判型策略：读取 `Space` 后尝试读 `kk`（成功则轴对称，失败/不足则球对称），再读 `omega`、`lambda`。
- 为球对称始终写 `lambda=0`（可写占位 `kk=0` 以完全统一）。
- 合并数据输出工具：`dataout.cpp`, `sphdataout.cpp`, `initdataout.cpp` 合并为单一工具（如 `export_data.cpp`），根据判型与模式参数输出。

#### 2.4 一次性旧→新转换脚本

- 编写 `tools/convert_old_to_new.cpp`：读取旧格式（无 lambda/无 kk 的球对称或旧轴对称），按新规范重写，生成新文件名或覆盖。
- 在 README/注释中给出使用说明，避免误覆盖原始数据。

### 阶段三：重命名文件


| 原文件名 | 新文件名 | 说明 |

|---------|---------|------|

| `sph.cpp` | `boson_star_spherical.cpp` | 球对称玻色星求解器 |

| `rbs.cpp` | `boson_star_rotating.cpp` | 轴对称旋转玻色星求解器 |

| `msol.cpp` | `boson_star_scan.cpp` | 参数扫描求解器 |

| `reader.cpp` | `compute_quantities.cpp` | 计算物理量工具 |

| `dataout.cpp` | `export_metric_fields.cpp` | 输出度量场 |

| `sphdataout.cpp` | `export_spherical_data.cpp` | 输出球对称数据 |

| `initdataout.cpp` | `export_initial_data.cpp` | 输出初始数据 |

| `in.cpp` | `generate_initial_phi.cpp` | 生成初始场 |



### 阶段四：代码优化

1. **添加命令行参数解析**：使用统一的参数解析库
2. **改进错误处理**：统一的错误处理机制
3. **添加日志系统**：统一的日志输出
4. **代码注释**：为每个函数添加详细注释
5. **配置文件支持**：支持从配置文件读取参数

### 阶段五：文档完善

1. **README.md**：说明目录结构和各文件用途
2. **各文件头部注释**：说明文件功能、使用方法、参数
3. **使用示例**：提供典型使用案例

## 实施步骤

### 步骤 1：创建新目录结构

- 创建 `solvers/`, `tools/`, `utils/`, `archive/` 目录
- 创建子目录 `spherical/`, `axisymmetric/`, `analysis/`, `output/`, `initialization/`

### 步骤 2：提取公共代码

- 创建 `utils/mpi_utils.hpp`
- 创建 `utils/gpu_utils.hpp`
- 创建 `utils/space_utils.hpp`

### 步骤 3：统一 I/O 接口并提供转换脚本

- 创建 `include/io_commons.hpp` 或 `src/utils/io_commons.hpp`，实现 `save_solution / load_solution`，统一写 `lambda`，判型 `kk`。
- 编写 `tools/convert_old_to_new.cpp`，将旧格式（无 lambda/无 kk）转为新格式（`lambda` 始终写出，球对称可补 `kk=0`）。

### 步骤 4：重构求解器

- 重构 `sph.cpp` → `solvers/spherical/boson_star_spherical.cpp`，输出时写入 `omega` 与 `lambda=0`（可写 `kk=0` 占位）。
- 重构 `rbs.cpp` → `solvers/axisymmetric/boson_star_rotating.cpp`（含 kk, lambda）。
- 重构 `msol.cpp` → `solvers/axisymmetric/boson_star_scan.cpp`，统一用新 I/O。

### 步骤 5：重构工具程序

- 重构 `reader.cpp` → `tools/analysis/compute_quantities.cpp`，使用统一 `load_solution` 判型。
- 合并并重构数据输出工具 → `tools/output/export_data.cpp`（支持球/轴对称，模式参数）。
- 重构 `in.cpp` → `tools/initialization/generate_initial_phi.cpp`。

### 步骤 6：归档备份文件

- 移动 `rbscopy.cpp` → `archive/`
- 移动 `msolcopy.cpp` → `archive/`
- 删除 `test.cpp` 或移动到 `archive/`

### 步骤 7：更新 CMakeLists.txt

- 更新可执行文件路径
- 添加新的工具库链接

### 步骤 8：编写文档

- 创建 `src/README.md`
- 为每个文件添加头部注释

## 预期效果

1. **目录结构清晰**：按功能分类，易于查找
2. **代码复用性提高**：公共代码提取到 utils
3. **命名规范统一**：文件名直观反映功能
4. **维护性提升**：代码组织合理，易于修改和扩展
5. **文档完善**：新用户能快速理解代码结构

## 风险评估

- **低风险**：目录重组、文件重命名（git 可以追踪）
- **中风险**：代码重构（需要充分测试）
- **建议**：先在 git 中创建新分支，逐步迁移，保留原文件直到新代码测试通过

