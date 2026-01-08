# Kadath Boson Star Solvers

基于 Kadath 的玻色星求解器与工具集，包含轴对称与球对称求解、数据转换与分析导出。

## 目录结构
- `src/solvers/axisymmetric/`
  - `rbs.cpp`：轴对称旋转玻色星求解（含 lambda）
  - `msol.cpp`：轴对称参数扫描（内部配置输入文件、tar=0 扫 omega，tar=1 扫 lambda）
- `src/solvers/spherical/`
  - `sph.cpp`：球对称玻色星求解（始终写出 lambda）
- `src/tools/analysis/reader.cpp`：读取解并计算 ADM/Komar/J（自动判型轴/球）
- `src/tools/output/dataout.cpp`：轴对称场 2D 导出
- `src/tools/output/sphdataout.cpp`：球对称导出
- `src/tools/convert/convert_old_to_new.cpp`：旧格式转新格式（补写 lambda）
- `src/utils/io_commons.hpp`：统一读写/判型 I/O 辅助
- `src/archive/`：备份文件（rbscopy, msolcopy）

## 构建
1. 确保环境变量 `HOME_KADATH` 指向 Kadath 源码路径，且其 `build/` 已编译。
2. 在项目根执行：
   ```bash
   cmake -S . -B out/build
   cmake --build out/build -j
   ```
3. 可执行文件输出到 `out/build/bin/`。

## 统一数据格式
- 轴对称：`Space_polar` → `kk:int` → `omega:double` → `lambda:double` → `nu, incA, incB, incbt, phi`
- 球对称：`Space_polar` → `omega:double` → `lambda:double` → `psi, nu, phi`
- 均为大端读写（Kadath `save`/`fwrite_be`）。

## 工具与用法
- 求解轴对称初解：`out/build/bin/rbs`（内部参数见源码）
- 扫描轴对称：`out/build/bin/msol`（在源码顶部配置 `input_file/tar/step/number`，无需命令行参数）
- 求解球对称：`out/build/bin/sph`
- 转换旧数据：
  ```bash
  out/build/bin/convert_old_to_new <old.dat> [output.dat] [lambda_override]
  ```
  默认输出到 `src/tools/converted/<basename>.dat`，若旧文件无 lambda 则补 0 或使用覆盖值。
- 读取与物理量：
  ```bash
  out/build/bin/reader <solution.dat>
  ```
  自动判定轴/球，轴对称输出 ADM/Komar/Js/Jv。
- 导出轴对称切片：`out/build/bin/dataout <solution.dat>` 生成 `metricfields_2d.txt`
- 导出球对称切片：`out/build/bin/sphdataout <solution.dat>` 生成 `sphresult.txt`

## 轴对称扫描（msol.cpp）配置
在 `src/solvers/axisymmetric/msol.cpp` 顶部配置：
```cpp
const char* input_file = "bosinit.dat"; // 初始解
int tar = 0;   // 0: 扫 omega; 1: 扫 lambda
double step = 0.003; // 步长（tar=1 默认可设 1.0）
int number = 8;      // 步数（tar=1 默认可设 2）
```
输出文件名：`bos_<kk>_<omega>_<lambda>.dat`。

## 数据转换注意
- 旧球对称文件若无 lambda，转换后将写入 `lambda=0`。
- 旧轴对称文件若无 lambda，转换后补写 0（或命令行覆盖）。

## 常见问题
- 找不到 `mpi.h`：确认 MPI 及 Kadath 编译环境正确，CMake 会自动添加 MPI include。
- 读取崩溃：确保文件符合统一格式；旧文件请先运行转换工具。

## 许可证
原 Kadath 库遵循其许可证；本仓库代码遵循同一许可证（见 `Kadath/COPYING`）。
