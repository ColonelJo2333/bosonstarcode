# Kadath Boson Star Solvers

基于 Kadath 的玻色星求解器与工具集，包含轴对称与球对称求解、数据转换与分析导出。

## 目录结构
- `src/solvers/axisymmetric/`
  - `rbs.cpp`：轴对称旋转玻色星求解（含 lambda）
  - `msol.cpp`：轴对称参数扫描（内部配置输入文件、tar=0 扫 omega，tar=1 扫 lambda）
- `src/solvers/spherical/`
  - `sph.cpp`：球对称玻色星求解（始终写出 lambda）
- `src/tools/analysis/reader.cpp`：读取解并计算/导出（自动判型轴/球，命令行传入模式与路径）
- `src/tools/convert/convert_old_to_new.cpp`：旧格式转新格式（补写 lambda，输入/输出路径在源码顶部配置）
- `src/utils/io_commons.hpp`：统一读写/判型 I/O 辅助
- `src/rbscopy.cpp` / `src/msolcopy.cpp`：备份文件
- `src/plan.md`：开发记录与规划

## 构建
1. 确保 Kadath 源码与其 `build/` 已编译。
2. 在项目根执行（本仓库的 `CMakeLists.txt` 仅本地使用）：
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
  - 在 `src/tools/convert/convert_old_to_new.cpp` 顶部修改 `input` / `output` / `lambda_override` 配置，重新编译后直接运行 `out/build/bin/convert_old_to_new`；若旧文件无 lambda 则补 0 或使用覆盖值。
- 读取与物理量：
  - `out/build/bin/reader <solution.dat> 0`：计算模式（自动判定轴/球；轴对称输出 ADM/Komar/Js/Jv，球对称输出 ADM/Komar）
  - `out/build/bin/reader <solution.dat> 1 [output.txt]`：导出模式（自动判定轴/球；导出真实度规场与标量场）
  - 导出文件首行注释含 `omega/lambda`，第二行注释为列名，后续为逐点数据

## 轴对称扫描（msol.cpp）配置
在 `src/solvers/axisymmetric/msol.cpp` 顶部配置：
```cpp
const char* input_file = "bosinit.dat"; // 初始解
int tar = 0;   // 0: 扫 omega; 1: 扫 lambda
double step = 0.003; // 步长（tar=1 默认可设 1.0）
int number = 8;      // 步数（tar=1 默认可设 2）
```
输出文件名：`bos_<kk>_<omega>_<lambda>.dat`。

## reader 导出字段说明
- 球对称：`Psi=exp(psi)`，`N=exp(nu)`，导出 `Psi N phi`
- 轴对称：`ap=exp(nu)`，`A=exp(incA-nu)`，`B=(incB.div_rsint()+1)/ap`，`bt=incbt.div_rsint()`，导出 `ap A B bt phi`

## 数据转换注意
- 旧球对称文件若无 lambda，转换后将写入 `lambda=0`。
- 旧轴对称文件若无 lambda，转换后补写 0（或命令行覆盖）。

## 常见问题
- 找不到 `mpi.h`：确认 MPI 及 Kadath 编译环境正确，CMake 会自动添加 MPI include。
- 读取崩溃：确保文件符合统一格式；旧文件请先运行转换工具。

## 许可证
原 Kadath 库遵循其许可证；本仓库代码遵循同一许可证。
