# Kadath 库核心架构中文指南

## 目录

1. [概述](#概述)
2. [核心类解析](#核心类解析)
3. [内存管理](#内存管理)
4. [实战求解流程](#实战求解流程)
5. [常见问题与注意事项](#常见问题与注意事项)

---

## 概述

Kadath 是一个基于谱方法的数值求解库，主要用于求解偏微分方程（PDE）。其核心思想是将计算区域分解为多个**域（Domain）**，每个域内使用谱展开（Chebyshev 或 Legendre 多项式）进行离散化。

### 核心概念

- **配置点（Colocation Points）**：物理空间中的离散点，用于计算场值
- **系数空间（Coefficient Space）**：谱展开的系数，用于存储谱表示
- **域（Domain）**：计算区域的一个子区域，具有特定的几何形状和坐标映射
- **空间（Space）**：由多个域组成的完整计算区域
- **系统（System_of_eqs）**：包含未知量、方程和边界条件的求解系统

---

## 核心类解析

### 1. Space 类及其派生类

`Space` 是抽象基类，表示整个计算区域。它管理多个 `Domain` 对象，定义了计算空间的拓扑结构。

#### Space 基类

```cpp
class Space {
protected:
    int nbr_domains;      // 域的数量
    int ndim;             // 维度数
    int type_base;        // 基类型（CHEB_TYPE 或 LEG_TYPE）
    Domain** domains;     // 指向各个域的指针数组
};
```

**关键方法**：
- `get_domain(int i)`：获取第 i 个域的指针
- `get_nbr_domains()`：返回域的数量
- `get_ndim()`：返回维度数

#### 主要派生类

##### Space_spheric（球坐标系空间）

用于三维球对称或轴对称问题。

**构造函数**：
```cpp
Space_spheric(int ttype, const Point& cr, const Dim_array& nbr, 
              const Array<double>& bounds, bool withzec=true)
```

**参数说明**：
- `ttype`：基类型，`CHEB_TYPE`（1）或 `LEG_TYPE`（2）
- `cr`：球坐标中心点（`Point` 对象，包含 X, Y, Z 坐标）
- `nbr`：每个域的配置点数（`Dim_array`，如 `Dim_array(3, 17, 9, 9)` 表示 3 维，径向 17 点，θ 方向 9 点，φ 方向 9 点）
- `bounds`：各域边界的半径数组（`Array<double>`），例如 `{0.5, 1.0, 2.0}` 表示三个域
- `withzec`：是否包含紧致化域（用于处理无穷远边界）

**示例**：
```cpp
#include "Kadath/include/kadath.hpp"
using namespace Kadath;

// 定义中心点
Point center(3);
center.set(1) = 0.0;  // X 坐标
center.set(2) = 0.0;  // Y 坐标
center.set(3) = 0.0;  // Z 坐标

// 定义配置点数：3维，径向17点，θ方向9点，φ方向9点
Dim_array res(3, 17, 9, 9);

// 定义域边界：核域半径0.5，壳层1.0，外域2.0
Array<double> bounds(3);
bounds.set(0) = 0.5;
bounds.set(1) = 1.0;
bounds.set(2) = 2.0;

// 创建球坐标系空间（使用 Chebyshev 基）
Space_spheric space(CHEB_TYPE, center, res, bounds, true);
```

##### Space_polar（极坐标空间）

用于二维轴对称问题（关于 φ 对称）。

**构造函数**：
```cpp
Space_polar(int ttype, const Point& cr, const Dim_array& nbr, 
            const Array<double>& bounds)
```

**参数说明**：
- `ttype`：基类型
- `cr`：中心点（2维：X, Z）
- `nbr`：配置点数（2维：径向和 θ 方向）
- `bounds`：域边界半径数组

**示例**：
```cpp
Point center(2);
center.set(1) = 0.0;  // X 坐标
center.set(2) = 0.0;  // Z 坐标

Dim_array res(2, 17, 9);  // 2维，径向17点，θ方向9点

Array<double> bounds(2);
bounds.set(0) = 0.5;
bounds.set(1) = 1.0;

Space_polar space(CHEB_TYPE, center, res, bounds);
```

##### Space_oned（一维空间）

用于一维球对称问题。

**构造函数**：
```cpp
Space_oned(int ttype, const Dim_array& nbr, const Array<double>& bounds)
```

**示例**：
```cpp
Dim_array res(1, 17);  // 1维，17个配置点

Array<double> bounds(2);
bounds.set(0) = 0.5;
bounds.set(1) = 1.0;

Space_oned space(CHEB_TYPE, res, bounds);
```

---

### 2. Domain 类及其派生类

`Domain` 是抽象基类，表示单个计算域。每个域定义了：
- 数值坐标到物理坐标的映射
- 配置点的位置
- 谱基的选择

#### Domain 基类关键成员

```cpp
class Domain {
protected:
    int num_dom;                    // 域编号
    int ndim;                       // 维度
    Dim_array nbr_points;           // 配置点数量
    Dim_array nbr_coefs;            // 系数数量
    int type_base;                  // 基类型
    Memory_mapped_array<Array<double>*> coloc;  // 配置点
    mutable Memory_mapped_array<Val_domain*> cart;  // 笛卡尔坐标
    mutable Val_domain* radius;      // 广义半径
};
```

#### 主要派生类

##### Domain_nucleus（核域）

包含原点的球域，用于处理奇点。

**构造函数**：
```cpp
Domain_nucleus(int num, int ttype, double radius, 
               const Point& cr, const Dim_array& nbr)
```

**参数说明**：
- `num`：域编号
- `ttype`：基类型
- `radius`：域的外半径
- `cr`：中心点
- `nbr`：配置点数

**坐标映射**：
- 数值坐标：`0 ≤ x ≤ 1`, `0 ≤ θ* ≤ π/2`, `0 ≤ φ* < 2π`
- 物理坐标：`r = αx`（其中 α 是半径）

##### Domain_shell（壳层域）

两个同心球之间的壳层。

**构造函数**：
```cpp
Domain_shell(int num, int ttype, double r_int, double r_ext, 
             const Point& cr, const Dim_array& nbr)
```

**参数说明**：
- `r_int`：内半径
- `r_ext`：外半径

**坐标映射**：
- 数值坐标：`-1 ≤ x ≤ 1`
- 物理坐标：`r = αx + β`（线性映射）

##### Domain_compact（紧致化域）

用于处理无穷远边界，通过坐标紧致化将无穷远映射到有限区间。

**构造函数**：
```cpp
Domain_compact(int num, int ttype, double r_int, 
               const Point& cr, const Dim_array& nbr)
```

**坐标映射**：
- 数值坐标：`-1 ≤ x ≤ 1`
- 物理坐标：`r = 1/(αx - 1)`（将无穷远映射到 x=1）

---

### 3. System_of_eqs 类

`System_of_eqs` 是求解系统的核心类，负责管理未知量、方程和求解过程。

#### 构造函数

```cpp
// 在整个空间求解
System_of_eqs(const Space& so)

// 在指定域范围求解
System_of_eqs(const Space& so, int dom_min, int dom_max)

// 在单个域求解
System_of_eqs(const Space& so, int dom)
```

#### 关键成员

```cpp
class System_of_eqs {
protected:
    const Space& espace;           // 关联的空间
    int dom_min, dom_max;          // 求解域范围
    int nvar;                      // 未知场数量
    MMPtr_array<Tensor> var;       // 未知场指针数组
    MMPtr_array<char> names_var;   // 未知场名称
    int nvar_double;               // 未知标量数量
    MMPtr_array<double> var_double;// 未知标量指针
    int neq;                       // 方程数量
    MMPtr_array<Equation> eq;     // 方程指针数组
};
```

#### 添加未知量

```cpp
// 添加未知场（Tensor）
void add_var(const char* name, Tensor& var)

// 添加未知标量（double）
void add_var(const char* name, double& var)

// 添加常数场
void add_cst(const char* name, const Tensor& cst)

// 添加常数标量
void add_cst(const char* name, double cst)
```

**示例**：
```cpp
Space_spheric space(...);
System_of_eqs syst(space);

// 创建未知标量场
Scalar phi(space);
syst.add_var("phi", phi);

// 创建未知标量
double omega = 0.5;
syst.add_var("omega", omega);

// 添加常数
syst.add_cst("pi", 3.14159);
```

#### 添加方程

```cpp
// 域内方程（假设二阶）
void add_eq_inside(int dom, const char* eq, 
                   int n_cmp=-1, Array<int>** p_cmp=nullptr)

// 边界条件
void add_eq_bc(int dom, int bb, const char* eq, ...)

// 匹配条件（域间连续性）
void add_eq_matching(int dom, int bb, const char* eq, ...)

// 指定阶数的方程
void add_eq_order(int dom, int order, const char* eq, ...)
```

**方程字符串语法**：
- `"lap(phi) = 0"`：拉普拉斯方程
- `"phi = 1"`：边界条件
- `"phi = phi"`：匹配条件（左右两侧的场值相等）
- `"dr(phi) = dr(phi)"`：匹配条件（法向导数连续）

**示例**：
```cpp
// 在所有域添加拉普拉斯方程
for (int d=0; d<space.get_nbr_domains(); d++) {
    syst.add_eq_inside(d, "lap(phi) = 0");
}

// 在原点添加边界条件
syst.add_eq_bc(0, INNER_BC, "phi = 1");

// 在域间添加匹配条件
syst.add_eq_matching(0, OUTER_BC, "phi = phi");
syst.add_eq_matching(0, OUTER_BC, "dr(phi) = dr(phi)");
```

#### 求解方法

```cpp
// 执行一步牛顿迭代
template<Computational_model computational_model = default_computational_model>
bool do_newton(double prec, double &error)

// 带线搜索的牛顿迭代
bool do_newton_with_linesearch(double precision, double& error, 
                               int ntrymax=10, double stepmax=1.0)
```

**示例**：
```cpp
double error = 1.0;
double precision = 1e-10;
int max_iter = 100;

for (int iter=0; iter<max_iter; iter++) {
    bool converged = syst.do_newton(precision, error);
    std::cout << "迭代 " << iter << ": 误差 = " << error << std::endl;
    
    if (converged) {
        std::cout << "收敛！" << std::endl;
        break;
    }
}
```

---

## 内存管理

Kadath 的内存管理有其特殊性，需要特别注意。

### 需要 new 的对象

以下对象**必须**使用 `new` 创建，并由用户管理生命周期：

1. **Space 对象**：通常使用栈分配，但如果需要动态创建，应使用 `new`
2. **Tensor/Scalar 对象**：作为未知量添加到系统时，必须保证在系统生命周期内有效

**正确示例**：
```cpp
// Space 通常在栈上创建
Space_spheric space(CHEB_TYPE, center, res, bounds);

// Scalar 作为未知量，必须在系统生命周期内有效
Scalar phi(space);
System_of_eqs syst(space);
syst.add_var("phi", phi);  // phi 必须保持有效

// 求解后可以安全使用 phi
syst.do_newton(1e-10, error);
// phi 现在包含解
```

### 自动管理的对象

以下对象由 Kadath 内部管理，**不需要**手动 `delete`：

1. **Domain 对象**：由 `Space` 管理，在 `Space` 析构时自动释放
2. **Val_domain 对象**：由 `Scalar`/`Tensor` 管理
3. **Term_eq 对象**：由 `System_of_eqs` 内部管理
4. **Equation 对象**：由 `System_of_eqs` 管理

### 反直觉的内存模型

#### 1. 配置空间与系数空间的延迟计算

`Val_domain` 对象在配置空间和系数空间之间自动转换，但转换是**延迟的**：

```cpp
Val_domain val(domain);
val.set(Index(1,1,1)) = 1.0;  // 在配置空间设置值
val.coef();                    // 显式转换到系数空间
// 此时配置空间的值可能被清除
val.coef_i();                  // 转换回配置空间
```

**注意事项**：
- 修改配置空间的值会清除系数空间
- 修改系数空间的值会清除配置空间
- 需要同时访问两者时，要小心操作顺序

#### 2. 引用语义 vs 值语义

`Tensor` 和 `Scalar` 对象使用**引用语义**：

```cpp
Scalar phi1(space);
Scalar phi2 = phi1;  // 这是浅拷贝！phi1 和 phi2 共享数据

// 如果需要深拷贝
Scalar phi3(phi1, true);  // 第二个参数 true 表示深拷贝
```

#### 3. 系统对变量的引用

`System_of_eqs::add_var()` 存储的是**引用**，不是拷贝：

```cpp
Scalar phi(space);
syst.add_var("phi", phi);

// phi 必须在 syst 的整个生命周期内有效
// 不能提前销毁 phi
```

#### 4. 域的坐标缓存

`Domain` 的坐标（如 `cart`、`radius`）是**延迟计算**和**缓存**的：

```cpp
const Domain* dom = space.get_domain(0);
const Val_domain& x = dom->get_cart(1);  // 第一次调用时计算并缓存
const Val_domain& y = dom->get_cart(1);  // 直接返回缓存值
```

---

## 实战求解流程

### 标准求解流程

#### 步骤 1：定义空间

```cpp
#include "Kadath/include/kadath.hpp"
using namespace Kadath;

// 1. 定义中心点
Point center(3);
center.set(1) = 0.0;
center.set(2) = 0.0;
center.set(3) = 0.0;

// 2. 定义配置点数
Dim_array res(3, 17, 9, 9);  // 3维，径向17点，θ方向9点，φ方向9点

// 3. 定义域边界
Array<double> bounds(2);
bounds.set(0) = 1.0;  // 第一个域半径
bounds.set(1) = 2.0;  // 第二个域半径

// 4. 创建空间
Space_spheric space(CHEB_TYPE, center, res, bounds, true);
```

#### 步骤 2：创建未知量

```cpp
// 创建未知标量场
Scalar phi(space);
phi = 1.0;  // 设置初始猜测值

// 创建未知标量（如果需要）
double omega = 0.5;
```

#### 步骤 3：定义系统

```cpp
// 创建求解系统
System_of_eqs syst(space);

// 添加未知量
syst.add_var("phi", phi);
syst.add_var("omega", omega);

// 添加常数（如果需要）
syst.add_cst("pi", 3.14159);
```

#### 步骤 4：添加方程

```cpp
int n_dom = space.get_nbr_domains();

// 在所有域添加域内方程
for (int d=0; d<n_dom; d++) {
    // 例如：拉普拉斯方程 lap(phi) = -omega^2 * phi
    syst.add_eq_inside(d, "lap(phi) + omega^2 * phi = 0");
}

// 在原点添加边界条件
syst.add_eq_bc(0, INNER_BC, "phi = 1");

// 在无穷远添加边界条件（如果使用紧致化域）
int last_dom = n_dom - 1;
syst.add_eq_bc(last_dom, OUTER_BC, "phi = 0");

// 在域间添加匹配条件
for (int d=0; d<n_dom-1; d++) {
    // 场值连续
    syst.add_eq_matching(d, OUTER_BC, "phi = phi");
    // 法向导数连续
    syst.add_eq_matching(d, OUTER_BC, "dr(phi) = dr(phi)");
}
```

#### 步骤 5：求解

```cpp
double precision = 1e-10;
double error = 1.0;
int max_iter = 100;

std::cout << "开始牛顿迭代..." << std::endl;

for (int iter=0; iter<max_iter; iter++) {
    bool converged = syst.do_newton(precision, error);
    
    std::cout << "迭代 " << iter << ": 误差 = " << error << std::endl;
    
    if (converged) {
        std::cout << "收敛达到精度要求！" << std::endl;
        break;
    }
    
    if (error > 1e10) {
        std::cerr << "误差过大，可能发散！" << std::endl;
        break;
    }
}

std::cout << "求解完成。最终误差 = " << error << std::endl;
```

#### 步骤 6：提取结果

```cpp
// phi 现在包含解
// 可以访问各个域的值
for (int d=0; d<space.get_nbr_domains(); d++) {
    const Val_domain& val = phi(d);
    // 使用 val 进行后续处理
}

// omega 也被更新
std::cout << "omega = " << omega << std::endl;
```

### 完整示例：求解球对称拉普拉斯方程

```cpp
#include "Kadath/include/kadath.hpp"
#include <iostream>
using namespace Kadath;

int main() {
    // 步骤 1：定义空间
    Point center(3);
    center.set(1) = 0.0;
    center.set(2) = 0.0;
    center.set(3) = 0.0;
    
    Dim_array res(3, 17, 9, 9);
    Array<double> bounds(1);
    bounds.set(0) = 1.0;
    
    Space_spheric space(CHEB_TYPE, center, res, bounds, true);
    
    // 步骤 2：创建未知量
    Scalar phi(space);
    phi = 1.0;  // 初始猜测
    
    // 步骤 3：定义系统
    System_of_eqs syst(space);
    syst.add_var("phi", phi);
    
    // 步骤 4：添加方程
    // 域内方程：拉普拉斯方程
    syst.add_eq_inside(0, "lap(phi) = 0");
    
    // 边界条件：原点 phi = 1
    syst.add_eq_bc(0, INNER_BC, "phi = 1");
    
    // 边界条件：无穷远 phi = 0
    syst.add_eq_bc(1, OUTER_BC, "phi = 0");
    
    // 匹配条件
    syst.add_eq_matching(0, OUTER_BC, "phi = phi");
    syst.add_eq_matching(0, OUTER_BC, "dr(phi) = dr(phi)");
    
    // 步骤 5：求解
    double precision = 1e-10;
    double error = 1.0;
    
    for (int iter=0; iter<50; iter++) {
        bool converged = syst.do_newton(precision, error);
        std::cout << "迭代 " << iter << ": 误差 = " << error << std::endl;
        if (converged) break;
    }
    
    std::cout << "求解完成！" << std::endl;
    
    return 0;
}
```

---

## 常见问题与注意事项

### 1. 边界条件常量

Kadath 定义了以下边界常量：

```cpp
#define OUTER_BC 1      // 外边界
#define INNER_BC 2      // 内边界
#define CHI_ONE_BC 3    // χ=1 边界
#define ETA_PLUS_BC 4   // η=+1 边界
#define ETA_MINUS_BC 5  // η=-1 边界
#define TIME_INIT 6     // 时间初始边界
```

### 2. 基类型常量

```cpp
#define CHEB_TYPE 1     // Chebyshev 多项式
#define LEG_TYPE 2      // Legendre 多项式
```

### 3. 方程字符串中的运算符

- `+`, `-`, `*`, `/`：基本运算
- `^`：幂运算（如 `phi^2`）
- `lap()`：拉普拉斯算子
- `dr()`：径向导数
- `dt()`：θ 方向导数
- `dp()`：φ 方向导数
- `partial_i()`：对第 i 个坐标的偏导数

### 4. 常见错误

#### 错误 1：变量生命周期问题

```cpp
// 错误！
{
    Scalar phi(space);
    syst.add_var("phi", phi);
}  // phi 在这里被销毁，但 syst 仍在使用它

// 正确：确保 phi 在 syst 的整个生命周期内有效
Scalar phi(space);
syst.add_var("phi", phi);
// ... 使用 syst ...
```

#### 错误 2：忘记设置初始猜测

```cpp
// 错误：未设置初始值可能导致求解失败
Scalar phi(space);
syst.add_var("phi", phi);

// 正确：设置合理的初始猜测
Scalar phi(space);
phi = 1.0;  // 或使用其他初始值
syst.add_var("phi", phi);
```

#### 错误 3：边界条件不完整

```cpp
// 错误：缺少边界条件
syst.add_eq_inside(0, "lap(phi) = 0");
// 缺少边界条件，系统欠定

// 正确：添加足够的边界条件
syst.add_eq_inside(0, "lap(phi) = 0");
syst.add_eq_bc(0, INNER_BC, "phi = 1");
syst.add_eq_bc(0, OUTER_BC, "phi = 0");
```

### 5. 性能优化建议

1. **合理选择配置点数**：太多会增加计算量，太少会降低精度
2. **使用合适的基类型**：Chebyshev 通常更快，Legendre 在某些情况下更精确
3. **初始猜测要合理**：好的初始猜测可以显著减少迭代次数
4. **域的数量要适中**：太多域会增加匹配条件的数量

### 6. 调试技巧

1. **检查方程数量**：
```cpp
int n_unknowns = syst.get_nbr_unknowns();
int n_conditions = syst.get_nbr_conditions();
std::cout << "未知量数: " << n_unknowns << std::endl;
std::cout << "方程数: " << n_conditions << std::endl;
// 应该满足 n_unknowns == n_conditions
```

2. **检查残差**：
```cpp
Array<double> errors = syst.check_equations();
std::cout << "各方程残差: " << errors << std::endl;
```

3. **逐步添加方程**：先添加简单的方程，确认系统工作正常，再逐步添加复杂方程

---

## 总结

Kadath 的核心架构围绕三个主要类展开：

1. **Space**：定义计算区域的几何结构
2. **Domain**：处理单个域的坐标映射和谱展开
3. **System_of_eqs**：管理未知量、方程和求解过程

理解这些类的构造函数参数、内存管理特性和标准求解流程，是有效使用 Kadath 的关键。在实际使用中，要注意变量的生命周期、边界条件的完整性，以及合理的初始猜测设置。

---

*本指南基于 Kadath 库的头文件分析编写，如有疑问请参考官方文档或源代码。*
