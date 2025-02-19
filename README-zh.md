# DiffusionX.jl
[English](README.md) | 简体中文

> [!NOTE]
> 开发进行中。DiffusionX.jl 是一个基于 Julia 的随机过程模拟库。Rust/Python 版本也在同步开发中，可见 [DiffusionX](https://github.com/tangxiangong/diffusionx)。

## 使用示例
```julia
using DiffusionX

# 模拟标准布朗运动路径
B = Bm()              
Bₜ = B(10)            
t, x = simulate(Bₜ; τ=0.01)

# 计算阶矩
T = collect(10:10:100)
@. 𝔼(B(T))
@. 𝔼((B(T))^2)

# 模拟首次通过时间
fpt = FPT((-1, 1), B)
simulate(fpt) 
𝔼(fpt)

# 计算时间平均均方位移.
𝔼(δ̄²(B, T=100, Δ=0.1))
```

## 可扩展性

DiffusionX.jl 的设计充分利用了 Julia 的类型系统和多重分派特性，使得扩展新的随机过程变得简单直观：

1. **自定义过程**: 通过继承 `StochasticProcess` 抽象类型，可以轻松定义新的随机过程：
```julia
struct MyProcess <: StochasticProcess
    # 定义过程参数
    x₀::Float64  # 例如：初始位置
end
```

2. **实现模拟方法**: 只需要为新过程实现 `simulate` 方法：
```julia
function simulate(p::MyProcess, T::Union{Int,Float64}; τ::Float64=0.01)
    # 实现具体的模拟逻辑
    # 返回时间序列 t 和对应的轨迹 x
    t, x
end
```

3. **泛函计算**: 可以直接使用内置的泛函计算工具，如：
```julia
using DiffusionX
# 首次通过时间
fpt = FPT((-1, 1), process)  # 计算到达区间 [-1,1] 边界的时间
simulate(fpt)
# 占据时间
ot = OccupationTime(T, (a, b), process)  # 计算在区间 [a,b] 中停留的时间
simulate(ot)

# 时间平均
ta = δ̄²(process, T, Δ)  # 计算时间平均均方位移
```

4. **统计分析**: 所有过程自动支持矩的计算和期望分析：
```julia
# 计算期望
𝔼(process(T))  # T 时刻的期望
𝔼(process(T)^2)  # T 时刻的二阶矩
𝔼(fpt)  # 首次通过时间的期望
𝔼(ot)  # 占据时间的期望
𝔼(δ̄²(process, 100, 0.1))  # 时间平均均方位移
```

所有扩展的功能都能无缝集成到现有框架中，并自动支持轨迹模拟、泛函计算和统计分析等核心功能。

## 许可证

本项目采用 [MIT 许可证](https://opensource.org/licenses/MIT)。

