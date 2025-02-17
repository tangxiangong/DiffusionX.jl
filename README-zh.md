# DiffusionX.jl
[English](README.md) | 简体中文

> [!NOTE]
> 开发进行中。DiffusionX.jl 是一个基于 Julia 的随机过程模拟库。Rust/Python 版本也在同步开发中，可见 [DiffusionX](https://github.com/tangxiangong/diffusionx)。

## 使用示例
```julia
using DiffusionX.Simulation

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
𝔼(δ̄²(B; T=100, Δ=0.1))
```

## 许可证

本项目采用 [MIT 许可证](https://opensource.org/licenses/MIT)。

