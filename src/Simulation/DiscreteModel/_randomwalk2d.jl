# """
#     RandomWalk2D(step_num, start_point, step_size, prob)

# 二维随机游走
# """
# struct RandomWalk2D
#     step_num::Integer
#     start_point::NTuple{2, Real}
#     step_size::Real
#     prob::NTuple{4, Real}
#     function RandomWalk2D(step_num::Integer, start_point::NTuple{2, Real}, step_size::Real, prob::NTuple{4, Real})
#         step_num > 0 || throw(DomainError("需为正整数"))
#         step_size > 0 || throw(DomainError("需为正数"))
#         all(p->0<=p<=1, prob) || throw(DomainError("概率需在 0 与 1 之间"))
#         sum(prob) == 1 || throw(DomainError("概率之和需为 1"))
#         return new(step_num, start_point, step_size, prob)
#     end
# end

# """
# 外置构造函数
# """
# RandomWalk2D(step_num::Integer; start_point::NTuple{2, Real}=(0, 0), step_size::Real=1, prob::NTuple{4, Real}=(0.25, 0.25, 0.25, 0.25)) = RandomWalk2D(step_num, start_point, step_size, prob)

"""
输出函数
"""
Base.show(io::IO, ::MIME"text/plain", model::RandomWalk2D) = print(io, "二维随机游走, 模拟步数 $(model.step_num), 起点 $(model.start_point), 步长 $(model.step_size), (上, 下, 左, 右)转概率 $(model.prob).")

"""
路径模拟
"""
function simulate(model::RandomWalk2D)
    x = Vector{promote_type(typeof(model.step_size), typeof(model.start_point[1]))}(undef, model.step_num+1)
    y = Vector{promote_type(typeof(model.step_size), typeof(model.start_point[2]))}(undef, model.step_num+1)
    x[1] = model.start_point[1]
    y[1] = model.start_point[2]
    prob  = cumsum(model.prob)
    @inbounds for k in 2:model.step_num+1
        rnd = rand()
        if rnd < prob[1]
            x[k] = x[k-1]
            y[k] = y[k-1] + model.step_size
        elseif model.prob[1] <= rnd < prob[2]
            x[k] = x[k-1]
            y[k] = y[k] - model.step_size
        elseif prob[2] <= rnd < prob[3]
            x[k] = x[k-1] - model.step_size
            y[k] = y[k-1]
        else
            x[k] = x[k-1] - model.step_size
            y[k] = y[k-1]
        end
    end
    return [x y]
end