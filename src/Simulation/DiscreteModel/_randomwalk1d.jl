# """
#     RandomWalk1D(step_num, step_size, prob)

# 1维随机游走
# """
# struct RandomWalk1D
#     step_num::Integer
#     start_point::Real
#     step_size::Real
#     prob::Real
#     function RandomWalk1D(step_num::Integer, start_point::Real, step_size::Real, prob::Real)
#         step_num > 0 || throw(DomainError("需为正整数"))
#         step_size > 1 || throw(DomainError("需为正数"))
#         0 <= prob <= 1 || throw(DomainError("概率需在0与1之间"))
#         return new(step_num, start_point, step_size, prob)
#     end
# end

# """
# 外置构造函数
# """
# RandomWalk1D(step_num::Integer; start_point::Real=0, step_size::Real=1, prob::Real=0.5) = RandomWalk1D(step_num, start_point, step_size, prob)

"""
输出函数
"""
Base.show(io::IO, ::MIME"text/plain", model::RandomWalk1D) = print(io, "一维随机游走, 模拟步数 $(model.step_num), 起点 $(model.start_point), 步长 $(model.step_size), (左, 右)转概率 $(model.prob).")

"""
路径模拟
"""
function simulate(model::RandomWalk1D)
    x = Vector{promote_type(typeof(model.step_size), typeof(model.start_point))}(undef, model.step_num+1)
    x[1] = model.start_point
    for k in 2:model.step_num+1
        if rand() < model.prob[1]
            direction = -1
        else
            direction = 1
        end
        @inbounds x[k] = x[k-1] + direction * model.step_size
    end
    return x
end
