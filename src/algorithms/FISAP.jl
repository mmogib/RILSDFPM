
"""
Algorithm 2.1 (FISAP) in
N. Zhang, J.K. Liu, L.Q. Zhang, Z.L. Lu,
A fast inertial self-adaptive projection based algorithm for solving large-scale nonlinear monotone equations,
Journal of Computational and Applied Mathematics,
Volume 426,
2023,
115087,
ISSN 0377-0427,
https://doi.org/10.1016/j.cam.2023.115087.

"""
function FISAP(problem::TestProblem;
    options::AlgorithmOptions=AlgorithmOptions(),
    params::FISAPParameters=FISAPParameters()
)

    result, elapsed_time = @timed FISAP(problem.f, problem.x0, options, params)
    if isa(result, SuccessfullResult)
        NumericalResult(problem, length(problem.x0), options, params, SuccessfullResult(result, elapsed_time))
    else
        NumericalResult(problem, length(problem.x0), options, params, result)
    end

end
function FISAP(F::Function, x0::Vector{Float64}, options::AlgorithmOptions, params::FISAPParameters)
    σ, ρ, r, ν, τ, ρk = params.σ, params.ρ, params.r, params.ν, params.τ, params.ρk
    ϵ, maxiters = options.tol, options.maxiters
    evals = 0
    dkparams = nothing
    d0 = nothing
    Fx0 = F(x0)
    w1 = w0 = x1 = x0
    evals += 1
    for i in 1:maxiters
        normF0 = norm(Fx0)
        if normF0 <= ϵ
            return SuccessfullResult(i, evals, x0, normF0)
        end
        w0, w1 = if i == 1
            x0, x0
        else
            w1, x0 + ρk * (x1 - x0)
        end
        Fw0 = F(w0)
        evals += 1
        normFw0 = norm(Fw0)
        if normFw0 <= ϵ
            return SuccessfullResult(i, evals, w0, normFw0)
        end
        Fw1 = F(w1)
        evals += 1
        d1 = findFISAPDk(Fw0, Fw1, d0, dkparams)
        αk, nevals = findFISAPαk(F, w1, d1, ρ, σ)
        evals += nevals
        if isnothing(αk)
            return FailedResult("Failed to find αk at iteration=$i")
        end
        z1 = w1 + αk * d1
        Fz1 = F(z1)
        evals += 1
        normFz1 = norm(Fz1)
        if normFz1 <= ϵ
            return SuccessfullResult(i, evals, z1, normFz1)
        end
        μk = dot(Fz1, w1 - z1) / normFz1^2
        x0, x1, d0 = x1, w1 - ν * μk * Fz1, d1
        dkparams = (τ, r)
    end
    SuccessfullResult(maxiters, evals, x0, norm(Fx0), :maxiter_reached)
end
function findFISAPαk(F, wk, dk, ρ, σ)
    maxitrs = 10_000
    for i in 0:maxitrs
        αk = ρ^i
        Fvalue = F(wk + αk * dk)
        lhs = -dot(Fvalue, dk)
        rhs = σ * αk * norm(Fvalue) * (norm(dk))^2
        if lhs >= rhs
            return αk, i
        end
    end
    return nothing, maxitrs
end
function findFISAPDk(Fw0::Vector{Float64}, Fw1::Vector{Float64}, d0::Union{Nothing,Vector{Float64}}, params::Union{Nothing,Tuple{Float64,Float64}})
    if isnothing(params)
        return -Fw0
    end
    τ, r = params
    normd0 = norm(d0)
    normd02 = normd0^2
    y0 = Fw1 - Fw0
    normy0 = norm(y0)
    λ0 = 1 + max(0, -(dot(d0, y0)) / (normd02))
    h0 = y0 + λ0 * d0
    η1 = max(normd02, dot(d0, h0), -dot(d0, Fw0), τ * normy0 * normd0)
    θ1 = 1 + ((dot(Fw1, d0) * normy0)^2) / (r * (2 * norm(Fw1) * η1)^2)
    β1 = (dot(Fw1, y0)) / η1
    return -θ1 * Fw1 + β1 * d0
end
