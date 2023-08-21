
"""
0)
Our Algorithm 

"""
function BBSTTDFP(problem::TestProblem, paramsv::Vector{BBSTTDFPParams};
    options::AlgorithmOptions=AlgorithmOptions(200, 1e-5)
)
    map(param -> BBSTTDFP(problem; options=options, params=param), paramsv)
end

function BBSTTDFP(problem::TestProblem;
    options::AlgorithmOptions=AlgorithmOptions(200, 1e-5),
    params::BBSTTDFPParams=BBSTTDFPParams(
        0.5, 0.5, 0.5, 1.5, 0.5, 0.8, 0.5, 0.5
    ),
    inertial::Bool=true
)

    result, elapsed_time = @timed BBSTTDFP(problem.f, problem.P, problem.Pcheck, problem.x0, options, params; inertial=inertial)
    if isa(result, SuccessfullResult)
        NumericalResult(problem, length(problem.x0), options, params, SuccessfullResult(result, elapsed_time))
    else
        NumericalResult(problem, length(problem.x0), options, params, result)
    end

end
function BBSTTDFP(F::Function, P::Function, Pcheck::Function, x1::Vector{Float64}, options::AlgorithmOptions,
    params::BBSTTDFPParams; inertial::Bool=true)
    if inertial
        BBSTTDFP(F, P, Pcheck, x1, options, params)
    else
        BBSTTDFPNOIN(F, P, Pcheck, x1, options, params)
    end
end
function BBSTTDFP(F::Function, P::Function, Pcheck::Function, x1::Vector{Float64}, options::AlgorithmOptions,
    params::BBSTTDFPParams)
    t, β, σ, γ, η, ξ = params.t, params.β, params.σ, params.γ, params.η, params.ξ
    αmin, αmax, r, ψ = params.αmin, params.αmax, params.r, params.ψ
    ϵ, maxiters, stop_criterion = options.tol, options.maxiters, options.stopping
    evals = 0
    u1 = x0 = x1
    Fu1 = F(u1)
    d1 = -Fu1
    for k in 1:maxiters
        τ1, tevals = findtk(F, d1, u1, σ, β, η, ξ)
        evals = evals + tevals
        if isnothing(τ1)
            return FailedResult("Failed to find tₖ at iteration=$k, σ = $σ, and β = $β")
        end

        z1 = u1 + τ1 * d1
        Fz1 = F(z1)
        evals = evals + 1
        Fz1_norm = norm(Fz1)
        if Pcheck(z1) && Fz1_norm <= ϵ
            return SuccessfullResult(
                k, evals, z1, Fz1_norm
            )
        end
        μ1 = dot(Fz1, u1 - z1) / (Fz1_norm^2)
        x1, x0, u0, Fu0, d0 = P(u1 - γ * μ1 * Fz1), x1, u1, Fu1, d1
        x1_x0 = x1 - x0
        x1_x0_norm = norm(x1_x0)
        t1 = if 10 * x1_x0_norm <= ϵ # this means x0 and x1 are almost the same
            t
        else
            min(t, 1 / (k^2 * x1_x0_norm))
        end
        u1 = x1 + t1 * (x1_x0)
        Fu1 = F(u1)
        evals = evals + 1
        Fu1_norm = norm(Fu1)
        if Fu1_norm <= ϵ
            return SuccessfullResult(
                k, evals, u1, Fu1_norm
            )
        end
        d1 = findDk(Fu0, Fu1, d0, u0, u1, r, ψ, αmin, αmax)
        if 10 * norm(d1) <= ϵ
            return SuccessfullResult(
                k, evals, u1, Fu1_norm
            )
        end
        if stop_criterion(Fu0, Fu1)
            return SuccessfullResult(
                k, evals, u1, Fu1_norm
            )
        end
        # x1, x0, u0, Fu0, d0 = if Fz1_norm <= ϵ
        #     z1, x1, u1, Fu1, d1
        # else
        #     μ1 = dot(Fz1, u1 - z1) / dot(Fz1, Fz1)
        #     P(u1 - γ * μ1 * Fz1), x1, u1, Fu1, d1
        # end
    end
    SuccessfullResult(maxiters, evals, u1, norm(Fu1), :maxiter_reached)
end
function BBSTTDFPNOIN(F::Function, P::Function, Pcheck::Function, x1::Vector{Float64}, options::AlgorithmOptions,
    params::BBSTTDFPParams
)
    t, β, σ, γ, η, ξ = params.t, params.β, params.σ, params.γ, params.η, params.ξ
    αmin, αmax, r, ψ = params.αmin, params.αmax, params.r, params.ψ
    ϵ, maxiters, stop_criterion = options.tol, options.maxiters, options.stopping
    evals = 0
    x0 = x1
    Fx1 = F(x1)
    evals += 1
    d1 = -Fx1
    for k in 1:maxiters
        τ1, tevals = findtk(F, d1, x1, σ, β, η, ξ)
        evals = evals + tevals
        if isnothing(τ1)
            return FailedResult("Failed to find tₖ at iteration=$k, σ = $σ, and β = $β")
        end

        z1 = x1 + τ1 * d1
        Fz1 = F(z1)
        evals = evals + 1
        Fz1_norm = norm(Fz1)
        if Pcheck(z1) && Fz1_norm <= ϵ
            return SuccessfullResult(
                k, evals, z1, Fz1_norm
            )
        end
        μ1 = dot(Fz1, x1 - z1) / (Fz1_norm^2)
        x1, x0, Fx0 = P(x1 - γ * μ1 * Fz1), x1, Fx1

        Fx1 = F(x1)
        evals += 1
        Fx1_norm = norm(Fx1)
        d1 = findDk(Fx0, Fx1, d1, x0, x1, r, ψ, αmin, αmax)
        if 10 * norm(d1) <= ϵ
            return SuccessfullResult(
                k, evals, x1, Fx1_norm
            )
        end
        if stop_criterion(Fx0, Fx1)
            return SuccessfullResult(
                k, evals, x1, Fx1_norm
            )
        end

    end
    SuccessfullResult(maxiters, evals, x1, norm(Fx1), :maxiter_reached)
end
function findDk(Fu0, Fu1, d0, u0, u1, r, ψ, αmin, αmax)
    y0 = Fu1 - Fu0
    s0 = u1 - u0 + r * y0
    v1 = max(ψ * norm(d0) * norm(y0), dot(d0, y0), norm(Fu0)^2)
    β1 = dot(Fu1, y0) / v1
    α12 = dot(Fu1, d0) / v1
    α11 = min(αmax, max(αmin, dot(s0, y0) / dot(y0, y0)))
    return -α11 * Fu1 + β1 * d0 - α12 * y0
end
function findtk(F::Function, d::Vector{Float64}, u::Vector{Float64}, σ::Float64, β::Float64, _η::Nothing, _ξ::Nothing)
    # ϵ = 1e-6
    max_iters = 10_000
    for i in 0:max_iters
        tk = β^i
        arg = F(u + tk .* d)
        lhs = -dot(arg, d)
        rhs = σ * tk * norm(arg) * norm(d)^2
        # rhs = σ * tk * norm(d)^2
        if lhs >= rhs || tk <= 1e-5
            return tk, i
        end
    end
    nothing, max_iters
end

function findtk(F::Function, d::Vector{Float64}, u::Vector{Float64}, σ::Float64, β::Float64, η::Float64, ξ::Float64)
    # ϵ = 1e-6
    max_iters = 10_000
    for i in 0:max_iters
        tk = β^i
        arg = F(u + tk .* d)
        lhs = -dot(arg, d)
        χ = norm(arg)
        Pηξχ = min(ξ, max(χ, η))
        rhs = σ * tk * Pηξχ * norm(d)^2
        # rhs = σ * tk * norm(d)^2
        if lhs >= rhs || tk <= 1e-5
            return tk, i
        end
    end
    nothing, max_iters
end

