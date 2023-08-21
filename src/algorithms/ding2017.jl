"""
Algorithm 2.1 in 
Yanyun Ding, Yunhai Xiao & Jianwei Li (2017) A class of conjugate gradient
methods for convex constrained monotone equations, Optimization, 66:12, 2309-2328, DOI:
10.1080/02331934.2017.1372438
"""
function ding2017(problem::TestProblem;
    options::AlgorithmOptions=AlgorithmOptions(200, 1e-5),
    params::Ding2017Parameters=Ding2017Parameters(0.74, 1e-4, 0.1, 0.5, 1))
    F, P, x0 = problem.f, problem.P, problem.x0
    result, time_elapsed = @timed ding2017(F, P, x0, options, params)
    if isa(result, SuccessfullResult)
        return NumericalResult(problem, length(x0), options, params, SuccessfullResult(result, time_elapsed))
    else
        return NumericalResult(problem, length(x0), options, params, result)
    end
end
function ding2017(F::Function, P::Function, x0::Vector{Float64}, options::AlgorithmOptions, params::Ding2017Parameters)
    maxiters, tol = options.maxiters, options.tol
    ρ, σ, θ, η, ξ = params.ρ, params.σ, params.θ, params.η, params.ξ
    Fevals = 0
    F0 = F(x0)
    Fevals = Fevals + 1
    F0Norm = norm(F0, 2)
    if (F0Norm < tol)
        return SuccessfullResult(0, Fevals, x0, F0Norm)
    end
    d0 = -F0
    t0 = ding2017Findtk(F, x0, d0, ξ, ρ, σ)
    if (isnothing(t0))
        return FailedResult("Line search failed at initial step.")
    end
    z0 = x0 + t0 * d0
    Fz0 = F(z0)
    Fevals = Fevals + 1
    Fz0Norm = norm(Fz0, 2)
    if (Fz0Norm < tol)
        return SuccessfullResult(0, Fevals, z0, Fz0Norm)
    end
    α0 = dot(Fz0, x0 - z0) / Fz0Norm^2
    xk = P(x0 - α0 * Fz0)
    iterations = 0
    for k in 1:maxiters
        iterations = k
        Fk = F(xk)
        Fevals = Fevals + 1
        FkNorm = norm(Fk, 2)
        if (FkNorm < tol)
            return SuccessfullResult(iterations, Fevals, xk, FkNorm)
        end
        dk = ding2017Direction(Fk, F0, d0, t0, θ, η)
        tk = ding2017Findtk(F, xk, dk, ξ, ρ, σ)
        if (isnothing(tk))
            return FailedResult("Line search failed at iteration $k")
        end
        zk = xk + tk * dk
        Fzk = F(zk)
        Fevals = Fevals + 1
        FzkNorm = norm(Fzk, 2)
        if (FzkNorm < tol)
            return SuccessfullResult(iterations, Fevals, zk, FzkNorm)
        end
        αk = dot(Fzk, xk - zk) / FzkNorm^2
        x0 = xk
        F0 = F(x0)
        Fevals = Fevals + 1
        xk = P(xk - αk * Fzk)
    end
    # return SuccessfullResult(iterations, Fevals, xk, norm(F(xk), 2))
    return FailedResult("Maximum number ($maxiters) of iterations reached.")
end
function ding2017Direction(Fk, F0, d0, t0, θ, η)
    γ0 = Fk - F0
    s0 = t0 * d0
    F0Norm = norm(F0, 2)
    λ0 = 1 + (1 / F0Norm) * max(0, -dot(γ0, s0) / norm(s0, 2)^2)
    y0 = γ0 + λ0 * t0 * F0Norm * d0
    y0Norm = norm(y0, 2)
    s0Doty0 = dot(s0, y0)
    τ0A = y0Norm^2 / s0Doty0
    τ0B = s0Doty0 / norm(s0, 2)^2
    τ0 = θ * τ0A + (1 - θ) * τ0B
    βk = (dot(Fk, y0) / dot(d0, y0)) - (τ0 + τ0A - τ0B) * (dot(Fk, s0) / dot(d0, y0))
    βkplus = max(βk, η * dot(Fk, d0) / norm(d0, 2)^2)
    dk = -Fk + βkplus * d0
    dk
end

function ding2017Findtk(F, xk, dk, ξ, ρ, σ)
    tks = [ξ * ρ^i for i in 1:1000]
    for tk in tks
        lhs = -dot(F(xk + tk * dk), dk)
        rhs = σ * tk * norm(dk, 2)^2
        if lhs >= rhs
            return tk
        end
    end
    return nothing
end
