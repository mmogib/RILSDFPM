"""
2)
Algorithm 2.1 (CGDFP) in 
Li Zheng, Lei Yang, Yong Liang,
A conjugate gradient projection method for solving equations with convex constraints,
Journal of Computational and Applied Mathematics,
Volume 375,
2020,
112781,
ISSN 0377-0427,
https://doi.org/10.1016/j.cam.2020.112781.

"""
function CGDFP(problem::TestProblem;
    options::AlgorithmOptions=AlgorithmOptions(),
    params::CGDFPParameters=CGDFPParameters()
)

    result, elapsed_time = @timed CGDFP(problem.f, problem.P,
        problem.x0, options, params)
    if isa(result, SuccessfullResult)
        NumericalResult(problem, length(problem.x0), options, params, SuccessfullResult(result, elapsed_time))
    else
        NumericalResult(problem, length(problem.x0), options, params, result)
    end

end
function CGDFP(F::Function, P::Function, x0::Vector{Float64},
    options::AlgorithmOptions, params::CGDFPParameters)
    ρ, β, σ, σ1, σ2 = params.ρ, params.β, params.σ, params.σ1, params.σ2
    ϵ, maxiters = options.tol, options.maxiters
    evals = 0
    Fx0 = F(x0)
    evals += 1
    normF0 = norm(Fx0)
    if normF0 <= ϵ
        return SuccessfullResult(1, evals, x0, normF0)
    end
    d0 = -Fx0
    for i in 1:maxiters
        αk, nevals = findCGDFPαk(F, x0, d0, β, ρ, σ)
        evals += nevals
        if isnothing(αk)
            return FailedResult("Failed to find αk at iteration=$i")
        end
        z0 = x0 + αk * d0
        Fz0 = F(z0)
        evals += 1
        normFz0 = norm(Fz0)
        τ0 = (dot(x0 - z0, Fz0) / normFz0^2)
        x1 = P(x0 - τ0 * Fz0)
        Fx1 = F(x1)
        evals += 1
        normF1 = norm(Fx1)
        if normF1 <= ϵ
            return SuccessfullResult(i, evals, x1, normF1)
        end
        d0 = findCGDFPDk(Fx0, Fx1, d0, (σ1, σ2))
        x0, Fx0 = x1, Fx1
    end
    SuccessfullResult(maxiters, evals, x0, norm(Fx0), :maxiter_reached)

end
function findCGDFPDk(F0::Vector{Float64}, F1::Vector{Float64},
    d0::Vector{Float64},
    params::Tuple{Float64,Float64})

    σ1, σ2 = params
    y0 = F1 - F0
    λk = 1 + max(0, -dot(y0, d0) / (norm(d0)^2))
    w0 = y0 + λk * d0
    prductd0w0 = dot(w0, d0)
    ck = (dot(F1, d0)) / (prductd0w0)
    βk = (dot(F1, w0) - 2 * ck * norm(w0)^2) / (prductd0w0)
    return -σ1 * F1 + βk * d0 + σ2 * ck * w0
end
function findCGDFPαk(F, x0, dk, β, ρ, σ)
    maxitrs = 10_000
    for i in 0:maxitrs
        αk = β * ρ^i
        Fvalue = F(x0 + αk * dk)
        lhs = -dot(Fvalue, dk)
        rhs = σ * αk * norm(Fvalue) * (norm(dk))^2
        if lhs >= rhs
            return αk, i
        end
    end
    return nothing, maxitrs
end