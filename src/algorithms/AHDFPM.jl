"""
3)
Algorithm 1: AHDFPM3
Liu P, Wu X, Shao H, Zhang Y, Cao S. Three adaptive hybrid derivative-free
projection methods for constrained monotone nonlinear equations and their applications. Numer Linear Algebra
Appl. 2023;30(2):e2471. https://doi.org/10.1002/nla.2471
"""
function AHDFPM(problem::TestProblem;
    options::AlgorithmOptions=AlgorithmOptions(),
    params::AHDFPMParameters=AHDFPMParameters()
)

    result, elapsed_time = @timed AHDFPM(problem.f, problem.P,
        problem.x0, options, params)
    if isa(result, SuccessfullResult)
        NumericalResult(problem, length(problem.x0), options, params, SuccessfullResult(result, elapsed_time))
    else
        NumericalResult(problem, length(problem.x0), options, params, result)
    end

end
function AHDFPM(F::Function, P::Function, x0::Vector{Float64},
    options::AlgorithmOptions, params::AHDFPMParameters)
    ρ, γ, σ, μ, ν, η1, η2, η3, ξ = params.ρ, params.γ, params.σ, params.μ, params.ν, params.η1, params.η2, params.η3, params.ξ
    ϵ, maxiters = options.tol, options.maxiters
    evals = 0
    x0 = P(x0)
    Fx0 = F(x0)
    evals += 1
    normF0 = norm(Fx0)
    if normF0 <= ϵ
        return SuccessfullResult(1, evals, x0, normF0)
    end
    d0 = -Fx0
    for i in 1:maxiters
        α0, nevals = findAHDFPMαk(F, x0, d0, ξ, ρ, σ, μ, ν)
        if isnothing(α0)
            return FailedResult("Failed to find αk at iteration=$i")
        end
        evals += nevals
        z0 = x0 + α0 * d0
        Fz0 = F(z0)
        evals += 1
        normFz0 = norm(Fz0)
        ξ0 = dot(Fz0, x0 - z0) / normFz0^2
        x1 = P(x0 - γ * ξ0 * Fz0)
        Fx1 = F(x1)
        evals += 1
        d1 = findAHDFPMdk(Fx0, Fx1, d0, η1, η2, η3)
        normF1 = norm(Fx1)
        if 10 * norm(d1) <= ϵ || normF1 <= ϵ
            return SuccessfullResult(i, evals, x1, normF1)
        end
        x0, d0, Fx0 = x1, d1, Fx1
    end
    SuccessfullResult(maxiters, evals, x0, norm(F(x0)), :maxiter_reached)
end
function findAHDFPMdk(Fx0, Fx1, d0, η1, η2, η3)
    y0 = Fx1 - Fx0
    normF0 = norm(Fx0)
    normd0 = norm(d0)
    ω0 = max(dot(d0, y0), normF0^2, normd0^2)
    dotFkd0 = dot(Fx1, d0)
    c0 = dotFkd0 / ω0
    part1 = -η1 * Fx1
    part2 = (1 / ω0) * (η1 * dot(Fx1, y0) - η2^2 * c0 * norm(y0)^2 - η3 * dotFkd0) * d0
    part3 = (2 * η2 - η1) * c0 * y0
    dk = part1 + part2 + part3
    dk
end
function LSPAHDFPM(μ, ν, χ)
    return if χ <= μ
        μ
    elseif χ >= ν
        ν
    else
        χ
    end
end
function findAHDFPMαk(F, x0, dk, ξ, ρ, σ, μ, ν)
    maxitrs = 10_000
    for i in 0:maxitrs
        αk = ξ * ρ^i
        Fvalue = F(x0 + αk * dk)
        lhs = -dot(Fvalue, dk)
        rhs = σ * αk * LSPAHDFPM(μ, ν, norm(Fvalue)) * (norm(dk))^2
        if lhs >= rhs
            return αk, i
        end
    end
    return nothing, maxitrs
end