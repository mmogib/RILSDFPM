""" 
1)
Algorithm 2.1 (MOPCG) in
Jamilu Sabi'u, Abdullah Shah, Predrag S. Stanimirović, Branislav Ivanov, Mohammed Yusuf Waziri, Modified optimal Perry conjugate gradient method 
for solving system of monotone equations with applications,  
Applied Numerical Mathematics, Volume 184,  2023, 
Pages 431-445, ISSN 0168-9274, 
https://doi.org/10.1016/j.apnum.2022.10.016

"""
function MOPCG(problem::TestProblem;
    options::AlgorithmOptions=AlgorithmOptions(),
    params::MOPCGParameters=MOPCGParameters()
)

    result, elapsed_time = @timed MOPCG(problem.f, problem.P, problem.Pcheck, problem.x0, options, params)
    if isa(result, SuccessfullResult)
        # npams = BBSTTDFPParams(params, result[2])
        NumericalResult(problem, length(problem.x0), options, params, SuccessfullResult(result, elapsed_time))
    else
        NumericalResult(problem, length(problem.x0), options, params, result)
    end

end
function MOPCG(F::Function, P::Function, Pckeck::Function, x0::Vector{Float64}, options::AlgorithmOptions, params::MOPCGParameters)
    ρ, η, ξ, λ = params.ρ, params.η, params.ξ, params.λ
    ϵ, maxiters = options.tol, options.maxiters
    evals = 0
    Fx0 = F(x0)
    evals += 1
    d0 = -Fx0
    for i in 1:maxiters
        normF0 = norm(Fx0)
        if normF0 <= ϵ
            return SuccessfullResult(i, evals, x0, normF0)
        end
        αk, nevals = findMOPCGαk(F, x0, d0, ξ, ρ, η)
        if isnothing(αk)
            return FailedResult("Failed to find αk at iteration=$i")
        end
        evals = evals + nevals
        μ0 = x0 + αk * d0
        Fμ0 = F(μ0)
        evals += 1
        normFμ0 = norm(Fμ0)
        if Pckeck(μ0) && normFμ0 <= ϵ
            return SuccessfullResult(i, evals, μ0, normFμ0)
        end
        b0 = (dot(Fμ0, x0 - μ0)) / (normFμ0^2)
        x1 = P(x0 - b0 * Fμ0)
        s0 = μ0 - x0
        Fx1 = F(x1)
        evals += 1
        y0 = Fx1 - Fx0 + λ * s0
        θ1 = (dot(s0, y0)) / (dot(s0, s0))
        β1 = (dot(y0 - θ1 * s0, Fx1)) / (dot(d0, y0))
        d0 = findDkMOPCG(Fx1, (β1, d0))
        x0, Fx0 = x1, Fx1
    end
    SuccessfullResult(maxiters, evals, x0, norm(Fx0), :maxiter_reached)
end
function findMOPCGαk(F, xk, dk, ξ, ρ, η)
    maxitrs = 10_000
    for i in 0:maxitrs
        αk = ρ * η^i
        lhs = -dot(F(xk + αk * dk), dk)
        rhs = ξ * αk * (norm(dk))^2
        if lhs >= rhs
            return αk, i
        end
    end
    return nothing, maxitrs
end
function findDkMOPCG(Fk, params::Tuple{Float64,Vector{Float64}})
    βk, dk = params
    return -(1 + βk * (dot(Fk, dk)) / (dot(Fk, Fk))) * Fk + βk * dk
end
