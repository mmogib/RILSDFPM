function RILSDFPM(problem::TestProblem, paramsv::Vector{RILSDFPMParameters}, options::AlgorithmOptions)
    map(param -> RILSDFPM(problem; options=options, params=param), paramsv)
end

function RILSDFPM(problem::TestProblem; options::AlgorithmOptions, params::RILSDFPMParameters)

    result, elapsed_time = @timed RILSDFPM(problem.f, problem.x0, params, options)
    if isa(result, SuccessfullResult)
        NumericalResult(problem, length(problem.x0), options, params, SuccessfullResult(result, elapsed_time))
    else
        NumericalResult(problem, length(problem.x0), options, params, result)
    end

end

# Relaxed Inertial Liu Storey derivative-free projection method with restart procedure (RILSDFPM)
function RILSDFPM(F::Function, x0::Vector{Float64}, params::RILSDFPMParameters, options::AlgorithmOptions)
    maxiters, ϵ, StopIt = options.maxiters, options.tol, options.stopping
    α, γ, ρ = params.α, params.γ, params.ρ
    function_evaluations = 0
    Fx0 = F(x0)
    function_evaluations += 1
    Fx0Norm = norm(Fx0)
    if (Fx0Norm <= ϵ)
        return SuccessfullResult(1, function_evaluations, x0, Fx0Norm, :solved_x0)
    end
    x1 = v1 = x0
    Fx1 = F(x1)
    d1 = -Fx1
    function_evaluations += 1

    t1, extra_evals = computeRILSDFPMtk(F, v1, d1, options, params)
    function_evaluations += extra_evals
    z1 = v1 + t1 * d1
    Fz1 = F(z1)
    function_evaluations += 1
    Fz1Norm = norm(Fz1)
    if (Fz1Norm <= ϵ)
        return SuccessfullResult(1, function_evaluations, z1, Fz1Norm, :solved_zk)
    end
    λ1 = dot(Fz1, v1 - z1) / Fz1Norm^2
    x0, x1 = x1, (1 - ρ) * v1 + ρ * (v1 - γ * λ1 * Fz1)
    for k in 1:maxiters
        v0, v1 = v1, x1 + α * (x1 - x0)
        Fv0, Fv1 = F(v0), F(v1)
        function_evaluations += 2
        Fv1Norm = norm(Fv1)
        if (Fv1Norm <= ϵ)
            return SuccessfullResult(k, function_evaluations, v1, Fv1Norm, :solved_vk)
        end
        p1 = d1#Fv1 - Fv0 # d1 ,vk−1
        d0 = d1
        d1 = computeRILSDFPMdk(Fv0, Fv1, p1, d0, params)
        if 1000 * norm(d1) <= ϵ
            return SuccessfullResult(k, function_evaluations, v1, Fv1Norm, :vanishing_direction_dk)
        end
        t1, extra_evals = computeRILSDFPMtk(F, v1, d1, options, params)
        function_evaluations += extra_evals
        z1 = v1 + t1 * d1
        Fz1 = F(z1)
        function_evaluations += 1
        Fz1Norm = norm(Fz1)
        if (Fz1Norm <= ϵ)
            return SuccessfullResult(k, function_evaluations, z1, Fz1Norm, :solved_zk)
        end
        λ1 = dot(Fz1, v1 - z1) / Fz1Norm^2
        x0, x1 = x1, (1 - ρ) * v1 + ρ * (v1 - γ * λ1 * Fz1)
        Fx1 = F(x1)
        function_evaluations += 1
        Fx1Norm = norm(Fx1)
        if (Fx1Norm <= ϵ)
            return SuccessfullResult(k, function_evaluations, x1, Fx1Norm, :solved_xk)
        end
    end
    return SuccessfullResult(maxiters, function_evaluations + 1, x1, norm(F(x1)), :maxiter_reached)
end

function computeRILSDFPMtk(F, vk, dk, options, params)
    bt_maxiters = options.bt_maxiters
    ζ, σ, β, ρ = params.ζ, params.σ, params.β, params.ρ
    nrmdk2 = σ * norm(dk)^2
    for i in 1:bt_maxiters
        tk = ζ * β^i
        Fvk = F(vk + tk * dk)
        lhs = dot(-Fvk, dk)
        rhs = tk * norm(Fvk) * nrmdk2
        if (lhs >= rhs)
            return tk, i
        end
    end
    return ζ * ρ^bt_maxiters, bt_maxiters
end

function computeRILSDFPMdk(Fv0, Fv1, p1, d0, params)
    ξ, μ = params.ξ, params.μ
    y0 = Fv1 - Fv0
    βk = abs(-dot(Fv1, y0) / dot(d0, Fv0))
    check = βk * dot(Fv1, d0) < μ * norm(Fv1)^2
    dk = -Fv1 + (ξ * (dot(Fv1, p1) / (norm(p1)^2)) * p1)
    if (check)
        return dk + βk * d0
    else
        return dk
    end
end
