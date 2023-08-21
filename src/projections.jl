function projectOnBoxCheck(x::Vector{<:Real}; bounds::Union{Nothing,Tuple{<:Real,<:Real}}=nothing)
    if isnothing(bounds)
        true
    else
        l, u = bounds
        all(l .<= x .<= u)
    end
end
function projectOnBox(x::Vector{<:Real}; bounds::Union{Nothing,Tuple{<:Real,<:Real}}=nothing)
    if isnothing(bounds)
        x
    else
        l, u = bounds
        x .|> t -> min(t, u) .|> t -> max(t, l)
    end
end

function projectOnHalfSpaceCheck(x; β::Union{Nothing,Real}=nothing)
    n = length(x)
    b = isnothing(β) ? n : β
    sum(x) <= b
end

function projectOnHalfSpace(x; β::Union{Nothing,Real}=nothing)
    # % project an n-dimensional vector x to the the halfspace Sᵦ
    # % Sᵦ = { x : x ∈ Rⁿ, sum(x) <= β}
    # % Coded by
    n = length(x)
    b = isnothing(β) ? n : β
    if (sum(x) <= b)
        return x
    end

    bget = false
    s = sort(x, rev=true)
    tmpsum = 0
    for ii = 1:n-1
        tmpsum = tmpsum + s[ii]
        tmax = (tmpsum - b) / ii
        if tmax >= s[ii+1]
            bget = true
            break
        end
    end
    tmpsum = tmpsum + s[n]
    if ~bget
        tmax = (tmpsum - b) / n
    end

    P = if (tmpsum <= b)
        x
    else
        x .- tmax .|> r -> max(r, 0)
    end
    P
end
function projectOnTriangleCheck(x; β::Union{Nothing,Real}=nothing, lower=0)
    projectOnBoxCheck(x, bounds=(lower, Inf64)) && projectOnHalfSpaceCheck(x, β=β)
end
function projectOnTriangle(x; β::Union{Nothing,Real}=nothing, lower=0)
    # % project an n-dim vector x to the convex set Cn
    # % Cn = { x : x n-dim, x >= lower, sum(x) <= β}
    P = projectOnBox(x, bounds=(lower, Inf64)) |> p -> projectOnHalfSpace(p, β=β)
    P
end

function projectOnRn(x)
    x
end


# void condaty(double* y, double* x, const unsigned int length, const double a) {
# 	double*	aux = (x==y ? (double*)malloc(length*sizeof(double)) : x);
# function projectionOnSimplex(x::Vector{T}, a::T) where {T<:Real}
#     aux = similar(x)
#     aux0 = aux
#     tau = aux[1] = x[1] - a
#     auxlength = 1
#     auxlengthold = -1
#     for i in 2:length(y)
#         if x[i] > tau
#             if (tau += (aux[auxlength+1] = x[i]) - tau / (auxlength - auxlengthold)) <= x[i] - a
#                 tau = x[i] - a
#                 auxlengthold = auxlength
#             end
#             auxlength += 1
#         end
#     end

#     if auxlengthold >= 0
#         auxlength -= auxlengthold + 1
#         aux = view(aux, (auxlengthold+1):auxlength)
#         while auxlengthold >= 0
#             if aux0[auxlengthold+1] > tau
#                 tau += ((aux[length(aux)+1] = aux0[auxlengthold+1]) - tau) / (auxlength + 1)
#                 pop!(aux)
#             end
#             auxlengthold -= 1
#         end
#     end

#     while true
#         auxlengthold = length(aux) - 1
#         for i in 0:auxlengthold
#             if aux[i+1] > tau
#                 push!(aux, aux[i+1])
#             else
#                 tau += (tau - aux[i+1]) / (auxlengthold - i + length(aux))
#             end
#         end
#         if length(aux) > auxlengthold
#             break
#         end
#     end
#     println(tau)
#     P = similar(x)
#     for i in eachindex(P)
#         P[i] = (x[i] > tau ? x[i] - tau : zero(T))
#     end

#     P
# end

# function projectionOnSimplex1(x; codim::Union{Real,Nothing}=nothing, lower::Real=0)
#     n = length(x)
#     a = isnothing(codim) ? n : codim
#     vt = Vector{Union{Nothing,eltype(x)}}(nothing, n)
#     v = Vector{Union{Nothing,eltype(x)}}(nothing, n)
#     v[1] = x[1]
#     ρ = x[1] - a
#     vlength = 1
#     for i in 2:n
#         if x[i] > ρ
#             ρ = ρ + (x[i] - ρ) / (vlength + 1)
#             if ρ > (x[i] - a)
#                 vlength = vlength + 1
#                 v[vlength] = x[i]
#             else
#                 vt = [v...]
#                 v .= nothing
#                 v[1] = x[i]
#                 vlength = 1
#                 ρ = x[i] - a
#             end
#         end
#     end
#     if !all(isnothing.(vt))
#         for w in filter(r -> !isnothing(r), vt)
#             if w > ρ
#                 vlength = vlength + 1
#                 v[vlength] = w
#                 ρ = ρ + (y - ρ) / vlength
#             end
#         end
#     end
#     if !all(isnothing.(v))
#         for w in filter(r -> !isnothing(r), v)
#             if w <= ρ
#                 vlength = vlength - 1
#                 ρ = ρ + (w - ρ) / vlength
#             end
#         end
#     end
#     x .|> r -> max(r - ρ, lower)
# end

# function Pj(y, l, u)
#     # P=max(l,min.(y,u));
#     # P
#     y .|> t -> min(t, u) .|> t -> max(t, l)
# end

# function Pr(y, _, _)
#     y
# end

# function Pj2(x::Vector{<:Real}, l::Real)
#     Pj2(x, undef, undef)
# end

# function Pj2(x, _, _)
#     # % project an n-dim vector x to the convex set Cn
#     # % Cn = { x : x n-dim, x > -1, sum(x) <= n}
#     # % Coded by
#     # % Hassan Mohammad hmuhd.mth@buk.edu.ng
#     # % October, 2018
#     n = length(x)
#     bget = false
#     s = sort(x, rev=true)
#     tmpsum = 0

#     for ii in 1:n-1
#         tmpsum = tmpsum + s[ii]
#         tmax = (tmpsum - n) / ii
#         if tmax >= s[ii+1]
#             bget = true
#             break
#         end
#     end
#     tmpsum = tmpsum + s[n]
#     if ~bget
#         tmax = (tmpsum - n) / n
#     end

#     P = if (tmpsum <= n && min(x...) > -1)
#         x
#     else
#         x .- tmax .|> r -> max(r, 0)
#     end
#     P
# end


# function Pj2(x, _, _)
#     # % project an n-dim vector x to the convex set Cn
#     # % Cn = { x : x n-dim, x > -1, sum(x) <= n}
#     # % Coded by
#     # % Hassan Mohammad hmuhd.mth@buk.edu.ng
#     # % October, 2018
#     n = length(x)
#     bget = false
#     s = sort(x, rev=true)
#     tmpsum = 0

#     for ii = 1:n-1
#         tmpsum = tmpsum + s[ii]
#         tmax = (tmpsum - n) / ii
#         if tmax >= s[ii+1]
#             bget = true
#             break
#         end
#     end
#     tmpsum = tmpsum + s[n]
#     if ~bget
#         tmax = (tmpsum - n) / n
#     end

#     P = if (tmpsum <= n && min(x...) > -1)
#         x
#     else
#         x .- tmax .|> r -> max(r, 0)
#     end
#     P
# end
# function Pj3(x)
#     Pj3(x, 0, 0)
# end

# function Pj3(x, _, _)
#     # % project an n-dim vector x to the convex set Cn
#     # % Cn = { x : x n-dim, x >= 0, sum(x) <= n}
#     # % Coded by
#     # % Hassan Mohammad hmuhd.mth@buk.edu.ng
#     # % October, 2018
#     n = length(x)
#     bget = false
#     s = sort(x, rev=true)
#     tmpsum = 0

#     for ii = 1:n-1
#         tmpsum = tmpsum + s[ii]
#         tmax = (tmpsum - n) / ii
#         if tmax >= s[ii+1]
#             bget = true
#             break
#         end
#     end
#     tmpsum = tmpsum + s[n]
#     if ~bget
#         tmax = (tmpsum - n) / n
#     end

#     P = if (tmpsum <= n && min(x...) >= 0)
#         x
#     else
#         x .- tmax .|> r -> max(r, 0)
#     end
#     P
# end

# function Pj9(x, _, _)
#     # % project an n-dim vector x to the convex set Cn
#     # % Cn = { x : x n-dim, x >= 0, sum(x) = n-1}
#     # % coded by
#     # % Hassan Mohammad hmuhd.mth@buk.edu.ng
#     # % October, 2018
#     n = length(x)
#     if (sum(x) == n - 1 && min(x...) >= 0)
#         return x
#     end

#     bget = false
#     s = sort(x, rev=true)
#     tmpsum = 0
#     for ii = 1:n-1
#         tmpsum = tmpsum + s[ii]
#         tmax = (tmpsum - (n - 1)) / ii
#         if tmax >= s[ii+1]
#             bget = true
#             break
#         end
#     end
#     tmpsum = tmpsum + s[n]
#     if ~bget
#         tmax = (tmpsum - (n - 1)) / n
#     end

#     P = if (tmpsum <= n - 1 && min(x...) >= 0)
#         x
#     else
#         x .- tmax .|> r -> max(r, 0)
#     end
#     P
# end


# function Pj9n(x, _, _)
#     # % project an n-dim vector x to the convex set Cn
#     # % Cn = { x : x n-dim, x >= 0, sum(x) = n}
#     # % coded by
#     # % Hassan Mohammad hmuhd.mth@buk.edu.ng
#     # % October, 2018
#     n = length(x)
#     if (sum(x) == n && min(x...) >= 0)
#         return x
#     end

#     bget = false
#     s = sort(x, rev=true)
#     tmpsum = 0
#     for ii = 1:n-1
#         tmpsum = tmpsum + s[ii]
#         tmax = (tmpsum - n) / ii
#         if tmax >= s[ii+1]
#             bget = true
#             break
#         end
#     end
#     tmpsum = tmpsum + s[n]
#     if ~bget
#         tmax = (tmpsum - n) / n
#     end

#     P = if (tmpsum <= n && min(x...) >= 0)
#         x
#     else
#         x .- tmax .|> r -> max(r, 0)
#     end
#     P
# end


# function Pj15(x, _, _)
#     # % project an n-dim vector x to the convex set Cn
#     # % Cn = { x : x n-dim, x >= -1, sum(x) <= n}
#     # % Coded by
#     # % Hassan Mohammad hmuhd.mth@buk.edu.ng
#     # % October, 2018
#     n = length(x)
#     if (sum(x) <= n && min(x...) >= -1)
#         return x
#     end

#     bget = false
#     s = sort(x, rev=true)
#     tmpsum = 0
#     for ii = 1:n-1
#         tmpsum = tmpsum + s[ii]
#         tmax = (tmpsum - n) / ii
#         if tmax >= s[ii+1]
#             bget = true
#             break
#         end
#     end
#     tmpsum = tmpsum + s[n]
#     if ~bget
#         tmax = (tmpsum - n) / n
#     end

#     P = if (tmpsum <= n && min(x...) >= -1)
#         x
#     else
#         x .- tmax .|> r -> max(r, 0)
#     end
#     P
# end
