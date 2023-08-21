
"""
PSEUDOMONOTONE FUNCTION 
"""

"""
Source:
Problem 5.9 in 
journalArticle2022JinbaoJianYinChunmingTangDaolanHan
Jinbao Jian, Yin, J., Chunming Tang, & Daolan Han. (2022). 
A family of inertial derivative-free projection methods for constrained nonlinear pseudo-monotone equations with applications. 
Computational & Applied Mathematics, 41(7). https://doi.org/10.1007/s40314-022-02019-6




Projected on [-10, 10]

Original Source:
- id: journalArticle2019ShehuDongJiang
  accessed:
    - year: 2023
      month: 8
      day: 17
  author:
    - family: Shehu
      given: Yekini
    - family: Dong
      given: Qiao-Li
    - family: Jiang
      given: Dan
  citation-key: journalArticle2019ShehuDongJiang
  container-title: Optimization
  DOI: 10.1080/02331934.2018.1522636
  ISSN: 0233-1934
  issue: '1'
  issued:
    - year: 2019
      month: 1
      day: 2
  page: 385-409
  publisher: Taylor & Francis
  source: Taylor and Francis+NEJM
  title: >-
    Single projection method for pseudo-monotone variational inequality in
    Hilbert spaces
  type: article-journal
  URL: https://doi.org/10.1080/02331934.2018.1522636
  volume: '68'

"""

function PseudoMonotone1(x)
    vcat((x[1] + (x[2] - 1)^2) * (1 + x[2]),
        -x[1]^3 - x[1] * (x[2] - 1)^2
    )

end

"""
Source:
Problem 5.10 in 
journalArticle2022JinbaoJianYinChunmingTangDaolanHan
Jinbao Jian, Yin, J., Chunming Tang, & Daolan Han. (2022). 
A family of inertial derivative-free projection methods for constrained nonlinear pseudo-monotone equations with applications. 
Computational & Applied Mathematics, 41(7). https://doi.org/10.1007/s40314-022-02019-6



Projected on [-5, 5]
"""
function PseudoMonotone2(x)
    M = [
        1 0 -1
        0 1.5 0
        -1 0 2
    ]
    (exp(-norm(x)^2) + 2) * M * x
end
#### OLD ONES ------------------



"""
Source:
Problem 4.6 in
Jamilu Sabi'u, Abdullah Shah, Predrag S. Stanimirović, Branislav Ivanov, Mohammed Yusuf Waziri, 
Modified optimal Perry conjugate gradient method for solving system of monotone equations with applications,Applied Numerical Mathematics,
Volume 184, 2023, Pages 431-445, ISSN 0168-9274, https://doi.org/10.1016/j.apnum.2022.10.016.

Projected on [0, ∞] or [-1, ∞]

"""
function PolynomialSineCosine(x)
    N = 1 / length(x)
    s = N * sum(x)
    x .* cos.(x .- N) .* (
        sin.(x) .- 1 .- (x .- 1) .^ 2 .- s
    )
end

"""
Source: 
Problem 5.2
Zhou, Weijun, and Donghui Li. “LIMITED MEMORY BFGS METHOD FOR NONLINEAR MONOTONE EQUATIONS.” 
Journal of Computational Mathematics, vol. 25, no. 1, 2007, pp. 89–96. https://www.global-sci.org/intro/article_detail/jcm/8675.html#
Example: 4.1 in
Wang, C., Wang, Y. & Xu, C. A projection method for a system of nonlinear monotone equations with convex constraints. 
Math Meth Oper Res 66, 33–46 (2007). https://doi.org/10.1007/s00186-006-0140-y
"""
function ExponetialI(x)
    exp(x) - 1
end



"""
Source:
Problem 3 in Appendix
William La Cruz & Marcos Raydan (2003) Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems, 
Optimization Methods and Software, 18:5, 583-599, DOI: 10.1080/10556780310001610493

"""
function ExponetialIII(x)
    X2 = x .^ 2
    Ii = [i / 10 for i in 1:length(x)]
    Y = vcat(
        X2[1:end-1] .+ exp.(-X2[1:end-1]),
        exp(-X2[end])
    )
    Ii .* (1 .- Y)
end


"""
Source:
Problem 9 
J. Sabi’u, A. Shah, M.Y. Waziri, M.K. Dauda, A new hybrid approach for solving large-scale monotone nonlinear equations, 
J. Math. Fund. Sci. 52 (2020) 17–26 https://doi.org/10.5614/j.math.fund.sci.2020.52.1.2


"""
function PolynomialI(x)
    vcat(
        4 .* x[1:end-1] .+ (x[2:end] .- 2 .* x[1:end-1]) .- ((x[2:end] .^ 2) / 3),
        4 * x[end] + (x[end-1] - 2 * x[end]) - x[end-1]^2 / 3
    )
end



"""
Source:
Problem 2 in
Zhou, Weijun, and Donghui Li. “LIMITED MEMORY BFGS METHOD FOR NONLINEAR MONOTONE EQUATIONS.” 
Journal of Computational Mathematics, vol. 25, no. 1, 2007, pp. 89–96.

Project to Rn₊

Modified Problem 1. in
W.J. Zhou, D.H. Li, A globally convergent BFGS method for nonlinear monotone equations without any merit functions, 
Math. Comput. 77 (264) (2008) 2231–2240.
Projected on [-2, ∞]


"""
function SmoothSine(x)
    v = try

        2 * x - sin(x)
    catch
        println(max(x...), min(x...))
        throw("problem with this z")
    end
    v
end



"""
Source:
Problem 1 in
Zhou, Weijun, and Donghui Li. “LIMITED MEMORY BFGS METHOD FOR NONLINEAR MONOTONE EQUATIONS.” 
Journal of Computational Mathematics, vol. 25, no. 1, 2007, pp. 89–96.tional and Applied Mathematics, Volume 196, Issue 2, 2006, Pages 478-484, ISSN 0377-0427, https://doi.org/10.1016/j.cam.2005.10.002.

Project to Rn₊
"""
function NonsmoothSine(x)
    2 * x - sin(abs(x))
end


"""
Source:
Problem 2 in
Zhensheng Yu, Ji Lin, Jing Sun, Yunhai Xiao, Liying Liu, Zhanhui Li, Spectral gradient projection method for monotone nonlinear equations with convex constraints,
Applied Numerical Mathematics, Volume 59, Issue 10, 2009, Pages 2416-2423, ISSN 0168-9274, https://doi.org/10.1016/j.apnum.2009.04.004.

with Projection on 
    C = {x ∈ Rⁿ | ∑xᵢ ≤ n, xᵢ≥ -1}

Modified from Problem 1 in
Li Zhang, Weijun Zhou, Spectral gradient projection method for solving nonlinear monotone equations, 
Journal of Computational and Applied Mathematics, Volume 196, Issue 2, 2006, Pages 478-484, ISSN 0377-0427, https://doi.org/10.1016/j.cam.2005.10.002.

"""
function ModifiedNonsmoothSine(x)
    x - sin(abs(x - 1))
end

"""
Source:
Problem 2
Gao PT, He CJ (2018) An efficient three-term conjugate gradient method for nonlinear monotone equations
with convex constraints. Calcolo. https://doi.org/10.1007/s10092-018-0291-2


with Projection on 
    C = {x ∈ Rⁿ | ∑xᵢ ≤ n, xᵢ≥ -1}

Modified from Problem 1 in
Li Zhang, Weijun Zhou, Spectral gradient projection method for solving nonlinear monotone equations, 
Journal of Computational and Applied Mathematics, Volume 196, Issue 2, 2006, Pages 478-484, ISSN 0377-0427, https://doi.org/10.1016/j.cam.2005.10.002.

"""
function ModifiedNonsmoothSine2(x)
    x - sin(abs(x) - 1)
end


"""
Source:
Problem 1
Gao, P., He, C. An efficient three-term conjugate gradient method for nonlinear monotone equations with convex constraints. 
Calcolo 55, 53 (2018). https://doi.org/10.1007/s10092-018-0291-2

with Projection on 
    Rn+

"""
function ExponetialSineCosine(x)
    (exp.(x)) .^ 2 + 3 * sin.(x) * cos.(x) - 1
end

"""
Source: 
Problem 4.2 in
Li Zheng, Lei Yang, Yong Liang, A conjugate gradient projection method for solving equations with convex constraints,
Journal of Computational and Applied Mathematics, Volume 375, 2020, 112781, ISSN 0377-0427,
https://doi.org/10.1016/j.cam.2020.112781.
Modified from Problem 3
Zhou, Weijun, and Donghui Li. “LIMITED MEMORY BFGS METHOD FOR NONLINEAR MONOTONE EQUATIONS.” 
Journal of Computational Mathematics, vol. 25, no. 1, 2007, pp. 89–96.
"""
function ModifiedTrigI(x)
    vcat(
        x[1] + sin(x[1]) - 1,
        -x[1:end-2] + 2 * x[2:end-1] + sin.(x[2:end-1]) .- 1,
        x[end] + sin(x[end]) - 1
    )
end

"""
Source:
Problem 10 in 
Y. Bing, G. Lin, An efficient implementation of Merrill’s method for sparse or partially separable systems of nonlinear equations, 
SIAM. J. Optim.1 (2) (1991) 206–221.
"""
function Tridiagonal(x)
    Ii = [i for i in 1:length(x)]
    x0 = vcat(0, x)
    x1 = vcat(x, 0)
    x .- exp.(cos.(Ii .* (x0[1:end-1] .+ x .+ x1[2:end])))
end


"""
Source:
Problem 4.4 in
Li Zheng, Lei Yang, Yong Liang, A conjugate gradient projection method for solving equations with convex constraints,
Journal of Computational and Applied Mathematics, Volume 375, 2020, 112781, ISSN 0377-0427,
https://doi.org/10.1016/j.cam.2020.112781.
Modified from Problem 10 in 
Y. Bing, G. Lin, An efficient implementation of Merrill’s method for sparse or partially separable systems of nonlinear equations, 
SIAM. J. Optim.1 (2) (1991) 206–221.
"""
function ModifiedTridiagonal(x)
    n = length(x)
    vcat(
        x[1] - exp(cos((x[1] + x[2]) / (n + 1))),
        -x[2:end-1] - exp.(cos.((x[1:end-2] .+ x[2:end-1] .+ x[3:end]) ./ (n + 1))),
        x[end] - exp(cos((x[end-1] + x[end]) / (n + 1)))
    )
end


"""
Source:
Problem 10 in Appendix
William La Cruz & Marcos Raydan (2003) Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems, 
Optimization Methods and Software, 18:5, 583-599, DOI: 10.1080/10556780310001610493


Problem 4.5 in
Li Zheng, Lei Yang, Yong Liang, A conjugate gradient projection method for solving equations with convex constraints,
Journal of Computational and Applied Mathematics, Volume 375, 2020, 112781, ISSN 0377-0427, https://doi.org/10.1016/j.cam.2020.112781

Projected on [-1, ∞]

"""
function Logarithmic(x)
    log.(x .+ 1) .- (x / length(x))
end


"""
Source:
Problem 10 in Appendix
William La Cruz & Marcos Raydan (2003) Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems, 
Optimization Methods and Software, 18:5, 583-599, DOI: 10.1080/10556780310001610493


Modified in Problem 4.2 in
Jamilu Sabi'u, Abdullah Shah, Predrag S. Stanimirović, Branislav Ivanov, Mohammed Yusuf Waziri, 
Modified optimal Perry conjugate gradient method for solving system of monotone equations with applications,Applied Numerical Mathematics,
Volume 184, 2023, Pages 431-445, ISSN 0168-9274, https://doi.org/10.1016/j.apnum.2022.10.016.

Projected on [0, ∞] or [-1, ∞]

"""
function NonmoothLogarithmic(x)
    log.(abs.(x) .+ 1) .- (x / length(x))
end

function ding2017FunII(x)
    x - sin(abs(x - 1))
end

function ding2017Diagonal6Fun(x)
    cp5(x)
end

# see https://www.cuter.rl.ac.uk//Problems/classification.shtml for classification

# NAME: ARWHEAD         see (https://bitbucket.org/optrove/sif/raw/HEAD/ARWHEAD.SIF)
# see also (https://vanderbei.princeton.edu/ampl/nlmodels/cute/arwhead.mod)
#  classification OUR2-AN-V-0
# O: none of the aboves
# U: the problem is unconstrained
# R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# 2: the degree of the highest derivatives provided analytically within the problem description.
# A: academic problem
# N: the problem description does not contain any explicit internal variables
# V: the number of variables in the problem can be chosen by the user
# 0: a nonnegative integer giving the actual (fixed) number of problem constraints.
function ARWHEADFun(x::Vector{<:Real})
    # sum {i in 1..N-1} (-4*x[i]+3.0) + sum {i in 1..N-1} (x[i]^2+x[N]^2)^2
    n = length(x)
    sum((-4 * x[i] + 3.0) for i in 1:n-1) + sum((x[i]^2 + x[n]^2)^2 for i in 1:n-1)
end
function ARWHEADGrad(x::Vector{<:Real})
    vcat(-4 .+ 4 * x[1:end-1] .* (x[1:end-1] .^ 2 .+ x[end]^2),
        4 * x[end] * sum(x[1:end-1] .^ 2 .+ x[end]^2))
end


# name: PENALTY1 see (https://bitbucket.org/optrove/sif/raw/HEAD/PENALTY1.SIF) 
# see also (https://vanderbei.princeton.edu/ampl/nlmodels/cute/penalty1.mod)
# classification SUR2-AN-V-0
# S: the objective function is a sum of squares
# U: the problem is unconstrained
# R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# 2: the degree of the highest derivatives provided analytically within the problem description.
# A: the problem is academic, that is, has been constructed specifically by researchers to test one or more algorithms,
# N: the problem description does not contain any explicit internal variables.
# V: the number of variables in the problem can be chosen by the user
# 0: a nonnegative integer giving the actual (fixed) number of problem constraints.

function PENALTY1Fun(x::Vector{<:Real})
    a = 1e-5
    N = length(x)
    # sum {i in 1..N} a*(x[i]-1)^2 + ( sum {j in 1..N} x[j]^2 - 1/4 )^2;
    sum(a * (x[i] - 1)^2 for i in 1:N) + (sum(x[i]^2 for i in 1:N) - 0.25)^2
end
function PENALTY1Grad(x::Vector{<:Real})
    t = sum(x .^ 2)
    c = 1e-5
    2 * c .* (x .- 1) + 4 * (t - 0.25) .* x
end

# NAME: DIXON3DQ see (https://bitbucket.org/optrove/sif/raw/HEAD/DIXON3DQ.SIF)
# SEE ALSO (https://vanderbei.princeton.edu/ampl/nlmodels/cute/dixon3dq.mod)
# classification QUR2-AN-V-0
# Q: the objective function is quadratic,
# U: the problem is unconstrained
# R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# 2: the degree of the highest derivatives provided analytically within the problem description.
# A: academic problem
# N: the problem description does not contain any explicit internal variables
# V: the number of variables in the problem can be chosen by the user
# 0: a nonnegative integer giving the actual (fixed) number of problem constraints.
function DIXON3DQFun(x::Vector{<:Real})
    # (x[1]-1.0)^2 + sum {j in 2..n-1} (x[j]-x[j+1])^2 + (x[n]-1.0)^2
    n = length(x)
    (x[1] - 1.0)^2 + sum((x[j] - x[j+1])^2 for j in 2:n-1) + (x[n] - 1.0)^2
end
function DIXON3DQGrad(x::Vector{<:Real})
    n = length(x)
    return vcat(2 * (x[1] - 1),
        2 * (x[2] - x[3]),
        2 * (2 * x[3:n-1] - x[2:n-2] - x[4:n]),
        2 * (2 * x[n] - x[n-1] - 1)
    )
end

# NAME: GENHUMPS, see (https://bitbucket.org/optrove/sif/raw/HEAD/GENHUMPS.SIF)
# SEE ALSO (https://vanderbei.princeton.edu/ampl/nlmodels/cute/genhumps.mod)
# classification OUR2-AN-V-0
# O: none of the aboves
# U: the problem is unconstrained
# R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# 2: the degree of the highest derivatives provided analytically within the problem description.
# A: academic problem
# N: the problem description does not contain any explicit internal variables
# V: the number of variables in the problem can be chosen by the user
# 0: a nonnegative integer giving the actual (fixed) number of problem constraints.
function GENHUMPSFun(x::Vector{<:Real}; ζ=20)
    # sum {i in 1..N-1} ( sin (zeta*x[i])^2*sin(zeta*x[i+1])^2+0.05*(x[i]^2+x[i+1]^2) );
    N = length(x)
    sum((sin(ζ * x[i]) * sin(ζ * x[i+1]))^2 + 0.05 * (x[i]^2 + x[i+1]^2) for i in 1:(N-1))
end
function GENHUMPSGrad(x::Vector{<:Real}; ζ=20)
    TNTHX = 0.1 * x
    ZX = ζ * x
    ZX2 = 2 * ZX
    SZX = (sin.(ZX)) .^ 2
    SZX2 = sin.(ZX2)
    ZSZX2 = ζ * SZX2
    vcat(
        ZSZX2[1] * SZX[2] + TNTHX[1],
        ZSZX2[2:end-1] .* (SZX[1:end-2] .+ SZX[3:end]) .+ (2 * TNTHX[2:end-1]),
        ZSZX2[end] * SZX[end-1] + TNTHX[end],
    )
end

# NAME: ENGVAL1 see (https://bitbucket.org/optrove/sif/raw/HEAD/ENGVAL1.SIF)
# SEE ALSO (https://vanderbei.princeton.edu/ampl/nlmodels/cute/engval1.mod)
# classification OUR2-AN-V-0
# O: none of the aboves
# U: the problem is unconstrained
# R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# 2: the degree of the highest derivatives provided analytically within the problem description.
# A: academic problem
# N: the problem description does not contain any explicit internal variables
# V: the number of variables in the problem can be chosen by the user
# 0: a nonnegative integer giving the actual (fixed) number of problem constraints.
function ENGVAL1Fun(x::Vector{<:Real})
    # sum {i in 1..N-1} (x[i]^2+x[i+1]^2)^2 + sum {i in 1..N-1} (-4*x[i]+3.0);
    N = length(x)
    sum((x[i]^2 + x[i+1]^2)^2 for i in 1:(N-1)) + sum(-4 * x[i] + 3.0 for i in 1:(N-1))
end
function ENGVAL1Grad(x::Vector{<:Real})
    X2 = x .^ 2
    vcat(
        4 * (x[1] * (X2[1] + X2[2]) - 1),
        4 * (x[2:end-1] .* (X2[1:end-2] + 2 * X2[2:end-1] + X2[3:end]) .- 1),
        4 * x[end] * (X2[end-1] + X2[end])
    )
end
# NAME: DIXMAANH see (https://bitbucket.org/optrove/sif/raw/HEAD/DIXMAANH.SIF)
# SEE ALSO (https://vanderbei.princeton.edu/ampl/nlmodels/cute/dixmaanh.mod)
# classification OUR2-AN-V-0
# O: none of the aboves
# U: the problem is unconstrained
# R: the problem is regular, that is its first and second derivatives exist and are continuous everywhere
# 2: the degree of the highest derivatives provided analytically within the problem description.
# A: academic problem
# N: the problem description does not contain any explicit internal variables
# V: the number of variables in the problem can be chosen by the user
# 0: a nonnegative integer giving the actual (fixed) number of problem constraints.
function DIXMAANHFun(x::Vector{<:Real}; α=1.0, β=0.26, γ=0.26, δ=0.26, K=[1 0 0 1])
    N = length(x)
    if N % 3 != 0
        throw("The length of `x0` must be divisible by 3")
    end
    M = fld(N, 3)
    1.0 + α * sum(x[i]^2 * (i / N)^K[1] for i in 1:N) +
    β * sum(x[i]^2 * (x[i+1] + x[i+1]^2)^2 for i in 1:N-1) +
    γ * sum(x[i]^2 * x[i+M]^4 for i in 1:2*M) +
    δ * sum(x[i] * x[i+2*M] * (i / N)^K[4] for i in 1:M)
end
function DIXMAANHGrad(x::Vector{<:Real}; α=1.0, β=0.26, γ=0.26, δ=0.26, K=[1 0 0 1])
    N = length(x)
    if N % 3 != 0
        throw("The length of `x0` must be divisible by 3")
    end
    M = fld(N, 3)
    X2 = x .^ 2
    X4 = X2 .* X2
    IOverN = [(i / N)^K[1] for i in 1:N]
    IOverM = [(i / N)^K[4] for i in 1:M]
    TwoAlphaXOverN = 2 * α * IOverN .* x
    OnePlusX = 1 .+ x
    OnePlus2X = OnePlusX .+ x
    OnePlusX2 = OnePlusX .^ 2
    G1 = TwoAlphaXOverN[1] + 2 * β * x[1] * X2[2] * OnePlusX2[2] + 2 * γ * x[1] * X4[1+M] + δ * IOverM[1] * x[1+2*M]
    G2M = TwoAlphaXOverN[2:M] .+ 2 * β * (X2[1:M-1] .* x[2:M] .* OnePlusX[2:M] .* OnePlus2X[2:M] .+ x[2:M] .* X2[3:M+1] .* OnePlusX2[3:M+1]) .+
          2 * γ * x[2:M] .* X4[2+M:2*M] .+ δ * IOverM[2:M] .* x[2+2*M:3*M]

    # check
    GM12M = TwoAlphaXOverN[M+1:2*M] .+
            2 * β * (
                X2[M:2*M-1] .* x[M+1:2*M] .* OnePlusX[M+1:2*M] .* OnePlus2X[M+1:2*M] .+ x[M+1:2*M] .* X2[M+2:2*M+1] .* OnePlusX2[M+2:2*M+1]
            ) .+
            γ * (2 * x[M+1:2*M] .* X4[2*M+1:3*M] .+ 4 * X2[1:M] .* x[M+1:2*M] .* X2[M+1:2*M])


    # check
    G2M1Nm1 = TwoAlphaXOverN[2*M+1:N-1] .+
              2 * β * (
                  X2[2*M:N-2] .* x[2*M+1:N-1] .* OnePlusX[2*M+1:N-1] .* OnePlus2X[2*M+1:N-1] .+ x[2*M+1:N-1] .* X2[2*M+2:N] .* OnePlusX2[2*M+2:N]
              ) .+
              (γ * 4 * X2[M+1:2*M-1] .* x[2*M+1:N-1] .* X2[2*M+1:N-1]) .+
              δ * IOverM[1:M-1] .* x[1:M-1]

    GN = TwoAlphaXOverN[N] + 2 * β * (X2[N-1] * x[N] * OnePlusX[N] * OnePlus2X[N]) +
         γ * (4 * X2[2*M] * x[N] * X2[N]) + δ * IOverM[M] .* x[M]
    vcat(
        G1,
        G2M,
        GM12M,
        G2M1Nm1,
        GN
    )

end

function DIXMAANIFun(x::Vector{<:Real})
    DIXMAANHFun(x; α=1.0, β=0.0, γ=0.125, δ=0.125, K=[2 0 0 2])
end
function DIXMAANIGrad(x::Vector{<:Real})
    DIXMAANHGrad(x; α=1.0, β=0.0, γ=0.125, δ=0.125, K=[2 0 0 2])
end




# function cp1(x)
#     # % Evaluate the modified exponential function 2 of La Cruz et. al. 2004
#     # % for convex constraints
#     # % x0=(1/(n^2),1/(n^2),...)' 
#     # % to call the function evaluation and projection type [P,F]=feval(f,x)
#     # % to call function evaluation only type [~,F]=feval(f,x)
#     # % last update 30/03/2018

#     F = [i == 1 ? exp(x[i]) - 1 : exp(x[i]) + x[i] - 1 for i in eachindex(x)]
#     F
#     # %P=max(x,0); % projection
# end

# cp2(x) = log(x + 1) - x / length(x);


# function cp3(x)
#     # % Evaluate a nonsmooth and monotone function
#     # %
#     n = length(x)
#     # %F = zeros(n,1);
#     # %h=1/(n+1);
#     2 * x - sin(abs(x))

# end

# function cp4(x)
#     # % Evaluate the  Mononote problem La Cruz 2017 for convex constraint
#     # % x0=(1,1,...) x0=(1,1/2,1/3,...,1/n)
#     # % to call the function evaluation and projection type [P,F]=feval(f,x)
#     # % to call function evaluation only type [~,F]=feval(f,x)
#     # % last update 30/03/2018
#     n = length(x)
#     i = 1:n
#     # %F(i)=min(min(abs(x(i)),x(i).^2),max(abs(x(i)),x(i).^3));
#     min(min(x..., [x .^ 2]...), max(x..., [x .^ 3]...))

# end
# function cp5(x)
#     # % Evaluate the Strictly convex function 1 La Cruz et. al. 2004
#     # % x0=(1/n,2/n,...,1)' (1, 1/2, 1/3,.., 1/n)'
#     # % to call the function evaluation and projection type [P,F]=feval(f,x)
#     # % to call function evaluation only type [~,F]=feval(f,x)
#     # % last update 26/03/2018
#     exp(x) - 1
# end
# function cp51(x)
#     # % Evaluate the Strictly convex function 1 La Cruz et. al. 2004
#     # % x0=(1/n,2/n,...,1)' (1, 1/2, 1/3,.., 1/n)'
#     # % to call the function evaluation and projection type [P,F]=feval(f,x)
#     # % to call function evaluation only type [~,F]=feval(f,x)
#     # % last update 26/03/2018
#     (exp.(x)) .^ 2 + 3 * sin.(x) .* cos.(x) .- 1

# end
# function cp6(x)
#     # % Evaluate the Strictly convex function 2 La Cruz et. al. 2004
#     # % x0=(1/n,2/n,...,1)' (1, 1/2, 1/3,.., 1/n)'
#     # % to call the function evaluation and projection type [P,F]=feval(f,x)
#     # % to call function evaluation only type [~,F]=feval(f,x)
#     # % last update 03/11/2018
#     n = length(x)
#     [(i / n) * exp(x[i]) - 1 for i in 1:n]

# end

# function cp7(x)
#     # % Evaluate the tridiagonal exponential function 
#     # % x0=(1.5,1.5,...)
#     n = length(x)
#     h = 1 / (n + 1)
#     vcat(x[1] - exp(cos(h * (x[1] + x[2]))),
#         x[2:end-1] .- exp.(cos.(h .* (x[1:end-2] + x[2:end-1] + x[3:end]))),
#         x[end] - exp(cos(h * (x[end-1] + x[end])))
#     )
# end

# function cp8(x)
#     # % Evaluate the monotone and nonsmooth (at x=1) function
#     # % Yu, Niu & Ma 2013 JIMO
#     x - sin(abs(x - 1))
# end


# function cp13(x, Extra)
#     # # % Evaluate the discretize of nondiffrentiable Dirichlet problem
#     # # % n must be perfect square
#     # n = length(x);
#     # r=sqrt(n);
#     # h=1/(r+1);
#     # B=(gallery('tridiag',r,-1,4,-1));
#     # A=kron(-eye(r),B);
#     # F(i)=A*x-h^2*((max(x(i)-1,0.5*x(i)-0.5))+ones(n,1));
#     @error("cp13 not implemeted")
# end
# function cp14(x)
#     # % Evaluate the Trigexp function  prob 52 Luksan & VIcek 2003
#     # % solution x*=(1,1,...,1)'
#     # % last update 31/10/2018
#     vcat(3 * x[1]^3 + 2 * x[2] - 5 + sin(x[1] - x[2]) * sin(x[1] + x[2]),
#         3 * x[2:end-1] .^ 3 + 2 * x[3:end] .- 5 + sin.(x[2:end-1] .- x[3:end]) .* sin.(x[2:end-1] + x[3:end]) + 4 * x[2:end-1] .- x[1:end-2] .* exp.(x[1:end-2] .- x[2:end-1]) .- 3,
#         -x[end-1] * exp(x[end-1] - x[end]) + 4 * x[end] - 3
#     )

# end

# function cp141(x)
#     # % Pursuit-Evasion problem 21/1/19
#     √8 * x - 1
# end

# function cp15(x)
#     # % Evaluate the semismooth function Yamashita and Fukushima 1997 Modified
#     # % Newton method for solving semismooth... Math. Program. 76
#     # % length of x is 4, that is F is from R^4 to R^4
#     # % x*=(2,0,1,0)'
#     # % last modified 30/10/2018
#     n = length(x)
#     # %F = zeros(n,1);
#     A = [1 0 0 0; 0 1 -1 0; 0 1 1 0; 0 0 0 0]
#     b = [-10; 1; -3; 0]
#     y = [x[1]^3; x[2]^3; 2 * x[3]^3; 2 * x[4]^3]
#     # %F=zeros(n,1);
#     F = A * x + y + b
#     # .+y .+ b;
#     F
# end

# function cp16(x)
#     # % Evaluate Problem Penalty 1
#     # % x0=(100,100,...)
#     c = 1e-5
#     2 * c * (x[1:end] .- 1) .+ 4 * (sum(x[1:end]) - 0.25) * x[1:end]

# end

