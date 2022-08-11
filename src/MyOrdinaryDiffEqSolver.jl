module MyOrdinaryDiffEqSolver

using UnPack
using Reexport
@reexport using DiffEqBase
using SciMLBase
@reexport using OrdinaryDiffEq
import OrdinaryDiffEq as ODE

############################################################################################
# Algorithm
############################################################################################
struct MyAlg <: ODE.OrdinaryDiffEqAdaptiveAlgorithm end
export MyAlg
ODE.alg_order(alg::MyAlg) = 1
ODE.isfsal(::MyAlg) = false
ODE.isadaptive(::MyAlg) = false

############################################################################################
# Cache
############################################################################################
mutable struct MyCache{uT} <: ODE.OrdinaryDiffEqCache
    du::uT
    tmp::uT  # without `tmp` it errors
end
function ODE.alg_cache(
    alg::MyAlg,
    u,
    rate_prototype,
    ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits},
    ::Type{tTypeNoUnits},
    uprev,
    uprev2,
    f,
    t,
    dt,
    reltol,
    p,
    calck,
    ::Val{IIP},
) where {IIP,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    @debug "alg_cache"
    du = similar(u)
    tmp = similar(u)
    return MyCache(du, tmp)
end

############################################################################################
# Actual solver: Initialization and step
############################################################################################
function ODE.initialize!(integ, cache::MyCache)
    @debug "initialize!"
end
function ODE.perform_step!(integ, cache::MyCache, repeat_step = false)
    @unpack du, u, f, p, t, dt = integ
    @unpack du = cache
    f(du, u, p, t)
    integ.destats.nf += 1
    @. u = u + dt * du
    if integ.opts.adaptive
        # compute some error estimate here
        # integ.EEst =
    end
end

############################################################################################
# Customizing other calls (optional)
############################################################################################
function ODE.postamble!(integ::ODE.ODEIntegrator{<:MyAlg})
    ODE._postamble!(integ)
end
function DiffEqBase.savevalues!(
    integ::ODE.ODEIntegrator{<:MyAlg},
    force_save = false,
    reduce_size = true,
)
    out = ODE._savevalues!(integ, force_save, reduce_size)
end

############################################################################################
# Solution and Interpolation
############################################################################################
struct MyODESolution{T,N,uType,uType2,DType,tType,rateType,P,A,IType,DE} <:
       SciMLBase.AbstractODESolution{T,N,uType}
    u::uType
    u_analytic::uType2
    errors::DType
    t::tType
    k::rateType
    prob::P
    alg::A
    interp::IType
    dense::Bool
    tslocation::Int
    destats::DE
    retcode::Symbol
end
MyODESolution{T,N}(
    u, u_analytic, errors, t, k, prob, alg, interp, dense, tslocation, destats, retcode,
) where {T,N} = MyODESolution{
    T,
    N,
    typeof(u),
    typeof(u_analytic),
    typeof(errors),
    typeof(t),
    typeof(k),
    typeof(prob),
    typeof(alg),
    typeof(interp),
    typeof(destats),
}(
    u,
    u_analytic,
    errors,
    t,
    k,
    prob,
    alg,
    interp,
    dense,
    tslocation,
    destats,
    retcode,
)

function SciMLBase.build_solution(
    prob::SciMLBase.AbstractODEProblem,
    alg::MyAlg,
    t,
    u;
    timeseries_errors = length(u) > 2,
    dense = false,
    dense_errors = dense,
    calculate_error = true,
    k = nothing,
    interp = nothing,
    retcode = :Default,
    destats = nothing,
    kwargs...,
)
    @debug "build_solution"
    T = eltype(eltype(u))

    if prob.u0 === nothing
        N = 2
    else
        N = length((size(prob.u0)..., length(u)))
    end

    if typeof(prob.f) <: Tuple
        f = prob.f[1]
    else
        f = prob.f
    end

    interp = MyInterpolation(t, u)

    if SciMLBase.has_analytic(f)
        u_analytic = Vector{typeof(prob.u0)}()
        errors = Dict{Symbol,real(eltype(prob.u0))}()
        sol = MyODESolution{T,N}(
            u,
            u_analytic,
            errors,
            t,
            k,
            prob,
            alg,
            interp,
            dense,
            0,
            destats,
            retcode,
        )
        if calculate_error
            calculate_solution_errors!(
                sol;
                timeseries_errors = timeseries_errors,
                dense_errors = dense_errors,
            )
        end
        return sol
    else
        return MyODESolution{T,N}(
            u,
            nothing,
            nothing,
            t,
            k,
            prob,
            alg,
            interp,
            dense,
            0,
            destats,
            retcode,
        )
    end
end

function SciMLBase.solution_new_retcode(sol::MyODESolution{T,N}, retcode) where {T,N}
    return MyODESolution{T,N}(
        sol.u, sol.u_analytic, sol.errors, sol.t, sol.k, sol.prob, sol.alg, sol.interp,
        sol.dense, sol.tslocation, sol.destats, retcode,
    )
end

struct MyInterpolation <: DiffEqBase.AbstractDiffEqInterpolation
    ts::Any
    timeseries::Any
end
DiffEqBase.interp_summary(interp::MyInterpolation) = "My custom interpolation"

(interp::MyInterpolation)(tval::Number, idxs, deriv, p, continuity::Symbol = :left) = begin
    @debug "called MyInterpolation(::Number)"
    @unpack ts, timeseries = interp
    tdir = sign(ts[end] - ts[begin])
    t = tval
    iprev = max(firstindex(ts), ODE._searchsortedlast(ts, t, 2, tdir > 0))
    inext = t > ts[end] ? nothing : min(lastindex(ts), iprev + 1)
    return timeseries[iprev]
end

(interp::MyInterpolation)(tvals, idxs, deriv, p, continuity::Symbol = :left) = begin
    @debug "called MyInterpolation()"
    @unpack ts = interp
    tdir = sign(ts[end] - ts[begin])
    idx = sortperm(tvals, rev = tdir < 0)

    vals = map(idx) do j
        t = tvals[j]
        return interp(t, idxs, deriv, p, continuity)
    end
    invpermute!(vals, idx)
    return SciMLBase.DiffEqArray(vals, tvals)
end

end
