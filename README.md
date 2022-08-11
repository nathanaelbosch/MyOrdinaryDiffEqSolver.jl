# How to build a solver on top of OrdinaryDiffEq.jl

Goal of this repo is to show a minimal example for how to build a custom ODE solver on top of OrdinaryDiffEq.jl.
If you want to contribute to OrdinaryDiffEq.jl by adding a solver, this is not the place to look - for that check out the corresponding section in the [SciML Developer Documentation](https://devdocs.sciml.ai/dev/contributing/adding_algorithms/#Adding-new-algorithms-to-OrdinaryDiffEq-1).
Here, we care about ODE solvers that don't quite fit the framework in OrdinaryDiffEq.jl, e.g. because they need custom solution objects (such as [ProbNumDiffEq.jl](https://github.com/nathanaelbosch/ProbNumDiffEq.jl)).

In a nutshell, define:
- `MyAlgorithm` + traits
- `MyCache` + `OrdinaryDiffEq.alg_cache`
- `OrdinaryDiffEq.initialize!` + `OrdinaryDiffEq.perform_step!`
- `MyODESolution` with multiple constructors and make it callable
- `MyInterpolation` to handle the interpolation

## Example
```julia
using MyOrdinaryDiffEqSolver, Plots

# Define ODEProblem
function f(du, u, p, t)
    a, b, c = p
    du[1] = c * (u[1] - u[1]^3 / 3 + u[2])
    du[2] = -(1 / c) * (u[1] - a - b * u[2])
end
u0 = [-1.0, 1.0]
tspan = (0.0, 20.0)
p = (0.2, 0.2, 3.0)
prob = ODEProblem(f, u0, tspan, p)

# Initializing the integrator works as usual
integ = init(prob, MyAlg(), adaptive = false, dt = 1e-2, dense = false);
# Solving as well
sol = solve(prob, MyAlg(), adaptive = false, dt = 1e-2, dense = false);

# `MyODESolution` subtypes `AbstractODESolution` and can do the things you expect it to do
sol.t, sol.u
sol[2:3]
sol(0.1)
plot(sol)
```
