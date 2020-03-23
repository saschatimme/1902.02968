using Pkg;
Pkg.activate(joinpath(@__DIR__, "..",))

include(joinpath(@__DIR__, "../implementation/PathTracking.jl"))
using .PathTracking

include(joinpath(@__DIR__, "utils.jl"))

@var x t q
H = Homotopy([x^2 - (t - 1/2)^2 - q^2], [x], t, [q])

for k in 1:7
    p = exp10(-k)
    S = [[sqrt(0.25 + p^2)], [-sqrt(0.25 + p^2)]]
    tracker = Tracker(ModelKitHomotopy(H, [p]))
    out = norm.(solution.(track.(tracker, S, 1.0, 0.0)) .- S)
    println("p = $p : ", sum(out) < 1e-12 ? "ok" : "failed")
end
