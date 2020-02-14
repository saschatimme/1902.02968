using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."))

if !(@isdefined PathTracking)
    include(joinpath(@__DIR__, "../implementation/PathTracking.jl"))
    using .PathTracking
end

using LinearAlgebra, DataFrames, CSV
import HomotopyContinuation
const HC = HomotopyContinuation

result_dir = joinpath(@__DIR__, "../results/katsura/")
mkpath(result_dir)

N = 20
df = DataFrame(n = Int[], t = Float64[], nsuccess = Int[], njumps= Int[])
for n = 8:N
    x = [Variable(Symbol(:x, i)) for i in 0:n]
    K = [
        (sum(x[abs(l)+1] * x[abs(m - l)+1] for l = -n:n if abs(m - l) <= n) -
         x[m+1] for m = 0:n-1)...,
        x[1] + 2 * sum(x[i+1] for i = 1:n) - 1,
    ]

    F = System(K, x)

    H, starts = total_degree_homotopy(F)
    S = collect(starts)
    tracker = Tracker(H)

    # warmup
    track.(tracker, S[1:1], 1, 0)

    t = @elapsed (res = track.(tracker, S, 1, 0))
    nsuccess = count(is_success, res)
    njumps = begin
        mults = HC.multiplicities(solution, res)
        isempty(mults) ? 0 : sum(length, mults)
    end
    @show n
    @show t, nsuccess, njumps
    push!(df, (n, t, nsuccess, njumps))

    # write dataframe to file
    open(joinpath(result_dir, "benchmark.txt"), "w") do io
        CSV.write(io, df; delim = " ")
    end
end
