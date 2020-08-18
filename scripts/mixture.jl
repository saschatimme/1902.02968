using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."))

if !(@isdefined PathTracking)
    include(joinpath(@__DIR__, "../implementation/PathTracking.jl"))
    using .PathTracking
end

include(joinpath(@__DIR__, "utils.jl"))

using LinearAlgebra, DataFrames, CSV, Statistics, Printf, PrettyTables

mixture_system = let
    lambda = Variable(:lambda)
    mu1 = Variable(:mu1)
    mu2 = Variable(:mu2)
    s1 = Variable(:s1)
    s2 = Variable(:s2)
    m1 = Variable(:m1)
    m2 = Variable(:m2)
    m3 = Variable(:m3)
    m4 = Variable(:m4)
    m5 = Variable(:m5)

    eqs = [
        -m1 + lambda * mu1 + (1 - lambda) * mu2,
        -m2 + lambda * (s1 + (1 / 2) * mu1^2) + (1 - lambda) * (s2 + (1 / 2) * mu2^2),
        -m3 + lambda * (s1 * mu1 + (1 / 6) * mu1^3) + (1 - lambda) * (s2 * mu2 + (1 / 6) * mu2^3),
        -m4 +
        lambda * ((1 / 2) * s1 * mu1^2 + (1 / 2) * s1^2 + (1 / 24) * mu1^4) +
        (1 - lambda) * ((1 / 2) * s2 * mu2^2 + (1 / 2) * s2^2 + (1 / 24) * mu2^4),
        -m5 +
        lambda * ((1 / 6) * s1 * mu1^3 + (1 / 2) * s1^2 * mu1 + (1 / 120) * mu1^5) +
        (1 - lambda) * ((1 / 6) * s2 * mu2^3 + (1 / 2) * s2^2 * mu2 + (1 / 120) * mu2^5),
    ]

    System(eqs, [lambda, mu1, s1, mu2, s2], [m1, m2, m3, m4, m5])
end

data_dir = joinpath(@__DIR__, "../data/mixture/")
S = read_solution_file(data_dir * "start_solutions", 5)
p = read_parameters_file(data_dir * "start_parameters")


df = DataFrame(accepted_steps = Int[],
            rejected_steps = Int[], nΔs1 = Int[], nΔs2 = Int[])
for i in 1:50
    q = read_parameters_file(data_dir * "target_parameters_$i")
    tracker = Tracker(ParameterHomotopy(mixture_system, p, q))
    for s in S
        r = track(tracker, s, 1.0, 0.0)
        push!(df, (r.accepted_steps, r.rejected_steps, r.nΔs1, r.nΔs2))
    end
end

table_data = let
    reshape(
        [
            @sprintf("%.2f", mean(df.accepted_steps + df.rejected_steps)),
            @sprintf("%.2f", mean(df.accepted_steps)),
            @sprintf("%.2f", mean(df.rejected_steps)),
            @sprintf("%.2f", mean(df.accepted_steps ./ (df.accepted_steps .+ df.rejected_steps))),
            @sprintf("%.2f", 100 .* mean(df.nΔs1 ./ (df.accepted_steps .+ df.rejected_steps))),
            @sprintf("%.2f", 100 - 100 .* mean(df.nΔs1 ./ (df.accepted_steps .+ df.rejected_steps))),
        ],
        1,
       6,
    )
end


header = [
    "avg. steps" "avg. steps acc." "avg. steps rej." "% acc." "% Δt₁" "% Δt₂"
]

pretty_table(table_data, header)
