using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."))

if !(@isdefined PathTracking)
    include(joinpath(@__DIR__, "../implementation/PathTracking.jl"))
    using .PathTracking
end

include(joinpath(@__DIR__, "utils.jl"))

using LinearAlgebra, DataFrames, CSV, Statistics, Printf, PrettyTables

mixture_system = let
    λ = Variable(:λ)
    μ₁ = Variable(:μ₁)
    μ₂ = Variable(:μ₂)
    s₁ = Variable(:s₁)
    s₂ = Variable(:s₂)
    m₁ = Variable(:m₁)
    m₂ = Variable(:m₂)
    m₃ = Variable(:m₃)
    m₄ = Variable(:m₄)
    m₅ = Variable(:m₅)

    eqs = [
        -m₁ + λ * μ₁ + (1 - λ) * μ₂,
        -m₂ + λ * (s₁ + (1 / 2) * μ₁^2) + (1 - λ) * (s₂ + (1 / 2) * μ₂^2),
        -m₃ + λ * (s₁ * μ₁ + (1 / 6) * μ₁^3) + (1 - λ) * (s₂ * μ₂ + (1 / 6) * μ₂^3),
        -m₄ +
        λ * ((1 / 2) * s₁ * μ₁^2 + (1 / 2) * s₁^2 + (1 / 24) * μ₁^4) +
        (1 - λ) * ((1 / 2) * s₂ * μ₂^2 + (1 / 2) * s₂^2 + (1 / 24) * μ₂^4),
        -m₅ +
        λ * ((1 / 6) * s₁ * μ₁^3 + (1 / 2) * s₁^2 * μ₁ + (1 / 120) * μ₁^5) +
        (1 - λ) * ((1 / 6) * s₂ * μ₂^3 + (1 / 2) * s₂^2 * μ₂ + (1 / 120) * μ₂^5),
    ]

    System(eqs, [λ, μ₁, s₁, μ₂, s₂], [m₁, m₂, m₃, m₄, m₅])
end

data_dir = joinpath(@__DIR__, "../data/mixture/")
S = read_solution_file(data_dir * "start_solutions", 5)
p = read_parameters_file(data_dir * "start_parameters")


df = DataFrame(accepted_steps = Int[],
            rejected_steps = Int[], nΔs₁ = Int[], nΔs₂ = Int[])
for i in 1:50
    q = read_parameters_file(data_dir * "target_parameters_$i")
    tracker = Tracker(ParameterHomotopy(mixture_system, p, q))
    for s in S
        r = track(tracker, s, 1.0, 0.0)
        push!(df, (r.accepted_steps, r.rejected_steps, r.nΔs₁, r.nΔs₂))
    end
end


table_data = let
    reshape(
        [
            @sprintf("%.2f", mean(df.accepted_steps + df.rejected_steps)),
            @sprintf("%.2f", mean(df.accepted_steps)),
            @sprintf("%.2f", mean(df.rejected_steps)),
            @sprintf("%.2f", mean(df.accepted_steps ./ (df.accepted_steps .+ df.rejected_steps))),
            @sprintf("%.2f", mean(df.nΔs₁ ./ (df.accepted_steps .+ df.rejected_steps))),
            @sprintf("%.2f", mean(df.nΔs₂ ./ (df.accepted_steps .+ df.rejected_steps))),
        ],
        1,
       6,
    )
end


header = [
    "avg. steps" "avg. steps acc." "avg. steps rej." "% acc." "% 1" "% 2"
]

pretty_table( table_data, header)
