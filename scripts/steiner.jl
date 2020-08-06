using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."))

if !(@isdefined PathTracking)
    include(joinpath(@__DIR__, "../implementation/PathTracking.jl"))
    using .PathTracking
end

using LinearAlgebra, DataFrames, CSV, Statistics, Printf, PrettyTables

RUN_EXPERIMENTS = true

include(joinpath(@__DIR__, "utils.jl"))

function steiner_system()
    x = [Variable(Symbol(:x, i)) for i = 1:2]
    a = [Variable(Symbol(:a, i)) for i = 1:5]
    c = [Variable(Symbol(:c, i)) for i = 1:6]
    y = [Variable(Symbol(:y, i, j)) for i = 1:2, j = 1:5]
    v = [Variable(Symbol(:v, i, j)) for i = 1:6, j = 1:5]

    #tangential conics
    f = a[1] * x[1]^2 + a[2] * x[1] * x[2] + a[3] * x[2]^2 + a[4] * x[1] + a[5] * x[2] + 1
    ∇ = differentiate(f, x)
    #5 conics
    g =
        c[1] * x[1]^2 +
        c[2] * x[1] * x[2] +
        c[3] * x[2]^2 +
        c[4] * x[1] +
        c[5] * x[2] +
        c[6]
    ∇_2 = differentiate(g, x)
    #the general system
    #f_a_0 is tangent to g_b₀ at x₀
    function Incidence(f, a₀, g, b₀, x₀)
        fᵢ = f(x => x₀, a => a₀)
        ∇ᵢ = [∇ᵢ(x => x₀, a => a₀) for ∇ᵢ in ∇]
        Cᵢ = g(x => x₀, c => b₀)
        ∇_Cᵢ = [∇ⱼ(x => x₀, c => b₀) for ∇ⱼ in ∇_2]

        [fᵢ; Cᵢ; det([∇ᵢ ∇_Cᵢ])]
    end
    F = vcat(map(i -> Incidence(f, a, g, v[:, i], y[:, i]), 1:5)...)

    System(F, [a; vec(y)], vec(v))
end

data_dir = joinpath(@__DIR__, "../data/steiner/")
result_dir = joinpath(@__DIR__, "../results/steiner/")
mkpath(result_dir)

N_J = 5

function run_hc_steiner()
    df = DataFrame(t = Float64[], count = Int[], duplicates = Int[])
    F = steiner_system()

    #warm up
    S = read_solution_file(data_dir * "solutions_1", 15)
    p = read_parameters_file(data_dir * "start_parameters_1")
    q = read_parameters_file(data_dir * "target_parameters_1_1")
    tracker = Tracker(ParameterHomotopy(F, p, q))
    track.(tracker, S[1:1], 1.0, 0.0)

    for i = 1:10
        S = read_solution_file(data_dir * "solutions_$i", 15)
        p = read_parameters_file(data_dir * "start_parameters_$i")
        for j = 1:N_J
            q = read_parameters_file(data_dir * "target_parameters_$(i)_$(j)")

            tracker = Tracker(ParameterHomotopy(F, p, q))
            t = @elapsed (res = track.(tracker, S, 1.0, 0.0))
            c = count(is_success, res)
            dups = nduplicates(res)
            # nsteps = steps.(res)
            @show t, c, dups
            push!(df, (t, c, dups))

            # write dataframe to file
            open(result_dir * "hc", "w") do io
                CSV.write(io, df; delim = " ")
            end
        end
    end
    return df
end

function run_bertini_steiner(;
    MPTYPE::Int = throw(UndefKeywordError(:MPTYPE)),
    track_tol = nothing,
)
    df = DataFrame(t = Float64[], count = Int[], duplicates = Int[])
    F = steiner_system()
    if MPTYPE == 0
        result_name = "bertini_dp"
    elseif MPTYPE == 2 && isnothing(track_tol)
        result_name = "bertini_ap"
    elseif MPTYPE == 2
        result_name = "bertini_ap_track_tol_8"
    end

    # create only one directory
    file_path = mktempdir()

    for i = 1:10
        S = read_solution_file(data_dir * "solutions_$i", 15)
        p = read_parameters_file(data_dir * "start_parameters_$i")
        for j = 1:N_J
            q = read_parameters_file(data_dir * "target_parameters_$(i)_$(j)")
            if !isnothing(track_tol)
                res = bertini(
                    F,
                    S,
                    start_parameters = p,
                    final_parameters = q,
                    MPTYPE = MPTYPE,
                    file_path = file_path,
                    TRACKTOLBEFOREEG = track_tol,
                    TRACKTOLDURINGEG = track_tol,
                )
            else
                res = bertini(
                    F,
                    S,
                    start_parameters = p,
                    final_parameters = q,
                    MPTYPE = MPTYPE,
                    file_path = file_path
                )
            end
            t = res[:runtime]
            c = length(res[:finite_solutions])
            dups = 0
            @show t, c, dups
            push!(df, (t, c, dups))
            # write dataframe to file
            open(result_dir * result_name, "w") do io
                CSV.write(io, df; delim = " ")
            end
        end
    end
    return df
end

if RUN_EXPERIMENTS
    run_hc_steiner()
    run_bertini_steiner(MPTYPE = 0)
    run_bertini_steiner(MPTYPE = 2)
    run_bertini_steiner(MPTYPE = 2, track_tol = 1e-8)
end

hc_df = CSV.read(result_dir * "hc")
bertini_dp_df = CSV.read(result_dir * "bertini_dp")
bertini_ap_df = CSV.read(result_dir * "bertini_ap")
bertini_ap_track_tol_df = CSV.read(result_dir * "bertini_ap_track_tol")

function table_data(df)
    reshape(
        [
            @sprintf("%.2f", mean(df.t)),
            @sprintf("%.2f", median(df.t)),
            @sprintf("%.2f", minimum(df.t)),
            @sprintf("%.2f", maximum(df.t)),
            "",
            @sprintf("%.1f", mean(df.count)),
            @sprintf("%.0f", median(df.count)),
            @sprintf("%.0f", minimum(df.count)),
            @sprintf("%.0f", maximum(df.count)),
        ],
        1,
        9,
    )
end


header = [
    "" "t" "t" "t" "t" "" "count" "count" "count" "count"
    "" "mean" "median" "min" "max" "" "mean" "median" "min" "max"
]

hc = table_data(hc_df)
bertini_dp = table_data(bertini_dp_df)
bertini_ap = table_data(bertini_ap_df)
bertini_ap_track_tol  = table_data(bertini_ap_track_tol_df)
cat = ["\\texttt{HC.jl}", "\\texttt{Bertini DP}", "\\texttt{Bertini AP}", "\\texttt{Bertini AP (tol)}"]

data = [cat [hc; bertini_dp; bertini_ap, bertini_ap_track_tol]]


open(result_dir * "table.tex", "w") do io
    pretty_table(io, data, header, backend = :latex)
end

open(result_dir * "table.txt", "w") do io
    pretty_table(io, data, header, backend = :text)
end
