using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."))

if !(@isdefined PathTracking)
    include(joinpath(@__DIR__, "../implementation/PathTracking.jl"))
    using .PathTracking
end

using LinearAlgebra, DataFrames, CSV, Statistics, Printf, PrettyTables

RUN_EXPERIMENTS = true

data_dir = joinpath(@__DIR__, "../data/fourbar/")
result_dir = joinpath(@__DIR__, "../results/fourbar/")
mkpath(result_dir)

include(joinpath(@__DIR__, "utils.jl"))

function fourbar_system()
    @var x a y b x_hat a_hat y_hat b_hat
    gamma = [Variable(Symbol(:gamma, i)) for i = 1:8]
    delta = [Variable(Symbol(:delta, i)) for i = 1:8]
    gamma_hat = [Variable(Symbol(:gamma_hat, i)) for i = 1:8]
    delta_hat = [Variable(Symbol(:delta_hat, i)) for i = 1:8]
    #system of polynomials
    D1 = [
        (a_hat * x - delta_hat[i] * x) * gamma[i] +
        (a * x_hat - delta[i] * x_hat) * gamma_hat[i] +
        (a_hat - x_hat) * delta[i] +
        (a - x) * delta_hat[i] - delta[i] * delta_hat[i] for i = 1:8
    ]
    D2 = [
        (b_hat * y - delta_hat[i] * y) * gamma[i] +
        (b * y_hat - delta[i] * y_hat) * gamma_hat[i] +
        (b_hat - y_hat) * delta[i] +
        (b - y) * delta_hat[i] - delta[i] * delta_hat[i] for i = 1:8
    ]
    D3 = [gamma[i] * gamma_hat[i] + gamma[i] + gamma_hat[i] for i = 1:8]
    System(
        [D1; D2; D3],
        [x; a; y; b; x_hat; a_hat; y_hat; b_hat; gamma; gamma_hat],
        [delta; delta_hat],
    )
end



function run_hc_fourbar()
    df = DataFrame(t = Float64[], count = Int[], duplicates = Int[])
    F = fourbar_system()

    #warm up
    S = read_solution_file(data_dir * "solutions_p", 24)
    p = read_parameters_file(data_dir * "p")
    q = read_parameters_file(data_dir * "q_1")
    tracker = Tracker(ParameterHomotopy(F, p, q))
    track.(tracker, S[1:1], 1.0, 0.0)

    for i = 1:10
        q = read_parameters_file(data_dir * "q_$i")
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
    return df
end

function run_bertini_fourbar(;
    MPTYPE::Int = throw(UndefKeywordError(:MPTYPE)),
    track_tol = nothing,
)
    df = DataFrame(t = Float64[], count = Int[], duplicates = Int[])
    F = fourbar_system()
    if MPTYPE == 0
        result_name = "bertini_dp"
    elseif MPTYPE == 2 && isnothing(track_tol)
        result_name = "bertini_ap"
    elseif MPTYPE == 2
        result_name = "bertini_ap_track_tol"
    end
    # create only one directory
    file_path = mktempdir()

    S = read_solution_file(data_dir * "solutions_p", 24)
    p = read_parameters_file(data_dir * "p")
    for i = 1:10
        q = read_parameters_file(data_dir * "q_$i")
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
                MAXNORM = 1e16,
            )
        else
            res = bertini(
                F,
                S,
                start_parameters = p,
                final_parameters = q,
                MPTYPE = MPTYPE,
                file_path = file_path,
                MAXNORM = 1e16,
            )
        end
        t = res[:runtime]
        c = length(res[:finite_solutions])
        dups = 0
        @show t, c, dups
        push!(df, (t, c, dups))
        # write dataframe to file
        open(result_dir * result_name, "w+") do io
            CSV.write(io, df; delim = " ")
        end
    end
    return df
end

if RUN_EXPERIMENTS
    run_hc_fourbar()
    run_bertini_fourbar(MPTYPE = 2)
    run_bertini_fourbar(MPTYPE = 2, track_tol = 1e-12)
else
    @warn "Experiments no run"
end



hc_df = CSV.read(result_dir * "hc")
bertini_ap_df = CSV.read(result_dir * "bertini_ap")
bertini_ap_df_track_tol = CSV.read(result_dir * "bertini_ap_track_tol")

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
bertini_ap = table_data(bertini_ap_df)
bertini_ap_track_tol = table_data(bertini_ap_df_track_tol)
cat = ["\\texttt{HC.jl}", "\\texttt{Bertini AP}", "\\texttt{Bertini AP (tol)}"]

data = [cat [hc; bertini_ap; bertini_ap_track_tol]]


open(result_dir * "table.tex", "w") do io
    pretty_table(io, data, header, backend = :latex)
end

open(result_dir * "table.txt", "w") do io
    pretty_table(io, data, header, backend = :text)
end
