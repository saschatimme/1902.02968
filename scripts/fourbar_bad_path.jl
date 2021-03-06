using Pkg;
Pkg.activate(joinpath(@__DIR__))

include(joinpath(@__DIR__, "utils.jl"))

if !(@isdefined PathTracking)
    include(joinpath(@__DIR__, "../implementation/PathTracking.jl"))
    using .PathTracking
end

using LinearAlgebra

F = let
    @var x a y b x_hat a_hat y_hat b_hat
    gamma = [Variable(Symbol(:gamma, i)) for i = 1:8]
    delta = [Variable(Symbol(:delta, i)) for i = 1:8]
    gamma_hat = [Variable(Symbol(:gamma_hat, i)) for i = 1:8]
    delta_hat = [Variable(Symbol(:delta_hat, i)) for i = 1:8]
    #system of polynomials
    D1 = [(a_hat * x - delta_hat[i] * x) * gamma[i] +
          (a * x_hat - delta[i] * x_hat) * gamma_hat[i] +
          (a_hat - x_hat) * delta[i] + (a - x) * delta_hat[i] -
          delta[i] * delta_hat[i] for i = 1:8]
    D2 = [(b_hat * y - delta_hat[i] * y) * gamma[i] +
          (b * y_hat - delta[i] * y_hat) * gamma_hat[i] +
          (b_hat - y_hat) * delta[i] + (b - y) * delta_hat[i] -
          delta[i] * delta_hat[i] for i = 1:8]
    D3 = [gamma[i] * gamma_hat[i] + gamma[i] + gamma_hat[i] for i = 1:8]
    System(
        [D1; D2; D3],
        [x; a; y; b; x_hat; a_hat; y_hat; b_hat; gamma; gamma_hat],
        [delta; delta_hat],
    )
end

s = Complex{Float64}[
    -2.6370887+0.35329306im,
    -3.7223306+2.5267087im,
    -0.70272787-0.74162847im,
    -0.52112188-1.5249im,
    -0.29444549+0.2698171im,
    0.23877446+1.1593521im,
    -0.2226372+0.24523739im,
    -0.23666974+0.87507397im,
    2.5277031-6.3353057im,
    -0.88010494+1.3133313im,
    -0.73606124+0.15544107im,
    0.80744564-3.0613903im,
    -0.82122921+0.80545901im,
    0.68534751-2.5000172im,
    0.65872999-2.3719828im,
    -0.87767826-0.35435132im,
    -0.93290889+0.12048708im,
    -0.93106365-0.75512927im,
    1.8130785-1.6567022im,
    -0.85699423+0.24221833im,
    -0.73738109-1.1832401im,
    -0.81460307+0.2750148im,
    -0.80200623+0.28313097im,
    -0.12955283+2.5215805im,
]

s_q = Complex{Float64}[
    -0.11052551340465448+0.1689556312781767im,
    0.919958012855403+0.24366762989216786im,
    0.901696352903974-1.6064926647211917im,
    -8.338280112808492-8.919160579975916im,
    0.006804245219856744+0.6666626834147511im,
    -0.010616472136481038+0.6647920903228585im,
    0.0054175505582474415+0.6594332738427415im,
    -0.06615618428899468+0.7160098329660299im,
    -0.6653465961148477+0.05690565132133756im,
    -0.8365994540764602+1.2293643774139256im,
    -1.165781038353456-0.40227703934042947im,
    -0.636262867891756-7.511530491963335im,
    -0.1555564468886374-4.36664898084851im,
    2.1006494673301965-5.4340054725677485im,
    -0.03216014143666483-1.1684119692446884im,
    -0.7479475298436062-0.1751684440725321im,
    1.9041920748392966-0.4938391173146806im,
    -0.8937602945251837-0.7993076684015589im,
    -1.8757112160341456+2.124962654367882im,
    -0.9935684803917256+0.13281722261129017im,
    -0.9573097412376701+0.22075291383283657im,
    -0.9207854885450023+0.1388264275883713im,
    -0.5795474872001224+0.5075857788948817im,
    1.675306384751507+1.859252783928112im,
]

p = Complex{Float64}[
    7.297148065259444-3.7969299471707454im,
    0.06125202773860465+0.9251492821389162im,
    1.4586512320807385+0.641527879044937im,
    3.3162292880980537-2.979942348702623im,
    -0.6613217458488018+0.8321586519581641im,
    2.6212360706406703-2.6766427578876666im,
    2.4522388848243555-2.5938408784246976im,
    -1.7265529790073868-0.29337325907998707im,
    0.24011427062306998+1.375412535460542im,
    -0.19336068638852033+0.5669754468190339im,
    0.147943349839693-0.2784406507316879im,
    -0.01641513182688631+1.5963735195205522im,
    -0.33880136827427765+0.2772103200265199im,
    -0.19788506286822416+1.6632361622942196im,
    -0.26097545384785903+1.6764055030939693im,
    0.44012684131549085+1.0212717350758722im,
]
q = Complex{Float64}[
    -0.35263257737537096+1.2598825875492277im,
    -0.5353479512983561-0.5963278468670689im,
    1.4096839876427854-1.2227632666172732im,
    -0.2550909822371234-0.641270846350799im,
    0.24180730168363096-0.46611650665985527im,
    0.29953332004980354-0.8469219548644632im,
    0.35386303195558433-0.3683866291866925im,
    1.2021910424897007+1.1148794002014533im,
    -0.35263257737537096-1.2598825875492277im,
    -0.5353479512983561+0.5963278468670689im,
    1.4096839876427854+1.2227632666172732im,
    -0.2550909822371234+0.641270846350799im,
    0.24180730168363096+0.46611650665985527im,
    0.29953332004980354+0.8469219548644632im,
    0.35386303195558433+0.3683866291866925im,
    1.2021910424897007-1.1148794002014533im,
]

tracker = Tracker(ParameterHomotopy(F, p, q))
r_q = track(tracker, s, 1, 0)

t_s = []
ω_s = []
μ_s = []
kappa_s = []
kappa_x_s = []
cond_s = []
cond_x_s = []

A = tracker.state.jacobian.J.A
x = tracker.state.x

PathTracking.init!(tracker, r_q, 0, 1)
while PathTracking.is_tracking(tracker.state.condition)

    push!(kappa_s, cond(A))
    push!(kappa_x_s, cond(A * diagm(0 => abs.(x))))
    push!(cond_s, condskeel(A))
    # push!(cond_x_s, condskeel(A * diagm(0 => abs.(x))))
    # There is a bug in Julia , 1.5 where condskeel(A,x) is actually cond(A, x) * norm(x, Inf)
    push!(cond_x_s, condskeel(A, x) / norm(x, Inf))

    push!(t_s, abs(tracker.state.t))
    push!(ω_s, tracker.state.ω)
    push!(μ_s, tracker.state.μ)

    PathTracking.step!(tracker)
end

@show maximum(kappa_s)
@show maximum(cond_s)
@show maximum(cond_x_s)

using Plots
pgfplotsx()

k = findfirst(t -> t > 0.5, t_s)
k = 1

ω_plot = begin
    plot(xlabel = "\$t\$", grid = false)
    plot!(
        t_s,
        ω_s,
        yscale = :log10,
        linestyle = :solid,
        label = "\\omega",
        color = :black,
        # legend = :topleft,
    )
    plot!(
        t_s,
        ω_s .* μ_s,
        yscale = :log10,
        linestyle = :dot,
        label = "\\omega \\mu",
        color = :DODGERBLUE,
    )

    plot!(
        t_s,
        μ_s,
        yscale = :log10,
        linestyle = :dash,
        label = "\\mu",
        color = :black,
    )

end
cond_plot = begin
    plot(xlabel = "\$t\$", grid = false)

    plot!(
        t_s[k:end],
        kappa_s[k:end],
        yscale = :log10,
        label = "\\kappa (J)",
        linestyle = :solid,
        # legend = :topleft,
        color = :black,
    )
    plot!(
        t_s[k:end],
        cond_s[k:end],
        yscale = :log10,
        label = "cond(J)",
        linestyle = :dash,
        color = :black,
    )
    plot!(
        t_s[k:end],
        cond_x_s[k:end],
        yscale = :log10,
        label = "cond(J,x)",
        linestyle = :dot,
        color = :DODGERBLUE,
    )

    plot!(
        [t_s[k], 1],
        [maximum(cond_s), maximum(cond_s)],
        linestyle = :dash,
        label = "",
        color = :BLACK,
    )
    plot!(
        [t_s[k], 1],
        [maximum(cond_x_s), maximum(cond_x_s)],
        linestyle = :dot,
        label = "",
        color = :DODGERBLUE,
    )
    plot!(
        [t_s[k], 1],
        [maximum(kappa_s), maximum(kappa_s)],
        linestyle = :solid,
        label = "",
        # alpha = 0.5,
        color = :BLACK,
    )
end

final_plt = plot(ω_plot, cond_plot, layout = (2, 1), size = (800, 540), dpi=300)
savefig(final_plt, joinpath(@__DIR__, "../../latex/figures/fourbar-bad-path.png"))

# run bertini

bertini(F, [s], start_parameters = p, final_parameters = q)
