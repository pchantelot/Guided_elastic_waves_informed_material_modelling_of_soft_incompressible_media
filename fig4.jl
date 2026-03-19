# %% Load packages and files
using CairoMakie, ColorSchemes, LaTeXStrings, MathTeXEngine
using MAT, FileIO, UnPack, NaturalSort
using LsqFit
using POnG
include("preamble.jl")
include("helper_func.jl")

#%% Load parameters
data = load(joinpath(@__DIR__, "parameters.jld2"))
@unpack ρs, VL, τ, n, μA, λs, λA, μP, λP = data
# get filenames for my data
A0path = joinpath(@__DIR__, "data_pierre/")
A0files = filter(x -> occursin("def", x), readdir(A0path))
sort!(A0files, lt=natural)
# get Alex data
datapar = matread(joinpath(@__DIR__, "data_Alex/Alex_EML_par_disprel.mat"))
dataperp = matread(joinpath(@__DIR__, "data_Alex/Alex_EML_perp_disprel.mat"))
# Select frequencies in plane and out-of-plane
fip = 170
fop = 50
# Velocities of the A0 mode in parallel and perpendicular directions
VApar, VAperp = getA0velocities(A0path, A0files, fop)
# Velocities of the SH0 and S0 mode in parallel and perpendicular directions
VSHpar, VSpar, VSHperp, VSperp = getSHSvelocities(datapar, dataperp, fip)

# %% Define engineering stress for the W considered in the small to moderate strain regime
# Load fit parameters
data = load(joinpath(@__DIR__, "fitparameters.jld2"))
@unpack MRfit, GTfit, Cfit, βMR, βGT, βC = data
# This figure requires to compute the maximum relative error
# Load traction test data
data = matread(joinpath(@__DIR__, "data_Alex/Ecoflex30_RubanMax.mat"))
lbd = vec(data["extension"] / 200 .+ 1)
# engineering stress
sig = vec(data["load"] / (40e-3 * 3e-3))
# Mooney-plot
MP = sig ./ (2 * (lbd .- lbd .^ (-2)))
st = 2000
# Mooney-Rivlin Fit with varying end point
lbdmax = 1.25:0.01:2
errMR = Vector{Float64}()
for l in lbdmax
    fit = curve_fit(MR, 1 ./ lbd[lbd.<l][st:end], MP[lbd.<l][st:end], 1 ./ MP[lbd.<l][st:end] .^ 2, [10e3, 1e3])
    # define error as in Destrade
    append!(errMR, maximum(abs.(MR(1 ./ lbd[1.2 .< lbd .< l], fit.param) ./ MP[1.2 .< lbd .< l] .- 1)))
end
# Gent-Thomas Fit
errGT = Vector{Float64}()
for l in lbdmax
    fit = curve_fit(GT, 1 ./ (2 * lbd[lbd.<l][st:end] .^ 2 .+ 1 ./ lbd[lbd.<l][st:end]), MP[lbd.<l][st:end], 1 ./ MP[lbd.<l][st:end] .^ 2, [10e3, 1e3])
    append!(errGT, maximum(abs.(GT(1 ./ (2 * lbd[1.2 .< lbd .< l] .^ 2 .+ 1 ./ lbd[1.2 .< lbd .< l]),
        fit.param) ./ MP[1.2 .< lbd .< l] .- 1)))
end
# Carroll Fit
errC = Vector{Float64}()
for l in lbdmax
    fit = curve_fit(C, 1 ./ (lbd[lbd.<l][st:end] .* sqrt.(2 * lbd[lbd.<l][st:end] .+ 1 ./ lbd[lbd.<l][st:end] .^ 2)), MP[lbd.<l][st:end], 1 ./ MP[lbd.<l][st:end] .^ 2, [10e3, 1e3])
    append!(errC, maximum(abs.(C(1 ./ (lbd[1.2 .< lbd .< l] .* sqrt.(2 * lbd[1.2 .< lbd .< l] .+ 1 ./ lbd[1.2 .< lbd .< l] .^ 2)), fit.param) ./ MP[1.2 .< lbd .< l] .- 1)))
end

#%% Get velocity for the A0 perp mode which is not well described by a long-wavelength approximation
simu = filter(x -> occursin("MR", x), readdir(joinpath(@__DIR__, "fig4")))
VSMR = Vector{Float64}()
λS = Vector{Float32}()
for name in simu
    local data = load(joinpath(@__DIR__, "fig4", name))
    @unpack k, f = data
    _, idx = findmin(abs.(f .- fop))
    append!(VSMR, real(2 * pi * f[idx] / k[idx]))
    append!(λS, parse(Float32, name[end-7:end-5]))
end
simu = filter(x -> occursin("GT", x), readdir(joinpath(@__DIR__, "fig4")))
VSGT = Vector{Float64}()
for name in simu
    local data = load(joinpath(@__DIR__, "fig4", name))
    @unpack k, f = data
    _, idx = findmin(abs.(f .- fop))
    append!(VSGT, real(2 * pi * f[idx] / k[idx]))
end
simu = filter(x -> occursin("C", x), readdir(joinpath(@__DIR__, "fig4")))
VSC = Vector{Float64}()
for name in simu
    local data = load(joinpath(@__DIR__, "fig4", name))
    @unpack k, f = data
    _, idx = findmin(abs.(f .- fop))
    append!(VSC, real(2 * pi * f[idx] / k[idx]))
end

#%% Plot figure
λ = 1:0.02:2.5
λc = 1.5
with_theme(My_theme, palette=(color=ColorSchemes.Set1_3, marker=[:circle])) do

    fig = Figure(size=(2 * 402, 4 / (3 * 1.2) * 402), figure_padding=(2, 12, 2, 2))

    ax1 = Axis(fig[1, 1])
    ax1.xlabel = L"1 / λ"
    ax1.ylabel = L"$\mathcal{M}$ (kPa)"
    ax1.limits = (0.3, 1, 10, 12)
    band!(ax1, [0, 1 / λc], [10, 10], [12, 12], color=tuple.(:grey, 0.2))
    scatter!(ax1, 1 ./ lbd[st:200:end], MP[st:200:end] / 1e3, color=:black, markersize=10)
    lines!(ax1, 1 ./ λ, MR(1 ./ λ, MRfit) / 1e3, label="Mooney-Rivlin")
    lines!(ax1, 1 ./ λ, GT(1 ./ (2 * λ .^ 2 .+ 1 ./ λ), GTfit) / 1e3, label="Gent-Thomas")
    lines!(ax1, 1 ./ λ, C(1 ./ (λ .* sqrt.(2 * λ .+ 1 ./ λ .^ 2)), Cfit) / 1e3, label="Carroll")

    ax21 = Axis(fig[2, 1])
    ax21.xlabel = L"λ"
    ax21.ylabel = L"$\mathcal{E}$ (%)"
    ax21.limits = (1.25, 1.7, 0, 1)
    band!(ax21, [λc, 3], [0, 0], [4, 4], color=tuple.(:grey, 0.2))
    lines!(ax21, lbdmax, errMR * 100, label="Mooney-Rivlin")
    lines!(ax21, lbdmax, errGT * 100, label="Gent-Thomas")
    lines!(ax21, lbdmax, errC * 100, label="Carroll")
    axislegend(ax21, position=:lt)

    ip = fig[1, 2] = GridLayout()
    Label(ip[1, 1:2, Top()], L"$S_0$ and $SH_0$, $f = %$(fip)$ Hz", valign=:bottom,
        padding=(0, 0, 10, 0), fontsize=20)
    ax2 = Axis(ip[1, 1])
    ax2.xlabel = L"λ"
    ax2.ylabel = L"V_∥ / V_T"
    ax2.limits = (1, 2.5, 1, 4)
    band!(ax2, [λc, 3], [0, 0], [4, 4], color=tuple.(:grey, 0.2))
    scatter!(ax2, λA, VSHpar ./ sqrt(μA / ρs), color=:black)
    scatter!(ax2, λA, VSpar ./ sqrt(μA / ρs), color=:black, marker=:utriangle)

    ax3 = Axis(ip[1, 2])
    ax3.xlabel = L"λ"
    ax3.ylabel = L"V_⟂ / V_T"
    ax3.limits = (1, 2.5, 1, 3)
    band!(ax3, [λc, 3], [0, 0], [4, 4], color=tuple.(:grey, 0.2))
    scatter!(ax3, λA, VSHperp ./ sqrt(μA / ρs), color=:black)
    scatter!(ax3, λA, VSperp ./ sqrt(μA / ρs), color=:black, marker=:utriangle)

    op = fig[2, 2] = GridLayout()
    Label(op[1, 1:2, Top()], L"$A_0$, $f = %$(fop)$ Hz", valign=:bottom,
        padding=(0, 0, 10, 0), fontsize=20)
    ax22 = Axis(op[1, 1])
    ax22.xlabel = L"λ"
    ax22.ylabel = L"V_∥ / V_T"
    ax22.limits = (1, 2.5, 0, 2.2)
    band!(ax22, [λc, 3], [0, 0], [4, 4], color=tuple.(:grey, 0.2))
    scatter!(ax22, λP, VApar / sqrt(μP / ρs), color=:black)

    ax23 = Axis(op[1, 2])
    ax23.xlabel = L"λ"
    ax23.ylabel = L"V_⟂ / V_T"
    ax23.limits = (1, 2.5, 0, 1)
    band!(ax23, [λc, 3], [0, 0], [4, 4], color=tuple.(:grey, 0.2))
    scatter!(ax23, λP, VAperp / sqrt(μP / ρs), color=:black)

    # add predictions
    CMRip = C_MR.(λ, λ .^ (-0.41), ρs, VL, MRfit[1], MRfit[2], fip, τ, n, βMR)
    VSHpart, VSHperpt, VSpart, VSperpt = lw_mode_velocity(λ, ρs, CMRip)
    lines!(ax2, λ, real(VSHpart) ./ sqrt(sum(MRfit) / ρs), color=ColorSchemes.Set1_3[1])
    lines!(ax2, λ, real(VSpart) ./ sqrt(sum(MRfit) / ρs), color=ColorSchemes.Set1_3[1])
    lines!(ax3, λ, real(VSHperpt) ./ sqrt(sum(MRfit) / ρs), color=ColorSchemes.Set1_3[1])
    lines!(ax3, λ, real(VSperpt) ./ sqrt(sum(MRfit) / ρs), color=ColorSchemes.Set1_3[1])
    CGTip = C_GT.(λ, λ .^ (-0.41), ρs, VL, GTfit[1], GTfit[2], fip, τ, n, βGT)
    VSHpart, VSHperpt, VSpart, VSperpt = lw_mode_velocity(λ, ρs, CGTip)
    lines!(ax2, λ, real(VSHpart) ./ sqrt(sum(GTfit) / ρs), color=ColorSchemes.Set1_3[2])
    lines!(ax2, λ, real(VSpart) ./ sqrt(sum(GTfit) / ρs), color=ColorSchemes.Set1_3[2])
    lines!(ax3, λ, real(VSHperpt) ./ sqrt(sum(GTfit) / ρs), color=ColorSchemes.Set1_3[2])
    lines!(ax3, λ, real(VSperpt) ./ sqrt(sum(GTfit) / ρs), color=ColorSchemes.Set1_3[2])
    CCip = C_C.(λ, λ .^ (-0.41), ρs, VL, Cfit[1], Cfit[2], fip, τ, n, βC)
    VSHpart, VSHperpt, VSpart, VSperpt = lw_mode_velocity(λ, ρs, CCip)
    lines!(ax2, λ, real(VSHpart) ./ sqrt(sum(Cfit) / ρs), color=ColorSchemes.Set1_3[3])
    lines!(ax2, λ, real(VSpart) ./ sqrt(sum(Cfit) / ρs), color=ColorSchemes.Set1_3[3])
    lines!(ax3, λ, real(VSHperpt) ./ sqrt(sum(Cfit) / ρs), color=ColorSchemes.Set1_3[3])
    lines!(ax3, λ, real(VSperpt) ./ sqrt(sum(Cfit) / ρs), color=ColorSchemes.Set1_3[3])
    # CNHip = C_MR.(λ, λ .^ (-0.41), ρs, VL, MRfit[1] + MRfit[2], 0, fip, τ, n, βMR)
    # VSHpart, VSHperpt, VSpart, VSperpt = lw_mode_velocity(λ, ρs, CNHip)
    # lines!(ax2, λ, real(VSHpart) ./ sqrt(sum(MRfit) / ρs), color=:black)
    # lines!(ax2, λ, real(VSpart) ./ sqrt(sum(MRfit) / ρs), color=:black)
    # lines!(ax3, λ, real(VSHperpt) ./ sqrt(sum(MRfit) / ρs), color=:black)
    # lines!(ax3, λ, real(VSperpt) ./ sqrt(sum(MRfit) / ρs), color=:black)

    CMRop = C_MR.(λ, λ .^ (-0.39), ρs, VL, MRfit[1], MRfit[2], fop, τ, n, βMR)
    _, _, _, _, VApart, VAperpt = lw_mode_velocity(λ, ρs, CMRop)
    lines!(ax22, λ, real(VApart) ./ sqrt(sum(MRfit) / ρs), color=ColorSchemes.Set1_3[1])
    # lines!(ax23, λ, real(VAperpt) ./ sqrt(sum(MRfit.param) / ρs), color = ColorSchemes.Set1_3[1], linestyle = :dash)
    lines!(ax23, λS, VSMR ./ sqrt(sum(MRfit) / ρs), color=ColorSchemes.Set1_3[1])
    CGTop = C_GT.(λ, λ .^ (-0.39), ρs, VL, GTfit[1], GTfit[2], fop, τ, n, βGT)
    _, _, _, _, VApart, VAperpt = lw_mode_velocity(λ, ρs, CGTop)
    lines!(ax22, λ, real(VApart) ./ sqrt(sum(GTfit) / ρs), color=ColorSchemes.Set1_3[2])
    lines!(ax23, λS, VSGT ./ sqrt(sum(GTfit) / ρs), color=ColorSchemes.Set1_3[2])
    # lines!(ax23, λ, real(VAperpt) ./ sqrt(sum(GTfit.param) / ρs), color = ColorSchemes.Set1_3[2], linestyle = :dash)
    CCop = C_C.(λ, λ .^ (-0.39), ρs, VL, Cfit[1], Cfit[2], fop, τ, n, βC)
    _, _, _, _, VApart, VAperpt = lw_mode_velocity(λ, ρs, CCop)
    lines!(ax22, λ, real(VApart) ./ sqrt(sum(Cfit) / ρs), color=ColorSchemes.Set1_3[3])
    lines!(ax23, λS, VSC ./ sqrt(sum(Cfit) / ρs), color=ColorSchemes.Set1_3[3])
    # lines!(ax23, λ, real(VAperpt) ./ sqrt(sum(Cfit.param) / ρs), color = ColorSchemes.Set1_3[3], linestyle = :dash)
    # CNHop = C_MR.(λ, λ .^ (-0.39), ρs, VL, MRfit[1] + MRfit[2], 0, fop, τ, n, βMR)
    # _, _, _, _, VApart, VAperpt = lw_mode_velocity(λ, ρs, CNHop)
    # lines!(ax22, λ, real(VApart) ./ sqrt(sum(MRfit) / ρs), color=:black)
    # lines!(ax23, λ, real(VAperpt) ./ sqrt(sum(MRfit) / ρs), color=:black, linestyle=:dash)

    # add labels
    Label(fig[1, 1, TopLeft()], "(a)", padding=(0, 30, -30, 0), fontsize=16)
    Label(fig[2, 1, TopLeft()], "(b)", padding=(0, 30, -30, 0), fontsize=16)
    Label(ip[1, 1, TopLeft()], "(c)", padding=(0, 30, -30, 0), fontsize=16)
    Label(ip[1, 2, TopLeft()], "(d)", padding=(0, 30, -30, 0), fontsize=16)
    Label(op[1, 1, TopLeft()], "(e)", padding=(0, 30, -30, 0), fontsize=16)
    Label(op[1, 2, TopLeft()], "(f)", padding=(0, 30, -30, 0), fontsize=16)

    rowgap!(fig.layout, 0)
    colsize!(fig.layout, 2, Relative(7 / 10))
    display(fig)
    save(joinpath(@__DIR__, "figure4.pdf"), fig; pt_per_unit=1.0)

end


