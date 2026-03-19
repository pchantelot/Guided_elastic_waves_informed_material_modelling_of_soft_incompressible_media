# %% Load packages and files
using CairoMakie, ColorSchemes, LaTeXStrings, MathTeXEngine
using MAT, FileIO, UnPack, NaturalSort
using LsqFit
using POnG
include("preamble.jl")
include("helper_func.jl")

#%% Load parameters
data = load(joinpath(@__DIR__, "parameters.jld2"))
@unpack ρs, VL, β, τ, n, μA, λs, λA, μP, λP = data
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

# %% 
# Load fit parameters
data = load(joinpath(@__DIR__, "fitparameters.jld2"))
@unpack MRfit, βMR, MRSHfit, βMRSH, GTSHfit, βGTSH, CSHfit, βCSH = data
# This figure requires to compute the maximum relative error
# Load traction test data
data = matread(joinpath(@__DIR__, "data_Alex/Ecoflex30_RubanMax.mat"))
lbd = vec(data["extension"] / 200 .+ 1)
# engineering stress
sig = vec(data["load"] / (40e-3 * 3e-3))
# Mooney-plot
MP = sig ./ (2 * (lbd .- lbd .^ (-2)))
st = 2000

# %% Get simulation results
simu = filter(x -> occursin("MRSH1", x), readdir(joinpath(@__DIR__, "fig5")))
VSMRSH1 = Vector{Float64}()
λS = Vector{Float32}()
for name in simu
    data = load(joinpath(@__DIR__, "fig5", name))
    @unpack k, f = data
    _, idx = findmin(abs.(f .- fop))
    append!(VSMRSH1, real(2 * pi * f[idx] / k[idx]))
    append!(λS, parse(Float32, name[end-7:end-5]))
end

simu = filter(x -> occursin("GTSH1", x), readdir(joinpath(@__DIR__, "fig5")))
VSGTSH = Vector{Float64}()
for name in simu
    data = load(joinpath(@__DIR__, "fig5", name))
    @unpack k, f = data
    _, idx = findmin(abs.(f .- fop))
    append!(VSGTSH, real(2 * pi * f[idx] / k[idx]))
end

simu = filter(x -> occursin("CSH", x), readdir(joinpath(@__DIR__, "fig5")))
VSCSH = Vector{Float64}()
for name in simu
    data = load(joinpath(@__DIR__, "fig5", name))
    @unpack k, f = data
    _, idx = findmin(abs.(f .- fop))
    append!(VSCSH, real(2 * pi * f[idx] / k[idx]))
end

simu = filter(x -> occursin("MR", x), readdir(joinpath(@__DIR__, "fig4")))
VSMR = Vector{Float64}()
for name in simu
    data = load(joinpath(@__DIR__, "fig4", name))
    @unpack k, f = data
    _, idx = findmin(abs.(f .- fop))
    append!(VSMR, real(2 * pi * f[idx] / k[idx]))
end
# %% Figure
λ = 1:0.02:2.5
λc = 2.5
with_theme(My_theme, palette=(color=ColorSchemes.Set1_3, marker=[:circle])) do
    fig = Figure(size=(2 * 402, 4 / (3 * 1.2) * 402), figure_padding=(2, 12, 2, 2))

    ax1 = Axis(fig[1, 1])
    ax1.xlabel = L"1 / λ"
    ax1.ylabel = L"$\mathcal{M}$ (kPa)"
    ax1.limits = (0.3, 1, 10, 12)
    band!(ax1, [0, 1 / λc], [10, 10], [12, 12], color=tuple.(:grey, 0.2))
    scatter!(ax1, 1 ./ lbd[st:200:end], MP[st:200:end] / 1e3, color=:black, markersize=10)
    lines!(ax1, 1 ./ λ, MR(λ, MRfit), label="Mooney-Rivlin",
        color=:black, linestyle=:dash)
    lines!(ax1, 1 ./ λ, MRSH(λ, MRSHfit) / 1e3, label=L"Mooney-Rivlin with $I_1^N$ term")
    lines!(ax1, 1 ./ λ, GTSH(λ, GTSHfit) / 1e3, label=L"Gent-Thomas with $I_1^N$ term")
    lines!(ax1, 1 ./ λ, CSH(λ, CSHfit) / 1e3, label=L"Carroll with $I_1^N$ term")
    # lines!(ax1, 1 ./ λ, MRSH2(λ, MRSH2fit) / 1e3, linestyle=:dash, label=L"Mooney-Rivlin with $I_2^N$ term")
    # lines!(ax1, 1 ./ λ, GTSH2(λ, GTSH2fit) / 1e3, color=ColorSchemes.Set1_3[2], linestyle=:dash, label=L"Gent-Thomas with $I_2^N$ term")

    fig[2, 1] = Legend(fig, ax1, labelsize=16, tellwidth=false)

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
    VSHpart, VSHperpt, VSpart, VSperpt, VApart, VAperpt = lw_mode_velocity(λ, ρs, CMRip)
    lines!(ax2, λ, real(VSHpart) ./ sqrt(sum(MRfit) / ρs), color=:black, linestyle=:dash)
    lines!(ax2, λ, real(VSpart) ./ sqrt(sum(MRfit) / ρs), color=:black, linestyle=:dash)
    lines!(ax3, λ, real(VSHperpt) ./ sqrt(sum(MRfit) / ρs), color=:black, linestyle=:dash)
    lines!(ax3, λ, real(VSperpt) ./ sqrt(sum(MRfit) / ρs), color=:black, linestyle=:dash)
    CMRSHip = C_MRSH.(λ, λ .^ (-0.41), ρs, VL, MRSHfit[1], MRSHfit[2], MRSHfit[3], MRSHfit[4], fip, τ, n, βMRSH)
    VSHpart, VSHperpt, VSpart, VSperpt, VApart, VAperpt = lw_mode_velocity(λ, ρs, CMRSHip)
    lines!(ax2, λ, real(VSHpart) ./ sqrt(sum(MRSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[1])
    lines!(ax2, λ, real(VSpart) ./ sqrt(sum(MRSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[1])
    lines!(ax3, λ, real(VSHperpt) ./ sqrt(sum(MRSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[1])
    lines!(ax3, λ, real(VSperpt) ./ sqrt(sum(MRSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[1])
    CGTSHip = C_GTSH.(λ, λ .^ (-0.41), ρs, VL, GTSHfit[1], GTSHfit[2], GTSHfit[3], GTSHfit[4], fip, τ, n, βGTSH)
    VSHpart, VSHperpt, VSpart, VSperpt, VApart, VAperpt = lw_mode_velocity(λ, ρs, CGTSHip)
    lines!(ax2, λ, real(VSHpart) ./ sqrt(sum(GTSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[2])
    lines!(ax2, λ, real(VSpart) ./ sqrt(sum(GTSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[2])
    lines!(ax3, λ, real(VSHperpt) ./ sqrt(sum(GTSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[2])
    lines!(ax3, λ, real(VSperpt) ./ sqrt(sum(GTSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[2])
    CCSHip = C_CSH.(λ, λ .^ (-0.41), ρs, VL, CSHfit[1], CSHfit[2], CSHfit[3], CSHfit[4], fip, τ, n, βCSH)
    VSHpart, VSHperpt, VSpart, VSperpt, VApart, VAperpt = lw_mode_velocity(λ, ρs, CCSHip)
    lines!(ax2, λ, real(VSHpart) ./ sqrt(sum(CSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[3])
    lines!(ax2, λ, real(VSpart) ./ sqrt(sum(CSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[3])
    lines!(ax3, λ, real(VSHperpt) ./ sqrt(sum(CSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[3])
    lines!(ax3, λ, real(VSperpt) ./ sqrt(sum(CSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[3])
    # CMRSH2ip = C_MRSH2.(λ, λ .^ (-0.41), ρs, VL, MRSH2fit[1], MRSH2fit[2], MRSH2fit[3], MRSH2fit[4], fip, τ, n, βMRSH2)
    # VSHpart, VSHperpt, VSpart, VSperpt, VApart, VAperpt = lw_mode_velocity(λ, ρs, CMRSH2ip)
    # lines!(ax2, λ, real(VSHpart) ./ sqrt(sum(MRSH2fit[1:3]) / ρs), color=ColorSchemes.Set1_3[1], linestyle=:dash)
    # lines!(ax2, λ, real(VSpart) ./ sqrt(sum(MRSH2fit[1:3]) / ρs), color=ColorSchemes.Set1_3[1], linestyle=:dash)
    # lines!(ax3, λ, real(VSHperpt) ./ sqrt(sum(MRSH2fit[1:3]) / ρs), color=ColorSchemes.Set1_3[1], linestyle=:dash)
    # lines!(ax3, λ, real(VSperpt) ./ sqrt(sum(MRSH2fit[1:3]) / ρs), color=ColorSchemes.Set1_3[1], linestyle=:dash)
    # CGTSH2ip = C_GTSH2.(λ, λ .^ (-0.41), ρs, VL, GTSH2fit[1], GTSH2fit[2], GTSH2fit[3], GTSH2fit[4], fip, τ, n, βMRSH2)
    # VSHpart, VSHperpt, VSpart, VSperpt, VApart, VAperpt = lw_mode_velocity(λ, ρs, CGTSH2ip)
    # lines!(ax2, λ, real(VSHpart) ./ sqrt(sum(GTSH2fit[1:3]) / ρs), color=ColorSchemes.Set1_3[2], linestyle=:dash)
    # lines!(ax2, λ, real(VSpart) ./ sqrt(sum(GTSH2fit[1:3]) / ρs), color=ColorSchemes.Set1_3[2], linestyle=:dash)
    # lines!(ax3, λ, real(VSHperpt) ./ sqrt(sum(GTSH2fit[1:3]) / ρs), color=ColorSchemes.Set1_3[2], linestyle=:dash)
    # lines!(ax3, λ, real(VSperpt) ./ sqrt(sum(GTSH2fit[1:3]) / ρs), color=ColorSchemes.Set1_3[2], linestyle=:dash)

    CMRop = C_MR.(λ, λ .^ (-0.39), ρs, VL, MRfit[1], MRfit[2], fop, τ, n, βMR)
    VSHpart, VSHperpt, VSpart, VSperpt, VApart, VAperpt = lw_mode_velocity(λ, ρs, CMRop)
    lines!(ax22, λ, real(VApart) ./ sqrt(sum(MRfit) / ρs), color=:black, linestyle=:dash)
    # lines!(ax23, λ, real(VAperpt) ./ sqrt(sum(MRfit) / ρs), color=:black, linestyle=:dash)
    lines!(ax23, λS, real(VSMR) ./ sqrt(sum(MRfit) / ρs), color=:black, linestyle=:dash)
    CMRSHop = C_MRSH.(λ, λ .^ (-0.39), ρs, VL, MRSHfit[1], MRSHfit[2], MRSHfit[3], MRSHfit[4], fop, τ, n, βMRSH)
    VSHpart, VSHperpt, VSpart, VSperpt, VApart, VAperpt = lw_mode_velocity(λ, ρs, CMRSHop)
    lines!(ax22, λ, real(VApart) ./ sqrt(sum(MRSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[1])
    # lines!(ax23, λ, real(VAperpt) ./ sqrt(sum(MRSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[1], linestyle=:dash)
    lines!(ax23, λS, VSMRSH1 ./ sqrt(sum(MRSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[1])
    CGTSHop = C_GTSH.(λ, λ .^ (-0.39), ρs, VL, GTSHfit[1], GTSHfit[2], GTSHfit[3], GTSHfit[4], fop, τ, n, βGTSH)
    VSHpart, VSHperpt, VSpart, VSperpt, VApart, VAperpt = lw_mode_velocity(λ, ρs, CGTSHop)
    lines!(ax22, λ, real(VApart) ./ sqrt(sum(GTSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[2])
    lines!(ax23, λS, VSGTSH ./ sqrt(sum(GTSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[2])
    CCSHop = C_CSH.(λ, λ .^ (-0.39), ρs, VL, CSHfit[1], CSHfit[2], CSHfit[3], CSHfit[4], fop, τ, n, βCSH)
    VSHpart, VSHperpt, VSpart, VSperpt, VApart, VAperpt = lw_mode_velocity(λ, ρs, CCSHop)
    lines!(ax22, λ, real(VApart) ./ sqrt(sum(CSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[3])
    lines!(ax23, λS, VSCSH ./ sqrt(sum(CSHfit[1:3]) / ρs), color=ColorSchemes.Set1_3[3])
    # CGTSH2op = C_GTSH2.(λ, λ .^ (-0.39), ρs, VL, GTSH2fit[1], GTSH2fit[2], GTSH2fit[3], GTSH2fit[4], fop, τ, n, βGTSH2)
    # VSHpart, VSHperpt, VSpart, VSperpt, VApart, VAperpt = lw_mode_velocity(λ, ρs, CGTSH2op)
    # lines!(ax22, λ, real(VApart) ./ sqrt(sum(GTSH2fit[1:3]) / ρs), color=ColorSchemes.Set1_3[2], linestyle=:dash)
    # # lines!(ax23, λ, real(VAperpt) ./ sqrt(sum(MRSH2fit[1:3]) / ρs), color=ColorSchemes.Set1_3[2], linestyle=:dash)
    # lines!(ax23, λS, VSGTSH2 ./ sqrt(sum(GTSH2fit[1:3]) / ρs), color=ColorSchemes.Set1_3[2], linestyle=:dash)

    # add labels
    Label(fig[1, 1, TopLeft()], "(a)", padding=(0, 30, -30, 0), fontsize=16)
    Label(ip[1, 1, TopLeft()], "(b)", padding=(0, 30, -30, 0), fontsize=16)
    Label(ip[1, 2, TopLeft()], "(c)", padding=(0, 30, -30, 0), fontsize=16)
    Label(op[1, 1, TopLeft()], "(d)", padding=(0, 30, -30, 0), fontsize=16)
    Label(op[1, 2, TopLeft()], "(e)", padding=(0, 30, -30, 0), fontsize=16)

    colsize!(fig.layout, 2, Relative(7 / 10))
    rowgap!(fig.layout, 0)
    display(fig)
    save(joinpath(@__DIR__, "figure5.pdf"), fig; pt_per_unit=1.0)

end
