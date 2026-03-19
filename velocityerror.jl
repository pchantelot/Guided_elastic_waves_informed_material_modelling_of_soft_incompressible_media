# %% Load packages and files
using CairoMakie, ColorSchemes, LaTeXStrings, MathTeXEngine
using MAT, FileIO, UnPack, NaturalSort
using LsqFit, Dierckx
using POnG
include("preamble.jl")
include("helper_func.jl")

#%% Load parameters
data = load(joinpath(@__DIR__, "parameters.jld2"))
@unpack ρs, VL, τ, n, μA, λs, λA, μP, λP = data
# Load fit parameters
data = load(joinpath(@__DIR__, "fitparameters.jld2"))
@unpack MRfit, βMR, GTfit, βGT, Cfit, βC, MRSHfit, βMRSH, GTSHfit, βGTSH, CSHfit, βCSH,
GMRfit, βGMR, GGfit, βGG, DCMRfit, βDCMR, DCGTfit, βDCGT = data
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
VApar, VAperp = getA0velocities(A0path, A0files, fop) ./ sqrt(μP / ρs)
# Velocities of the SH0 and S0 mode in parallel and perpendicular directions
VSHpar, VSpar, VSHperp, VSperp = getSHSvelocities(datapar, dataperp, fip) ./ sqrt(μA / ρs)

# %% Estimate relative error on velocity in the small to moderate deformation regime
λc = 1.5
idxA = λA .< λc
idxP = λP .< λc
VAparc = VApar[idxP]
while length(VAparc) < length(VSpar[idxA])
    append!(VAparc, 1)
end
VAperc = VAperp[idxP]
while length(VAperc) < length(VSpar[idxA])
    append!(VAperc, 1)
end
target1 = hcat(VSHpar[idxA], VSpar[idxA], VSHperp[idxA], VSperp[idxA], VAparc, VAperc)
#%% Mooney-Rivlin
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
splMR = Spline1D(λS, VSMR; k=1)
function fMR(λ1, λ2, p, splMR)
    CMRip = C_MR.(λ1, λ1 .^ (-0.41), ρs, VL, MRfit[1], MRfit[2], fip, τ, n, p[1])
    VSHpart, VSHperpt, VSpart, VSperpt = lw_mode_velocity(λ1, ρs, CMRip)
    CMRop = C_MR.(λ2, λ2 .^ (-0.39), ρs, VL, MRfit[1], MRfit[2], fop, τ, n, p[1])
    _, _, _, _, VApart = lw_mode_velocity(λ2, ρs, CMRop)
    while length(VApart) < length(VSpar[idxA])
        append!(VApart, sqrt(sum(MRfit) / ρs))
    end
    VAperpt = evaluate(splMR, λP[idxP])
    while length(VAperpt) < length(VSpar[idxA])
        append!(VAperpt, sqrt(sum(MRfit) / ρs))
    end
    return hcat(real(VSHpart), real(VSpart), real(VSHperpt), real(VSperpt), real(VApart), VAperpt) ./ sqrt(sum(MRfit) / ρs)
end
errorMR = findmax(abs.(fMR(λA[idxA], λP[idxP], βMR[1], splMR) ./ target1 .- 1) * 100, dims=1)

#%% Gent-Thomas
simu = filter(x -> occursin("GT", x), readdir(joinpath(@__DIR__, "fig4")))
VSGT = Vector{Float64}()
for name in simu
    local data = load(joinpath(@__DIR__, "fig4", name))
    @unpack k, f = data
    _, idx = findmin(abs.(f .- fop))
    append!(VSGT, real(2 * pi * f[idx] / k[idx]))
end
splGT = Spline1D(λS, VSGT; k=1)
function fGT(λ1, λ2, p, splGT)
    CGTip = C_GT.(λ1, λ1 .^ (-0.41), ρs, VL, GTfit[1], GTfit[2], fip, τ, n, p[1])
    VSHpart, VSHperpt, VSpart, VSperpt = lw_mode_velocity(λ1, ρs, CGTip)
    CGTop = C_GT.(λ2, λ2 .^ (-0.39), ρs, VL, GTfit[1], GTfit[2], fop, τ, n, p[1])
    _, _, _, _, VApart = lw_mode_velocity(λ2, ρs, CGTop)
    while length(VApart) < length(VSpar[idxA])
        append!(VApart, sqrt(sum(GTfit) / ρs))
    end
    VAperpt = evaluate(splGT, λP[idxP])
    while length(VAperpt) < length(VSpar[idxA])
        append!(VAperpt, sqrt(sum(GTfit) / ρs))
    end
    return hcat(real(VSHpart), real(VSpart), real(VSHperpt), real(VSperpt), real(VApart), VAperpt) ./ sqrt(sum(GTfit) / ρs)
end
errorGT = findmax(abs.(fGT(λA[idxA], λP[idxP], βGT[1], splGT) ./ target1 .- 1) * 100, dims=1)

#%% Carroll
simu = filter(x -> occursin("C", x), readdir(joinpath(@__DIR__, "fig4")))
VSC = Vector{Float64}()
for name in simu
    local data = load(joinpath(@__DIR__, "fig4", name))
    @unpack k, f = data
    _, idx = findmin(abs.(f .- fop))
    append!(VSC, real(2 * pi * f[idx] / k[idx]))
end
splC = Spline1D(λS, VSC; k=1)
function fC(λ1, λ2, p, spl)
    CCip = C_C.(λ1, λ1 .^ (-0.41), ρs, VL, Cfit[1], Cfit[2], fip, τ, n, p[1])
    VSHpart, VSHperpt, VSpart, VSperpt = lw_mode_velocity(λ1, ρs, CCip)
    CCop = C_C.(λ2, λ2 .^ (-0.39), ρs, VL, Cfit[1], Cfit[2], fop, τ, n, p[1])
    _, _, _, _, VApart = lw_mode_velocity(λ2, ρs, CCop)
    while length(VApart) < length(VSpar[idxA])
        append!(VApart, sqrt(sum(Cfit) / ρs))
    end
    VAperpt = evaluate(spl, λP[idxP])
    while length(VAperpt) < length(VSpar[idxA])
        append!(VAperpt, sqrt(sum(Cfit) / ρs))
    end
    return hcat(real(VSHpart), real(VSpart), real(VSHperpt), real(VSperpt), real(VApart), VAperpt) ./ sqrt(sum(Cfit) / ρs)
end
errorC = findmax(abs.(fC(λA[idxA], λP[idxP], βC[1], splC) ./ target1 .- 1) * 100, dims=1)

#%% Regime strain hardening
while length(VApar) < length(VSpar)
    append!(VApar, 1)
end
while length(VAperp) < length(VSpar)
    append!(VAperp, 1)
end
target2 = hcat(VSHpar, VSpar, VSHperp, VSperp, VApar, VAperp)
#%% MRSH I1
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
splMRSH = Spline1D(λS, VSMRSH1; k=1)
function fMRSH(λ1, λ2, p, spl)
    CMRSHip = C_MRSH.(λ1, λ1 .^ (-0.41), ρs, VL, MRSHfit[1], MRSHfit[2], MRSHfit[3], MRSHfit[4], fip, τ, n, p[1])
    VSHpart, VSHperpt, VSpart, VSperpt, _, _ = lw_mode_velocity(λ1, ρs, CMRSHip)
    CMRSHop = C_MRSH.(λ2, λ2 .^ (-0.39), ρs, VL, MRSHfit[1], MRSHfit[2], MRSHfit[3], MRSHfit[4], fop, τ, n, p[1])
    _, _, _, _, VApart, = lw_mode_velocity(λ2, ρs, CMRSHop)
    while length(VApart) < length(VSpar)
        append!(VApart, sqrt(sum(MRSHfit[1:3]) / ρs))
    end
    VAperpt = evaluate(spl, λP)
    while length(VAperpt) < length(VSpar)
        append!(VAperpt, sqrt(sum(MRSHfit[1:3]) / ρs))
    end
    return hcat(real(VSHpart), real(VSpart), real(VSHperpt), real(VSperpt), real(VApart), VAperpt) ./ sqrt(sum(MRSHfit[1:3]) / ρs)
end
errorMRSH = findmax(abs.(fMRSH(λA, λP, βMRSH[1], splMRSH) ./ target2 .- 1) * 100, dims=1)

#%% GTSH
simu = filter(x -> occursin("GTSH1", x), readdir(joinpath(@__DIR__, "fig5")))
VSGTSH = Vector{Float64}()
for name in simu
    data = load(joinpath(@__DIR__, "fig5", name))
    @unpack k, f = data
    _, idx = findmin(abs.(f .- fop))
    append!(VSGTSH, real(2 * pi * f[idx] / k[idx]))
end
splGTSH = Spline1D(λS, VSGTSH; k=1)
function fGTSH(λ1, λ2, p, spl)
    CGTSHip = C_GTSH.(λ1, λ1 .^ (-0.41), ρs, VL, GTSHfit[1], GTSHfit[2], GTSHfit[3], GTSHfit[4], fip, τ, n, p[1])
    VSHpart, VSHperpt, VSpart, VSperpt, _, _ = lw_mode_velocity(λ1, ρs, CGTSHip)
    CGTSHop = C_GTSH.(λ2, λ2 .^ (-0.39), ρs, VL, GTSHfit[1], GTSHfit[2], GTSHfit[3], GTSHfit[4], fop, τ, n, p[1])
    _, _, _, _, VApart, = lw_mode_velocity(λ2, ρs, CGTSHop)
    while length(VApart) < length(VSpar)
        append!(VApart, sqrt(sum(GTSHfit[1:3]) / ρs))
    end
    VAperpt = evaluate(spl, λP)
    while length(VAperpt) < length(VSpar)
        append!(VAperpt, sqrt(sum(GTSHfit[1:3]) / ρs))
    end
    return hcat(real(VSHpart), real(VSpart), real(VSHperpt), real(VSperpt), real(VApart), VAperpt) ./ sqrt(sum(GTSHfit[1:3]) / ρs)
end
errorGTSH = findmax(abs.(fGTSH(λA, λP, βGTSH[1], splGTSH) ./ target2 .- 1) * 100, dims=1)

#%% CSH
simu = filter(x -> occursin("CSH", x), readdir(joinpath(@__DIR__, "fig5")))
VSCSH = Vector{Float64}()
for name in simu
    data = load(joinpath(@__DIR__, "fig5", name))
    @unpack k, f = data
    _, idx = findmin(abs.(f .- fop))
    append!(VSCSH, real(2 * pi * f[idx] / k[idx]))
end
splCSH = Spline1D(λS, VSCSH; k=1)
function fCSH(λ1, λ2, p, spl)
    CCSHip = C_CSH.(λ1, λ1 .^ (-0.41), ρs, VL, CSHfit[1], CSHfit[2], CSHfit[3], CSHfit[4], fip, τ, n, p[1])
    VSHpart, VSHperpt, VSpart, VSperpt, _, _ = lw_mode_velocity(λ1, ρs, CCSHip)
    CCSHop = C_CSH.(λ2, λ2 .^ (-0.39), ρs, VL, CSHfit[1], CSHfit[2], CSHfit[3], CSHfit[4], fop, τ, n, p[1])
    _, _, _, _, VApart, = lw_mode_velocity(λ2, ρs, CCSHop)
    while length(VApart) < length(VSpar)
        append!(VApart, sqrt(sum(CSHfit[1:3]) / ρs))
    end
    VAperpt = evaluate(spl, λP)
    while length(VAperpt) < length(VSpar)
        append!(VAperpt, sqrt(sum(CSHfit[1:3]) / ρs))
    end
    return hcat(real(VSHpart), real(VSpart), real(VSHperpt), real(VSperpt), real(VApart), VAperpt) ./ sqrt(sum(CSHfit[1:3]) / ρs)
end
errorCSH = findmax(abs.(fCSH(λA, λP, βCSH[1], splCSH) ./ target2 .- 1) * 100, dims=1)

#%% relative error
with_theme(My_theme) do
    fig = Figure(size=(1.2 * 360, 360))
    ax = Axis(fig[1, 1])
    ax.ylabel = L"$\mathcal{E}$ (%)"
    ax.xticks = (1:6, [L"V_{SH_0,∥}", L"V_{S_0,∥}", L"V_{SH_0,⟂}", L"V_{S_0,⟂}", L"V_{A_0,∥}", L"V_{A_0,⟂}"])
    lines!(ax, vec(errorMR[1]))
    lines!(ax, vec(errorGT[1]))
    lines!(ax, vec(errorC[1]))
    display(fig)
end
with_theme(My_theme) do
    fig = Figure(size=(1.2 * 360, 360))
    ax = Axis(fig[1, 1])
    ax.ylabel = L"$\mathcal{E}$ (%)"
    ax.xticks = (1:6, [L"V_{SH_0,∥}", L"V_{S_0,∥}", L"V_{SH_0,⟂}", L"V_{S_0,⟂}", L"V_{A_0,∥}", L"V_{A_0,⟂}"])
    lines!(ax, vec(errorMRSH[1]))
    lines!(ax, vec(errorGTSH[1]))
    lines!(ax, vec(errorCSH[1]))
    display(fig)
end

