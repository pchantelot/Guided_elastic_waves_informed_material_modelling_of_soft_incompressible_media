#%% Save all parameters useful to plot figures in one place 
using MAT, FileIO, UnPack, NaturalSort
using LsqFit, CMAEvolutionStrategy
using POnG
include("helper_func.jl")

#%% Load traction test data
data = matread(joinpath(@__DIR__, "data_Alex/Ecoflex30_RubanMax.mat"))
lbd = vec(data["extension"] / 200 .+ 1)
# engineering stress
sig = vec(data["load"] / (40e-3 * 3e-3))
# Mooney-plot
MP = sig ./ (2 * (lbd .- lbd .^ (-2)))
st = 2000

#%% Load dynamic experiment data to determine β' (only using in-plane modes)
data = load(joinpath(@__DIR__, "parameters.jld2"))
@unpack ρs, VL, τ, n, μA, λA = data
datapar = matread(joinpath(@__DIR__, "data_Alex/Alex_EML_par_disprel.mat"))
dataperp = matread(joinpath(@__DIR__, "data_Alex/Alex_EML_perp_disprel.mat"))
# Velocities @f = 170 Hz
fip = 170
VSHpar, VSpar, VSHperp, VSperp = getSHSvelocities(datapar, dataperp, fip)

#%% Fits in the  small to moderate strain regime (hyperelastic model parameters + β')
λc = 1.5 # End of regime (determined quantitatively from error)
fittarget = vcat(VSHpar[λA.<λc], VSpar[λA.<λc], VSHperp[λA.<λc], VSperp[λA.<λc],) ./ sqrt(μA / ρs)
# Perform fits
MRfit = curve_fit(MR, 1 ./ lbd[lbd.<λc][st:end], MP[lbd.<λc][st:end], 1 ./ MP[lbd.<λc][st:end] .^ 2, [10e3, 10e3])
function fMR(λ, p)
    CMRip = C_MR.(λ, λ .^ (-0.41), ρs, VL, MRfit.param[1], MRfit.param[2], fip, τ, n, p[1])
    VSHpart, VSHperpt, VSpart, VSperpt = lw_mode_velocity(λ, ρs, CMRip)
    return vcat(real(VSHpart), real(VSpart), real(VSHperpt), real(VSperpt)) ./ sqrt(sum(MRfit.param) / ρs)
end
βMR = curve_fit(fMR, λA[λA.<λc], fittarget, 1 ./ fittarget .^ 2, [0.5], lower=[0.0])

GTfit = curve_fit(GT, 1 ./ (2 * lbd[lbd.<λc][st:end] .^ 2 .+ 1 ./ lbd[lbd.<λc][st:end]),
    MP[lbd.<λc][st:end], 1 ./ MP[lbd.<λc][st:end] .^ 2, [10e3, 1e3])
function fGT(λ, p)
    CGTip = C_GT.(λ, λ .^ (-0.41), ρs, VL, GTfit.param[1], GTfit.param[2], fip, τ, n, p[1])
    VSHpart, VSHperpt, VSpart, VSperpt = lw_mode_velocity(λ, ρs, CGTip)
    return vcat(real(VSHpart), real(VSpart), real(VSHperpt), real(VSperpt)) ./ sqrt(sum(GTfit.param) / ρs)
end
βGT = curve_fit(fGT, λA[λA.<λc], fittarget, 1 ./ fittarget .^ 2, [0.5], lower=[0.0])

Cfit = curve_fit(C, 1 ./ (lbd[lbd.<λc][st:end] .* sqrt.(2 * lbd[lbd.<λc][st:end] .+ 1 ./ lbd[lbd.<λc][st:end] .^ 2)), MP[lbd.<λc][st:end], 1 ./ MP[lbd.<λc][st:end] .^ 2, [10e3, 1e3])
function fC(λ, p)
    CCip = C_C.(λ, λ .^ (-0.41), ρs, VL, Cfit.param[1], Cfit.param[2], fip, τ, n, p[1])
    VSHpart, VSHperpt, VSpart, VSperpt = lw_mode_velocity(λ, ρs, CCip)
    return vcat(real(VSHpart), real(VSpart), real(VSHperpt), real(VSperpt)) ./ sqrt(sum(Cfit.param) / ρs)
end
βC = curve_fit(fC, λA[λA.<λc], fittarget, 1 ./ fittarget .^ 2, [0.5], lower=[0.0])

#%% Strain hardening regime
λc2 = 2.5 # qualitative input 
fittarget = vcat(VSHpar[λA.<λc2], VSpar[λA.<λc2], VSHperp[λA.<λc2], VSperp[λA.<λc2],) ./ sqrt(μA / ρs)
nrange = range(1.0, 3.0, step=0.01)
error = Vector{Float64}()
fitresults = Array{Float64,2}(undef, length(nrange), 3)
for i in eachindex(nrange)
    @. MRSHn(λ, p) = 1 / 2 * p[1] .+ 1 / λ * 1 / 2 * p[2] + p[3] * 3^(1 - nrange[i]) / 2 * (λ^2 + 2 / λ)^(nrange[i] - 1)
    fit = curve_fit(MRSHn, lbd[lbd.<λc2][st:end], MP[lbd.<λc2][st:end], 1 ./ MP[lbd.<λc2][st:end] .^ 2, [MRfit.param[1], MRfit.param[2], 1e3])
    append!(error, sum(fit.resid .^ 2))
    fitresults[i, :] = fit.param
end
MRSHfit = vcat(fitresults[argmin(error), :], nrange[argmin(error)])
function fMRSH(λ, p)
    CMRSHip = C_MRSH.(λ, λ .^ (-0.41), ρs, VL, MRSHfit[1], MRSHfit[2], MRSHfit[3], MRSHfit[4], fip, τ, n, p[1])
    VSHpart, VSHperpt, VSpart, VSperpt, _, _ = lw_mode_velocity(λ, ρs, CMRSHip)
    return vcat(real(VSHpart), real(VSpart), real(VSHperpt), real(VSperpt)) ./ sqrt(sum(MRSHfit[1:3]) / ρs)
end
βMRSH = curve_fit(fMRSH, λA[λA.<λc2], fittarget, 1 ./ fittarget .^ 2, [0.5], lower=[0.0])

error = Vector{Float64}()
fitresults = Array{Float64,2}(undef, length(nrange), 3)
for i in eachindex(nrange)
    @. GTSHn(λ, p) = p[1] / 2 + 1 / (2 * λ^2 + 1 / λ) * p[2] * 3 / 2 + p[3] * 3^(1 - nrange[i]) / 2 * (λ^2 + 2 / λ)^(nrange[i] - 1)
    fit = curve_fit(GTSHn, lbd[lbd.<λc2][st:end], MP[lbd.<λc2][st:end], 1 ./ MP[lbd.<λc2][st:end] .^ 2, [GTfit.param[1], GTfit.param[2], 1e3])
    append!(error, sum(fit.resid .^ 2))
    fitresults[i, :] = fit.param
end
GTSHfit = vcat(fitresults[argmin(error), :], nrange[argmin(error)])
function fGTSH(λ, p)
    CGTSHip = C_GTSH.(λ, λ .^ (-0.41), ρs, VL, GTSHfit[1], GTSHfit[2], GTSHfit[3], GTSHfit[4], fip, τ, n, p[1])
    VSHpart, VSHperpt, VSpart, VSperpt = lw_mode_velocity(λ, ρs, CGTSHip)
    return vcat(real(VSHpart), real(VSpart), real(VSHperpt), real(VSperpt)) ./ sqrt(sum(GTSHfit[1:3]) / ρs)
end
βGTSH = curve_fit(fGTSH, λA[λA.<λc2], fittarget, 1 ./ fittarget .^ 2, [0.5], lower=[0.0])

error = Vector{Float64}()
fitresults = Array{Float64,2}(undef, length(nrange), 3)
for i in eachindex(nrange)
    @. CSHn(λ, p) = p[1] / 2 + sqrt(3) / 2 * 1 / (λ * sqrt(2 * λ + 1 / λ^2)) * p[2] + p[3] * 3^(1 - nrange[i]) / 2 * (λ^2 + 2 / λ)^(nrange[i] - 1)
    fit = curve_fit(CSHn, lbd[lbd.<λc2][st:end], MP[lbd.<λc2][st:end], 1 ./ MP[lbd.<λc2][st:end] .^ 2, [Cfit.param[1], Cfit.param[2], 1e3])
    append!(error, sum(fit.resid .^ 2))
    fitresults[i, :] = fit.param
end
CSHfit = vcat(fitresults[argmin(error), :], nrange[argmin(error)])
function fCSH(λ, p)
    CCSHip = C_CSH.(λ, λ .^ (-0.41), ρs, VL, CSHfit[1], CSHfit[2], CSHfit[3], CSHfit[4], fip, τ, n, p[1])
    VSHpart, VSHperpt, VSpart, VSperpt = lw_mode_velocity(λ, ρs, CCSHip)
    return vcat(real(VSHpart), real(VSpart), real(VSHperpt), real(VSperpt)) ./ sqrt(sum(CSHfit[1:3]) / ρs)
end
βCSH = curve_fit(fCSH, λA[λA.<λc2], fittarget, 1 ./ fittarget .^ 2, [0.5], lower=[0.0])

# %% Limiting-chain regime
fittarget = vcat(VSHpar, VSpar, VSHperp, VSperp) ./ sqrt(μA / ρs)
fGMR(p) = sum((GMR(lbd[st:end], p) ./ MP[st:end] .- 1) .^ 2)
globalGMR = minimize(fGMR, [MRfit.param[1], MRfit.param[2], 10], 1e3, lower=[0.0, 0.0, 1.0], upper=[23e3, 23e3, 100])
function f2GMR(λ, p)
    CGMRip = C_GMR.(λ, λ .^ (-0.41), ρs, VL, population_mean(globalGMR)[1], population_mean(globalGMR)[2],
        population_mean(globalGMR)[3], fip, τ, n, p[1])
    VSHpart, VSHperpt, VSpart, VSperpt, _, _ = lw_mode_velocity(λ, ρs, CGMRip)
    return vcat(real(VSHpart), real(VSpart), real(VSHperpt), real(VSperpt)) ./ sqrt(sum(population_mean(globalGMR)[1:2]) / ρs)
end
βGMR = curve_fit(f2GMR, λA, fittarget, 1 ./ fittarget .^ 2, [0.5], lower=[0.0])

fGG(p) = sum((GG(lbd[st:end], p) ./ MP[st:end] .- 1) .^ 2)
globalGG = minimize(fGG, [GTfit.param[1], GTfit.param[2], 10], 1e3, lower=[0.0, 0.0, 1.0], upper=[23e3, 23e3, 100])
function f2GG(λ, p)
    CGGip = C_GG.(λ, λ .^ (-0.41), ρs, VL, population_mean(globalGG)[1], population_mean(globalGG)[2],
        population_mean(globalGG)[3], fip, τ, n, p[1])
    VSHpart, VSHperpt, VSpart, VSperpt, _, _ = lw_mode_velocity(λ, ρs, CGGip)
    return vcat(real(VSHpart), real(VSpart), real(VSHperpt), real(VSperpt)) ./ sqrt(sum(population_mean(globalGG)[1:2]) / ρs)
end
βGG = curve_fit(f2GG, λA, fittarget, 1 ./ fittarget .^ 2, [0.5], lower=[0.0])

fDCMR(p) = sum((DCMR(lbd[st:end], p) ./ MP[st:end] .- 1) .^ 2)
globalDCMR = minimize(fDCMR, [MRfit.param[1], MRfit.param[2], 10], 1e3, lower=[0.0, 0.0, 1.0], upper=[23e3, 23e3, 100])
function f2DCMR(λ, p)
    CDCMRip = C_DCMR2.(λ, λ .^ (-0.41), ρs, VL, population_mean(globalDCMR)[1], population_mean(globalDCMR)[2],
        population_mean(globalDCMR)[3], fip, τ, n, p[1])
    VSHpart, VSHperpt, VSpart, VSperpt, _, _ = lw_mode_velocity(λ, ρs, CDCMRip)
    return vcat(real(VSHpart), real(VSpart), real(VSHperpt), real(VSperpt)) ./ sqrt(sum(population_mean(globalDCMR)[1:2]) / ρs)
end
βDCMR = curve_fit(f2DCMR, λA, fittarget, 1 ./ fittarget .^ 2, [0.5], lower=[0.0])

#%%
fDCGT(p) = sum((DCGT(lbd[st:end], p) ./ MP[st:end] .- 1) .^ 2)
globalDCGT = minimize(fDCGT, [GTfit.param[1], GTfit.param[2], 10], 1e3, lower=[0.0, 0.0, 1.0], upper=[23e3, 23e3, 100])
function f2DCGT(λ, p)
    CDCGTip = C_DCGT.(λ, λ .^ (-0.41), ρs, VL, population_mean(globalDCGT)[1], population_mean(globalDCGT)[2],
        population_mean(globalDCGT)[3], fip, τ, n, p[1])
    VSHpart, VSHperpt, VSpart, VSperpt, _, _ = lw_mode_velocity(λ, ρs, CDCGTip)
    return vcat(real(VSHpart), real(VSpart), real(VSHperpt), real(VSperpt)) ./ sqrt(sum(population_mean(globalDCGT)[1:2]) / ρs)
end
βDCGT = curve_fit(f2DCGT, λA, fittarget, 1 ./ fittarget .^ 2, [0.5], lower=[0.0])

# %%
save(joinpath(@__DIR__, "fitparameters.jld2"),
    "MRfit", MRfit.param, "βMR", βMR.param,
    "GTfit", GTfit.param, "βGT", βGT.param,
    "Cfit", Cfit.param, "βC", βC.param,
    "MRSHfit", MRSHfit, "βMRSH", βMRSH.param,
    "GTSHfit", GTSHfit, "βGTSH", βGTSH.param,
    "CSHfit", CSHfit, "βCSH", βCSH.param,
    "GGfit", population_mean(globalGG), "βGG", βGG.param,
    "GMRfit", population_mean(globalGMR), "βGMR", βGMR.param,
    "DCMRfit", population_mean(globalDCMR), "βDCMR", βDCMR.param,
    "DCGTfit", population_mean(globalDCGT), "βDCGT", βDCGT.param)
