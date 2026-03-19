#%% Save all parameters useful to plot figures in one place 
using MAT, FileIO, UnPack, NaturalSort
using LsqFit

#%% 0. For all data
# Fixed from previous studies
ρs = 1070
VL = 1000
# β value from Alex fit (λ < 1.8). We will obtain our own value.
β = 0.29
# rheology fit from exp with Samuel
data = load(joinpath(@__DIR__, "data_Alex/rheoOO30.jld2"))
@unpack freq, Gp, Gpp = data
# Not fitting the norm of G to have significant weight on imaginary part
f(x, p) = vcat(real(p[1] * (1 .+ (1im * x * p[2]) .^ (p[3]))), imag(p[1] * (1 .+ (1im * x * p[2]) .^ (p[3]))))
rheofit = curve_fit(f, 2 * pi * freq, vcat(Gp, Gpp), 1 ./ vcat(Gp, Gpp) .^ 2, [20e3, 1e-4, 0.3])
τ = rheofit.param[2]
n = rheofit.param[3]

#%% 1. For Alex data
# get shear modulus from Alex plate extension test assuming it is uniaxial
ls = [1.218, 1.357, 1.477, 1.607, 1.751]
sigs = [19.89, 30.08, 41.36, 54.38, 69.45]
ls3 = ls .^ (-0.5) # uniaxial deformation
MP = sigs ./ (2 * ls .* (ls .- 1 ./ (ls .^ 3 .* ls3 .^ 2)))
@. MR(x, p) = 1 / 2 * p[1] .+ 1 / 2 * x * p[2]
MRfit = curve_fit(MR, ls3[ls.<1.61] .^ 2, MP[ls.<1.61], 1 ./ MP[ls.<1.61] .^ 2, [10e3, 1e3])
μA = sum(MRfit.param) * 1e3
# Compute static deformation in the middle of the plate
λs = 1 + ρs * 9.81 * 60e-2 / (6 * μA)
# get Alex lambdas
data = matread(joinpath(@__DIR__, "data_Alex/Alex_EML.mat"))
@unpack lbd = data
λA = vec(lbd) * λs

#%% 2. For Pierre data
μP = 23.4e3 # NB: usual shear modulus obtained from Mooney-plots
# Deformation for my data
sc = 40e-2 / 680
d0th = 3e-2 / sc
dpar = [56.45, 57.37, 63.84, 71.29, 76.87, 83.27, 89.03, 94.78, 100.51, 106.75, 113.26]
dperp = [46.18, 45.85, 44.19, 42.77, 41.35, 40.5, 39.25, 38, 37.27, 36.68, 36.7]
λP = dpar / d0th;

#%%
save(joinpath(@__DIR__, "parameters.jld2"), "ρs", ρs, "VL", VL, "β", β, "τ", τ, "n", n,
    "μA", μA, "λs", λs, "λA", λA,
    "μP", μP, "λP", λP)
