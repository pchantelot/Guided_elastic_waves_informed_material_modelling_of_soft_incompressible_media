using CairoMakie, ColorSchemes, LaTeXStrings, MathTeXEngine
using MAT, FileIO, UnPack, NaturalSort
include("preamble.jl")

# Load traction test data
data = matread(joinpath(@__DIR__, "data_Alex/Ecoflex30_RubanMax.mat"))
lbd = vec(data["extension"] / 200 .+ 1)
# engineering stress
sig = vec(data["load"] / (40e-3 * 3e-3))
# Cauchy stress
sigstar = sig .* lbd #NB: assumes uniaxial extension
# Define the Mooney space and start of fit
MP = sig ./ (2 * (lbd .- lbd .^ (-2)))
st = 2000

#%%
with_theme(My_theme) do
    fig = Figure(size=(2 * 402 * 2 / 3, 2 / (3 * 1.2) * 402), figure_padding=(2, 15, 2, 2))

    ax1 = Axis(fig[1, 1])
    Label(fig[1, 1, TopLeft()], "(a)", padding=(0, 30, -5, 0), fontsize=16)
    ax1.xlabel = L"λ"
    ax1.ylabel = L"$σ^e$ (kPa)"
    ax1.limits = (1, 3, 0, 70)
    scatter!(ax1, lbd[1:200:end], sig[1:200:end] / 1e3, color=:black, markersize=10)

    ax2 = Axis(fig[1, 2])
    Label(fig[1, 2, TopLeft()], "(b)", padding=(0, 30, -5, 0), fontsize=16)
    ax2.xlabel = L"1 / λ"
    ax2.ylabel = L"$\mathcal{M}$ (kPa)"
    ax2.limits = (0.3, 1, 10, 12)
    scatter!(ax2, 1 ./ lbd[st:200:end], MP[st:200:end] / 1e3, color=:black, markersize=10)

    display(fig)
    save(joinpath(@__DIR__, "figure2.pdf"), fig; pt_per_unit=1.0)
end
