# %%
using ColorSchemes, MathTeXEngine, LaTeXStrings, GLMakie
using FileIO, UnPack, MAT, NaturalSort
# include(joinpath(homedir(), "Documents/GitHub/JuliaToolBox/general.jl"))
include("preamble.jl")
include("helper_func.jl")


# Visualize fit of uniaxial test
# Load traction test data
data = matread(joinpath(@__DIR__, "data_Alex/Ecoflex30_RubanMax.mat"))
lbd = vec(data["extension"] / 200 .+ 1)
# engineering stress
sig = vec(data["load"] / (40e-3 * 3e-3))
# Mooney-plot
MP = sig ./ (2 * (lbd .- lbd .^ (-2)))
st = 2000 # exclude early spurious points

# Visualize dynamic measurements
# Load parameters
data = load(joinpath(@__DIR__, "parameters.jld2"))
@unpack ρs, VL, β, τ, n, μA, λs, λA, μP, λP = data
# A0 data
A0files = filter(x -> occursin("def", x), readdir(joinpath(@__DIR__, "data_pierre/")))
sort!(A0files, lt=natural)
# SH0 and S0, Alex data
dataApar = matread(joinpath(@__DIR__, "data_Alex/Alex_EML_par_disprel.mat"))
dataAperp = matread(joinpath(@__DIR__, "data_Alex/Alex_EML_perp_disprel.mat"))

# Load functions to lift the velocity of each mode at the chosen frequency and to lift predictions
include("visualize_func.jl")

# predictions
λ = 1:0.05:3
λS = 1:0.1:2.5 # same λ range for all simulations
data = load(joinpath(@__DIR__, "fitparameters.jld2"))
@unpack MRfit, βMR, GTfit, βGT, Cfit, βC, MRSHfit, βMRSH, GTSHfit, βGTSH, CSHfit, βCSH,
GMRfit, βGMR, GGfit, βGG, DCMRfit, βDCMR, DCGTfit, βDCGT = data

# %%
with_theme(My_theme, palette=(color=ColorSchemes.tab10, marker=[:circle])) do

    fig = Figure(size=(2 * 1.2 * 360, 3 * 360))

    # setup axes
    ax0 = Axis(fig[0, 1])
    ax0.xlabel = L"λ"
    ax0.ylabel = L"$σ^e$ (kPa)"
    ax0.limits = (1, 3, 0, 70)
    ax00 = Axis(fig[0, 2])
    ax00.xlabel = L"1/λ"
    ax00.ylabel = L"$\mathcal{M}$ (kPa)"
    ax00.limits = (0.3, 1, 10, 12)
    ax1 = Axis(fig[1, 1])
    ax1.xlabel = L"λ"
    ax1.ylabel = L"V_∥/V_T"
    ax1.limits = (1, 2.5, 1, 4)
    ax2 = Axis(fig[1, 2])
    ax2.xlabel = L"λ"
    ax2.ylabel = L"V_⟂/V_T"
    ax2.limits = (1, 2.5, 1, 3)
    ax3 = Axis(fig[2, 1])
    ax3.xlabel = L"λ"
    ax3.ylabel = L"V_∥/V_T"
    ax3.limits = (1, 2.5, 0, 2.2)
    ax4 = Axis(fig[2, 2])
    ax4.xlabel = L"λ"
    ax4.ylabel = L"V_⟂/V_T"
    ax4.limits = (1, 2.5, 0, 1)

    Label(fig[0, 1:2, Top()], "Uniaxial extension test data", valign=:bottom,
        padding=(0, 0, 10, 0), fontsize=20)

    # setup sliders for frequency choice
    sg = SliderGrid(fig[3, 1:2], (label=L"f_{ip}", range=50:5:300, startvalue=170, format="{:.0d}Hz"),
        (label=L"f_{op}", range=10:5:200, startvalue=50, format="{:.0d}Hz"))
    sliderobservables = [s.value for s in sg.sliders]
    global fs = lift(sliderobservables...) do slvalues...
        [slvalues...]
    end
    # Setup corresponding labels
    Label(fig[1, 1:2, Top()], @lift("In-plane modes, f = $($(fs)[1]) Hz"), valign=:bottom,
        padding=(0, 0, 10, 0), fontsize=20)
    Label(fig[2, 1:2, Top()], @lift("Out-of-plane mode, f = $($(fs)[2]) Hz"), valign=:bottom,
        padding=(0, 0, 10, 0), fontsize=20)

    # setup some toggles for hyperelastic model choice
    tgrid = GridLayout(fig[:, 3])
    toggle_colors = resample(ColorSchemes.tab10, 10, (alpha) -> 0.75)
    toggleMR = Toggle(tgrid[1, 1], active=true, framecolor_active=toggle_colors[1], buttoncolor=:black)
    Label(tgrid[1, 2], "Mooney-Rivlin", halign=:right)
    toggleGT = Toggle(tgrid[2, 1], active=false, framecolor_active=toggle_colors[2], buttoncolor=:black)
    Label(tgrid[2, 2], "Gent-Thomas", halign=:right)
    toggleC = Toggle(tgrid[3, 1], active=false, framecolor_active=toggle_colors[3], buttoncolor=:black)
    Label(tgrid[3, 2], "Carroll", halign=:right)
    toggleMRSH = Toggle(tgrid[4, 1], active=false, framecolor_active=toggle_colors[4], buttoncolor=:black)
    Label(tgrid[4, 2], L"MR + $I_1^N$", halign=:right)
    toggleGTSH = Toggle(tgrid[5, 1], active=false, framecolor_active=toggle_colors[5], buttoncolor=:black)
    Label(tgrid[5, 2], L"GT + $I_1^N$", halign=:right)
    toggleCSH = Toggle(tgrid[6, 1], active=false, framecolor_active=toggle_colors[6], buttoncolor=:black)
    Label(tgrid[6, 2], L"C + $I_1^N$", halign=:right)
    toggleGMR = Toggle(tgrid[7, 1], active=false, framecolor_active=toggle_colors[7], buttoncolor=:black)
    Label(tgrid[7, 2], "Gent-MR", halign=:right)
    toggleGG = Toggle(tgrid[8, 1], active=false, framecolor_active=toggle_colors[8], buttoncolor=:black)
    Label(tgrid[8, 2], "Gent-Gent", halign=:right)
    toggleDCMR = Toggle(tgrid[9, 1], active=false, framecolor_active=toggle_colors[9], buttoncolor=:black)
    Label(tgrid[9, 2], "DC-MR", halign=:right)
    toggleDCGT = Toggle(tgrid[10, 1], active=false, framecolor_active=toggle_colors[10], buttoncolor=:black)
    Label(tgrid[10, 2], "DC-GT", halign=:right)


    # Static plots
    scatter!(ax00, 1 ./ lbd[st:200:end], MP[st:200:end] / 1e3, color=:black, markersize=10)
    scatter!(ax0, lbd[st:200:end], sig[st:200:end] / 1e3, color=:black, markersize=10)

    lines!(ax0, λ, 2 * (λ .- λ .^ (-2)) .* MR(1 ./ λ, MRfit) / 1e3, alpha=@lift($(toggleMR.active) ? 1.0 : 0.0))
    lines!(ax00, 1 ./ λ, MR(1 ./ λ, MRfit) / 1e3, alpha=@lift($(toggleMR.active) ? 1.0 : 0.0))
    lines!(ax0, λ, 2 * (λ .- λ .^ (-2)) .* GT(1 ./ (2 * λ .^ 2 .+ 1 ./ λ), GTfit) / 1e3, alpha=@lift($(toggleGT.active) ? 1.0 : 0.0))
    lines!(ax00, 1 ./ λ, GT(1 ./ (2 * λ .^ 2 .+ 1 ./ λ), GTfit) / 1e3, alpha=@lift($(toggleGT.active) ? 1.0 : 0.0))
    lines!(ax0, λ, 2 * (λ .- λ .^ (-2)) .* C(1 ./ (λ .* sqrt.(2 * λ .+ 1 ./ λ .^ 2)), Cfit) / 1e3, alpha=@lift($(toggleC.active) ? 1.0 : 0.0))
    lines!(ax00, 1 ./ λ, C(1 ./ (λ .* sqrt.(2 * λ .+ 1 ./ λ .^ 2)), Cfit) / 1e3, alpha=@lift($(toggleC.active) ? 1.0 : 0.0))
    lines!(ax0, λ, 2 * (λ .- λ .^ (-2)) .* MRSH(λ, MRSHfit) / 1e3, alpha=@lift($(toggleMRSH.active) ? 1.0 : 0.0))
    lines!(ax00, 1 ./ λ, MRSH(λ, MRSHfit) / 1e3, alpha=@lift($(toggleMRSH.active) ? 1.0 : 0.0))
    lines!(ax0, λ, 2 * (λ .- λ .^ (-2)) .* GTSH(λ, GTSHfit) / 1e3, alpha=@lift($(toggleGTSH.active) ? 1.0 : 0.0))
    lines!(ax00, 1 ./ λ, GTSH(λ, GTSHfit) / 1e3, alpha=@lift($(toggleGTSH.active) ? 1.0 : 0.0))
    lines!(ax0, λ, 2 * (λ .- λ .^ (-2)) .* CSH(λ, CSHfit) / 1e3, alpha=@lift($(toggleCSH.active) ? 1.0 : 0.0))
    lines!(ax00, 1 ./ λ, CSH(λ, CSHfit) / 1e3, alpha=@lift($(toggleCSH.active) ? 1.0 : 0.0))
    lines!(ax0, λ, 2 * (λ .- λ .^ (-2)) .* GMR(λ, GMRfit) / 1e3, alpha=@lift($(toggleGMR.active) ? 1.0 : 0.0))
    lines!(ax00, 1 ./ λ, GMR(λ, GMRfit) / 1e3, alpha=@lift($(toggleGMR.active) ? 1.0 : 0.0))
    lines!(ax0, λ, 2 * (λ .- λ .^ (-2)) .* GG(λ, GGfit) / 1e3, alpha=@lift($(toggleGG.active) ? 1.0 : 0.0))
    lines!(ax00, 1 ./ λ, GG(λ, GGfit) / 1e3, alpha=@lift($(toggleGG.active) ? 1.0 : 0.0))
    lines!(ax0, λ, 2 * (λ .- λ .^ (-2)) .* DCMR(λ, DCMRfit) / 1e3, alpha=@lift($(toggleDCMR.active) ? 1.0 : 0.0))
    lines!(ax00, 1 ./ λ, DCMR(λ, DCMRfit) / 1e3, alpha=@lift($(toggleDCMR.active) ? 1.0 : 0.0))
    lines!(ax0, λ, 2 * (λ .- λ .^ (-2)) .* DCGT(λ, DCGTfit) / 1e3, alpha=@lift($(toggleDCGT.active) ? 1.0 : 0.0))
    lines!(ax00, 1 ./ λ, DCGT(λ, DCGTfit) / 1e3, alpha=@lift($(toggleDCGT.active) ? 1.0 : 0.0))

    # Dynamic plots
    # Plot the current phase velocities
    scatter!(ax1, λA, lift(getVSHpar, fs), color=:black)
    scatter!(ax1, λA, lift(getVSpar, fs), color=:black, marker=:utriangle)
    scatter!(ax2, λA, lift(getVSHperp, fs), color=:black)
    scatter!(ax2, λA, lift(getVSperp, fs), color=:black, marker=:utriangle)
    scatter!(ax3, λP, lift(getVApar, fs), color=:black)
    scatter!(ax4, λP, lift(getVAperp, fs), color=:black)

    # prediction
    # Mooney-Rivlin model
    lines!(ax1, λ, lift(predicVSHparMR, fs), alpha=@lift($(toggleMR.active) ? 1.0 : 0.0), color=toggle_colors[1])
    lines!(ax1, λ, lift(predicVSparMR, fs), alpha=@lift($(toggleMR.active) ? 1.0 : 0.0), color=toggle_colors[1])
    lines!(ax2, λ, lift(predicVSHperpMR, fs), alpha=@lift($(toggleMR.active) ? 1.0 : 0.0), color=toggle_colors[1])
    lines!(ax2, λ, lift(predicVSperpMR, fs), alpha=@lift($(toggleMR.active) ? 1.0 : 0.0), color=toggle_colors[1])
    lines!(ax3, λ, lift(predicVAparMR, fs), alpha=@lift($(toggleMR.active) ? 1.0 : 0.0), color=toggle_colors[1])
    lines!(ax4, λS, lift(predicVAperpMR, fs), alpha=@lift($(toggleMR.active) ? 1.0 : 0.0), color=toggle_colors[1])
    # Gent-Thomas model
    lines!(ax1, λ, lift(predicVSHparGT, fs), alpha=@lift($(toggleGT.active) ? 1.0 : 0.0), color=toggle_colors[2])
    lines!(ax1, λ, lift(predicVSparGT, fs), alpha=@lift($(toggleGT.active) ? 1.0 : 0.0), color=toggle_colors[2])
    lines!(ax2, λ, lift(predicVSHperpGT, fs), alpha=@lift($(toggleGT.active) ? 1.0 : 0.0), color=toggle_colors[2])
    lines!(ax2, λ, lift(predicVSperpGT, fs), alpha=@lift($(toggleGT.active) ? 1.0 : 0.0), color=toggle_colors[2])
    lines!(ax3, λ, lift(predicVAparGT, fs), alpha=@lift($(toggleGT.active) ? 1.0 : 0.0), color=toggle_colors[2])
    lines!(ax4, λS, lift(predicVAperpGT, fs), alpha=@lift($(toggleGT.active) ? 1.0 : 0.0), color=toggle_colors[2])
    # Carroll model
    lines!(ax1, λ, lift(predicVSHparC, fs), alpha=@lift($(toggleC.active) ? 1.0 : 0.0), color=toggle_colors[3])
    lines!(ax1, λ, lift(predicVSparC, fs), alpha=@lift($(toggleC.active) ? 1.0 : 0.0), color=toggle_colors[3])
    lines!(ax2, λ, lift(predicVSHperpC, fs), alpha=@lift($(toggleC.active) ? 1.0 : 0.0), color=toggle_colors[3])
    lines!(ax2, λ, lift(predicVSperpC, fs), alpha=@lift($(toggleC.active) ? 1.0 : 0.0), color=toggle_colors[3])
    lines!(ax3, λ, lift(predicVAparC, fs), alpha=@lift($(toggleC.active) ? 1.0 : 0.0), color=toggle_colors[3])
    lines!(ax4, λS, lift(predicVAperpC, fs), alpha=@lift($(toggleC.active) ? 1.0 : 0.0), color=toggle_colors[3])
    #  MRSH model
    lines!(ax1, λ, lift(predicVSHparMRSH, fs), alpha=@lift($(toggleMRSH.active) ? 1.0 : 0.0), color=toggle_colors[4])
    lines!(ax1, λ, lift(predicVSparMRSH, fs), alpha=@lift($(toggleMRSH.active) ? 1.0 : 0.0), color=toggle_colors[4])
    lines!(ax2, λ, lift(predicVSHperpMRSH, fs), alpha=@lift($(toggleMRSH.active) ? 1.0 : 0.0), color=toggle_colors[4])
    lines!(ax2, λ, lift(predicVSperpMRSH, fs), alpha=@lift($(toggleMRSH.active) ? 1.0 : 0.0), color=toggle_colors[4])
    lines!(ax3, λ, lift(predicVAparMRSH, fs), alpha=@lift($(toggleMRSH.active) ? 1.0 : 0.0), color=toggle_colors[4])
    lines!(ax4, λS, lift(predicVAperpMRSH, fs), alpha=@lift($(toggleMRSH.active) ? 1.0 : 0.0), color=toggle_colors[4])
    #  GTSH model
    lines!(ax1, λ, lift(predicVSHparGTSH, fs), alpha=@lift($(toggleGTSH.active) ? 1.0 : 0.0), color=toggle_colors[5])
    lines!(ax1, λ, lift(predicVSparGTSH, fs), alpha=@lift($(toggleGTSH.active) ? 1.0 : 0.0), color=toggle_colors[5])
    lines!(ax2, λ, lift(predicVSHperpGTSH, fs), alpha=@lift($(toggleGTSH.active) ? 1.0 : 0.0), color=toggle_colors[5])
    lines!(ax2, λ, lift(predicVSperpGTSH, fs), alpha=@lift($(toggleGTSH.active) ? 1.0 : 0.0), color=toggle_colors[5])
    lines!(ax3, λ, lift(predicVAparGTSH, fs), alpha=@lift($(toggleGTSH.active) ? 1.0 : 0.0), color=toggle_colors[5])
    lines!(ax4, λS, lift(predicVAperpGTSH, fs), alpha=@lift($(toggleGTSH.active) ? 1.0 : 0.0), color=toggle_colors[5])
    #  CSH model
    lines!(ax1, λ, lift(predicVSHparCSH, fs), alpha=@lift($(toggleCSH.active) ? 1.0 : 0.0), color=toggle_colors[6])
    lines!(ax1, λ, lift(predicVSparCSH, fs), alpha=@lift($(toggleCSH.active) ? 1.0 : 0.0), color=toggle_colors[6])
    lines!(ax2, λ, lift(predicVSHperpCSH, fs), alpha=@lift($(toggleCSH.active) ? 1.0 : 0.0), color=toggle_colors[6])
    lines!(ax2, λ, lift(predicVSperpCSH, fs), alpha=@lift($(toggleCSH.active) ? 1.0 : 0.0), color=toggle_colors[6])
    lines!(ax3, λ, lift(predicVAparCSH, fs), alpha=@lift($(toggleCSH.active) ? 1.0 : 0.0), color=toggle_colors[6])
    lines!(ax4, λS, lift(predicVAperpCSH, fs), alpha=@lift($(toggleCSH.active) ? 1.0 : 0.0), color=toggle_colors[6])
    #  GMR model
    lines!(ax1, λ, lift(predicVSHparGMR, fs), alpha=@lift($(toggleGMR.active) ? 1.0 : 0.0), color=toggle_colors[7])
    lines!(ax1, λ, lift(predicVSparGMR, fs), alpha=@lift($(toggleGMR.active) ? 1.0 : 0.0), color=toggle_colors[7])
    lines!(ax2, λ, lift(predicVSHperpGMR, fs), alpha=@lift($(toggleGMR.active) ? 1.0 : 0.0), color=toggle_colors[7])
    lines!(ax2, λ, lift(predicVSperpGMR, fs), alpha=@lift($(toggleGMR.active) ? 1.0 : 0.0), color=toggle_colors[7])
    lines!(ax3, λ, lift(predicVAparGMR, fs), alpha=@lift($(toggleGMR.active) ? 1.0 : 0.0), color=toggle_colors[7])
    lines!(ax4, λS, lift(predicVAperpGMR, fs), alpha=@lift($(toggleGMR.active) ? 1.0 : 0.0), color=toggle_colors[7])
    #  GG model
    lines!(ax1, λ, lift(predicVSHparGG, fs), alpha=@lift($(toggleGG.active) ? 1.0 : 0.0), color=toggle_colors[8])
    lines!(ax1, λ, lift(predicVSparGG, fs), alpha=@lift($(toggleGG.active) ? 1.0 : 0.0), color=toggle_colors[8])
    lines!(ax2, λ, lift(predicVSHperpGG, fs), alpha=@lift($(toggleGG.active) ? 1.0 : 0.0), color=toggle_colors[8])
    lines!(ax2, λ, lift(predicVSperpGG, fs), alpha=@lift($(toggleGG.active) ? 1.0 : 0.0), color=toggle_colors[8])
    lines!(ax3, λ, lift(predicVAparGG, fs), alpha=@lift($(toggleGG.active) ? 1.0 : 0.0), color=toggle_colors[8])
    lines!(ax4, λS, lift(predicVAperpGG, fs), alpha=@lift($(toggleGG.active) ? 1.0 : 0.0), color=toggle_colors[8])
    #  DCMR model
    lines!(ax1, λ, lift(predicVSHparDCMR, fs), alpha=@lift($(toggleDCMR.active) ? 1.0 : 0.0), color=toggle_colors[9])
    lines!(ax1, λ, lift(predicVSparDCMR, fs), alpha=@lift($(toggleDCMR.active) ? 1.0 : 0.0), color=toggle_colors[9])
    lines!(ax2, λ, lift(predicVSHperpDCMR, fs), alpha=@lift($(toggleDCMR.active) ? 1.0 : 0.0), color=toggle_colors[9])
    lines!(ax2, λ, lift(predicVSperpDCMR, fs), alpha=@lift($(toggleDCMR.active) ? 1.0 : 0.0), color=toggle_colors[9])
    lines!(ax3, λ, lift(predicVAparDCMR, fs), alpha=@lift($(toggleDCMR.active) ? 1.0 : 0.0), color=toggle_colors[9])
    lines!(ax4, λS, lift(predicVAperpDCMR, fs), alpha=@lift($(toggleDCMR.active) ? 1.0 : 0.0), color=toggle_colors[9])
    #  DCGT model
    lines!(ax1, λ, lift(predicVSHparDCGT, fs), alpha=@lift($(toggleDCGT.active) ? 1.0 : 0.0), color=toggle_colors[10])
    lines!(ax1, λ, lift(predicVSparDCGT, fs), alpha=@lift($(toggleDCGT.active) ? 1.0 : 0.0), color=toggle_colors[10])
    lines!(ax2, λ, lift(predicVSHperpDCGT, fs), alpha=@lift($(toggleDCGT.active) ? 1.0 : 0.0), color=toggle_colors[10])
    lines!(ax2, λ, lift(predicVSperpDCGT, fs), alpha=@lift($(toggleDCGT.active) ? 1.0 : 0.0), color=toggle_colors[10])
    lines!(ax3, λ, lift(predicVAparDCGT, fs), alpha=@lift($(toggleDCGT.active) ? 1.0 : 0.0), color=toggle_colors[10])
    lines!(ax4, λS, lift(predicVAperpDCGT, fs), alpha=@lift($(toggleDCGT.active) ? 1.0 : 0.0), color=toggle_colors[10])

    display(fig)
end
