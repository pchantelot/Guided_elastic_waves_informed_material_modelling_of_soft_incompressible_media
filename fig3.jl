using CairoMakie, ColorSchemes, LaTeXStrings, MathTeXEngine
using MAT, FileIO, UnPack, NaturalSort
using LsqFit
include("preamble.jl")
include("helper_func.jl")

# Load parameters
data = load(joinpath(@__DIR__, "parameters.jld2"))
@unpack ρs, VL, β, τ, n, μA, λs, λA, μP, λP = data
# get filenames for my data
A0files = filter(x -> occursin("def", x), readdir(joinpath(@__DIR__, "data_pierre/")))
sort!(A0files, lt=natural)

with_theme(My_theme) do

    fig = Figure(size=(2 * 402, 4 / (3 * 1.2) * 402), figure_padding=(2, 15, 2, 2))

    # Define figure layout
    ax1 = Axis(fig[1, 1])
    ax2 = Axis(fig[1, 2])
    ax3 = Axis(fig[1, 3])
    #ax4 = Axis(fig[2,1])
    ax5 = Axis(fig[2, 2])
    ax6 = Axis(fig[2, 3])
    # add labels
    Label(fig[1, 1, TopLeft()], "(a)", padding=(0, 30, -5, 0), fontsize=16)
    Label(fig[1, 2, TopLeft()], "(b)", padding=(0, 30, -5, 0), fontsize=16)
    Label(fig[1, 3, TopLeft()], "(c)", padding=(0, 30, -5, 0), fontsize=16)
    Label(fig[2, 2, TopLeft()], "(d)", padding=(0, 30, -5, 0), fontsize=16)
    Label(fig[2, 3, TopLeft()], "(e)", padding=(0, 30, -5, 0), fontsize=16)

    ax1.xlabel = L"$f$ (Hz)"
    ax1.ylabel = L"$\mu', \; \mu''$ (kPa)"
    ax1.limits = (0.05, 150, 0.5, 100)
    ax1.yticks = [1, 5, 10, 50]
    ax1.xticks = [0.1, 1, 10, 100]
    ax1.xscale = log10
    ax1.yscale = log10
    # rheology fit from exp with Samuel
    data = load(joinpath(@__DIR__, "data_Alex/rheoOO30.jld2"))
    @unpack freq, Gp, Gpp = data
    f(x, p) = vcat(real(p[1] * (1 .+ (1im * x * p[2]) .^ (p[3]))), imag(p[1] * (1 .+ (1im * x * p[2]) .^ (p[3]))))
    rheofit = curve_fit(f, 2 * pi * freq, vcat(Gp, Gpp), 1 ./ vcat(Gp, Gpp) .^ 2, [20e3, 1e-4, 0.3])
    scatter!(ax1, vec(freq), vec(Gp) / 1e3, label=L"μ'")
    scatter!(ax1, vec(freq), vec(Gpp) / 1e3, label=L"μ''")
    x = 0.01:0.1:200
    # Plot μ
    @. μ(x, p) = p[1] * (1 + (1im * x * p[2])^(p[3]))
    lines!(ax1, x, real(μ(2 * pi * x, rheofit.param)) / 1e3, color=:black, linestyle=:dash)
    lines!(ax1, x, imag(μ(2 * pi * x, rheofit.param)) / 1e3, color=:black, linestyle=:dash)
    axislegend(ax1, position=:rb, labelsize=16)

    ax2.xlabel = L"$k_\parallel$ (1/m)"
    ax2.ylabel = L"$f$ (Hz)"
    ax2.limits = (0, 300, 0, 300)
    data = matread(joinpath(@__DIR__, "data_Alex/Alex_EML_par_disprel.mat"))
    @unpack kpar, fpar = data
    kSH = reverse!(kpar[1:2:end, :], dims=1)
    fSH = reverse!(fpar[1:2:end, :], dims=1)
    kS = reverse!(kpar[2:2:end, :], dims=1)
    fS = reverse!(fpar[2:2:end, :], dims=1)
    for i in [1]
        scatter!(ax2, kSH[i, :], fSH[i, :],
            color=λA[end-1], colormap=ColorSchemes.Greens, colorrange=(1, 2.5),
            label=L"λ = %$(round(λA[i], sigdigits = 3))")
        scatter!(ax2, kS[i, :], fS[i, :], marker=:utriangle,
            color=λA[end-1], colormap=ColorSchemes.PuRd, colorrange=(1, 2.5),
            label=L"λ = %$(round(λA[i], sigdigits = 3))")
    end

    ax3.xlabel = L"$k_\perp$ (1/m)"
    ax3.ylabel = L"$f$ (Hz)"
    ax3.limits = (0, 300, 0, 300)
    data = matread(joinpath(@__DIR__, "data_Alex/Alex_EML_perp_disprel.mat"))
    @unpack kperp, fperp = data
    kSH = reverse!(kperp[1:2:end, :], dims=1)
    fSH = reverse!(fperp[1:2:end, :], dims=1)
    kS = reverse!(kperp[2:2:end, :], dims=1)
    fS = reverse!(fperp[2:2:end, :], dims=1)
    for i in [1]
        scatter!(ax3, kSH[i, :], fSH[i, :],
            color=λA[end-1], colormap=ColorSchemes.Greens, colorrange=(1, 2.5))
        scatter!(ax3, kS[i, :], fS[i, :], marker=:utriangle,
            color=λA[end-1], colormap=ColorSchemes.PuRd, colorrange=(1, 2.5))
    end

    ax5.xlabel = L"$k_\parallel$ (1/m)"
    ax5.ylabel = L"$f$ (Hz)"
    ax5.limits = (0, 400, 0, 200)
    for i in [1]
        data = load(joinpath(@__DIR__, "data_pierre/", A0files[i]))
        @unpack kpar, fpar = data
        scatter!(ax5, kpar, fpar,
            color=λP[end-1], colormap=ColorSchemes.Blues, colorrange=(1, 2.5),
            label=L"λ = %$(round(λP[i], sigdigits = 3))")
    end

    ax6.xlabel = L"$k_\perp$ (1/m)"
    ax6.ylabel = L"$f$ (Hz)"
    ax6.limits = (0, 400, 0, 200)
    for i in [1]
        data = load(joinpath(@__DIR__, "data_pierre/", A0files[i]))
        @unpack kperp, fperp = data
        scatter!(ax6, kperp, fperp,
            color=λP[end-1], colormap=ColorSchemes.Blues, colorrange=(1, 2.5))
    end

    # Long-wavelength predictions
    f = 0:300
    α = 0.28 # NB: from Alex
    CMRip = C_MR.(λA[1], λA[1] .^ (-0.41), ρs, VL, μA * (1 - α), μA * α, f, τ, n, β)
    VSHpar, VSHperp, VSpar, VSperp, VApar, VAperp = lw_mode_velocity(f, ρs, CMRip)
    lines!(ax2, 2 * pi * f ./ real(VSHpar), f, linestyle=:solid, color=λA[end-1],
        colormap=ColorSchemes.Greens, colorrange=(1, 2.5))
    lines!(ax2, 2 * pi * f ./ real(VSpar), f, linestyle=:solid, color=λA[end-1],
        colormap=ColorSchemes.PuRd, colorrange=(1, 2.5))
    lines!(ax3, 2 * pi * f ./ real(VSHperp), f, linestyle=:solid, color=λA[end-1],
        colormap=ColorSchemes.Greens, colorrange=(1, 2.5))
    lines!(ax3, 2 * pi * f ./ real(VSperp), f, linestyle=:solid, color=λA[end-1],
        colormap=ColorSchemes.PuRd, colorrange=(1, 2.5))

    CMRop = C_MR.(λP[1], λP[1] .^ (-0.39), ρs, VL, μP * (1 - α), μP * α, f, τ, n, β)
    VSHpar, VSHperp, VSpar, VSperp, VApar, VAperp = lw_mode_velocity(f, ρs, CMRop)
    lines!(ax5, 2 * pi * f ./ real(VApar), f, linestyle=:solid, color=λP[end-1],
        colormap=ColorSchemes.Blues, colorrange=(1, 2.5))
    lines!(ax6, 2 * pi * f ./ real(VAperp), f, linestyle=:solid, color=λP[end-1],
        colormap=ColorSchemes.Blues, colorrange=(1, 2.5))

    simu = load(joinpath(@__DIR__, "fig3/par_1.11.jld2"))
    @unpack k, f = simu
    lines!(ax5, real(k), f, color=:black, linestyle=:solid)

    simu = load(joinpath(@__DIR__, "fig3/perp_1.11.jld2"))
    @unpack k, f = simu
    lines!(ax6, real(k), f, color=:black, linestyle=:solid)

    # Making a Legend by hand
    marker_SH0 = MarkerElement(color=λA[end-1], colormap=ColorSchemes.Greens, colorrange=(1, 2.5),
        marker=:circle, markersize=16, strokecolor=:black, strokewidth=1)
    marker_S0 = MarkerElement(color=λA[end-1], colormap=ColorSchemes.PuRd, colorrange=(1, 2.5),
        marker=:circle, markersize=16, strokecolor=:black, strokewidth=1)
    marker_A0 = MarkerElement(color=λP[end-1], colormap=ColorSchemes.Blues, colorrange=(1, 2.5),
        marker=:circle, markersize=16, strokecolor=:black, strokewidth=1)
    lines_lw1 = LineElement(color=λA[end-1], colormap=ColorSchemes.Greens, colorrange=(1, 2.5),
        linewidth=3, points=Point2f[(0, 1), (1, 1)])
    lines_lw2 = LineElement(color=λA[end-1], colormap=ColorSchemes.PuRd, colorrange=(1, 2.5),
        linewidth=3, points=Point2f[(0, 0.5), (1, 0.5)])
    lines_lw3 = LineElement(color=λP[end-1], colormap=ColorSchemes.Blues, colorrange=(1, 2.5),
        linewidth=3, points=Point2f[(0, 0), (1, 0)])
    lines_simu = LineElement(color=:black, linewidth=3, linestyle=:solid)

    Legend(fig[2, 1], [marker_SH0, marker_S0, marker_A0, [lines_lw1, lines_lw2, lines_lw3], lines_simu],
        [L"$SH_0$, $λ = %$(round(λA[1], sigdigits = 3))$", L"$S_0$, $λ = %$(round(λA[1], sigdigits = 3))$",
            L"$A_0$, $λ = %$(round(λP[1], sigdigits = 3))$", "Long-wavelength \napproximations",
            "Theoretical predictions"],
        labelsize=16)

    colsize!(fig.layout, 1, Relative(1 / 3))
    rowgap!(fig.layout, 0)
    display(fig)
    save(joinpath(@__DIR__, "figure3.pdf"), fig; pt_per_unit=1.0)

end
