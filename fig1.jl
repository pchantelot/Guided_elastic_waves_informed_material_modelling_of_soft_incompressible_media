using CairoMakie, ColorSchemes, LaTeXStrings, MathTeXEngine
using MAT, FileIO, UnPack, NaturalSort
include("preamble.jl")

# Load parameters
data = load(joinpath(@__DIR__, "parameters.jld2"))
@unpack ρs, VL, β, τ, n, μA, λs, λA, μP, λP = data
# get filenames for my data
A0files = filter(x -> occursin("def", x), readdir(joinpath(@__DIR__, "data_pierre")))
sort!(A0files, lt=natural)
# get Alex data
datapar = matread(joinpath(@__DIR__, "data_Alex/Alex_EML_par_disprel.mat"))
dataperp = matread(joinpath(@__DIR__, "data_Alex/Alex_EML_perp_disprel.mat"))

with_theme(My_theme) do

    fig = Figure(size=(2 * 402, 4 / (3 * 1.2) * 402), figure_padding=(2, 2, 2, 2))

    # Define the figure layout
    p1 = fig[1, 1] = GridLayout()
    ax1 = Axis(p1[1, 1])
    hidedecorations!(ax1)
    hidespines!(ax1)
    ax2 = Axis(p1[1, 2])
    ax3 = Axis(p1[1, 3])
    p2 = fig[2, 1] = GridLayout()
    ax4 = Axis(p2[1, 1])
    hidedecorations!(ax4)
    hidespines!(ax4)
    ax5 = Axis(p2[1, 2])
    ax6 = Axis(p2[1, 3])
    cb1 = fig[1, 2] = GridLayout()
    cb2 = fig[2, 2] = GridLayout()
    # Add labels
    Label(p1[1, 1, TopLeft()], "(a)", padding=(0, 0, -5, 0), fontsize=16)
    Label(p1[1, 2, TopLeft()], "(b)", padding=(0, 30, -5, 0), fontsize=16)
    Label(p1[1, 3, TopLeft()], "(c)", padding=(0, 30, -5, 0), fontsize=16)
    Label(p2[1, 1, TopLeft()], "(d)", padding=(0, 0, -5, 0), fontsize=16)
    Label(p2[1, 2, TopLeft()], "(e)", padding=(0, 30, -5, 0), fontsize=16)
    Label(p2[1, 3, TopLeft()], "(f)", padding=(0, 30, -5, 0), fontsize=16)

    ax2.xlabel = L"$k_\parallel$ (1/m)"
    ax2.ylabel = L"$f$ (Hz)"
    ax2.limits = (0, 300, 0, 300)
    @unpack kpar, fpar = datapar
    kSH = reverse!(kpar[1:2:end, :], dims=1)
    fSH = reverse!(fpar[1:2:end, :], dims=1)
    kS = reverse!(kpar[2:2:end, :], dims=1)
    fS = reverse!(fpar[2:2:end, :], dims=1)
    for i in axes(kSH, 1)
        scatter!(ax2, kSH[i, :], fSH[i, :], markersize=10,
            color=λA[i], colormap=ColorSchemes.Greens, colorrange=(1, 2.5))
        scatter!(ax2, kS[i, :], fS[i, :], markersize=10, marker=:utriangle,
            color=λA[i], colormap=ColorSchemes.PuRd, colorrange=(1, 2.5))
    end
    text!(ax2, 30, 225, text=L"S_0", align=(:center, :center),
        fontsize=20, color=ColorSchemes.PuRd[end])
    text!(ax2, 180, 135, text=L"SH_0", align=(:center, :center),
        fontsize=20, color=ColorSchemes.Greens[end])

    ax3.xlabel = L"$k_\perp$ (1/m)"
    ax3.ylabel = L"$f$ (Hz)"
    ax3.limits = (0, 300, 0, 300)
    @unpack kperp, fperp = dataperp
    kSH = reverse!(kperp[1:2:end, :], dims=1)
    fSH = reverse!(fperp[1:2:end, :], dims=1)
    kS = reverse!(kperp[2:2:end, :], dims=1)
    fS = reverse!(fperp[2:2:end, :], dims=1)
    for i in axes(kSH, 1)
        scatter!(ax3, kSH[i, :], fSH[i, :], markersize=10,
            color=λA[i], colormap=ColorSchemes.Greens, colorrange=(1, 2.5))
        scatter!(ax3, kS[i, :], fS[i, :], markersize=10, marker=:utriangle,
            color=λA[i], colormap=ColorSchemes.PuRd, colorrange=(1, 2.5))
    end
    text!(ax3, 60, 225, text=L"S_0", align=(:center, :center),
        fontsize=20, color=ColorSchemes.PuRd[end])
    text!(ax3, 200, 135, text=L"SH_0", align=(:center, :center),
        fontsize=20, color=ColorSchemes.Greens[end])
    Colorbar(cb1[1, 1], limits=(1, 2.5), colormap=ColorSchemes.Greens,
        ticklabelsvisible=false, ticksvisible=false)
    Colorbar(cb1[1, 2], limits=(1, 2.5), colormap=ColorSchemes.PuRd,
        label="λ", labelrotation=0, flipaxis=true, labelsize=20)
    colgap!(cb1, 1, 0)

    ax5.xlabel = L"$k_\parallel$ (1/m)"
    ax5.ylabel = L"$f$ (Hz)"
    ax5.limits = (0, 400, 0, 200)
    for i in axes(A0files, 1)
        data = load(joinpath(@__DIR__, "data_pierre", A0files[i]))
        @unpack kpar, fpar = data
        scatter!(ax5, kpar, fpar, markersize=10,
            color=λP[i], colormap=ColorSchemes.Blues, colorrange=(1, 2.5))
    end
    text!(ax5, 40, 150, text=L"A_0", align=(:center, :center),
        fontsize=20, color=ColorSchemes.Blues[end])

    ax6.xlabel = L"$k_\perp$ (1/m)"
    ax6.ylabel = L"$f$ (Hz)"
    ax6.limits = (0, 400, 0, 200)
    for i in axes(A0files, 1)
        data = load(joinpath(@__DIR__, "data_pierre", A0files[i]))
        @unpack kperp, fperp = data
        scatter!(ax6, kperp, fperp, markersize=10,
            color=λP[i], colormap=ColorSchemes.Blues, colorrange=(1, 2.5))
    end
    text!(ax6, 250, 150, text=L"A_0", align=(:center, :center),
        fontsize=20, color=ColorSchemes.Blues[end])
    Colorbar(cb2[1, 1], limits=(1, 2.5), colormap=ColorSchemes.Blues,
        label="λ", labelrotation=0, flipaxis=true, labelsize=20)

    # Alex experiment
    ax1.aspect = DataAspect()
    ax1.limits = (-7.5, 6, -2.8, 9.5)
    poly!(ax1, Rect(-3, 0, 6, 6), color=(:grey25, 0),
        strokecolor=:grey25, strokewidth=2, linestyle=(:dash, :dense))
    ext = 1.5
    poly!(ax1, Rect(-3 * 1 / sqrt(ext), 0, 6 * 1 / sqrt(ext), 6 * ext), color=:hotpink2)
    poly!(ax1, Rect(-3 * 1 / sqrt(ext), 6 * ext, 6 * 1 / sqrt(ext), 0.6), color=:grey25)
    poly!(ax1, Rect(-3 * 1 / sqrt(ext), 0, 6 * 1 / sqrt(ext), -0.6), color=:grey25)
    # DIC points
    scat = rand(180, 2)
    scatter!(ax1, -(6 * 1 / sqrt(ext) / 2 - 0.5) .+ scat[:, 1] .* (6 * 1 / sqrt(ext) - 1), 0.5 .+ scat[:, 2] .* (6 * ext - 1),
        color=:black, markersize=2)
    # measurement zone
    # poly!(ax1, Rect(-2, 6 * ext / 2 - 2, 4, 4), color = (:red, 0), 
    #     strokewidth = 2, strokecolor = :white, linestyle = (:dash, :dense))
    # text!(ax1, -2, 6 * ext / 2 + 2, text = "ROI", align = (:left,:bottom), color = :white)
    # distances
    arrows2d!(ax1, [0, 0], [-1.15, -1.15], [-3, 3], [0, 0], tipwidth=10, tiplength=12)
    text!(ax1, 0, -1.5, text=L"$L_0 = 60$ cm", align=(:center, :top))
    arrows2d!(ax1, [3.75, 3.75], 6 * ext / 2 * [1, 1], [0, 0], 6 * ext / 2 * [-1, 1], tipwidth=10, tiplength=12)
    text!(ax1, 4.75, 6 * ext / 2, text=L"$L = λ L_0$", align=(:center, :center), rotation=π / 2)
    # source
    poly!(ax1, Rect(-1, 6 * ext / 2 - 0.25, 2, 0.5), color=:grey25)
    arrows2d!(ax1, [0, 0], 6 * ext / 2 * [1, 1], 1.5 * cos(π / 4) * [1, -1], 1.5 * sin(π / 4) * [1, -1],
        tipwidth=10, tiplength=12, color=:white)
    arrows2d!(ax1, [-3.5], [6 * ext / 2], [2.5], [0], tipwidth=10, tiplength=12)
    text!(ax1, -3.7, 6 * ext / 2, text="shaker", align=(:right, :center))
    # axis 
    arrows2d!(ax1, [-7, -7], [6.5, 6.5], [2.5, 0], [0, 2.5], tipwidth=10, tiplength=12)
    text!(ax1, -4.75, 6.75, text=L"\mathbf{e}_3", align=(:left, :top))
    text!(ax1, -6.75, 8.5, text=L"\mathbf{e}_1", align=(:left, :bottom))

    # My experiment
    ax4.aspect = DataAspect()
    ax4.limits = (-7.5, 6, -2.8, 9.5)
    poly!(ax4, Rect(-3, 0, 6, 6), color=(:grey25, 0),
        strokecolor=:grey25, strokewidth=2, linestyle=(:dash, :dense))
    ext = 1.5
    poly!(ax4, Rect(-3 * 1 / sqrt(ext), 0, 6 * 1 / sqrt(ext), 6 * ext), color=:hotpink2)
    poly!(ax4, Rect(-3 * 1 / sqrt(ext), 6 * ext, 6 * 1 / sqrt(ext), 0.6), color=:grey25)
    poly!(ax4, Rect(-3 * 1 / sqrt(ext), 0, 6 * 1 / sqrt(ext), -0.6), color=:grey25)
    # measurement zone
    poly!(ax4, Rect(-2, 6 * ext / 2 - 2, 4, 4), color=(:red, 0),
        strokewidth=2, strokecolor=:white, linestyle=(:dash, :dense))
    text!(ax4, 2, 6 * ext / 2 + 2, text="ROI", color=:white, align=(:right, :bottom))
    # source
    poly!(ax4, Rect(-2.25, 6 * ext / 2 - 1, 0.5, 2), color=:grey25)
    arrows2d!(ax4, [-3.5], [6 * ext / 2], [1.5], [0], tipwidth=10, tiplength=12)
    text!(ax4, -3.7, 6 * ext / 2, text="shaker", align=(:right, :center))
    poly!(ax4, Circle(Point2f(-1, 6 * ext / 2 + 1.25), 0.5),
        strokecolor=:white, strokewidth=2, color=(:white, 0))
    poly!(ax4, Circle(Point2f(-1, 6 * ext / 2 + 1.25), 0.15), color=:white)
    # Laser line
    x = 0:0.01:4
    poly!(ax4, Point2f[(-2, 6 * ext / 2), (-2 + 4, 6 * ext / 2), (0.2, -1), (-0.2, -1)], color=(:lime, 0.2))
    poly!(ax4, Rect(-0.3, -2.5, 0.6, 1.5), color=:black)
    lines!(ax4, -2 .+ x, 6 * ext / 2 .+ 0.5 .* cos.(5 * (x .- x[1])) .* exp.(-0.5 .* (x .- x[1])), color=:lime)
    arrows2d!(ax4, [-3.5], [-1.75], [3], [0], tipwidth=10, tiplength=12)
    text!(ax4, -3.7, -1.75, text="laser\nsheet", align=(:right, :center), lineheight=0.75)

    rowgap!(fig.layout, 0)
    display(fig)
    save(joinpath(@__DIR__, "figure1.pdf"), fig; pt_per_unit=1.0)
end
