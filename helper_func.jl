using POnG
using UnPack, FileIO

@. MR(x, p) = 1 / 2 * p[1] .+ 1 / 2 * x * p[2]
@. GT(x, p) = p[1] / 2 + x * p[2] * 3 / 2
@. C(x, p) = p[1] / 2 + sqrt(3) / 2 * x * p[2]
@. MRSH(λ, p) = 1 / 2 * p[1] .+ 1 / λ * 1 / 2 * p[2] + p[3] * 3^(1 - p[4]) / 2 * (λ^2 + 2 / λ)^(p[4] - 1)
@. MRSH2(λ, p) = 1 / 2 * p[1] + 1 / λ * (1 / 2 * p[2] + p[3] * 3^(1 - p[4]) / 2 * (2 * λ + 1 / λ^2)^(p[4] - 1))
@. GTSH(λ, p) = p[1] / 2 + 1 / (2 * λ^2 + 1 / λ) * p[2] * 3 / 2 + p[3] * 3^(1 - p[4]) / 2 * (λ^2 + 2 / λ)^(p[4] - 1)
@. GTSH2(λ, p) = p[1] / 2 + 1 / (2 * λ^2 + 1 / λ) * p[2] * 3 / 2 + 1 / λ * p[3] * 3^(1 - p[4]) / 2 * (2 * λ + 1 / λ^2)^(p[4] - 1)
@. CSH(λ, p) = p[1] / 2 + sqrt(3) / 2 * 1 / (λ * sqrt(2 * λ + 1 / λ^2)) * p[2] + p[3] * 3^(1 - p[4]) / 2 * (λ^2 + 2 / λ)^(p[4] - 1)
@. GMR(λ, p) = 1 / 2 * p[1] * 1 / (1 - (λ^2 + 2 / λ - 3) / p[3]) + 1 / 2 * p[2] * 1 / λ
@. GG(λ, p) = 1 / 2 * p[1] * 1 / (1 - (λ^2 + 2 / λ - 3) / p[3]) + 3 / 2 * p[2] * 1 / λ * 1 / (2 * λ + 1 / λ^2)
@. GC(λ, p) = 1 / 2 * p[1] * 1 / (1 - (λ^2 + 2 / λ - 3) / p[3]) + sqrt(3) / 2 * 1 ./ (λ .* sqrt.(2 * λ .+ 1 ./ λ .^ 2)) * p[2]
@. DCMR(λ, p) = p[1] / 6 * (1 + 2 * (1 - (λ^2 + 2 / λ - 3) / p[3])^(-2)) + p[2] * 1 / λ * 1 / 2
@. DCGT(λ, p) = p[1] / 6 * (1 + 2 * (1 - (λ^2 + 2 / λ - 3) / p[3])^(-2)) + 1 / (2 * λ^2 + 1 / λ) * p[2] * 3 / 2

function extract_plot(simu, kmax, α)

    @unpack result = simu

    idx = modedirection(result)

    result.k[real(result.k).>kmax] .= NaN
    shading = 1 .- abs.(imag.(result.k)) / (α * maximum(abs.(imag.(filter(!isnan, result.k)))))
    result.k[shading.<0] .= NaN

    kplot = result.k[:]
    fplot = repeat(result.ω, 1, Int(2 * length(result.udof) * result.N))[:] / (2 * pi)
    idxplot = idx[:]
    idxplot = convert(Vector{Float64}, idxplot)
    splot = shading[:]

    fplot[isnan.(kplot)] .= NaN
    idxplot[isnan.(kplot)] .= NaN
    splot[isnan.(kplot)] .= NaN

    filter!(.!isnan, kplot)
    filter!(.!isnan, fplot)
    filter!(.!isnan, idxplot)
    filter!(.!isnan, splot)

    p = sortperm(real(kplot))
    kplot = kplot[p]
    fplot = fplot[p]
    idxplot = idxplot[p]
    splot = splot[p]

    return kplot, fplot, idxplot, splot

end

function lw_mode_velocity(λ, ρs, C)

    VSHpar = Vector{ComplexF64}()
    VSHperp = Vector{ComplexF64}()
    VSpar = Vector{ComplexF64}()
    VSperp = Vector{ComplexF64}()
    VApar = Vector{ComplexF64}()
    VAperp = Vector{ComplexF64}()

    for i in eachindex(λ)

        append!(VSHpar, sqrt(C[i][1, 3, 3, 1] / ρs))
        append!(VSHperp, sqrt(C[i][3, 1, 1, 3] / ρs))
        append!(VSpar, sqrt((C[i][1, 1, 1, 1] - C[i][1, 1, 2, 2]^2 / C[i][2, 2, 2, 2]) / ρs))
        append!(VSperp, sqrt((C[i][3, 3, 3, 3] - C[i][3, 3, 2, 2]^2 / C[i][2, 2, 2, 2]) / ρs))
        append!(VApar, sqrt((C[i][1, 2, 2, 1] - C[i][2, 1, 1, 2]) / ρs))
        append!(VAperp, sqrt((C[i][3, 2, 2, 3] - C[i][2, 3, 3, 2]) / ρs))

    end

    return VSHpar, VSHperp, VSpar, VSperp, VApar, VAperp

end
"""
    getA0velocities(path, filenames, f)
Get the velocity of the A0 mode in the parallel and perpendicular direction at frequency `f`.
"""
function getA0velocities(path, filenames, f)

    VApar = Vector{Float64}()
    VAperp = Vector{Float64}()

    for name in filenames

        local data = load(joinpath(path, name))
        @unpack kperp, fperp, kpar, fpar = data
        _, idxpar = findmin(abs.(fpar .- f))
        _, idxperp = findmin(abs.(fperp .- f))
        append!(VApar, fpar[idxpar] * 2 * pi / kpar[idxpar])
        append!(VAperp, fperp[idxperp] * 2 * pi / kperp[idxperp])

    end

    return VApar, VAperp

end

function getSHSvelocities(datapar, dataperp, f)

    VSHpar = Vector{Float64}()
    VSpar = Vector{Float64}()
    VSHperp = Vector{Float64}()
    VSperp = Vector{Float64}()

    @unpack kpar, fpar = datapar
    kSH = reverse!(kpar[1:2:end, :], dims=1)
    kSH[isnan.(kSH)] .= 0
    fSH = reverse!(fpar[1:2:end, :], dims=1)
    fSH[isnan.(fSH)] .= 0
    kS = reverse!(kpar[2:2:end, :], dims=1)
    kS[isnan.(kS)] .= 0
    fS = reverse!(fpar[2:2:end, :], dims=1)
    fS[isnan.(fS)] .= 0
    _, idxSH = findmin(abs.(f .- fSH), dims=2)
    _, idxS = findmin(abs.(f .- fS), dims=2)
    for i in axes(kSH, 1)
        append!(VSHpar, fSH[idxSH[i]] * 2 * pi / kSH[idxSH[i]])
        append!(VSpar, fS[idxS[i]] * 2 * pi / kS[idxS[i]])
    end

    @unpack kperp, fperp = dataperp
    kSH = reverse!(kperp[1:2:end, :], dims=1)
    kSH[isnan.(kSH)] .= 0
    fSH = reverse!(fperp[1:2:end, :], dims=1)
    fSH[isnan.(fSH)] .= 0
    kS = reverse!(kperp[2:2:end, :], dims=1)
    kS[isnan.(kS)] .= 0
    fS = reverse!(fperp[2:2:end, :], dims=1)
    fS[isnan.(fS)] .= 0
    _, idxSH = findmin(abs.(f .- fSH), dims=2)
    _, idxS = findmin(abs.(f .- fS), dims=2)
    for i in axes(kSH, 1)
        append!(VSHperp, fSH[idxSH[i]] * 2 * pi / kSH[idxSH[i]])
        append!(VSperp, fS[idxS[i]] * 2 * pi / kS[idxS[i]])
    end

    return VSHpar, VSpar, VSHperp, VSperp
end
