# Define functions to lift the velocity of each mode
function getVApar(f)
    VApar = Vector{Float64}()
    for name in A0files
        data = load(joinpath(@__DIR__, "data_pierre/", name))
        @unpack kpar, fpar = data
        _, idxpar = findmin(abs.(fpar .- to_value(f)[2]))
        append!(VApar, fpar[idxpar] * 2 * pi / kpar[idxpar] / sqrt(μP / ρs))
    end
    return VApar
end

function getVAperp(f)
    VAperp = Vector{Float64}()
    for name in A0files
        data = load(joinpath(@__DIR__, "data_pierre", name))
        @unpack kperp, fperp = data
        _, idxperp = findmin(abs.(fperp .- to_value(f)[2]))
        append!(VAperp, fperp[idxperp] * 2 * pi / kperp[idxperp] / sqrt(μP / ρs))
    end
    return VAperp
end

function getVSHpar(f)
    VSHpar = Vector{Float64}()
    @unpack kpar, fpar = dataApar
    kSH = reverse!(kpar[1:2:end, :], dims=1)
    kSH[isnan.(kSH)] .= 0
    fSH = reverse!(fpar[1:2:end, :], dims=1)
    fSH[isnan.(fSH)] .= 0
    _, idxSH = findmin(abs.(to_value(f)[1] .- fSH), dims=2)
    for i in axes(kSH, 1)
        append!(VSHpar, fSH[idxSH[i]] * 2 * pi / kSH[idxSH[i]] / sqrt(μA / ρs))
    end
    return VSHpar
end

function getVSHperp(f)
    VSHperp = Vector{Float64}()
    @unpack kperp, fperp = dataAperp
    kSH = reverse!(kperp[1:2:end, :], dims=1)
    kSH[isnan.(kSH)] .= 0
    fSH = reverse!(fperp[1:2:end, :], dims=1)
    fSH[isnan.(fSH)] .= 0
    _, idxSH = findmin(abs.(to_value(f)[1] .- fSH), dims=2)
    for i in axes(kSH, 1)
        append!(VSHperp, fSH[idxSH[i]] * 2 * pi / kSH[idxSH[i]] / sqrt(μA / ρs))
    end
    return VSHperp
end

function getVSpar(f)
    VSpar = Vector{Float64}()
    @unpack kpar, fpar = dataApar
    kS = reverse!(kpar[2:2:end, :], dims=1)
    kS[isnan.(kS)] .= 0
    fS = reverse!(fpar[2:2:end, :], dims=1)
    fS[isnan.(fS)] .= 0
    _, idxS = findmin(abs.(to_value(f)[1] .- fS), dims=2)
    for i in axes(kS, 1)
        append!(VSpar, fS[idxS[i]] * 2 * pi / kS[idxS[i]] / sqrt(μA / ρs))
    end
    return VSpar
end

function getVSperp(f)
    VSperp = Vector{Float64}()
    @unpack kperp, fperp = dataAperp
    kS = reverse!(kperp[2:2:end, :], dims=1)
    kS[isnan.(kS)] .= 0
    fS = reverse!(fperp[2:2:end, :], dims=1)
    fS[isnan.(fS)] .= 0
    _, idxS = findmin(abs.(to_value(f)[1] .- fS), dims=2)
    for i in axes(kS, 1)
        append!(VSperp, fS[idxS[i]] * 2 * pi / kS[idxS[i]] / sqrt(μA / ρs))
    end
    return VSperp
end

# Define functions to lift the predicted velocity for each hyperelastic model
# Mooney-Rivlin model
function predicVSparMR(f)
    CMRip = C_MR.(λ, λ .^ (-0.41), ρs, VL, MRfit[1], MRfit[2], to_value(f)[1], τ, n, βMR[1])
    _, _, VSpart, _ = lw_mode_velocity(λ, ρs, CMRip)
    return real(VSpart) / sqrt(μP / ρs)
end

function predicVSperpMR(f)
    CMRip = C_MR.(λ, λ .^ (-0.41), ρs, VL, MRfit[1], MRfit[2], to_value(f)[1], τ, n, βMR[1])
    _, _, _, VSperpt = lw_mode_velocity(λ, ρs, CMRip)
    return real(VSperpt) / sqrt(μP / ρs)
end

function predicVSHparMR(f)
    CMRip = C_MR.(λ, λ .^ (-0.41), ρs, VL, MRfit[1], MRfit[2], to_value(f)[1], τ, n, βMR[1])
    VSHpart, _ = lw_mode_velocity(λ, ρs, CMRip)
    return real(VSHpart) / sqrt(μP / ρs)
end

function predicVSHperpMR(f)
    CMRip = C_MR.(λ, λ .^ (-0.41), ρs, VL, MRfit[1], MRfit[2], to_value(f)[1], τ, n, βMR[1])
    _, VSHperpt = lw_mode_velocity(λ, ρs, CMRip)
    return real(VSHperpt) / sqrt(μP / ρs)
end

function predicVAparMR(f)
    CMRop = C_MR.(λ, λ .^ (-0.39), ρs, VL, MRfit[1], MRfit[2], to_value(f)[2], τ, n, βMR[1])
    _, _, _, _, VApart = lw_mode_velocity(λ, ρs, CMRop)
    return real(VApart) / sqrt(μP / ρs)
end

function predicVAperpMR(fs)
    simu = filter(x -> occursin("MR", x), readdir(joinpath(@__DIR__, "fig4")))
    VAperpt = Vector{Float64}()
    for name in simu
        local data = load(joinpath(@__DIR__, "fig4", name))
        @unpack k, f = data
        idx = argmin(abs.(f .- to_value(fs)[2]))
        append!(VAperpt, real(2 * pi * f[idx] / k[idx]))
    end
    return VAperpt / sqrt(μP / ρs)
end

# Gent-Thomas model
function predicVSparGT(f)
    CGTip = C_GT.(λ, λ .^ (-0.41), ρs, VL, GTfit[1], GTfit[2], to_value(f)[1], τ, n, βGT[1])
    _, _, VSpart, _ = lw_mode_velocity(λ, ρs, CGTip)
    return real(VSpart) / sqrt(μP / ρs)
end

function predicVSperpGT(f)
    CGTip = C_GT.(λ, λ .^ (-0.41), ρs, VL, GTfit[1], GTfit[2], to_value(f)[1], τ, n, βGT[1])
    _, _, _, VSperpt = lw_mode_velocity(λ, ρs, CGTip)
    return real(VSperpt) / sqrt(μP / ρs)
end

function predicVSHparGT(f)
    CGTip = C_GT.(λ, λ .^ (-0.41), ρs, VL, GTfit[1], GTfit[2], to_value(f)[1], τ, n, βGT[1])
    VSHpart, _ = lw_mode_velocity(λ, ρs, CGTip)
    return real(VSHpart) / sqrt(μP / ρs)
end

function predicVSHperpGT(f)
    CGTip = C_GT.(λ, λ .^ (-0.41), ρs, VL, GTfit[1], GTfit[2], to_value(f)[1], τ, n, βGT[1])
    _, VSHperpt = lw_mode_velocity(λ, ρs, CGTip)
    return real(VSHperpt) / sqrt(μP / ρs)
end

function predicVAparGT(f)
    CGTop = C_GT.(λ, λ .^ (-0.39), ρs, VL, GTfit[1], GTfit[2], to_value(f)[2], τ, n, βGT[1])
    _, _, _, _, VApart = lw_mode_velocity(λ, ρs, CGTop)
    return real(VApart) / sqrt(μP / ρs)
end

function predicVAperpGT(fs)
    simu = filter(x -> occursin("GT", x), readdir(joinpath(@__DIR__, "fig4")))
    VAperpt = Vector{Float64}()
    for name in simu
        local data = load(joinpath(@__DIR__, "fig4", name))
        @unpack k, f = data
        idx = argmin(abs.(f .- to_value(fs)[2]))
        append!(VAperpt, real(2 * pi * f[idx] / k[idx]))
    end
    return VAperpt / sqrt(μP / ρs)
end

# Carroll model
function predicVSparC(f)
    CCip = C_C.(λ, λ .^ (-0.41), ρs, VL, Cfit[1], Cfit[2], to_value(f)[1], τ, n, βC[1])
    _, _, VSpart, _ = lw_mode_velocity(λ, ρs, CCip)
    return real(VSpart) / sqrt(μP / ρs)
end

function predicVSperpC(f)
    CCip = C_C.(λ, λ .^ (-0.41), ρs, VL, Cfit[1], Cfit[2], to_value(f)[1], τ, n, βC[1])
    _, _, _, VSperpt = lw_mode_velocity(λ, ρs, CCip)
    return real(VSperpt) / sqrt(μP / ρs)
end

function predicVSHparC(f)
    CCip = C_C.(λ, λ .^ (-0.41), ρs, VL, Cfit[1], Cfit[2], to_value(f)[1], τ, n, βC[1])
    VSHpart, _ = lw_mode_velocity(λ, ρs, CCip)
    return real(VSHpart) / sqrt(μP / ρs)
end

function predicVSHperpC(f)
    CCip = C_C.(λ, λ .^ (-0.41), ρs, VL, Cfit[1], Cfit[2], to_value(f)[1], τ, n, βC[1])
    _, VSHperpt = lw_mode_velocity(λ, ρs, CCip)
    return real(VSHperpt) / sqrt(μP / ρs)
end

function predicVAparC(f)
    CCop = C_C.(λ, λ .^ (-0.39), ρs, VL, Cfit[1], Cfit[2], to_value(f)[2], τ, n, βC[1])
    _, _, _, _, VApart = lw_mode_velocity(λ, ρs, CCop)
    return real(VApart) / sqrt(μP / ρs)
end

function predicVAperpC(fs)
    simu = filter(x -> occursin("C", x), readdir(joinpath(@__DIR__, "fig4")))
    VAperpt = Vector{Float64}()
    for name in simu
        local data = load(joinpath(@__DIR__, "fig4", name))
        @unpack k, f = data
        idx = argmin(abs.(f .- to_value(fs)[2]))
        append!(VAperpt, real(2 * pi * f[idx] / k[idx]))
    end
    return VAperpt / sqrt(μP / ρs)
end

# MRSH model
function predicVSparMRSH(f)
    CMRSHip = C_MRSH.(λ, λ .^ (-0.41), ρs, VL, MRSHfit[1], MRSHfit[2], MRSHfit[3], MRSHfit[4], to_value(f)[1], τ, n, βMRSH[1])
    _, _, VSpart, _ = lw_mode_velocity(λ, ρs, CMRSHip)
    return real(VSpart) / sqrt(μP / ρs)
end

function predicVSperpMRSH(f)
    CMRSHip = C_MRSH.(λ, λ .^ (-0.41), ρs, VL, MRSHfit[1], MRSHfit[2], MRSHfit[3], MRSHfit[4], to_value(f)[1], τ, n, βMRSH[1])
    _, _, _, VSperpt = lw_mode_velocity(λ, ρs, CMRSHip)
    return real(VSperpt) / sqrt(μP / ρs)
end

function predicVSHparMRSH(f)
    CMRSHip = C_MRSH.(λ, λ .^ (-0.41), ρs, VL, MRSHfit[1], MRSHfit[2], MRSHfit[3], MRSHfit[4], to_value(f)[1], τ, n, βMRSH[1])
    VSHpart, _ = lw_mode_velocity(λ, ρs, CMRSHip)
    return real(VSHpart) / sqrt(μP / ρs)
end

function predicVSHperpMRSH(f)
    CMRSHip = C_MRSH.(λ, λ .^ (-0.41), ρs, VL, MRSHfit[1], MRSHfit[2], MRSHfit[3], MRSHfit[4], to_value(f)[1], τ, n, βMRSH[1])
    _, VSHperpt = lw_mode_velocity(λ, ρs, CMRSHip)
    return real(VSHperpt) / sqrt(μP / ρs)
end

function predicVAparMRSH(f)
    CMRSHop = C_MRSH.(λ, λ .^ (-0.39), ρs, VL, MRSHfit[1], MRSHfit[2], MRSHfit[3], MRSHfit[4], to_value(f)[2], τ, n, βMRSH[1])
    _, _, _, _, VApart = lw_mode_velocity(λ, ρs, CMRSHop)
    return real(VApart) / sqrt(μP / ρs)
end

function predicVAperpMRSH(fs)
    simu = filter(x -> occursin("MRSH1", x), readdir(joinpath(@__DIR__, "fig5")))
    VAperpt = Vector{Float64}()
    for name in simu
        local data = load(joinpath(@__DIR__, "fig5", name))
        @unpack k, f = data
        idx = argmin(abs.(f .- to_value(fs)[2]))
        append!(VAperpt, real(2 * pi * f[idx] / k[idx]))
    end
    return VAperpt / sqrt(μP / ρs)
end

# GTSH model
function predicVSparGTSH(f)
    CGTSHip = C_GTSH.(λ, λ .^ (-0.41), ρs, VL, GTSHfit[1], GTSHfit[2], GTSHfit[3], GTSHfit[4], to_value(f)[1], τ, n, βGTSH[1])
    _, _, VSpart, _ = lw_mode_velocity(λ, ρs, CGTSHip)
    return real(VSpart) / sqrt(μP / ρs)
end

function predicVSperpGTSH(f)
    CGTSHip = C_GTSH.(λ, λ .^ (-0.41), ρs, VL, GTSHfit[1], GTSHfit[2], GTSHfit[3], GTSHfit[4], to_value(f)[1], τ, n, βGTSH[1])
    _, _, _, VSperpt = lw_mode_velocity(λ, ρs, CGTSHip)
    return real(VSperpt) / sqrt(μP / ρs)
end

function predicVSHparGTSH(f)
    CGTSHip = C_GTSH.(λ, λ .^ (-0.41), ρs, VL, GTSHfit[1], GTSHfit[2], GTSHfit[3], GTSHfit[4], to_value(f)[1], τ, n, βGTSH[1])
    VSHpart, _ = lw_mode_velocity(λ, ρs, CGTSHip)
    return real(VSHpart) / sqrt(μP / ρs)
end

function predicVSHperpGTSH(f)
    CGTSHip = C_GTSH.(λ, λ .^ (-0.41), ρs, VL, GTSHfit[1], GTSHfit[2], GTSHfit[3], GTSHfit[4], to_value(f)[1], τ, n, βGTSH[1])
    _, VSHperpt = lw_mode_velocity(λ, ρs, CGTSHip)
    return real(VSHperpt) / sqrt(μP / ρs)
end

function predicVAparGTSH(f)
    CGTSHop = C_GTSH.(λ, λ .^ (-0.39), ρs, VL, GTSHfit[1], GTSHfit[2], GTSHfit[3], GTSHfit[4], to_value(f)[2], τ, n, βGTSH[1])
    _, _, _, _, VApart = lw_mode_velocity(λ, ρs, CGTSHop)
    return real(VApart) / sqrt(μP / ρs)
end

function predicVAperpGTSH(fs)
    simu = filter(x -> occursin("GTSH1", x), readdir(joinpath(@__DIR__, "fig5")))
    VAperpt = Vector{Float64}()
    for name in simu
        local data = load(joinpath(@__DIR__, "fig5", name))
        @unpack k, f = data
        idx = argmin(abs.(f .- to_value(fs)[2]))
        append!(VAperpt, real(2 * pi * f[idx] / k[idx]))
    end
    return VAperpt / sqrt(μP / ρs)
end

# CSH model
function predicVSparCSH(f)
    CCSHip = C_CSH.(λ, λ .^ (-0.41), ρs, VL, CSHfit[1], CSHfit[2], CSHfit[3], CSHfit[4], to_value(f)[1], τ, n, βCSH[1])
    _, _, VSpart, _ = lw_mode_velocity(λ, ρs, CCSHip)
    return real(VSpart) / sqrt(μP / ρs)
end

function predicVSperpCSH(f)
    CCSHip = C_CSH.(λ, λ .^ (-0.41), ρs, VL, CSHfit[1], CSHfit[2], CSHfit[3], CSHfit[4], to_value(f)[1], τ, n, βCSH[1])
    _, _, _, VSperpt = lw_mode_velocity(λ, ρs, CCSHip)
    return real(VSperpt) / sqrt(μP / ρs)
end

function predicVSHparCSH(f)
    CCSHip = C_CSH.(λ, λ .^ (-0.41), ρs, VL, CSHfit[1], CSHfit[2], CSHfit[3], CSHfit[4], to_value(f)[1], τ, n, βCSH[1])
    VSHpart, _ = lw_mode_velocity(λ, ρs, CCSHip)
    return real(VSHpart) / sqrt(μP / ρs)
end

function predicVSHperpCSH(f)
    CCSHip = C_CSH.(λ, λ .^ (-0.41), ρs, VL, CSHfit[1], CSHfit[2], CSHfit[3], CSHfit[4], to_value(f)[1], τ, n, βCSH[1])
    _, VSHperpt = lw_mode_velocity(λ, ρs, CCSHip)
    return real(VSHperpt) / sqrt(μP / ρs)
end

function predicVAparCSH(f)
    CCSHop = C_CSH.(λ, λ .^ (-0.39), ρs, VL, CSHfit[1], CSHfit[2], CSHfit[3], CSHfit[4], to_value(f)[2], τ, n, βCSH[1])
    _, _, _, _, VApart = lw_mode_velocity(λ, ρs, CCSHop)
    return real(VApart) / sqrt(μP / ρs)
end

function predicVAperpCSH(fs)
    simu = filter(x -> occursin("CSH", x), readdir(joinpath(@__DIR__, "fig5")))
    VAperpt = Vector{Float64}()
    for name in simu
        local data = load(joinpath(@__DIR__, "fig5", name))
        @unpack k, f = data
        idx = argmin(abs.(f .- to_value(fs)[2]))
        append!(VAperpt, real(2 * pi * f[idx] / k[idx]))
    end
    return VAperpt / sqrt(μP / ρs)
end

# GMR model
function predicVSparGMR(f)
    CGMRip = C_GMR.(λ, λ .^ (-0.41), ρs, VL, GMRfit[1], GMRfit[2], GMRfit[3], to_value(f)[1], τ, n, βGMR[1])
    _, _, VSpart, _ = lw_mode_velocity(λ, ρs, CGMRip)
    return real(VSpart) / sqrt(μP / ρs)
end

function predicVSperpGMR(f)
    CGMRip = C_GMR.(λ, λ .^ (-0.41), ρs, VL, GMRfit[1], GMRfit[2], GMRfit[3], to_value(f)[1], τ, n, βGMR[1])
    _, _, _, VSperpt = lw_mode_velocity(λ, ρs, CGMRip)
    return real(VSperpt) / sqrt(μP / ρs)
end

function predicVSHparGMR(f)
    CGMRip = C_GMR.(λ, λ .^ (-0.41), ρs, VL, GMRfit[1], GMRfit[2], GMRfit[3], to_value(f)[1], τ, n, βGMR[1])
    VSHpart, _ = lw_mode_velocity(λ, ρs, CGMRip)
    return real(VSHpart) / sqrt(μP / ρs)
end

function predicVSHperpGMR(f)
    CGMRip = C_GMR.(λ, λ .^ (-0.41), ρs, VL, GMRfit[1], GMRfit[2], GMRfit[3], to_value(f)[1], τ, n, βGMR[1])
    _, VSHperpt = lw_mode_velocity(λ, ρs, CGMRip)
    return real(VSHperpt) / sqrt(μP / ρs)
end

function predicVAparGMR(f)
    CGMRop = C_GMR.(λ, λ .^ (-0.39), ρs, VL, GMRfit[1], GMRfit[2], GMRfit[3], to_value(f)[2], τ, n, βGMR[1])
    _, _, _, _, VApart = lw_mode_velocity(λ, ρs, CGMRop)
    return real(VApart) / sqrt(μP / ρs)
end

function predicVAperpGMR(fs)
    simu = filter(x -> occursin("GMR", x), readdir(joinpath(@__DIR__, "fig6")))
    VAperpt = Vector{Float64}()
    for name in simu
        local data = load(joinpath(@__DIR__, "fig6", name))
        @unpack k, f = data
        idx = argmin(abs.(f .- to_value(fs)[2]))
        append!(VAperpt, real(2 * pi * f[idx] / k[idx]))
    end
    return VAperpt / sqrt(μP / ρs)
end

# GG model
function predicVSparGG(f)
    CGGip = C_GG.(λ, λ .^ (-0.41), ρs, VL, GGfit[1], GGfit[2], GGfit[3], to_value(f)[1], τ, n, βGG[1])
    _, _, VSpart, _ = lw_mode_velocity(λ, ρs, CGGip)
    return real(VSpart) / sqrt(μP / ρs)
end

function predicVSperpGG(f)
    CGGip = C_GG.(λ, λ .^ (-0.41), ρs, VL, GGfit[1], GGfit[2], GGfit[3], to_value(f)[1], τ, n, βGG[1])
    _, _, _, VSperpt = lw_mode_velocity(λ, ρs, CGGip)
    return real(VSperpt) / sqrt(μP / ρs)
end

function predicVSHparGG(f)
    CGGip = C_GG.(λ, λ .^ (-0.41), ρs, VL, GGfit[1], GGfit[2], GGfit[3], to_value(f)[1], τ, n, βGG[1])
    VSHpart, _ = lw_mode_velocity(λ, ρs, CGGip)
    return real(VSHpart) / sqrt(μP / ρs)
end

function predicVSHperpGG(f)
    CGGip = C_GG.(λ, λ .^ (-0.41), ρs, VL, GGfit[1], GGfit[2], GGfit[3], to_value(f)[1], τ, n, βGG[1])
    _, VSHperpt = lw_mode_velocity(λ, ρs, CGGip)
    return real(VSHperpt) / sqrt(μP / ρs)
end

function predicVAparGG(f)
    CGGop = C_GG.(λ, λ .^ (-0.39), ρs, VL, GGfit[1], GGfit[2], GGfit[3], to_value(f)[2], τ, n, βGG[1])
    _, _, _, _, VApart = lw_mode_velocity(λ, ρs, CGGop)
    return real(VApart) / sqrt(μP / ρs)
end

function predicVAperpGG(fs)
    simu = filter(x -> occursin("GG", x), readdir(joinpath(@__DIR__, "fig6")))
    VAperpt = Vector{Float64}()
    for name in simu
        local data = load(joinpath(@__DIR__, "fig6", name))
        @unpack k, f = data
        idx = argmin(abs.(f .- to_value(fs)[2]))
        append!(VAperpt, real(2 * pi * f[idx] / k[idx]))
    end
    return VAperpt / sqrt(μP / ρs)
end

# DCMR model
function predicVSparDCMR(f)
    CDCMRip = C_DCMR2.(λ, λ .^ (-0.41), ρs, VL, DCMRfit[1], DCMRfit[2], DCMRfit[3], to_value(f)[1], τ, n, βDCMR[1])
    _, _, VSpart, _ = lw_mode_velocity(λ, ρs, CDCMRip)
    return real(VSpart) / sqrt(μP / ρs)
end

function predicVSperpDCMR(f)
    CDCMRip = C_DCMR2.(λ, λ .^ (-0.41), ρs, VL, DCMRfit[1], DCMRfit[2], DCMRfit[3], to_value(f)[1], τ, n, βDCMR[1])
    _, _, _, VSperpt = lw_mode_velocity(λ, ρs, CDCMRip)
    return real(VSperpt) / sqrt(μP / ρs)
end

function predicVSHparDCMR(f)
    CDCMRip = C_DCMR2.(λ, λ .^ (-0.41), ρs, VL, DCMRfit[1], DCMRfit[2], DCMRfit[3], to_value(f)[1], τ, n, βDCMR[1])
    VSHpart, _ = lw_mode_velocity(λ, ρs, CDCMRip)
    return real(VSHpart) / sqrt(μP / ρs)
end

function predicVSHperpDCMR(f)
    CDCMRip = C_DCMR2.(λ, λ .^ (-0.41), ρs, VL, DCMRfit[1], DCMRfit[2], DCMRfit[3], to_value(f)[1], τ, n, βDCMR[1])
    _, VSHperpt = lw_mode_velocity(λ, ρs, CDCMRip)
    return real(VSHperpt) / sqrt(μP / ρs)
end

function predicVAparDCMR(f)
    CDCMRop = C_DCMR2.(λ, λ .^ (-0.39), ρs, VL, DCMRfit[1], DCMRfit[2], DCMRfit[3], to_value(f)[2], τ, n, βDCMR[1])
    _, _, _, _, VApart = lw_mode_velocity(λ, ρs, CDCMRop)
    return real(VApart) / sqrt(μP / ρs)
end

function predicVAperpDCMR(fs)
    simu = filter(x -> occursin("DCMR", x), readdir(joinpath(@__DIR__, "fig6")))
    VAperpt = Vector{Float64}()
    for name in simu
        local data = load(joinpath(@__DIR__, "fig6", name))
        @unpack k, f = data
        idx = argmin(abs.(f .- to_value(fs)[2]))
        append!(VAperpt, real(2 * pi * f[idx] / k[idx]))
    end
    return VAperpt / sqrt(μP / ρs)
end

# DCGT model
function predicVSparDCGT(f)
    CDCGTip = C_DCGT.(λ, λ .^ (-0.41), ρs, VL, DCGTfit[1], DCGTfit[2], DCGTfit[3], to_value(f)[1], τ, n, βDCGT[1])
    _, _, VSpart, _ = lw_mode_velocity(λ, ρs, CDCGTip)
    return real(VSpart) / sqrt(μP / ρs)
end

function predicVSperpDCGT(f)
    CDCGTip = C_DCGT.(λ, λ .^ (-0.41), ρs, VL, DCGTfit[1], DCGTfit[2], DCGTfit[3], to_value(f)[1], τ, n, βDCGT[1])
    _, _, _, VSperpt = lw_mode_velocity(λ, ρs, CDCGTip)
    return real(VSperpt) / sqrt(μP / ρs)
end

function predicVSHparDCGT(f)
    CDCGTip = C_DCGT.(λ, λ .^ (-0.41), ρs, VL, DCGTfit[1], DCGTfit[2], DCGTfit[3], to_value(f)[1], τ, n, βDCGT[1])
    VSHpart, _ = lw_mode_velocity(λ, ρs, CDCGTip)
    return real(VSHpart) / sqrt(μP / ρs)
end

function predicVSHperpDCGT(f)
    CDCGTip = C_DCGT.(λ, λ .^ (-0.41), ρs, VL, DCGTfit[1], DCGTfit[2], DCGTfit[3], to_value(f)[1], τ, n, βDCGT[1])
    _, VSHperpt = lw_mode_velocity(λ, ρs, CDCGTip)
    return real(VSHperpt) / sqrt(μP / ρs)
end

function predicVAparDCGT(f)
    CDCGTop = C_DCGT.(λ, λ .^ (-0.39), ρs, VL, DCGTfit[1], DCGTfit[2], DCGTfit[3], to_value(f)[2], τ, n, βDCGT[1])
    _, _, _, _, VApart = lw_mode_velocity(λ, ρs, CDCGTop)
    return real(VApart) / sqrt(μP / ρs)
end

function predicVAperpDCGT(fs)
    simu = filter(x -> occursin("DCGT", x), readdir(joinpath(@__DIR__, "fig6")))
    VAperpt = Vector{Float64}()
    for name in simu
        local data = load(joinpath(@__DIR__, "fig6", name))
        @unpack k, f = data
        idx = argmin(abs.(f .- to_value(fs)[2]))
        append!(VAperpt, real(2 * pi * f[idx] / k[idx]))
    end
    return VAperpt / sqrt(μP / ρs)
end
