#%% Guided wave dispersion in a stretched viscoelastic plate using SCM 
# Load packages
# Easy identity matrix and kronecker tensor notation
using FileIO, UnPack, MAT
using LinearAlgebra, Kronecker
# POnG
using POnG
# Plotting 
using CairoMakie

## Problem definition: We look for Lamb waves propagating in direction 1 in a plate with thickness h
# 
#  e2        ---------------------------------- top
#  |_ e1      solid: μs, λs, ρs, h
#            ---------------------------------- bottom
#
# Choose degrees of freedom
udof = 1:3 # 1:2 -> Lamb waves only, 3 -> SH waves 

# Load parameters
data = load(joinpath(@__DIR__, "parameters.jld2"))
@unpack ρs, VL, τ, n = data
h = 2.5e-3
data = load(joinpath(@__DIR__, "fitparameters.jld2"))
@unpack MRfit, βMR, GTfit, βGT, Cfit, βC = data

# %% 
@time for λ1 in 1:0.1:2.5
    #λ1 = 1.10
    λ3 = λ1^(-0.5)
    # update geometry (λ1λ2λ3 = 1)
    hdef = h / (λ1 * λ3)

    # Discretization: We differentiate in direction 2 using SCM.
    N = 20 # number of nodes
    _, D = chebdif(N, 2)
    D2 = -2 * (hdef / hdef) * D[:, :, 1]
    D22 = 4 * (hdef / hdef)^2 * D[:, :, 2]

    ## Solve polynomial eigenvalue problem
    # Impose ω and get eigenvalues and eigenvectors
    ω = 2 * pi * range(1, 200, length=70)
    k = Array{ComplexF64,2}(undef, length(ω), 2 * length(udof) * N)
    u = Array{ComplexF64,3}(undef, length(ω), length(udof) * N, 2 * length(udof) * N)

    for i in axes(ω, 1)

        # C = C_MR(λ3, λ1, ρs, VL, MRfit[1], MRfit[2], ω[i] / (2 * pi), τ, n, βMR[1])
        # μs = sum(MRfit)
        C = C_GT(λ3, λ1, ρs, VL, GTfit[1], GTfit[2], ω[i] / (2 * pi), τ, n, βGT[1])
        μs = sum(GTfit)
        # C = C_C(λ3, λ1, ρs, VL, Cfit[1], Cfit[2], ω[i] / (2 * pi), τ, n, βC[1])
        # μs = sum(Cfit)
        C = C ./ μs

        # Material matrices used in the elastodynamic equation
        C11 = C[1, udof, udof, 1]
        C12 = C[1, udof, udof, 2]
        C21 = C[2, udof, udof, 1]
        C22 = C[2, udof, udof, 2]
        M = I(length(udof))

        # Discretized matrices to build the polynomial eigenvalue problem
        # (ik)^2 Lkk + (ik) Lk + L0 + ω^2 Md =0
        Lkk = C11 ⊗ I(N)
        Lkk = collect(Lkk)
        Lk = (C12 .+ C21) ⊗ D2
        Lk = collect(Lk)
        L0 = C22 ⊗ D22
        L0 = collect(L0)
        Md = M ⊗ I(N)
        Md = collect(Md)
        # Top traction operator : 1 for BC in top, N for BC bottom
        Bkt = C21 ⊗ I(N)[1, :]'
        Bkt = collect(Bkt)
        B0t = C22 ⊗ D2[1, :]'
        B0t = collect(B0t)
        # Bottom traction operator
        Bkb = C21 ⊗ I(N)[N, :]'
        Bkb = collect(Bkb)
        B0b = C22 ⊗ D2[N, :]'
        B0b = collect(B0b)

        # Impose BC: location of boundary vectors in discretized matrix
        bc = [(0:length(udof)-1) .* N .+ 1, (1:length(udof)) .* N]
        # top BC 
        Lkk[bc[1], :] .= 0
        Lk[bc[1], :] = Bkt
        L0[bc[1], :] = B0t
        Md[bc[1], :] .= 0
        # bottom BC 
        Lkk[bc[2], :] .= 0
        Lk[bc[2], :] = Bkb
        L0[bc[2], :] = B0b
        Md[bc[2], :] .= 0

        ## Solve polynomial eigenvalue problem in k to get complex wavenumber
        # create and solve polynomial eigenvalue problem
        temp = ω[i] ./ sqrt(μs / ρs) .* hdef
        problem = PEP([temp^2 .* Md .+ L0, Lk, Lkk])
        ktmp, utmp = polyeig(problem)
        k[i, :] = -1im * ktmp' / hdef
        u[i, :, :] = utmp

    end

    # Store result of a 1D calculation before (and for) postprocessing
    result = result1D(udof, N, ω, k, u)
    # Determine mode symmetry and principal displacement direction
    idx = modedirection(result)
    sym = modesymmetry(result)

    # remove negative real(k) and purely imaginary solutions
    k[(0.49*pi.<angle.(k)).|(angle.(k).<-0.49*pi)] .= NaN
    # remove large spurious k to improve shading
    k[real(k).>400] .= NaN
    shading = 1 .- abs.(imag.(k[:])) / (0.1 * maximum(abs.(imag.(filter(!isnan, k[:])))))

    kplot = k[:]
    idxplot = idx[:]
    idxplot = convert(Vector{Float64}, idxplot);
    symplot = sym[:]
    fplot = repeat(ω, 1, Int(2 * length(udof) * N))[:] / (2 * pi)
    splot = shading[:]

    fplot[isnan.(kplot)] .= NaN
    idxplot[isnan.(kplot)] .= NaN
    splot[isnan.(kplot)] .= NaN

    filter!(!isnan, fplot)
    filter!(!isnan, kplot)
    filter!(!isnan, idxplot)
    filter!(!isnan, splot)

    ksave = real(kplot[idxplot.==2 .&& splot .> 0])
    fsave = fplot[idxplot.==2 .&& splot .> 0]

    p = sortperm(ksave)
    ksave = ksave[p]
    fsave = fsave[p]

    with_theme(theme_latexfonts()) do 
        set_theme!(Theme(palette = (color = [:red, :blue, :green],),))

        fig = Figure(size = (1.2*360, 360))

            ax = Axis(fig[1,1])
            ax.limits = (0, 500, 0, 200)
            ax.xlabel = L"k"
            ax.ylabel = L"f"
            ax.xlabelsize = 20 
            ax.ylabelsize = 20
            ax.title = "Lamb modes"
            scatter!(ksave, fsave)
            colsize!(fig.layout, 1, Aspect(1, 1.2))

            resize_to_layout!(fig)
        display(fig)
    end

    save(joinpath(@__DIR__, string("GT_perp_", h, "m_", λ1, ".jld2")), "k", ksave, "f", fsave)

end

