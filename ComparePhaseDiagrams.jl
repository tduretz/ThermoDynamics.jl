using GLMakie, FileIO, ImageCore, Printf

function ReadBinFileVisualise(fname, nx, ny)
    data = zeros(Float64, nx, ny)
    ρ = read!(fname, data)
    return ρ 
end

function main()
    Tmin    = 298+1e-3+100;               # K
    Tmax    = 800+273;
    Pmin    = -101325;                  # Pa
    Pmax    = 50/10*1e9;
    P       = LinRange(Pmin, Pmax, 2500)
    T       = LinRange(Tmin, Tmax, 675)
    #   ρ        = ReadBinFileVisualise("SiO2.dat")
 
    # Read in from MD7
    path  = "/Users/tduretz/REPO/MDOODZ7.0/IMPORT/"
    ρ_nsm010 = ReadBinFileVisualise(path * "SiO2_nsm010.dat", 675, 2500)
    ρ_Gerya  = ReadBinFileVisualise(path * "SiO2_Gerya.dat", 675, 2500)
    # Read in from here
    ρ_old    = ReadBinFileVisualise("SiO2_old_lappy.dat", 675, 2500)
    ρ_cindy  = ReadBinFileVisualise("SiO2_revised.dat", 675, 2500)
    ρ_jl1    = ReadBinFileVisualise("SiO2_julia_revised_v1.dat", 675, 2500)
    ρ_jl2    = ReadBinFileVisualise("SiO2_julia_revised_v2.dat", 675, 2500)
    ρ_jl3    = ReadBinFileVisualise("SiO2_julia_revised_v3.dat", 675, 2500)
    ρ_HP98   = ReadBinFileVisualise("SiO2_julia_revised_v4_HP98.dat", 675, 2500)
    ρ_HP11   = ReadBinFileVisualise("SiO2_julia_revised_v4_HP11.dat", 675, 2500)


    # # Diffusion smoothing
    # K              = 1.0
    # dt             = min(dP,dT)^2/K/4.2
    # nt             = 10
    # ρ_Gerya_smooth = copy(ρ_Gerya)
    # for it=1:nt
    #     dρdT = diff(ρ_Gerya_smooth[:,2:end-1], dims=1) / dT
    #     d2ρdT2 = diff(dρdT, dims=1) / dT
    #     dρdP = diff(ρ_Gerya_smooth[2:end-1,:], dims=2) / dP
    #     d2ρdP2 = diff(dρdP, dims=2) / dP
    #     ρ_Gerya_smooth[2:end-1,2:end-1] .+= dt.*K.*(d2ρdT2 .+ d2ρdP2)
    # end

    # heatmap(T, P, ρ_nsm010' .- ρ_Gerya_smooth' ) 
    # heatmap(T, P, ρ_Gerya' .- ρ_Gerya_smooth' ) 
    # p = plot(xlabel="P [kbar]", ylabel="ρ [kg/m³]")
    # plot!( ρ_nsm010[ 400, :], label="HP")
    # plot!( ρ_Gerya_smooth[ 400, :], label="Gerya (2004)")

    f = Figure(resolution = (800, 600))
    ax = Axis(f[1, 1],     
        title = "@T = $(round(T[600], digits=2)) K",
        ylabel = "ρ [kg/m³]",
        xlabel = "P [kbar]",
    )
    # lines!(ax, P, ρ_nsm010[ 600, :], label="Holland & Powell (1998) nsm010")

    lines!(ax, P, ρ_old[ 600, :], label="Holland & Powell (1998) old", linestyle=:dash)

    # # # lines!(ax, P, ρ_jl1[ 600, :], label="Holland & Powell (1998)")
    # # lines!(ax, P, ρ_cindy[ 600, :], label="Holland & Powell (1998) cindy")
    # # lines!(ax, P, ρ_jl2[ 600, :], label="Holland & Powell (1998)", linestyle=:dash)

    # # lines!(ax, P, ρ_jl3[ 600, :], label="Holland & Powell (1998) jl3", linestyle=:dash)
    lines!(ax, P, ρ_HP98[ 600, :], label="Holland & Powell (1998)", linestyle=:dashdot)
    lines!(ax, P, ρ_HP11[ 600, :], label="Holland & Powell (2011)", linestyle=:dashdot)

    # lines!(ax, P, ρ_Gerya[ 600, :], label="Gerya (2004)")
    GLMakie.ylims!(ax, 2000, 3100)
    axislegend(framevisible = false, position = :lt)
    current_figure()
    DataInspector(f)
    display(f)

    # save("QuartzCoesiteHP98HP11.png", f, px_per_unit = 4)

end

main()