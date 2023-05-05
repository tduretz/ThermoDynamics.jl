using CairoMakie
import LinearAlgebra:norm

function main()

    Pmin = 1e5
    Pmax = 5e9
    P    = LinRange(Pmin, Pmax, 20)
    V0   = 1513.04 
    K0   = 172.6e9
    K01  = 4.0

    V    = abs.(V0.*10 .*rand(size(P)))

    BM   = 3/2*K0.*((V0./V).^(7/3) .- (V0./V).^(5/3)) #.* (1 .+ (3/4) .*(K01-4) .* ((V0./V).^(2/3) .- 1) )
    R    = P .- BM
    norm(BM)
    @show maximum(BM)
    @show ((V0./V).^(7/3) .- (V0./V).^(5/3))
end

main()