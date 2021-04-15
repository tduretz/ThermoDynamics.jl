using JuMP, GLPK
using Plots
using LaTeXStrings
#----------------------------------------------------------#
# This codes combines elements of the thermodynamics programming classes of Yury Podladchikov and from Annelore Bessat's Ph.D.
#----------------------------------------------------------#
@views function GibbsEnergy(T2,P2,a,b,c,d,Vref,Pref,Tref,a0,Kref,dSref,dHref,dGL,Tc0,Smax,Vmax)
    P2       = P2./1e5/1e3   # From Pa to kbar
    Cp       = a  .+ b.*T2 .+ c.*T2.^(-2) + d.*T2.^(-1/2);                                # J/K   - extensive
    Tc       = Tc0 .+ Vmax./Smax.*(P2.-Pref);
    if dGL == 1
        Cp      = Cp .+ (Smax.*T2) ./(2*sqrt.(Tc) .* sqrt.( (Tc-T2).* (T2.<=Tc) )) ;      # dG Landau for SiO2 (second order phase transitions)
    end
    V1T      = Vref.*(1 .+     a0.*(T2.-Tref) .- 20.0 .*a0 .*(sqrt.(T2) .- sqrt.(Tref))); # J/bar - extensive
    KT       = Kref.*(1 .- 1.5e-4.*(T2.-Tref));
    V1       = V1T .* (1 .- 4.0 .*P2./(KT .+ 4.0 .*P2)).^(1/4);                           # Murnaghan EOS - Holland - Powell (1998)
    int_Cp   = zeros(size(P2))
    int_Cp_T = zeros(size(P2))
    int_V    = zeros(size(P2))
    if dGL == 0
        @. int_Cp   = (2.0 *sqrt(T2).*d + T2.^2 *b/2 + T2.*a - c./T2) - (2.0 *sqrt(Tref).*d + Tref.^2 *b/2 + Tref.*a - c./Tref); 
        @. int_Cp_T = -0.5*T2.^(-2.0).*c - 2.0*1 /sqrt(T2).*d + 1.0*T2.^1.0.*b + 1.0*a.*log(T2) - (-0.5*Tref.^(-2.0).*c - 2.0*1 /sqrt(Tref).*d + 1.0*Tref.^1.0.*b + 1.0*a.*log(Tref)); 
    else
        @. int_Cp = (a*T2+(1.0/2.0)*b*T2.^2.0-c./T2+2*d*sqrt(T2)+(1.0/2.0)*Smax*((2.0/3.0)*(Tc-T2).^(3.0/2.0)-2.0*Tc.*sqrt(Tc-T2))./sqrt(Tc)) - (a*Tref+(1.0/2.0)*b*Tref.^2-c./Tref+2.0*d*sqrt(Tref)+(1.0/2.0)*Smax*((2.0/3.0)*(Tc-Tref).^(3.0/2.0)-2.0*Tc.*sqrt(Tc-Tref))./sqrt(Tc))
        @. int_Cp_T = (-2.0*d./sqrt(T2)-Smax*sqrt(Tc-T2)./sqrt(Tc)+b*T2-(1.0/2.0)*c./T2.^2.0+a.*log(T2)) - (-2.0*d./sqrt(Tref)-Smax*sqrt(Tc-Tref)./sqrt(Tc)+b*Tref-(1.0/2.0)*c./Tref.^2.0+a.*log(Tref))
    end
    @. int_V    = 1/3*(KT+4*P2).*V1T.*(KT./(KT+4*P2)).^(1.0/4.0) - 1.0/3.0*(KT+4*Pref).*V1T.*(KT./(KT+4*Pref)).^(1/4); 
    dH_T     = dHref .+ int_Cp;             
    dS_T     = dSref .+ int_Cp_T;           
    dG       = dH_T  .- T2.*dS_T  .+ int_V; 
    if dGL == 1 #  dG Landau for SiO2 and Al2SiO5 (second order phase transitions) -- validated with Hans
        Q_ref   = ( (1.0 .- Tref./Tc0).* (T2.<=Tc) ).^(1/4);
        Q       = ( (1.0 .-   T2./Tc ).* (T2.<=Tc) ).^(1/4);
        H_ref   = Smax.*Tc0.*(Q_ref.^2 .- (1/3).*Q_ref.^6);
        S_ref   = Smax.*Q_ref.^2.0;
        V_t     = Vmax.*Q_ref.^2.0 .* (1.0 .+ a0.*(T2 .- Tref) .- 20.0 .* a0 .* (sqrt.(T2) .-  sqrt.(Tref)));
        int_VdP = (1.0/3.0) .* V_t .* KT .* ((1.0 .+ (4.0 .* ( P2 .- Pref ))./KT).^(3.0/4.0) .- 1.0);               
        GL      = Smax .* ( (T2 .- Tc) .* Q.^2.0 + (1.0/3.0).*Tc.*Q.^6.0);
        Gex     = H_ref .- T2 .* S_ref + int_VdP + GL;
        dG      = dG + Gex .* (T2.<=Tc);
    end
    return dG
end
#----------------------------------------------------------#
@views function PhysicalProperties(T2,P2,dG,mmol,Pref,Tref,m0,m1,m2,k0,k1,k2,dP,dT)
    dG       = dG/mmol*1e3                                       # from kJ/mol to J/kg
    rho      =  1.0 ./ ( diff(dG, dims=2)/dP)
    rho_tcax = 0.5*(rho[1:end-1,:] + rho[2:end,:])
    rho_tcay = 0.5*(rho[:,1:end-1] + rho[:,2:end])
    alp      = -1.0 ./ rho_tcax.*diff(rho,dims=1)/dT         # thermodynamically consistant alpha at cst P
    bet      =  1.0 ./ rho_tcay.*diff(rho,dims=2)/dP         # thermodynamically consistant beta at cst T
    G        = m0 .+ m1.*(P2 .- Pref)./1e5 .+ m2.*(T2 .- Tref)
    K        = k0 .+ k1.*(P2 .- Pref)./1e5 .+ k2.*(T2 .- Tref)
    return K, G, alp, bet, rho 
end
#----------------------------------------------------------#
function AluminoSilicates()
    # Parameters - all thermodynamic parameters are taken from hp-ds55.txt database 
    nb_ph   = 2;
    Tref    = 298;                    # K
    Pref    = 1e-3;                   # kbar        
    R       = 8.3144;                 # J/K/mol
    dT      = 10*1
    dP      = 10*0.01*1e5*1e3
    T       = 373:dT:1273
    P       = 1e5:dP:1e9
    # Fake ngrid for P-T space  ;)
    T2      = repeat(T, 1, length(P))
    P2      = repeat(P, 1, length(T))'
    # ANDALUSITE
    # Cp - andalousite
    a       =  0.2773;     
    b       = -0.000006588;
    c       = -1914.1;
    d       = -2.2656;
    # V - andalousite
    Vref    = 5.153;      
    a0      = 0.0000411;      
    Kref    = 1334;      
    dSref   = 0.09270;   
    dHref   = -2588.77; 
    mmol    = 162.05e-3;
    dGL     = 0
    Tc0     = 0.0
    Smax    = 0.0
    Vmax    = 0.0
    dG_and  = GibbsEnergy(T2,P2,a,b,c,d,Vref,Pref,Tref,a0,Kref,dSref,dHref,dGL,Tc0,Smax,Vmax)
    K,G,alp_and,bet_and,rho_and = PhysicalProperties(T2,P2,dG_and,mmol,Pref,Tref,0,0,0,0,0,0,dP,dT)
    # ---------- Kyanite ---------- #
    # Cp
    a       =  0.2794;     
    b       = -0.000007124;
    c       = -2055.6;
    d       = -2.2894;
    # V
    Vref    = 4.4140;     
    a0      = 0.0000404;
    Kref    = 1590;       
    dSref   = 0.08350;     
    dHref   = -2593.11; 
    mmol    = 162.05e-3;
    dGL     = 0
    Tc0     = 0.0
    Smax    = 0.0
    Vmax    = 0.0
    # Gibbs energy
    dG_ky    = GibbsEnergy(T2,P2,a,b,c,d,Vref,Pref,Tref,a0,Kref,dSref,dHref,dGL,Tc0,Smax,Vmax)
    K,G,alp_ky,bet_ky,rho_ky = PhysicalProperties(T2,P2,dG_ky,mmol,Pref,Tref,0,0,0,0,0,0,dP,dT)
    # ---------- Sillimanite ---------- #
    # Cp
    a       =  0.2802;     
    b       = -0.000006900;
    c       = -1375.7;
    d       = -2.3994;
    # V
    Vref    = 4.9860;    
    a0      = 0.0000221;       
    Kref    = 1320.00;      
    dSref   = 0.09550;     
    dHref   = -2585.68;  
    mmol    = 162.05e-3;
    # dG Landau
    dGL     = 1
    Tc0     = 2200;
    Smax    = 0.00400;
    Vmax    = 0.0350;
    # Gibbs energy
    dG_sill  = GibbsEnergy(T2,P2,a,b,c,d,Vref,Pref,Tref,a0,Kref,dSref,dHref,dGL,Tc0,Smax,Vmax)
    K,G,alp_sill,bet_sill,rho_sill = PhysicalProperties(T2,P2,dG_sill,mmol,Pref,Tref,0,0,0,0,0,0,dP,dT)
    # Phase diagram
    Stab                = zeros(length(T),length(P));
    Stab[(dG_and .<dG_ky ) .& (dG_and  .< dG_sill)] .= 1;
    Stab[(dG_ky  .<dG_and) .& (dG_ky   .< dG_sill)] .= 2;
    Stab[(dG_sill.<dG_and) .& (dG_sill .< dG_ky  )] .= 3;
    # VERIFICATION: compare numerically obtained phase diagram with analytic phase diagram boundaries
    # and - ky --- HP 98
    Vand       = 5.153;  #342.25 * 1e-24 / 4 * 6.022*1e23 * 0.1; # conversion unit cell avogadro...
    Vky        = 4.4140; #293.41 * 1e-24 / 4 * 6.022*1e23 * 0.1;
    Vsill      = 4.9860; 
    Sand       = 92.70;
    Sky        = 83.50;
    Ssill      = 95.50;
    Gand       = -2440.97e3;
    Gky        = -2442.59e3;
    Gsill      = -2438.93e3;
    dV         = Vky-Vand;
    dS         = Sky-Sand;
    dG         = Gky-Gand;
    Tr         = (dG + dS*Tref) / dS;
    ord        = (-dS/dV*Tr - Pref*1000) /1000;
    P_and_ky   = (T.*dS./dV./1000 .+ ord) /10; # from kbar to GPa
    dV         = Vky-Vsill;
    dS         = Sky-Ssill;
    dG         = Gky-Gsill;
    Tr         = (dG + dS*Tref) / dS;
    ord        = (-dS/dV*Tr - Pref*1000) /1000;
    P_ky_sill  = (T.*dS./dV./1000 .+ ord) /10; # from kbar to GPa
    dV         = Vand-Vsill;
    dS         = Sand-Ssill;
    dG         = Gand-Gsill;
    Tr         = (dG + dS*Tref) / dS;
    ord        = (-dS/dV*Tr - Pref*1000) /1000;
    P_and_sill = (T.*dS./dV./1000 .+ ord) /10; # from kbar to GPa
    # Physical properties
    dG_and_p   = 0.50*( dG_and[:,1:end-1] .+ dG_and[:,2:end] ) 
    dG_sill_p  = 0.50*(dG_sill[:,1:end-1] .+ dG_sill[:,2:end])
    dG_ky_p    = 0.50*(  dG_ky[:,1:end-1] .+ dG_ky[:,2:end]  )
    dG_and_pt  = 0.25*(dG_and[2:end-0,1:end-1]  .+ dG_and[1:end-1,2:end-0]  .+ dG_and[1:end-1,1:end-1]  .+ dG_and[2:end,2:end])
    dG_sill_pt = 0.25*(dG_sill[2:end-0,1:end-1] .+ dG_sill[1:end-1,2:end-0] .+ dG_sill[1:end-1,1:end-1] .+ dG_sill[2:end,2:end])
    dG_ky_pt   = 0.25*(dG_ky[2:end-0,1:end-1]   .+ dG_ky[1:end-1,2:end-0]   .+ dG_ky[1:end-1,1:end-1]   .+ dG_ky[2:end,2:end])
    rho                                                         = rho_ky;
    rho[(dG_and_p  .< dG_ky_p) .& (dG_and_p  .< dG_sill_p)]    .= rho_and[(dG_and_p  .< dG_ky_p) .& (dG_and_p  .< dG_sill_p)] 
    rho[(dG_sill_p .< dG_ky_p) .& (dG_sill_p .< dG_and_p )]    .= rho_sill[(dG_sill_p .< dG_ky_p) .& (dG_sill_p .< dG_and_p )]
    alp                                                         = alp_ky;
    alp[(dG_and_pt .< dG_ky_pt) .& (dG_and_pt .< dG_sill_pt)]  .= alp_and[(dG_and_pt .< dG_ky_pt) .& (dG_and_pt .< dG_sill_pt)]
    alp[(dG_sill_pt.< dG_ky_pt) .& (dG_sill_pt.< dG_and_pt)]   .= alp_sill[(dG_sill_pt.< dG_ky_pt) .& (dG_sill_pt.< dG_and_pt)]
    dG_and_pi                                                   = dG_and[:,2:end-1]
    dG_ky_pi                                                    = dG_ky[:,2:end-1]
    dG_sill_pi                                                  = dG_sill[:,2:end-1]
    bet                                                         = bet_ky;
    bet[(dG_and_pi .< dG_ky_pi) .& (dG_and_pi  .< dG_sill_pi)] .= bet_and[(dG_and_pi .< dG_ky_pi) .& (dG_and_pi  .< dG_sill_pi)]
    bet[(dG_sill_pi.< dG_ky_pi) .& (dG_sill_pi .< dG_and_pi )] .= bet_sill[(dG_sill_pi.< dG_ky_pi) .& (dG_sill_pi .< dG_and_pi )]
    Pc                       =  0.5.*(  P[1:end-1] .+   P[2:end])
    Tc                       =  0.5.*(  T[1:end-1] .+   T[2:end])
    Pi                       = P[2:end-1]
    # Visualise
    p1 = heatmap( T.-273.15,  P/1e9, Stab', c=:inferno, xlabel=L"T_{ } [{\mathrm{C}}]", ylabel=L"P_{ } [\mathrm{GPa}]", title=L"\mathrm{AL2SiO5}" )
    p1 = plot!(T.-273.15,  P_and_ky, linecolor=:white)
    p1 = plot!(T.-273.15, P_ky_sill, linecolor=:white)
    p1 = plot!(T.-273.15,P_and_sill, linecolor=:white, xlims=(minimum(T.-273.15),maximum(T.-273.15)), ylims=(minimum(P./1e9),maximum(P./1e9)), legend = false)
    p2 = heatmap( T.-273.15, Pc/1e9,  rho', c=:inferno, xlabel=L"T_{ } [{\mathrm{C}}]", ylabel=L"P_{ } [\mathrm{GPa}]", title=L"\rho_{ } [\mathrm{kg.m}^{-3}]" )
    p3 = heatmap(Tc.-273.15, Pc/1e9,  log10.(alp)', c=:inferno, xlabel=L"T_{ } [{\mathrm{C}}]", ylabel=L"P_{ } [\mathrm{GPa}]", title=L"\log_{10} \alpha_{ } [\mathrm{K}^{-1}]" )
    p4 = heatmap( T.-273.15, Pi/1e9,  log10.(bet)', c=:inferno, xlabel=L"T_{ } [{\mathrm{C}}]", ylabel=L"P_{ } [\mathrm{GPa}]", title=L"\log_{10} \beta_{ } [\mathrm{Pa}^{-1}]" )
    display(plot(p1,p2,p3,p4, layout = 4))
end

@time AluminoSilicates()