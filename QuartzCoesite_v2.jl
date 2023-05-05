using Plots, LaTeXStrings, MAT
gr()
#----------------------------------------------------------#
# This codes combines elements of the thermodynamics programming classes of Yury Podladchikov and from Annelore Bessat's Ph.D.
#----------------------------------------------------------#
@views function GibbsEnergy(T2,P2,a,b,c,d,Vref,Pref,Tref,a0,Kref,dSref,dHref,dGL,Tc0,Smax,Vmax)
    # P2       = P2./1e5/1e3   # From Pa to kbar
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
    # if dGL == 0
        @. int_Cp   = (2.0 *sqrt(T2).*d + T2.^2 *b/2 + T2.*a - c./T2) - (2.0 *sqrt(Tref).*d + Tref.^2 *b/2 + Tref.*a - c./Tref); 
        @. int_Cp_T = -0.5*T2.^(-2.0).*c - 2.0*1 /sqrt(T2).*d + 1.0*T2.^1.0.*b + 1.0*a.*log(T2) - (-0.5*Tref.^(-2.0).*c - 2.0*1 /sqrt(Tref).*d + 1.0*Tref.^1.0.*b + 1.0*a.*log(Tref)); 
    # else
    #     @. int_Cp   = (a*T2+(1.0/2.0)*b*T2.^2.0-c./T2+2*d*sqrt(T2)+(1.0/2.0)*Smax*((2.0/3.0)*(Tc-T2).^(3.0/2.0)-2.0*Tc.*sqrt(Tc-T2))./sqrt(Tc)) - (a*Tref+(1.0/2.0)*b*Tref.^2-c./Tref+2.0*d*sqrt(Tref)+(1.0/2.0)*Smax*((2.0/3.0)*(Tc-Tref).^(3.0/2.0)-2.0*Tc.*sqrt(Tc-Tref))./sqrt(Tc))
    #     @. int_Cp_T = (-2.0*d./sqrt(T2)-Smax*sqrt(Tc-T2)./sqrt(Tc)+b*T2-(1.0/2.0)*c./T2.^2.0+a.*log(T2)) - (-2.0*d./sqrt(Tref)-Smax*sqrt(Tc-Tref)./sqrt(Tc)+b*Tref-(1.0/2.0)*c./Tref.^2.0+a.*log(Tref))
    # end
    @. int_V    = 1/3*(KT+4*P2).*V1T.*(KT./(KT+4*P2)).^(1/4) - 1/3*(KT+4*Pref).*V1T.*(KT./(KT+4*Pref)).^(1/4); 
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
    # dG       = dG/mmol*1e3                                       # from kJ/mol to J/kg
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
function QuartzCoesite()
    # Parameters - all thermodynamic parameters are taken from hp-ds55.txt database 
    nb_ph   = 2;
    Tref    = 298;                    # K
    Pref    = 1e-3                    # kbar        

    R       = 8.3144;                 # J/K/mol
    Tmin    = 673;               # K
    Tmax    = 873;
    Pmin    = 20;                  # kbar
    Pmax    = 35;

    Tmin    = 298+1e-3+100;               # K
    Tmax    = 800+273;
    Pmin    = -1e-3;                  # kbar
    Pmax    = 50;
    dP      = 0.02;                    # kbar
    dT      = 1;
    T       = Tmin:dT:Tmax;           # K
    P       = LinRange(Pmin, Pmax, 2500)
    T       = LinRange(Tmin, Tmax, 675)
    dP, dT  = P[2]-P[1], T[2]-T[1] 
    # Fake ngrid for P-T space  ;)
    T2      = repeat(T, 1, length(P))
    P2      = repeat(P, 1, length(T))'
    # COESITE
    a       = 0.0965;
    b       = -0.000000577;
    c       = -444.8;
    d       = -0.7982; 
    Vref    = 2.064;      
    a0      = 1.80e-5;      
    Kref    = 1000;      
    dSref   = 40.8*1e-3;   
    dHref   = -905.52; 
    mmol    = 60.08e-3;
    dGL     = 0
    Tc0     = 847;
    Smax    = 0.00495 ;
    Vmax    = 0.1188;
    dG_coe  = GibbsEnergy(T2,P2,a,b,c,d,Vref,Pref,Tref,a0,Kref,dSref,dHref,dGL,Tc0,Smax,Vmax)
    dG_coe  = dG_coe/mmol/1e5
    # Physical properties
    m0      = 616000; 
    m1      = 1.0541; 
    m2      = -29.095; 
    k0      = 974000; 
    k1      = 4.3; 
    k2      = -46.004;
    K_coe, G_coe, alp_coe, bet_coe, rho_coe = PhysicalProperties(T2,P2,dG_coe,mmol,Pref,Tref,m0,m1,m2,k0,k1,k2,dP,dT)
    # Qz
    a       = 0.1107;
    b       =-0.000005189;
    c       = 0.0;
    d       =-1.1283; 
    Vref    = 2.269;    
    a0      = 0.65e-5;       
    Kref    = 750;      
    dSref   = 41.5e-3;     
    dHref   = -910.84;  
    mmol    = 60.08e-3;
    dGL     = 1
    Tc0     = 847;
    Smax    = 0.00495 ;
    Vmax    = 0.1188;
    # Gibbs energy
    dG_q    = GibbsEnergy(T2,P2,a,b,c,d,Vref,Pref,Tref,a0,Kref,dSref,dHref,dGL,Tc0,Smax,Vmax)
    dG_q    = dG_q/mmol/1e5
    m0      = 448538.0025;
    m1      = 1.645295858;
    m2      = -17.5543599;
    k0      = 371254.6662;
    k1      = 8.16916871;
    k2      = -140.80;
    # Physical properties
    K_q, G_q, alp_q, bet_q, rho_q = PhysicalProperties(T2,P2,dG_q,mmol,Pref,Tref,m0,m1,m2,k0,k1,k2,dP,dT)
    # Phase diagram
    Stab                = zeros(length(T),length(P));
    Stab[dG_coe.<dG_q] .= 1;
    # Physical properties
    dG_q_p    =  0.5.*(  dG_q[:,1:end-1] .+   dG_q[:,2:end]) 
    dG_coe_p  =  0.5.*(dG_coe[:,1:end-1] .+ dG_coe[:,2:end]) 
    dG_q_pt   = 0.25.*(  dG_q[2:end-0,1:end-1] .+   dG_q[1:end-1,2:end-0] .+   dG_q[1:end-1,1:end-1] .+   dG_q[2:end,2:end]) 
    dG_coe_pt = 0.25.*(dG_coe[2:end-0,1:end-1] .+ dG_coe[1:end-1,2:end-0] .+ dG_coe[1:end-1,1:end-1] .+ dG_coe[2:end,2:end])
    dG_q_pi   = dG_q[:,2:end-1]
    dG_coe_pi = dG_coe[:,2:end-1]
    rho                      = rho_q
    rho[dG_coe_p .< dG_q_p ].= rho_coe[dG_coe_p .< dG_q_p]
    alp                      = alp_q
    alp[dG_coe_pt.< dG_q_pt].= alp_coe[dG_coe_pt.< dG_q_pt]
    bet                      = bet_q
    bet[dG_coe_pi.< dG_q_pi].= bet_coe[dG_coe_pi.< dG_q_pi]
    G                        = G_q
    G[dG_coe   .< dG_q   ]  .= G_coe[dG_coe   .< dG_q   ]
    K                        = K_q
    K[dG_coe   .< dG_q   ]  .= K_coe[dG_coe   .< dG_q   ]
    nu                       = (3.0 .* K .- 2.0 .*G) ./ 2.0 ./(3.0 .* K .+ G);
    Pc                       =  0.5.*(  P[1:end-1] .+   P[2:end])
    Tc                       =  0.5.*(  T[1:end-1] .+   T[2:end])
    Pi                       = P[2:end-1]
    # VERIFICATION: compare numerically obtained phase diagram with analytic phase diagram boundaries
    # q - coe --- HP 98
    Vcoe     = 2.064;  
    Vq       = 2.269;  
    Scoe     = 40.80;
    Sq       = 41.50;
    Gcoe     = -850.89e3;
    Gq       = -856.46e3;
    dV       = (Vq-Vcoe);
    dS       = (Sq-Scoe);
    dG       = (Gq-Gcoe);
    Tr       = (dG + dS*Tref) / dS;
    Pref    = 101325
    Pr       = (-dG+dV*Pref) / dV 
    ord      = (-dS/dV*Tr - Pref/1e5*1000) /1000
    P_coe_q  = (T.*dS./dV./1000 .+ ord) / 10 # from bkar to GPa
    # Visualise
    p1 = Plots.heatmap( T,  P./10, Stab', c=:inferno, xlabel=L"T_{ } [{\mathrm{K}}]", ylabel=L"P_{ } [\mathrm{GPa}]", title=L"\mathrm{Quartz/Coesite}" )
    p1 = Plots.plot!(   T,  P_coe_q, linecolor=:red, linewidth=3, linestyle=:dash, xlims=(minimum(T),maximum(T)), ylims=(minimum(P/10),maximum(P./10)), legend = false )
    p2 = Plots.heatmap( T, Pc./10,  rho', c=:inferno, xlabel=L"T_{ } [{\mathrm{K}}]", ylabel=L"P_{ } [\mathrm{GPa}]", title=L"\rho_{ } [\mathrm{kg.m}^{-3}]" )
    p3 = Plots.heatmap(Tc, Pc/1e9,  log10.(abs.(alp))', c=:inferno, xlabel=L"T_{ } [{\mathrm{K}}]", ylabel=L"P_{ } [\mathrm{GPa}]", title=L"\log_{10} \alpha_{ } [\mathrm{K}^{-1}]" )
    p4 = Plots.heatmap( T, Pi/1e9,  log10.(abs.(bet))', c=:inferno, xlabel=L"T_{ } [{\mathrm{K}}]", ylabel=L"P_{ } [\mathrm{GPa}]", title=L"\log_{10} \beta_{ } [\mathrm{Pa}^{-1}]" )
    display(Plots.plot(p1, p2, layout = 4)) #, p3, p4

end

@time QuartzCoesite()