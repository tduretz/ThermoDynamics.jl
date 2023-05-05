using Printf
import Statistics: mean

avP(a) = 0.5*(a[:,1:end-1].+a[:,2:end])

function main_qcoe_paper()

    TDdata = :HP98

    # Physics
    Tref    = 298+1e-3+100;           # K
    Pref    = -1e-3;                  # kbar
    dP      = 0.020;                  # kbar
    dT      = 1*1;
    T       = Tref:dT:800+273;        # K
    P       = Pref:dP:50;  
    nT = length(T)
    nP = length(P)
    # Fake ngrid for P-T space  ;)
    T2      = repeat(T, 1, length(P))
    P2      = repeat(P, 1, length(T))'
    anal_int        = 1;
    SAME_AS_PERPLEX = 1;
    land            = 1;
    ## COESITE
    if TDdata==:HP98
        # Cp - coesite
        a       = 0.0965;    # HP98 ORIGINAL ✓
        b       = -0.577e-6; # HP98 ORIGINAL ✓
        c       = -444.8;    # HP98 ORIGINAL ✓
        d       = -0.7982;   # HP98 ORIGINAL ✓
        # V - coesite
        Vref    = 2.064;     # HP98 ORIGINAL ✓
        a0      = 1.80e-5;   # HP98 ORIGINAL ✓    
        Kref    = 1000;      # HP98 ORIGINAL ✓ 
        dSref   = 40.8*1e-3; # HP98 ORIGINAL ✓ 
        dHref   = -905.52;   # HP98 ORIGINAL ✓
        mmol    = 60.08e-3;
    end

    if TDdata==:HP11
        # Cp - coesite
        a       = 0.1078;     # HP11 ORIGINAL ✓
        b       = -0.3279e-5; # HP11 ORIGINAL ✓
        c       = -190.3;     # HP11 ORIGINAL ✓
        d       = -1.0416;    # HP11 ORIGINAL ✓
        # V - coesite
        Vref    = 2.064      # HP11 ORIGINAL ✓
        dSref   = 39.6*1e-3  # HP11 ORIGINAL ✓ 
        dHref   = -907.02    # HP11 ORIGINAL ✓
        mmol    = 60.08e-3
        a0      = 1.23e-5    # HP11 ORIGINAL ✓    
        Kref    = 979.       # HP11 ORIGINAL ✓ 
        Kref1   = 4.19       # HP11 ORIGINAL ✓
        Kref2   = -.0043     # HP11 ORIGINAL ✓
    end

    # My try based on Yury's thermodynamics class
    Cp, V1T, KT, V1 = zero(T2), zero(T2), zero(T2),zero(T2)
    @. Cp      = a + b*T2 + c*T2.^(-2) + d*T2.^(-1/2);   
    if TDdata==:HP98                   # J/K   - extensive
        @. V1T     = Vref*(1 + a0*(T2-Tref) - 20*a0 .*(sqrt(T2) - sqrt(Tref))); # J/bar - extensive
        @. KT      = Kref*(1-1.5e-4*(T2-Tref));
        @. V1      = V1T .* (1 - 4 .*P2./(KT+4 .*P2)).^(1/4); # Murnaghan EOS - Holland - Powell (1998)
    elseif TDdata==:HP11
        a1 = (1 + Kref1) / (1 + Kref1 + Kref*Kref2)
        b1 = Kref1/Kref - Kref2/(1 + Kref1)
        c1 = (1 + Kref1 + Kref*Kref2) / (Kref1^2 + Kref1 - Kref*Kref2) 
        @. V1T  = Vref*(1 + a0*(T2-Tref) - 20*a0 .*(sqrt(T2) - sqrt(Tref))); # J/bar - extensive
        @. V1   = V1T *(1 - a1*(1 - (1 + b1*P2)^(-c1)) )
    end
    @printf("mean coe density = %2.2f\n", mmol/mean(V1[:])*1e5)
    @printf("mean coe heatcap = %2.2f\n", mean(Cp[:])/mmol*1000)

    int_Cp, int_Cp_T, int_V = zero(T2), zero(T2), zero(T2)
    dH_T, dS_T, dG_coe = zero(T2), zero(T2), zero(T2)
    if anal_int == 1 @. int_Cp   = (2.0*sqrt(T2).*d + T2.^2 .*b/2 + T2.*a - c./T2) - (2.0*sqrt(Tref).*d + Tref.^2 .*b/2 + Tref.*a - c./Tref); end
    if anal_int == 1 @. int_Cp_T = -0.5*T2.^(-2.0).*c - 2.0*1 ./sqrt(T2).*d + 1.0*T2.^1.0.*b + 1.0*a.*log(T2) - (-0.5*Tref.^(-2.0).*c - 2.0*1 ./sqrt(Tref).*d + 1.0*Tref.^1.0.*b + 1.0*a.*log(Tref)); end 
    if (anal_int == 1 && TDdata==:HP98) @. int_V    = 1/3 .*(KT+4*P2).*V1T.*(KT./(KT+4*P2)).^(1/4) - 1/3 .*(KT+4*Pref).*V1T.*(KT./(KT+4*Pref)).^(1/4); end
    if anal_int == 0 int_Cp   .= cumsum(Cp    , dims=1).*dT; end 
    if anal_int == 0 int_Cp_T .= cumsum(Cp./T2, dims=1).*dT; end
    if (anal_int == 0 || TDdata==:HP11) int_V    .= cumsum(V1    , dims=2).*dP; end
    @. dH_T    = dHref + int_Cp;            # dH_T = dH_ref + int dCp   dT
    @. dS_T    = dSref + int_Cp_T;          # dS_T = dS_ref + int dCp/T dT
    @. dG_coe  = dH_T  - T2.*dS_T  + int_V; # dG   = dH_T - T*dS_T + int dV dP 
    # Physical properties
    dG_coe   = dG_coe/mmol*1000;
    dP_Pa    = dP*1000*1e5; # from kbar to Pa
    rho_coe = zeros(nT, nP-1)
    rho_coe  .=  1 ./(          diff(dG_coe    , dims=2)./dP_Pa);
    # rho_tcax = (rho_coe[1:end-1,:]+rho_coe[2:end,:])/2;
    # rho_tcay = (rho_coe[:,1:end-1]+rho_coe[:,2:end])/2;
    # alp_coe  = -1./rho_tcax.*diff(rho_coe,1,1)/dT;               # thermodynamically consistant alpha at cst P
    # bet_coe  =  1./rho_tcay.*diff(rho_coe,1,2)/dP_Pa;            # thermodynamically consistant beta at cst T

    # m0= 616000; 
    # m1= 1.0541; 
    # m2= -29.095; 
    # k0= 974000; 
    # k1= 4.3; 
    # k2= -46.004;
    # G_coe = m0 + m1*(P2-Pref)*1e3 + m2*(T2-Tref);
    # K_coe = k0 + k1*(P2-Pref)*1e3 + k2*(T2-Tref);


    # ## Quartz
    # # q   1  1    1.00 10    2.00  0          -910.83   0.04150   2.2688
    # #          0.1107   -0.000005189       0.0   -1.1283         0.0000065    10.000   750.00     4.00   -0.11250    1    847   0.00495    0.1188

    if TDdata==:HP98
        # Cp - quartz
        a       = 0.1107;    # HP98 ORIGINAL ✓
        b       =-0.5189e-5; # HP98 ORIGINAL ✓
        c       = 0.0;       # HP98 ORIGINAL ✓
        d       =-1.1283;    # HP98 ORIGINAL ✓ 
        # V - quartz
        Vref    = 2.269;     # HP98 ORIGINAL ✓  
        a0      = 0.65e-5;   # HP98 ORIGINAL ✓      
        Kref    = 750;       # HP98 ORIGINAL ✓    
        dSref   = 41.5e-3;   # HP98 ORIGINAL ✓ 
        dHref   = -910.88;   # HP98 ORIGINAL ✓ 
        mmol    = 60.08e-3;
        # dG Landau
        Tc0     = 847;       # HP98 ORIGINAL ✓
        Smax    = 0.00495;   # HP98 ORIGINAL ✓
        Vmax    = 0.1188;    # HP98 ORIGINAL ✓
    end

    if TDdata==:HP11
        # Cp - quartz
        a       = 0.0929    # HP11 ORIGINAL ✓
        b       =-0.0642e-5 # HP11 ORIGINAL ✓
        c       = -714.9    # HP11 ORIGINAL ✓
        d       =-0.7161    # HP11 ORIGINAL ✓ 
        # V - quartz
        Vref    = 2.269      # HP11 ORIGINAL ✓     
        dSref   = 41.43e-3   # HP11 ORIGINAL ✓ 
        dHref   = -910.7     # HP11 ORIGINAL ✓ 
        mmol    = 60.08e-3
        a0      = 0.         # HP11 ORIGINAL ✓      
        Kref    = 730.       # HP11 ORIGINAL ✓ 
        Kref1   = 6.         # HP11 ORIGINAL ✓ 
        Kref2   = -0.0082    # HP11 ORIGINAL ✓
        # dG Landau
        Tc0     = 847        # HP11 ORIGINAL ✓
        Smax    = 0.00495    # HP11 ORIGINAL ✓
        Vmax    = 0.1188     # HP11 ORIGINAL ✓
    end

    Tc = zero(T2)
    @. Tc      = Tc0 + Vmax./Smax.*(P2-Pref);
    # My try based on Yury's thermodynamics class
    @. Cp      = a  + b*T2 + c*T2.^(-2) + d*T2.^(-1/2);                         # kJ/K   - extensive
    # @. Cp      = Cp + (Smax*T2) ./(2*sqrt(Tc) .* sqrt(Tc-T2)) .* (T2<=Tc);
    for i in eachindex(Cp)
        if T2[i]<=Tc[i]
            Cp[i]      +=  (Smax*T2[i]) ./(2*sqrt(Tc[i]) .* sqrt(Tc[i]-T2[i]))
            if isinf(Cp[i])
                @show sqrt(Tc[i]-T2[i])
                @show sqrt(Tc[i])
            end 
        end
    end
    if TDdata==:HP98
        @. V1T     = Vref*(1 +     a0*(T2-Tref) - 20*a0 .*(sqrt(T2) - sqrt(Tref))); # J/bar - extensive
        @. KT      = Kref*(1 - 1.5e-4*(T2-Tref));
        @. V1      = V1T .* (1 - 4 .*P2./(KT+4 .*P2)).^(1/4); # Murnaghan EOS - Holland - Powell (1998)
    elseif TDdata==:HP11
        a1 = (1 + Kref1) / (1 + Kref1 + Kref*Kref2)
        b1 = Kref1/Kref - Kref2/(1 + Kref1)
        c1 = (1 + Kref1 + Kref*Kref2) / (Kref1^2 + Kref1 - Kref*Kref2) 
        @. V1T  = Vref*(1 + a0*(T2-Tref) - 20*a0 .*(sqrt(T2) - sqrt(Tref))); # J/bar - extensive
        @. V1   = V1T *(1 - a1*(1 - (1 + b1*P2)^(-c1)) )
        @show minimum(V1)
        @show maximum(V1)
    end
    @printf("mean q density = %2.2f\n", mmol/mean(V1[:])*1e5)
    @printf("mean q heatcap = %2.2f\n", mean(Cp[:])/mmol*1000)

    if SAME_AS_PERPLEX ==1
        # This is the same as perplex
        if anal_int == 1 @. int_Cp   = (2.0*sqrt(T2).*d + T2.^2 .*b/2 + T2.*a - c./T2) - (2.0*sqrt(Tref).*d + Tref.^2 .*b/2 + Tref.*a - c./Tref); end
        if anal_int == 1 @. int_Cp_T = -0.5*T2.^(-2.0).*c - 2.0*1 ./sqrt(T2).*d + 1.0*T2.^1.0.*b + 1.0*a.*log(T2) - (-0.5*Tref.^(-2.0).*c - 2.0*1 ./sqrt(Tref).*d + 1.0*Tref.^1.0.*b + 1.0*a.*log(Tref)); end
    else
        # This is different than perplex
        if anal_int == 1 @. int_Cp   = (a*T2+(1/2)*b*T2.^2-c./T2+2*d*sqrt(T2)+(1/2)*(T2<Tc).*Smax.*((2/3).*(Tc-T2).^(3/2)-2*Tc.*sqrt(Tc-T2))./sqrt(Tc)) - (a*Tref+(1/2)*b*Tref.^2-c./Tref+2*d*sqrt(Tref)+(1/2)*(T2<Tc).*Smax.*((2/3).*(Tc-Tref).^(3/2)-2*Tc.*sqrt(Tc-Tref))./sqrt(Tc)); end
        if anal_int == 1 @. int_Cp_T = (-2*d./sqrt(T2)-(T2<Tc).*Smax.*sqrt(Tc-T2)./sqrt(Tc)+b*T2-(1/2)*c./T2.^2+a.*log(T2)) - (-2*d./sqrt(Tref)-(T2<Tc).*Smax.*sqrt(Tc-Tref)./sqrt(Tc)+b*Tref-(1/2)*c./Tref.^2+a.*log(Tref)); end
    end
    if (anal_int == 1 && TDdata==:HP98)  @. int_V    = 1/3 .*(KT+4*P2).*V1T.*(KT./(KT+4*P2)).^(1/4) - 1/3 .*(KT+4*Pref).*V1T.*(KT./(KT+4*Pref)).^(1/4); end
    if anal_int == 0 int_Cp   .= cumsum(Cp    , dims=1).*dT; end 
    if anal_int == 0 int_Cp_T .= cumsum(Cp./T2, dims=1).*dT; end
    if (anal_int == 0 || TDdata==:HP11) int_V    .= cumsum(V1    , dims=2).*dP; end
    dG_q = zero(T2)
    @. dH_T    = dHref + int_Cp;            # dH_T = dH_ref + int dCp   dT
    @. dS_T    = dSref + int_Cp_T;          # dS_T = dS_ref + int dCp/T dT
    @. dG_q    = dH_T  - T2.*dS_T  + int_V; # dG   = dH_T - T*dS_T + int dV dP 
    # ---> Hans - dG Landau
    Q_298, Q, h_298, s_298, v_t, intvdP, G_land, G_exc  = zero(T2), zero(T2), zero(T2), zero(T2),zero(T2), zero(T2),zero(T2) ,zero(T2)
    for i in eachindex(Cp)
        if T2[i]<=Tc[i]
            Q_298[i]   = (1-Tref./Tc0).^(1/4);
            Q[i]      = (1-  T2[i]./Tc[i] ).^(1/4);
        end
    end
    @. h_298   = Smax.*Tc0.*(Q_298.^2 - (1/3)*Q_298.^6);
    @. s_298   = Smax.*Q_298.^2;
    @. v_t     = Vmax.*Q_298.^2 .*(1 + a0.*(T2-Tref) - 20*a0 .* (sqrt(T2)-  sqrt(Tref)));
    if TDdata==:HP98
        @. intvdP  = (1/3)*v_t.*KT.*((1 + (4*(P2-Pref            ))./KT).^(3/4) - 1);               
    elseif TDdata==:HP11
        intvdP    = cumsum(v_t    , dims=2).*dP
    end
    @. G_land  = Smax .* ( (T2 - Tc) .*Q.^2 + (1/3).*Tc.*Q.^6);
    @. G_exc   = h_298 - T2.*s_298 + intvdP + G_land;

    dG_q_a  = dG_q/mmol*1000;
    dG_q_b  = (dG_q + land*G_exc)/mmol*1000;
    rho_q_a, rho_q_b = zeros(nT, nP-1), zeros(nT, nP-1)
    rho_q_a .=  1 ./(          diff(dG_q_a    , dims=2)./dP_Pa);
    rho_q_b .=  1 ./(          diff(dG_q_b    , dims=2)./dP_Pa);

    # Stab2                = zeros(length(T),length(P));
    # Stab2( dG_q_b  < dG_coe &  (T2<=Tc)) = 1;  # identify a-b qz transition
    # Stab2( dG_q_b  < dG_coe &  (T2 >Tc)) = 2;  # identify a-b qz transition

    rho2                     = rho_coe
    rho2[avP(dG_q_b)  .< avP(dG_coe) .&&  (avP(T2).<=avP(Tc)) ] .= rho_q_b[avP(dG_q_b)  .< avP(dG_coe) .&&  (avP(T2).<=avP(Tc)) ]
    rho2[avP(dG_q_b)  .< avP(dG_coe) .&&  (avP(T2).> avP(Tc)) ] .= rho_q_a[avP(dG_q_b)  .< avP(dG_coe) .&&  (avP(T2).> avP(Tc)) ]

    @show minimum(dG_coe)
    @show maximum(dG_coe)
    @show minimum(dG_q_a)
    @show maximum(dG_q_a)
    @show minimum(dG_q_b)
    @show maximum(dG_q_b)
    write("SiO2_julia_revised_v4_"*string(TDdata)*".dat", rho2)


    # dG_q    = dG_q + land*G_exc .* (T2<=Tc);
    # # Physical properties
    # dG_q     = dG_q/mmol*1000; # from kJ/mol to J/kg
    # dP_Pa    = dP*1000*1e5;            # from kbar to Pa
    # rho_q    =  1./(          diff(dG_q    ,1,2)/dP_Pa);
    # rho_tcax = (rho_q(1:end-1,:)+rho_q(2:end,:))/2;
    # rho_tcay = (rho_q(:,1:end-1)+rho_q(:,2:end))/2;
    # alp_q    = -1./rho_tcax.*diff(rho_q,1,1)/dT;               # thermodynamically consistant alpha at cst P
    # bet_q    =  1./rho_tcay.*diff(rho_q,1,2)/dP_Pa;            # thermodynamically consistant beta at cst T

    # m0= 448538.0025;
    # m1= 1.645295858;
    # m2= -17.5543599;
    # k0= 371254.6662;
    # k1= 8.16916871;
    # k2= -140.80;

    # G_q = m0 + m1*(P2-Pref)*1e3 + m2*(T2-Tref);
    # K_q = k0 + k1*(P2-Pref)*1e3 + k2*(T2-Tref);

    # ## Phase diagram
    # Stab                 = zeros(length(T),length(P));
    # Stab(dG_coe  < dG_q) = 1;
    # ## Thermodynamic properties
    # dG_q_p    = (dG_q(:,1:end-1)+dG_q(:,2:end))/2; 
    # dG_coe_p  = (dG_coe(:,1:end-1)+dG_coe(:,2:end))/2; 
    # dG_q_pt   = (dG_q(2:end-0,1:end-1)+dG_q(1:end-1,2:end-0)+dG_q(1:end-1,1:end-1)+dG_q(2:end,2:end))/4; 
    # dG_coe_pt = (dG_coe(2:end-0,1:end-1)+dG_coe(1:end-1,2:end-0)+dG_coe(1:end-1,1:end-1)+dG_coe(2:end,2:end))/4;
    # dG_q_pi   = dG_q(:,2:end-1);
    # dG_coe_pi = dG_coe(:,2:end-1);
    # rho                      = rho_q;
    # rho(dG_coe_p  < dG_q_p ) = rho_coe(dG_coe_p  < dG_q_p);
    # alp                      = alp_q;
    # alp(dG_coe_pt < dG_q_pt) = alp_coe(dG_coe_pt < dG_q_pt);
    # bet                      = bet_q;
    # bet(dG_coe_pi < dG_q_pi) = bet_coe(dG_coe_pi < dG_q_pi);
    # G                        = G_q;
    # G(dG_coe_pi < dG_q_pi)   = G_coe(dG_coe_pi < dG_q_pi);
    # K                        = K_q;
    # K(dG_coe_pi < dG_q_pi)   = K_coe(dG_coe_pi < dG_q_pi);
    # nu                       = (3.*K-2.*G)./2./(3.*K+G);

    # # ## VERIFICATION: compare numerically obtained phase diagram with analytic phase diagram boundaries
    # # q - coe --- HP 98
    # Vcoe     = 2.064;  
    # Vq       = 2.269;  
    # Scoe     = 40.80;
    # Sq       = 41.50;
    # Gcoe     = -850.89e3;
    # Gq       = -856.46e3;
    # dV       = (Vq-Vcoe);
    # dS       = (Sq-Scoe);
    # dG       = (Gq-Gcoe);
    # Tr       = (dG + dS*Tref) / dS;
    # Pr       = (-dG+dV*Pref) / dV ;
    # ord      = (-dS/dV*Tr - Pref*1000) /1000;
    # P_coe_q = T*dS/dV/1000 + ord-0*1.2;

end

main_qcoe_paper()