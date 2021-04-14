using MAT
using Base.Threads
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
        Cp      = Cp .+ (Smax.*T2) ./(2*sqrt.(Tc) .* sqrt.(Tc-T2)) .* (T2.<=Tc);          # dG Landau for SiO2 (second order phase transitions)
    end
    V1T      = Vref.*(1 .+     a0.*(T2.-Tref) .- 20.0 .*a0 .*(sqrt.(T2) .- sqrt.(Tref))); # J/bar - extensive
    KT       = Kref.*(1 .- 1.5e-4.*(T2.-Tref));
    V1       = V1T .* (1 .- 4.0 .*P2./(KT .+ 4.0 .*P2)).^(1/4);                           # Murnaghan EOS - Holland - Powell (1998)
    int_Cp   = zeros(size(P2))
    int_Cp_T = zeros(size(P2))
    int_V    = zeros(size(P2))
    @. int_Cp   = (2.0 *sqrt(T2).*d + T2.^2 *b/2 + T2.*a - c./T2) - (2.0 *sqrt(Tref).*d + Tref.^2 *b/2 + Tref.*a - c./Tref); 
    @. int_Cp_T = -0.5*T2.^(-2.0).*c - 2.0*1 /sqrt(T2).*d + 1.0*T2.^1.0.*b + 1.0*a.*log(T2) - (-0.5*Tref.^(-2.0).*c - 2.0*1 /sqrt(Tref).*d + 1.0*Tref.^1.0.*b + 1.0*a.*log(Tref)); 
    @. int_V    = 1/3*(KT+4*P2).*V1T.*(KT./(KT+4*P2)).^(1/4) - 1/3*(KT+4*Pref).*V1T.*(KT./(KT+4*Pref)).^(1/4); 
    dH_T     = dHref .+ int_Cp;             
    dS_T     = dSref .+ int_Cp_T;           
    dG       = dH_T  .- T2.*dS_T  .+ int_V; 
    return dG
end
#----------------------------------------------------------#
function CigarMelt()
    # Parameters - all thermodynamic parameters are taken from hp-ds55.txt database 
    nb_ph   = 4;
    Tref    = 298;                    # K
    Pref    = 1e-3;                   # kbar        
    R       = 8.3144;                 # J/K/mol
    dC      = 0.05
    dT      = 10
    dP      = 1e8
    C_sys   = 0:dC:1;
    T       = 1400:dT:2300
    P       = 1e5:dP:5e9
    # Fake ngrid for P-T space  ;)
    T2      = repeat(T, 1, length(P))
    P2      = repeat(P, 1, length(T))'
    # Forsterite
    a       = 0.2333;
    b       = 0.000001494;
    c       = -603.8;
    d       = -1.8697;
    Vref    = 4.3660;    
    a0      = 0.0000613;       
    Kref    = 1250.00;      
    dSref   = 0.09510;     
    dHref   = -2172.20;  
    mmol_fo = 140.69e-3; # kg/mol not in HP database
    dG_fo   = GibbsEnergy(T2,P2,a,b,c,d,Vref,Pref,Tref,a0,Kref,dSref,dHref,0,0,0,0)
    # Fayalite
    a       = 0.2011;
    b       = 0.000017330;
    c       = -1960.6;
    d       = -0.9009;
    Vref    = 4.6310;    
    a0      = 0.0000505;       
    Kref    = 1330.00;      
    dSref   = 0.15100;     
    dHref   = -1478.15;  
    mmol_fa = 207.8e-3; # kg/mol not in HP database
    dG_fa   = GibbsEnergy(T2,P2,a,b,c,d,Vref,Pref,Tref,a0,Kref,dSref,dHref,0,0,0,0)
    # Forsterite liquid
    a       = 0.2679;
    b       = 0.0;
    c       = 0.0;
    d       = 0.0;
    Vref    = 4.2430;    
    a0      = 0.0001450;       
    Kref    = 730.00;      
    dSref   = -0.05500;     
    dHref   = -2225.16;  
    mmol_foL= 140.69e-3; # kg/mol not in HP database
    dG_foL  = GibbsEnergy(T2,P2,a,b,c,d,Vref,Pref,Tref,a0,Kref,dSref,dHref,0,0,0,0)
    # Fayalite liquid
    a       = 0.2397;
    b       = 0.0;
    c       = 0.0;
    d       = 0.0;
    Vref    = 4.6950;    
    a0      = 0.0001690;       
    Kref    = 410.00;      
    dSref   = 0.10250;     
    dHref   = -1459.21;  
    mmol_faL= 207.8e-3; # kg/mol not in HP database
    dG_faL  = GibbsEnergy(T2,P2,a,b,c,d,Vref,Pref,Tref,a0,Kref,dSref,dHref,0,0,0,0)
    # Table of Gibbs ebergies for each phase and for each P and T condition
    Gibbs = zeros(length(T),length(P),nb_ph)
    Gibbs[:,:,1] = dG_fo*1e3
    Gibbs[:,:,2] = dG_fa*1e3
    Gibbs[:,:,3] = dG_foL*1e3
    Gibbs[:,:,4] = dG_faL*1e3
    Molar_mass   = [mmol_fo mmol_fa mmol_foL mmol_faL]
    # Volume: Pa*m^3/mol*Pa = m^3/mol
    Vol_c            = diff( Gibbs, dims=2)/dP;
    Vol              = zeros(size(Gibbs));
    Vol[:,2:end-1,:].= 0.5.*(Vol_c[:,1:end-1,:] .+ Vol_c[:,2:end,:]);
    Vol[:,  1,:]    .= Vol[:,    2,:];
    Vol[:,end,:]    .= Vol[:,end-1,:];
    # Entropy: J/mol*K
    S_c              = -diff(Gibbs, dims=1)/dT;
    S                = zeros(size(Gibbs));
    S[2:end-1,:,:]  .= 0.5.*(S_c[1:end-1,:,:] .+ S_c[2:end,:,:]);
    S[  1,:,:]      .= S[    2,:,:];
    S[end,:,:]      .= S[end-1,:,:];
    # Heat Capacity
    Cp_c             = T[2:end-1].*diff(S_c, dims=1)/dT;
    Cp               = zeros(size(Gibbs));
    Cp[2:end-1,:,:] .= Cp_c;
    Cp[  1,:,:]     .= Cp[2,:,:];
    Cp[end,:,:]     .= Cp[end-1,:,:];
    rho              = 1.0 ./Vol;
    # Enthalpy
    H                = Gibbs .+ T.*S;
    # Internal energy
    U                = Gibbs .+ T.*S .- P'.*Vol;
    TS               = T.*S;
    PV               = P'.*Vol;
    # Remove reference values
    for i= 1:4
        H[:,:,i]   .= H[:,:,i]  .- H[1,1,i];
        U[:,:,i]   .= U[:,:,i]  .- U[1,1,i];
        TS[:,:,i]  .= TS[:,:,i] .- TS[1,1,i];
        PV[:,:,i]  .= PV[:,:,i] .- PV[1,1,i];
    end
    #- PREPROCESSING ----------------------------------------------------------
    C                = zeros(length(C_sys),2)
    C[:,1]          .= C_sys;                    # X Mg
    C[:,2]          .= 1 .- C[:,1];              # X Fe
    Gibbs_id0        =  R.*sum(C.*log.(C .+ (C.==0)), dims=2);  
    C_all            = [C_sys;  C_sys];
    M_all            = zeros(size(C_all));
    V_all            = zeros(size(C_all));
    S_all            = zeros(size(C_all));
    Gibbs_all        = zeros(size(C_all));
    Cp_all           = zeros(size(C_all));
    i_sol            = 1:length(C_sys);                                       # index solid part
    i_mlt            = i_sol[end]+1:2*length(C_sys);                          # index melt part
    M_all[i_sol]    .= C * Molar_mass[1:2];    # (C*M_fo_sol + (1-C)*M_fa_sol)   ----   [kg/mol]
    M_all[i_mlt]    .= C * Molar_mass[3:4];    # (C*M_fo_lid + (1-C)*M_fa_liq)
    # For storing results
    alph_sol         = zeros(  length(T    ), length(P    ), length(C_sys) );
    alph_mlt         = zeros(  length(T    ), length(P    ), length(C_sys) );
    alph_stb         = zeros(2*length(C_sys), length(T    ), length(C_sys) );
    xMg_sol          = NaN*ones(  length(T    ), length(P    ), length(C_sys) );
    xMg_mlt          = NaN*ones(  length(T    ), length(P    ), length(C_sys) );
    Cp_sol           = zeros(  length(T    ), length(P    ), length(C_sys) );
    Cp_mlt           = zeros(  length(T    ), length(P    ), length(C_sys) );
    U_sol            = zeros(  length(T    ), length(P    ), length(C_sys) );
    U_mlt            = zeros(  length(T    ), length(P    ), length(C_sys) );
    rho_sol          = zeros(  length(T    ), length(P    ), length(C_sys) );
    rho_mlt          = zeros(  length(T    ), length(P    ), length(C_sys) );
    CFo_sol          = zeros(  length(T    ), length(P    ), length(C_sys) );
    CFa_sol          = zeros(  length(T    ), length(P    ), length(C_sys) );
    CFo_mlt          = zeros(  length(T    ), length(P    ), length(C_sys) );
    CFa_mlt          = zeros(  length(T    ), length(P    ), length(C_sys) );
    # Loop on T and P
    for iT =  1:length(T)
        for iP = 1:1#length(P)
            # Evaluate Gibbs energies for the solid and melt compositions
            V_all[i_sol]          = C *   Vol[iT,iP,1:2];
            V_all[i_mlt]          = C *   Vol[iT,iP,3:4];
            Cp_all[i_sol]         = C *    Cp[iT,iP,1:2];
            Cp_all[i_mlt]         = C *    Cp[iT,iP,3:4];            
            S_all[i_sol]          = C *     S[iT,iP,1:2] + Gibbs_id0;
            S_all[i_mlt]          = C *     S[iT,iP,3:4] + Gibbs_id0;
            Gibbs_all[i_sol]      = C * Gibbs[iT,iP,1:2] + T[iT]*Gibbs_id0;
            Gibbs_all[i_mlt]      = C * Gibbs[iT,iP,3:4] + T[iT]*Gibbs_id0;
            Aeq                   = [C_all ones(size(C_all))]';
            LB                    = zeros(1,length(C_all));
            UB                    =  ones(1,length(C_all));
            # Loop on system composition
            for iC = 1:length(C_sys)
                # modify equality constrain for system composition
                beq   = [C_sys[iC]; 1]; 
                # Use optimisation to find alpha such that:
                # min( Gibbs * alpha )
                # sum( alpha*C ) = C_sys
                # sum( alpha   ) = 1 
                local model = Model(GLPK.Optimizer)
                @variable(model, LB[i] <= x[i=1:length(LB)] <= UB[i])
                @objective(model, Min, Gibbs_all'x)
                @constraint(model, Aeq * x .== beq)
                optimize!(model)
                # Extract and store results
                alph                = value.(x)
                X_sol               = sum(alph[i_sol])
                X_mlt               = sum(alph[i_mlt])
                alph_stb[ :,iT,iC] .= alph;
                alph_sol[iT,iP,iC]  = X_sol;
                alph_mlt[iT,iP,iC]  = X_mlt; 
                # Extract fractionation of MgO in solid and melt
                local idx = findall(alph.>0.0)
                if length(idx)==2 # solid + melt are non-zero
                    xMg_sol[iT,iP,iC] = C_all[idx[1]]
                    xMg_mlt[iT,iP,iC] = C_all[idx[2]]
                end
                if length(idx)==1 && idx[1] <= length(C_sys) && alph_sol[iT,iP,iC] > 0.99 # solid only
                    xMg_sol[iT,iP,iC] = C_all[idx[1]]
                end
                if length(idx)==1 && idx[1]  > length(C_sys) && alph_mlt[iT,iP,iC] > 0.99 # melt only
                    xMg_mlt[iT,iP,iC] = C_all[idx[1]]
                end
                # println(alph[i_sol], Xsol)
                M_sol             = sum(    M_all[i_sol].*alph[i_sol] ) / X_sol            # mass_sol = sum( M_fa*X_fo_sol ) / Xsol
                M_mlt             = sum(    M_all[i_mlt].*alph[i_mlt] ) / X_mlt
                V_sol             = sum(    V_all[i_sol].*alph[i_sol] ) / X_sol / M_sol
                V_mlt             = sum(    V_all[i_mlt].*alph[i_mlt] ) / X_mlt / M_mlt
                S_sol             = sum(    S_all[i_sol].*alph[i_sol] ) / X_sol / M_sol
                S_mlt             = sum(    S_all[i_mlt].*alph[i_mlt] ) / X_mlt / M_mlt
                G_sol             = sum(Gibbs_all[i_sol].*alph[i_sol] ) / X_sol / M_sol
                G_mlt             = sum(Gibbs_all[i_mlt].*alph[i_mlt] ) / X_mlt / M_mlt
                # Heat capacity 
                Cp_sol[iT,iP,iC]  = sum(   Cp_all[i_sol].*alph[i_sol]') / X_sol / M_sol
                Cp_mlt[iT,iP,iC]  = sum(   Cp_all[i_mlt].*alph[i_mlt]') / X_mlt / M_mlt
                # Internal energy
                U_sol[iT,iP,iC]   = G_sol + T[iT]*S_sol - P[iP]*V_sol
                U_mlt[iT,iP,iC]   = G_mlt + T[iT]*S_mlt - P[iP]*V_mlt
                # Density
                rho_sol[iT,iP,iC] = 1.0 ./ V_sol
                rho_mlt[iT,iP,iC] = 1.0 ./ V_mlt
                # Molecular fractions 
                CFo_sol_mol       = sum( C_all[i_sol].*alph[i_sol]) ./ X_sol
                CFo_mlt_mol       = sum( C_all[i_mlt].*alph[i_mlt]) ./ X_mlt
                CFa_sol_mol       = 1.0 - CFo_sol_mol;
                CFa_mlt_mol       = 1.0 - CFo_mlt_mol;
                # Mass fractions [kg/kg]
                CFo_sol[iT,iP,iC] = CFo_sol_mol*Molar_mass[1] ./ ( CFo_sol_mol*Molar_mass[1] + CFa_sol_mol*Molar_mass[2] )
                CFo_mlt[iT,iP,iC] = CFo_mlt_mol*Molar_mass[3] ./ ( CFo_mlt_mol*Molar_mass[3] + CFa_mlt_mol*Molar_mass[4] )
                CFa_sol[iT,iP,iC] = CFa_sol_mol*Molar_mass[2] ./ ( CFo_sol_mol*Molar_mass[1] + CFa_sol_mol*Molar_mass[2] )
                CFa_mlt[iT,iP,iC] = CFa_mlt_mol*Molar_mass[4] ./ ( CFo_mlt_mol*Molar_mass[3] + CFa_mlt_mol*Molar_mass[4] )
            end
        end
    end
    # Visualise
    p1 = heatmap(C_sys, T, alph_sol[:,1,:], c=:inferno, xlabel=L"X_{\mathrm{Mg}}", ylabel=L"T_{ } [\mathrm{K}]", title=L"X^{\mathrm{sol}}" )
    p2 = heatmap(C_sys, T, alph_mlt[:,1,:], c=:inferno, xlabel=L"X_{\mathrm{Mg}}", ylabel=L"T_{ } [\mathrm{K}]", title=L"X^{\mathrm{mlt}}" )
    p3 = heatmap(C_sys, T,  xMg_sol[:,1,:], c=:inferno, xlabel=L"X_{\mathrm{Mg}}", ylabel=L"T_{ } [\mathrm{K}]", title=L"X^{\mathrm{sol}}_{\mathrm{Mg}}" )
    p4 = heatmap(C_sys, T,  xMg_mlt[:,1,:], c=:inferno, xlabel=L"X_{\mathrm{Mg}}", ylabel=L"T_{ } [\mathrm{K}]", title=L"X^{\mathrm{mlt}}_{\mathrm{Mg}}" )
    display(plot(p1,p2,p3,p4, layout = 4))
end
#----------------------------------------------------------#
@time CigarMelt()

