function main()
clear
% Data to disk
write2disk = 1;
% Physics
Tref    = 298+1e-3+100;               % K
Pref    = -1e-3;                  % kbar
dP      = 0.020;                    % kbar
dT      = 1*1;
T       = Tref:dT:800+273;           % K
P       = Pref:dP:50;            % kbars
[T2,P2] = ndgrid(T,P);            % make it 2D

%%%%%%%%%%%%%%%%% Compute dG reac as function of P and T
[G_reac, rho, stab, rho_ab] = Qcoe_thermo( T, P, T2, P2, dP, dT, Tref, Pref );

% figure(2), clf
% % subplot(211)
% % imagesc(T-273, P/10, rho'), axis xy, colorbar
% % caxis([2500 3200])
% % subplot(212)
% imagesc(T-273, P/10, rho_ab'), axis xy
% xlim([400 1300])
% ylim([0 5.5])
% colorbar, caxis([2500 3200])

%%%%%%%%%%%%%%%%% Pressure and temperature dependent kinetics

%% Experiments
My     = 1e6*365*24*3600;
R      = 8.314;
t0     = 0.0;
dt     = 1e12;
nt     = 1000;

% Mosenfelder and Bohlen 1997
Q_MB   = 243e3;   % Activation energy [J/mol] - Cahn equation
k0_MB  = 0.185;  % constant
% Q_MB  = 269e3;   % Activation energy [J/mol] - Cahn equation
% k0_MB = 26.5;  % constant
d      = 1e-4;    % Grain size [m]
S      =  3.35/d; % Grain boundary area [m-1]
tr     = 0.99;

gr_MB_TP = k0_MB.*T2  .*(exp(-Q_MB./R./T2).*(1-exp(-G_reac./R./T2))); % growth rate [m/s]
t_MB     = log(1-tr) * (-2*S.*gr_MB_TP).^-1;

%%%%%%%%%%%%%%%%% Growth
Tfix   = 400+273;
Pfix   = 2.7e9;
Gfix   = interpn(T2,P2,G_reac, Tfix, Pfix/1e9*10);
gr_MB  = k0_MB.*Tfix.*(exp(-Q_MB./R./Tfix).*(1-exp(-Gfix./R./Tfix))); % growth rate [m/s]

% Make time loop
X_MD  = zeros(nt,1);
X_MB  = zeros(nt,1);
ttab  = zeros(nt,1);
t     = 0;
X0    = 0;
Xeq   = 1.0;
tau   = log(1-0.666) * (-2*S.*gr_MB).^-1; % Kinetic time of MDoodz
for i=1:nt
    if i>1, X0 = X_MD(i-1); end
    t        = t + dt;
    X_MD(i)  = 1/(tau+dt)*(X0*tau + Xeq*dt);
    X_MB(i)  = 1-exp(-2*S.*gr_MB*t);
    ttab(i)  = t;  
end 

%%%%%%%%%%%%%%%%% Visu

figure(1), clf

subplot(221)
imagesc(T-273.15, P/10, rho'), axis xy, title('rho'), 
c=colorbar
xlim([300 500])
ylim([0 4])
caxis ([2500  3000]);
xlabel('Temperature [$^{o}$C]','interpreter','latex','Fontsize',14)
ylabel('Pressure [GPa]','interpreter','latex','Fontsize',14)
set(gca,'TickLabelInterpreter','latex','Fontsize',12)
title('Density','interpreter','latex','Fontsize',14)
set(c,'TickLabelInterpreter','latex','Fontsize',12)
box on

subplot(222), hold on
imagesc(T-273.15, P/10, G_reac'), axis xy, title('Greac'), 
c=colorbar
[C,h] = contour(T2-273.15, P2/10, G_reac, [0.0 0.0], 'w');
plot(Tfix-273.15, Pfix/1e9, 'wx')
hold off
xlim([300 500])
ylim([0 4])
xlabel('Temperature [$^{o}$C]','interpreter','latex','Fontsize',14)
ylabel('Pressure [GPa]','interpreter','latex','Fontsize',14)
set(gca,'TickLabelInterpreter','latex','Fontsize',12)
title('$\Delta G_{r}$','interpreter','latex','Fontsize',14)
set(c,'TickLabelInterpreter','latex','Fontsize',12)
t=text(Tfix-263.15, Pfix/1e9-0.2, ['T = ', num2str(Tfix-273, '%2.0f'), ' $^{o}$C'],'interpreter', 'latex','Fontsize', 12)
tt=text(Tfix-263.15, Pfix/1e9-0.4, ['P = ', num2str(Pfix/1e9, '%2.2f'), ' GPa'],'interpreter', 'latex','Fontsize',12)
set(t,'color','w')
set(tt,'color','w')

box on

subplot(223)
[C,h] = contour(T2-273.15, P2/10, t_MB/My, [0.1 10 1 100], 'k');
xlim([300 500])
ylim([0 4])
xlabel('Temperature [$^{o}$C]','interpreter','latex','Fontsize',14)
ylabel('Pressure [GPa]','interpreter','latex','Fontsize',14)
set(gca,'TickLabelInterpreter','latex','Fontsize',12)
title('Time required for 100 \% conversion [Ma]','interpreter','latex','Fontsize',14)
text(422, 1, '0.1','interpreter', 'latex','Fontsize', 12,'Rotation',80)
text(386, 1, '1','interpreter', 'latex','Fontsize', 12,'Rotation',80)
text(355, 1, '10','interpreter', 'latex','Fontsize', 12,'Rotation',80)
text(325, 1, '100','interpreter', 'latex','Fontsize', 12,'Rotation',80)
box on

subplot(224)
plot(ttab/My, X_MB*100, '-k', ttab(1:10:end)/My, X_MD(1:10:end)*100, 'ob')
l=legend('X$_{MB}$', 'X$_{MD}$','interpreter','latex','Fontsize',12,'Location','east', 'box', 'off')
xlabel('t'), ylabel('X')
xlabel('Time [Ma]','interpreter','latex','Fontsize',14)
ylabel('Amount of quartz [\%]','interpreter','latex','Fontsize',14)
set(gca,'TickLabelInterpreter','latex','Fontsize',12)
title('Volume fraction of product transformed','interpreter','latex','Fontsize',14)
text(24, 15, ['T = ', num2str(Tfix-273, '%2.0f'), ' $^{o}$C'],'interpreter', 'latex','Fontsize', 12)
text(24, 10, ['P = ', num2str(Pfix/1e9, '%2.2f'), ' GPa'],'interpreter', 'latex','Fontsize',12)
box on

if write2disk == 1
    fileID = fopen('SiO2_old_lappy.dat','w');
    fwrite(fileID, rho_ab, 'double');
    fclose(fileID);
    size(rho)
    
    fileID = fopen('dG_QuartzCoesite_old_lappy.dat','w');
    fwrite(fileID, G_reac, 'double');
    fclose(fileID);
end
size(rho_ab)



end

function [DeltaG, rho, Stab2, rho2] = Qcoe_thermo( T, P, T2, P2, dP, dT, Tref, Pref )

anal_int        = 1;
SAME_AS_PERPLEX = 1;
land            = 1;
%% COESITE
% Cp - coesite
a       = 0.0965;
b       = -0.000000577;
c       = -444.8;
d       = -0.7982; 
% V - coesite
Vref    = 2.064;      
a0      = 1.80e-5;      
Kref    = 1000;      
dSref   = 40.8*1e-3; % HP
%dSref   = 38.5*1e-3; % Robie 1978
dHref   = -905.52; % HP
%dHref   = -907.8; % Robie 1978
mmol    = 60.08e-3;
% My try based on Yury's thermodynamics class
Cp      = a + b*T2 + c*T2.^(-2) + d*T2.^(-1/2);                      % J/K   - extensive
V1T     = Vref*(1 + a0*(T2-Tref) - 20*a0 .*(sqrt(T2) - sqrt(Tref))); % J/bar - extensive
KT      = Kref*(1-1.5e-4*(T2-Tref));
V1      = V1T .* (1 - 4.*P2./(KT+4.*P2)).^(1/4); % Murnaghan EOS - Holland - Powell (1998)

fprintf('mean coe density = %2.2f\n', mmol/mean(V1(:))*1e5)
fprintf('mean coe heatcap = %2.2f\n', mean(Cp(:))/mmol*1000)

if anal_int == 1, int_Cp   = (2.0*sqrt(T2).*d + T2.^2.*b/2 + T2.*a - c./T2) - (2.0*sqrt(Tref).*d + Tref.^2.*b/2 + Tref.*a - c./Tref); end
if anal_int == 1, int_Cp_T = -0.5*T2.^(-2.0).*c - 2.0*1./sqrt(T2).*d + 1.0*T2.^1.0.*b + 1.0*a.*log(T2) - (-0.5*Tref.^(-2.0).*c - 2.0*1./sqrt(Tref).*d + 1.0*Tref.^1.0.*b + 1.0*a.*log(Tref)); end 
if anal_int == 1, int_V    = 1/3.*(KT+4*P2).*V1T.*(KT./(KT+4*P2)).^(1/4) - 1/3.*(KT+4*Pref).*V1T.*(KT./(KT+4*Pref)).^(1/4); end
if anal_int == 0, int_Cp   = cumsum(Cp    ,1)*dT; end 
if anal_int == 0, int_Cp_T = cumsum(Cp./T2,1)*dT; end
if anal_int == 0, int_V    = cumsum(V1    ,2)*dP; end
dH_T    = dHref + int_Cp;            % dH_T = dH_ref + int dCp   dT
dS_T    = dSref + int_Cp_T;          % dS_T = dS_ref + int dCp/T dT
dG_coe  = dH_T  - T2.*dS_T  + int_V; % dG   = dH_T - T*dS_T + int dV dP 
% Physical properties
dG_coe   = dG_coe/mmol*1000;
dP_Pa    = dP*1000*1e5; % from kbar to Pa
rho_coe  =  1./(          diff(dG_coe    ,1,2)/dP_Pa);
rho_tcax = (rho_coe(1:end-1,:)+rho_coe(2:end,:))/2;
rho_tcay = (rho_coe(:,1:end-1)+rho_coe(:,2:end))/2;
alp_coe  = -1./rho_tcax.*diff(rho_coe,1,1)/dT;               % thermodynamically consistant alpha at cst P
bet_coe  =  1./rho_tcay.*diff(rho_coe,1,2)/dP_Pa;            % thermodynamically consistant beta at cst T

m0= 616000; 
m1= 1.0541; 
m2= -29.095; 
k0= 974000; 
k1= 4.3; 
k2= -46.004;
G_coe = m0 + m1*(P2-Pref)*1e3 + m2*(T2-Tref);
K_coe = k0 + k1*(P2-Pref)*1e3 + k2*(T2-Tref);


%% Quartz
% q   1  1    1.00 10    2.00  0          -910.83   0.04150   2.2688
%          0.1107   -0.000005189       0.0   -1.1283         0.0000065    10.000   750.00     4.00   -0.11250    1    847   0.00495    0.1188

% Cp - quartz
a       = 0.1107;
b       =-0.000005189;
c       = 0.0;
d       =-1.1283; 
% V - quartz
Vref    = 2.269;    
a0      = 0.65e-5;       
Kref    = 750;      
dSref   = 41.5e-3;   % kJ/mol  
dHref   = -910.88; % HP
%dHref   = -910.7; % Robie 1978 
mmol    = 60.08e-3;
% dG Landau
Tc0     = 847;
Smax    = 0.00495 ;
Vmax    = 0.1188;
Tc      = Tc0 + Vmax./Smax.*(P2-Pref);
% My try based on Yury's thermodynamics class
Cp      = a  + b*T2 + c*T2.^(-2) + d*T2.^(-1/2);                         % kJ/K   - extensive
Cp      = Cp + (Smax*T2) ./(2*sqrt(Tc) .* sqrt(Tc-T2)) .* (T2<=Tc);
V1T     = Vref*(1 +     a0*(T2-Tref) - 20*a0 .*(sqrt(T2) - sqrt(Tref))); % J/bar - extensive
KT      = Kref*(1 - 1.5e-4*(T2-Tref));
V1      = V1T .* (1 - 4.*P2./(KT+4.*P2)).^(1/4); % Murnaghan EOS - Holland - Powell (1998)
fprintf('mean q density = %2.2f\n', mmol/mean(V1(:))*1e5)
fprintf('mean q heatcap = %2.2f\n', mean(Cp(:))/mmol*1000)

if SAME_AS_PERPLEX ==1
    % This is the same as perplex
    if anal_int == 1, int_Cp   = (2.0*sqrt(T2).*d + T2.^2.*b/2 + T2.*a - c./T2) - (2.0*sqrt(Tref).*d + Tref.^2.*b/2 + Tref.*a - c./Tref); end
    if anal_int == 1, int_Cp_T = -0.5*T2.^(-2.0).*c - 2.0*1./sqrt(T2).*d + 1.0*T2.^1.0.*b + 1.0*a.*log(T2) - (-0.5*Tref.^(-2.0).*c - 2.0*1./sqrt(Tref).*d + 1.0*Tref.^1.0.*b + 1.0*a.*log(Tref)); end
else
    % This is different than perplex
    if anal_int == 1, int_Cp   = (a*T2+(1/2)*b*T2.^2-c./T2+2*d*sqrt(T2)+(1/2)*(T2<Tc).*Smax.*((2/3).*(Tc-T2).^(3/2)-2*Tc.*sqrt(Tc-T2))./sqrt(Tc)) - (a*Tref+(1/2)*b*Tref.^2-c./Tref+2*d*sqrt(Tref)+(1/2)*(T2<Tc).*Smax.*((2/3).*(Tc-Tref).^(3/2)-2*Tc.*sqrt(Tc-Tref))./sqrt(Tc)); end
    if anal_int == 1, int_Cp_T = (-2*d./sqrt(T2)-(T2<Tc).*Smax.*sqrt(Tc-T2)./sqrt(Tc)+b*T2-(1/2)*c./T2.^2+a.*log(T2)) - (-2*d./sqrt(Tref)-(T2<Tc).*Smax.*sqrt(Tc-Tref)./sqrt(Tc)+b*Tref-(1/2)*c./Tref.^2+a.*log(Tref)); end
end
if anal_int == 1, int_V    = 1/3.*(KT+4*P2).*V1T.*(KT./(KT+4*P2)).^(1/4) - 1/3.*(KT+4*Pref).*V1T.*(KT./(KT+4*Pref)).^(1/4); end
if anal_int == 0, int_Cp   = cumsum(Cp    ,1)*dT; end 
if anal_int == 0, int_Cp_T = cumsum(Cp./T2,1)*dT; end
if anal_int == 0, int_V    = cumsum(V1    ,2)*dP; end
dH_T    = dHref + int_Cp;            % dH_T = dH_ref + int dCp   dT
dS_T    = dSref + int_Cp_T;          % dS_T = dS_ref + int dCp/T dT
dG_q    = dH_T  - T2.*dS_T  + int_V; % dG   = dH_T - T*dS_T + int dV dP 
% ---> Hans - dG Landau
Q_298   = (1-Tref./Tc0).^(1/4);
Q       = (1-  T2./Tc ).^(1/4);
h_298   = Smax.*Tc0.*(Q_298.^2 - (1/3)*Q_298.^6);
s_298   = Smax.*Q_298.^2;
v_t     = Vmax.*Q_298.^2.*(1 + a0.*(T2-Tref) - 20*a0 .* (sqrt(T2)-  sqrt(Tref)));
intvdP  = (1/3)*v_t.*KT.*((1 + (4*(P2-Pref            ))./KT).^(3/4) - 1);               
G_land  = Smax .* ( (T2 - Tc) .*Q.^2 + (1/3).*Tc.*Q.^6);
G_exc   = h_298 - T2.*s_298 + intvdP + G_land;

dG_q_a  = dG_q/mmol*1000;
dG_q_b  = (dG_q + land*G_exc)/mmol*1000;
rho_q_a =  1./(          diff(dG_q_a    ,1,2)/dP_Pa);
rho_q_b =  1./(          diff(dG_q_b    ,1,2)/dP_Pa);

Stab2                = zeros(length(T),length(P));
Stab2( dG_q_b  < dG_coe &  (T2<=Tc)) = 1;  % identify a-b qz transition
Stab2( dG_q_b  < dG_coe &  (T2 >Tc)) = 2;  % identify a-b qz transition

rho2                     = rho_coe;
rho2(dG_q_b  < dG_coe &  (T2<=Tc) ) = rho_q_b(dG_q_b  < dG_coe &  (T2<=Tc) );
rho2(dG_q_b  < dG_coe &  (T2> Tc) ) = rho_q_a(dG_q_b  < dG_coe &  (T2> Tc) );




dG_q    = dG_q + land*G_exc .* (T2<=Tc);
% Physical properties
dG_q     = dG_q/mmol*1000; % from kJ/mol to J/kg
dP_Pa    = dP*1000*1e5;            % from kbar to Pa
rho_q    =  1./(          diff(dG_q    ,1,2)/dP_Pa);
rho_tcax = (rho_q(1:end-1,:)+rho_q(2:end,:))/2;
rho_tcay = (rho_q(:,1:end-1)+rho_q(:,2:end))/2;
alp_q    = -1./rho_tcax.*diff(rho_q,1,1)/dT;               % thermodynamically consistant alpha at cst P
bet_q    =  1./rho_tcay.*diff(rho_q,1,2)/dP_Pa;            % thermodynamically consistant beta at cst T

m0= 448538.0025;
m1= 1.645295858;
m2= -17.5543599;
k0= 371254.6662;
k1= 8.16916871;
k2= -140.80;

G_q = m0 + m1*(P2-Pref)*1e3 + m2*(T2-Tref);
K_q = k0 + k1*(P2-Pref)*1e3 + k2*(T2-Tref);

%% Phase diagram
Stab                 = zeros(length(T),length(P));
Stab(dG_coe  < dG_q) = 1;
%% Thermodynamic properties
dG_q_p    = (dG_q(:,1:end-1)+dG_q(:,2:end))/2; 
dG_coe_p  = (dG_coe(:,1:end-1)+dG_coe(:,2:end))/2; 
dG_q_pt   = (dG_q(2:end-0,1:end-1)+dG_q(1:end-1,2:end-0)+dG_q(1:end-1,1:end-1)+dG_q(2:end,2:end))/4; 
dG_coe_pt = (dG_coe(2:end-0,1:end-1)+dG_coe(1:end-1,2:end-0)+dG_coe(1:end-1,1:end-1)+dG_coe(2:end,2:end))/4;
dG_q_pi   = dG_q(:,2:end-1);
dG_coe_pi = dG_coe(:,2:end-1);
rho                      = rho_q;
rho(dG_coe_p  < dG_q_p ) = rho_coe(dG_coe_p  < dG_q_p);
alp                      = alp_q;
alp(dG_coe_pt < dG_q_pt) = alp_coe(dG_coe_pt < dG_q_pt);
bet                      = bet_q;
bet(dG_coe_pi < dG_q_pi) = bet_coe(dG_coe_pi < dG_q_pi);
G                        = G_q;
G(dG_coe_pi < dG_q_pi)   = G_coe(dG_coe_pi < dG_q_pi);
K                        = K_q;
K(dG_coe_pi < dG_q_pi)   = K_coe(dG_coe_pi < dG_q_pi);
nu                       = (3.*K-2.*G)./2./(3.*K+G);

% %% VERIFICATION: compare numerically obtained phase diagram with analytic phase diagram boundaries
% q - coe --- HP 98
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
Pr       = (-dG+dV*Pref) / dV ;
ord      = (-dS/dV*Tr - Pref*1000) /1000;
P_coe_q = T*dS/dV/1000 + ord-0*1.2;

% Convert 
DeltaG = (dG_coe- dG_q)*mmol;

end