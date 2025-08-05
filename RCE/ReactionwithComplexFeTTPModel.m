clc;
clear;
% close all;

%#ok<*NUSED>
%#ok<*GVMIS>
%#ok<*INUSD>

rpm = 100;

Expdata = readtable('FeTTPReactionsH2OCO2.xlsx');
% Expdata.Properties.VariableNames
Evar = sprintf('x%dIRCorr_V_', rpm);
Ivar = sprintf('x%dCurrent_A_', rpm);
ExpE = Expdata.(Evar); %E in V
ExpI = Expdata.(Ivar); %I in A

xmesh = 200;
tmesh = xmesh;

% CV waveform
c.scan_rate = 0.2;      % V/s
c.E_start = 0.0;
c.E_end = -3.0;
c.half_cycle_time = abs(c.E_end - c.E_start)/c.scan_rate;

time = 4*c.half_cycle_time; %total time of experiment [s] to complete 2 full cycles

if mod(tmesh, 4) ~= 0
    error('tmesh must be divisible by 4');
end

t_half = linspace(0, c.half_cycle_time, floor(tmesh/4));

E1 = c.E_start - c.scan_rate * t_half;
E2 = c.E_end + c.scan_rate * t_half;
E3 = c.E_start - c.scan_rate * t_half;
E4 = c.E_end + c.scan_rate * t_half;
E = [E1, E2, E3, E4];

%Constants
delta = 0.001385/(rpm^0.59) % boundary layer thickness [m]

c.T = 298.0;
c.F = 96485.333; % in C/mol or A*s/mol
c.R = 8.314459848;
c.A = 0.0003; %Area of electrode in m2

c.C_Fe3_i = 0.001*1000.0;   %initial Fe(III) bulk concentration at t=0 in [mol/m3] units
c.C_Fe2_i = 0;              %initial Fe(II) bulk concentration at t=0 in [mol/m3] units
c.C_Fe1_i = 0;   %initial Fe(I) bulk concentration at t=0 in [mol/m3] units
c.C_Fe0_i = 0;   %initial Fe(0) bulk concentration at t=0 in [mol/m3] units
c.C_FeCO2_i = 0;
c.C_H2O_i = 0.06*1000;      %initial water bulk concentration [mol/m3]. water does not dissociate in MeCN
c.C_CO2_i = 0.3*1000;
c.C_CO_i = 0;
c.C_OH_i = 0;

%"new" diffusion coefficients
c.D0_Fe3 = 1.06e-9; %Diffusion coefficient of CO2 in water at 25C at infinite dilution [m2/s]                  e-10 if have adjusted diffusion coefficient. it is somewhere between e-10 and e-11
c.D0_Fe2 = 6.7e-9; %Diffusion coefficient of (CO3)2- in water at 25C at infinite dilution [m2/s]
c.D0_Fe1 = 4.6e-9; %Diffusion coefficient of HCO3- in water at 25C at infinite dilution [m2/s]
c.D0_Fe0 = 5.7e-9; %Diffusion coefficient of HCO3- in water at 25C at infinite dilution [m2/s]
c.D0_FeCO2 = 4e-9; %GUESSED PARAMETER THIS WILL PROBABLY NEED TO BE ADJUSTED
c.D0_H2O = 5.78e-9; %https://doi.org/10.1007/978-3-662-54089-3
c.D0_CO2 = 2.89e-9; %https://pubs.acs.org/doi/full/10.1021/acs.jpcc.3c03992
c.D0_CO = 6e-9;
c.D0_OH = 2.1e-9;

c.E0_3_2 = -0.2; %from data
c.E0_2_1 = -1.3;
c.E0_1_0 = -2.08;

function [rFeCO2, rCO2_CO] = HomoReaction(C_Fe0, C_FeCO2, C_H2O, C_CO2)
    %FeTTP(0) + CO2 = FeTTP(II)CO2
    kFeCO2 = 2*10^-1;
    rFeCO2 = kFeCO2.*C_Fe0.*C_CO2;
    %FeTTP(II)CO2 + H2O = FeTTP(II) + CO + 2HO-
    kco = 1*10^0;
    rCO2_CO = kco*C_FeCO2.*C_H2O;
end

function [r3_2, r2_1, r1_0] = ElecReactions(C, E, const)
    k0_3_2 = 0.00002; % rate constant Fe(III) to Fe(II) (m/s)                 %0.00002                     higher reaction rates mean steaper slopes
    k0_2_1 = 0.00002;% rate constant Fe(II) to Fe(I) (m/s)
    k0_1_0 = 0.00001;
    alpha = 0.5;

    r3_2 = k0_3_2*(C(1)*exp(-alpha*(E-const.E0_3_2)*const.F/const.R/const.T)-C(2)*exp((1-alpha)*(E-const.E0_3_2)*const.F/const.R/const.T)); %mol/s/m2
    r2_1 = k0_2_1*(C(2)*exp(-alpha*(E-const.E0_2_1)*const.F/const.R/const.T)-C(3)*exp((1-alpha)*(E-const.E0_2_1)*const.F/const.R/const.T));
    r1_0 = k0_1_0*(C(3)*exp(-alpha*(E-const.E0_1_0)*const.F/const.R/const.T)-C(4)*exp((1-alpha)*(E-const.E0_1_0)*const.F/const.R/const.T));
end

m = 0; 
xspan = linspace(0,delta,xmesh);
tspan = linspace(0,time,tmesh);

sol = pdepe(m, @(x,t,u,dudx) pde(x,t,u,dudx,c), @(x) pdeic(x,c), ...
            @(xl,ul,xr,ur,t) pdebc(xl,ul,xr,ur,t,c), xspan, tspan);
u1 = sol(:,:,1); %sol(t(i), x(j), component)
u2 = sol(:,:,2);
u3 = sol(:,:,3);
u4 = sol(:,:,4);

dx = xspan(xmesh) - xspan(xmesh-1);
current_Fe3 = -c.F*c.A*c.D0_Fe3*(sol(:,xmesh,1)-sol(:,xmesh-1,1))/dx; % reaction is at the right boundary
current_Fe2 = -c.F*c.A*c.D0_Fe2*(sol(:,xmesh,2)-sol(:,xmesh-1,2))/dx;
current_Fe1 = -c.F*c.A*c.D0_Fe1*(sol(:,xmesh,3)-sol(:,xmesh-1,3))/dx;
current_Fe0 = -c.F*c.A*c.D0_Fe0*(sol(:,xmesh,4)-sol(:,xmesh-1,4))/dx;
current_FeCO2 = -c.F*c.A*c.D0_FeCO2*(sol(:,xmesh,5)-sol(:,xmesh-1,5))/dx;
current_CO = -c.F*c.A*c.D0_CO*(sol(:,xmesh,8)-sol(:,xmesh-1,8))/dx;
global_currentOld = (current_Fe2+2*current_Fe1+3*current_Fe0+100*current_CO); %+current_FeCO2

% for i = 1:tmesh
%     t = tspan(i);
%     [r3_2, r2_1, r1_0] = ElecReactions(sol(i, xmesh-1, :), E(i), c);
%     current_Fe3(i) = (r3_2); %1*c.F*c.A*
%     current_Fe2(i) = (r2_1); %1*c.F*c.A*
%     current_Fe1(i) = (r1_0); %1*c.F*c.A*
%     current_Fe0(i) = (0); %1*c.F*c.A*
% end
% % current_Fe3
% global_current = (-current_Fe3-current_Fe2-current_Fe1)/1000; %-current_Fe0

% figure(1);
% surf(xspan,tspan,u1/1000.0,'edgecolor','none');
% xlim([0.0, delta]);
% %title('K^{+} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.008 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('Fe(III)[mol/L]');
% view(30,20);

figure(2);
surf(xspan,tspan,u2/1000.0,'edgecolor','none');
xlim([0.0, delta]);
%title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
xlabel('Distance x [m]');
ylabel('Time t [s]');
zlabel('Fe(II)[mol/L]');
view(30,20);

figure(3);
surf(xspan,tspan,u3/1000.0,'edgecolor','none');
xlim([0.0, delta]);
%title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
xlabel('Distance x [m]');
ylabel('Time t [s]');
zlabel('Fe(I)[mol/L]');
view(30,20);

figure(4);
surf(xspan,tspan,u4/1000.0,'edgecolor','none');
xlim([0.0, delta]);
%title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
xlabel('Distance x [m]');
ylabel('Time t [s]');
zlabel('Fe(0)[mol/L]');
view(30,20);

figure(5);
surf(xspan,tspan,sol(:,:,5)/1000.0,'edgecolor','none');
xlim([0.0, delta]);
%title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
xlabel('Distance x [m]');
ylabel('Time t [s]');
zlabel('FeCO2 [mol/L]');
view(30,20);

% figure(6);
% surf(xspan,tspan,sol(:,:,6)/1000.0,'edgecolor','none');
% xlim([0.0, delta]);
% %title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('CO2[mol/L]');
% view(30,20);

figure(7);
surf(xspan,tspan,sol(:,:,7)/1000.0,'edgecolor','none');
xlim([0.0, delta]);
%title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
xlabel('Distance x [m]');
ylabel('Time t [s]');
zlabel('CO[mol/L]');
view(30,20);

% figure(8)
% plot(tspan,current_Fe3,tspan,current_Fe2,tspan,current_Fe1,tspan,current_Fe0,'LineWidth',1.5); %,tspan,global_currentOld,
% ylim([-10, 10]);
% ylabel('reaction rates [mol/m2/s]');
% xlabel('Time t [s]');
% legend('Fe3', 'Fe2', 'Fe1', 'Fe0'); %, 'global'

figure(9)
plot(ExpE, ExpI, 'r-', E(floor(tmesh/2):end), global_currentOld(floor(tmesh/2):end), 'b--'); %(floor(tmesh/2):end) 
xlabel('E (V)');
ylabel('Current (A)');
title('Experimental vs Model Data');
legend('Experimental', 'Model');

% figure(10)
% plot(ExpE, ExpI, 'r-', E, global_current, 'b--'); %(floor(nmesh/2):end)
% xlabel('E (V)');
% ylabel('Current (A)');
% title('Experimental vs Model Data');
% legend('Experimental', 'Model');

function [c,f,s] = pde(x,t,u,dudx,const)
    [rFeCo2, rCO] = HomoReaction(u(4), u(5), u(6), u(7));
    c = [1; 1; 1; 1; 1; 1; 1; 1; 1];
    f = [const.D0_Fe3*dudx(1); const.D0_Fe2*dudx(2); const.D0_Fe1*dudx(3); const.D0_Fe0*dudx(4); ...
            const.D0_FeCO2*dudx(5); const.D0_H2O*dudx(6); const.D0_CO2*dudx(7); const.D0_CO*dudx(8); const.D0_OH*dudx(9)];
    s = [0; rCO; 0; -rFeCo2; rFeCo2-rCO; -rFeCo2; -rCO; rCO; rCO]; %we change this later so that the iron can react in the entire solution
end

function u0 = pdeic(x, c) 
    u0 = [c.C_Fe3_i; c.C_Fe2_i; c.C_Fe1_i; c.C_Fe0_i; c.C_FeCO2_i; c.C_H2O_i; c.C_CO2_i; c.C_CO_i; c.C_OH_i];
end

function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t,c)
    if t <= c.half_cycle_time
        E = c.E_start - c.scan_rate*t;
    elseif (t > c.half_cycle_time) && (t <= 2*c.half_cycle_time)
        E = c.E_end + c.scan_rate*(t-c.half_cycle_time);
    elseif (t > 2*c.half_cycle_time) && (t <= 3*c.half_cycle_time)
        E = c.E_start - c.scan_rate*(t-2*c.half_cycle_time);
    elseif (t > 3*c.half_cycle_time) && (t <= 4*c.half_cycle_time)
        E = c.E_end + c.scan_rate*(t-3*c.half_cycle_time);
    end

    [r3_2, r2_1, r1_0] = ElecReactions(ur, E, c); %formation of Fe2, formation of Fe1, formation of Fe0

    ql = [0; 0; 0; 0; 0; 0; 0; 0; 0];
    pl = [ul(1)-c.C_Fe3_i; ul(2)-c.C_Fe2_i; ul(3)-c.C_Fe1_i; ul(4)-c.C_Fe0_i; ul(5)-c.C_FeCO2_i; ul(6)-c.C_H2O_i; ul(7)-c.C_CO2_i; ul(8)-c.C_CO_i; ul(9)-c.C_OH_i];
    qr = [1; 1; 1; 1; 1; 1; 1; 1; 1];
    pr = [r3_2; (-r3_2+r2_1); (-r2_1+r1_0); (-r1_0); 0; 0; 0; 0; 0]; % f = -D*dC/dx = r     [m2/s] [mol/m3]/[m] = [mol/m2/s]; 
end