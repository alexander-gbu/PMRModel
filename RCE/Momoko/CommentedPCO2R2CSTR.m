clc;
clear;
% close all; %closes all the plots

%dC/dt = F/V*(Cin - C) - nFE/RT
%   where here C = [Cco2liq, Ccoliq, Cco2gas, Ccogas]

% SET PRESSURE
p.P = 5; %bar

% Experimental data
if p.P == 5
    tEXP    = [10, 30, 50, 70, 90, 110];
    CcoEXP  = [117, 638, 1364, 1596, 1725, 1779];
elseif p.P == 10
    tEXP    = [0, 30, 50, 70, 90, 110, 130];
    CcoEXP  = [0, 292, 596, 884, 1214, 1217, 1294];
elseif p.P == 20
    tEXP    = [10, 35, 60, 85, 110, 135];
    CcoEXP  = [74, 297, 633, 880, 1097, 1207];
elseif p.P == 25
    tEXP    = [10, 35, 60, 85, 110, 135, 160, 185];
    CcoEXP  = [70, 218, 500, 739, 908, 1051, 1155, 1250];
end

% Reactor constants [Can be adjusted]
T = 293; %K
Fstand = 50; %cm3/min standard ccm 
p.C0co = 0; %inlet CO concentration = 0
p.Vgas = 265; %cm3 gas volume in reactor approximated based on inserts
p.Vliq = 180; %cm3 liquid volume in reactor
hliq = 8; %cm liquid height
p.SAliq = p.Vliq/hliq; %cm2 liquid surface area for mass transfer
FE = 1.2; %faradaic efficiecy
p.R = 0.015*FE*60/(2*96485); %mol/min Reaction rate based on current
p.HCO = 9.5e-5; %mol/cm3 bar Henrys constant for CO
p.HCO2 = 3.3e-5; %mol/cm3 bar Henrys constant CO2

% = kG*kL/(H*kL+kG) 
p.K = 1; % cm/min mass transfer coeffienct for co in water (accross the liquid-gas interface)
    % kG = 0; % cm/min Gas mass transfer coefficient
    % kL = 0; % cm/min Liquid mass transfer coefficient

tend = 200; %minutes
tnum = 10000;
tspan = linspace(0, tend, tnum); %minutes

%moles [liq co2, liq co, gas co2, gas co]
p.F = Fstand/p.P; %actual volumetric flowrate cm3/min
p.C0co2 = p.P/(83.1446*T); %initial CO2 concentration in the gas phase (calculated from partial pressure)
n0 = [p.P*p.HCO2*p.Vliq; 0; p.C0co2*p.Vgas; 0]; %initialization of the moles present in both phases. (liquid phase based on henry's constant, gas phase based on partial pressure)

[t, n] = ode15s(@(t, n) ode(t, n, p), tspan, n0);

%plots outlet concentration of CO in ppm
figure();
hold on;
plot(t, n(:, 4)./(n(:, 3) + n(:, 4))*10^6, 'LineWidth', 2);
plot(tEXP, CcoEXP, 'o', 'LineWidth', 2);
xlabel('Time (min)');
ylabel('Concentration (ppm)');
title(['CO outlet concentration at ', num2str(p.P), ' bar'])
legend('Model','Experimental', 'Location', 'southeast');
hold off;

%plots concentrations in liquid in mol/cm3
figure();
hold on;
plot(t, n(:, 1)/p.Vliq, 'LineWidth', 2);
plot(t, n(:, 2)/p.Vliq, 'LineWidth', 2);
xlabel('Time (min)');
ylabel('Concentration (mol/cm3)');
title(['Liquid concentration at Pressure =', num2str(p.P), ' bar'])
legend('CO2','CO');
hold off;

%this is the ODE solvet that solves the ODEs
function dCdt = ode(t, n, p)
    % 2 film theory flux calculation
    JA = p.SAliq*p.K*(n(2)/p.Vliq - p.HCO*n(4)/(n(3)+n(4))*p.P);
    if JA < 0
        JA
        % JA = 0; %Can un-comment this if the simulation is breaking or graph is occilating too much
        'JA is negative. Maybe reduce stepsize'
    end

    %mass balance for gas phase
    %right now CO2 is just replenished at the same rate that
    dnco2liq = JA - p.R; %F*Cin is molar flow in 
    dncoliq = -JA + p.R;

    %mass balance for liquid phase
    dnco2gas = p.F*p.C0co2 - p.F*n(3)/p.Vgas - p.R; %assume that co2 in liquid is always in equilibrium with the gas
    dncogas = JA - p.F*n(4)/p.Vgas;

    dCdt = [dnco2liq; dncoliq; dnco2gas; dncogas];
end
