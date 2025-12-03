clc;
clear;
close all;

%#ok<*NUSED>
%#ok<*GVMIS>
%#ok<*INUSD>

%dC/dt = F/V*(Cin - C) - nFE/RT
%   where here C = [C1co2, C1co, C2co2, C2co]

% 5 bar experimental data
PEXP = [5, 10, 20, 25];
tEXP =  [10, 30, 50, 70, 90, 110, 0, 0;
         0, 30, 50, 70, 90, 110, 130, 0;
         10, 35, 60, 85, 110, 135, 0, 0;
         10, 35, 60, 85, 110, 135, 160, 185]; %min
CcoEXP =    [117, 638, 1364, 1596, 1725, 1779, 0, 0;
             0, 292, 596, 884, 1214, 1217, 1294, 0;
             74, 297, 633, 880, 1097, 1207, 0, 0;
             70, 218, 500, 739, 908, 1051, 1155, 1250]; %ppm
% ChEXP = [0, 36, 74, 80, 84] %ppm
T = 293; %K

Fstand = 50; %standard cm3/min
p.C0co = 0;
p.Vgas = 165; %cm3 approximated based on inserts
p.Vliq = 180; %cm3
p.SAliq = p.Vliq/6.66; %cm2
p.R = 0.015*1.2*60/(2*96485); %mol/min
kG = 0; % cm/min Gas mass transfer coefficient
kL = 0; % cm/min Liquid mass transfer coefficient
p.HCO = 9.5e-5; %mol/cm3 bar Henrys constant for CO
p.HCO2 = 3.3e-5; %mol/cm3 bar Henrys constant CO2
p.K = 0.25; % = kG*kL/(H*kL+kG) in cm/min mass transfer coeffienct for co in water 0.14 mm/s

tend = 200; %minutes
tnum = 100000;
tspan = linspace(0, tend, tnum); %minutes
% p.dt = tend/tnum;

%moles [liq co2, liq co, gas co2, gas co]
for i = [1:4]
    p.P = PEXP(i); %bar
    p.F = Fstand/p.P; %actual cm3/min
    p.C0co2 = p.P/(83.1446*T);
    n0 = [p.P*p.HCO2*p.Vliq; 0; p.C0co2*p.Vgas; 0]

    [t, n] = ode15s(@(t, n) ode(t, n, p), tspan, n0);

    figure();
    hold on;
    plot(t, n(:, 4)./(n(:, 3) + n(:, 4))*10^6, 'LineWidth', 2);
    plot(tEXP(i, :), CcoEXP(i, :), 'o', 'LineWidth', 2);
    xlabel('Time (min)');
    ylabel('Concentration (ppm)');
    title(['CO outlet concentration at ', num2str(p.P), ' bar'])
    legend('Model','Experimental', 'Location', 'southeast');
    hold off;

    figure();
    hold on;
    plot(t, n(:, 1)/p.Vliq, 'LineWidth', 2);
    plot(t, n(:, 2)/p.Vliq, 'LineWidth', 2);
    xlabel('Time (min)');
    ylabel('Liquid Concentrations (mol/cm3) at 5 bar');
    title(['Pressure =', num2str(p.P), ' bar'])
    legend('CO2','CO');
    hold off;

    % figure();
    % plot(t, n(:, 3)./(n(:, 3) + n(:, 4))*10^6, 'LineWidth', 2);
    % xlabel('Time (min)');
    % ylabel('Outlet Concentrations (ppm) at 5 bar');
    % legend('CO2');
end

function dCdt = ode(t, n, p)

    %equilibrium of CO2 in water
    % nco2eq = p.HCO2*n(3)/(n(3)+n(4))*p.P*p.Vliq;
    % nco2needed = nco2eq - n(1);
    % if nco2needed < 0
    %     nco2needed = 0;x
    % end
    % ndotinliq = nco2needed/p.dt; %average flowrate needed in liquid
    % if ndotinliq > p.F*p.C0co2
    %     ndotinliq = p.F*p.C0co2;
    % end

    % 2 film theory
    JA = p.SAliq*p.K*(n(2)/p.Vliq - p.HCO*n(4)/(n(3)+n(4))*p.P);
    if JA < 0
        JA = 0;
        'JA is negative. Maybe reduce stepsize'
    end

    %mass balance for gas phase
    dnco2liq = JA - p.R; %F*Cin is molar flow in
    dncoliq = -JA + p.R;

    %mass balance for liquid phase
    dnco2gas = p.F*p.C0co2 - p.F*n(3)/p.Vgas - p.R; %assume that co2 in liquid is always in equilibrium with the gas
    dncogas = JA - p.F*n(4)/p.Vgas;

    dCdt = [dnco2liq; dncoliq; dnco2gas; dncogas];
end

