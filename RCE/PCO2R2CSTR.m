clc;
clear;
close all;

%#ok<*NUSED>
%#ok<*GVMIS>
%#ok<*INUSD>

%dC/dt = F/V*(Cin - C) - nFE/RT
%   where here C = [C1co2, C1co, C2co2, C2co]

% 5 bar experimental data
tEXP = [10, 30, 50, 70, 90, 110]; %min
CcoEXP = [117, 638, 1364, 1596, 1725, 1779]; %ppm
% ChEXP = [0, 36, 74, 80, 84] %ppm
p.P = 5; %bar
T = 293; %K

D = 7; %cm
Fstand = 50; %standard cm3/min
p.F = Fstand/p.P; %actual cm3/min
p.C0co2 = p.P/(83.1446*T);
p.C0co = 0;
p.SAliq = pi*(D/2)^2;
Vtot = 360; %total volume of reactor cm3
GLratio = 4/6; %ratio of gas:liquid volume in cstr
p.Vgas = Vtot*(GLratio/(GLratio+1)); %cm3
p.Vliq = Vtot-p.Vgas; %cm3
p.R = 0.015*0.93*60/(2*96485) %mol/min
kG = 0; % cm/min Gas mass transfer coefficient
kL = 0; % cm/min Liquid mass transfer coefficient
p.HCO = 9.5e-5; %mol/cm3 bar Henrys constant for CO
p.HCO2 = 3.3e-5; %mol/cm3 bar Henrys constant CO2
p.K = 0.2; % = kG*kL/(H*kL+kG) in cm/min mass transfer coeffienct for co in water 0.14 mm/s

%moles [liq co2, liq co, gas co2, gas co]
n0 = [p.P*p.HCO2*p.Vliq; 0; p.C0co2*p.Vgas; 0]
tend = 100; %minutes
tnum = 10000;
tspan = linspace(0, tend, tnum); %minutes
p.dt = tend/tnum;

[t, n] = ode45(@(t, n) ode(t, n, p), tspan, n0);

function dCdt = ode(t, n, p)

    %equilibrium of CO2 in water
    nco2eq = p.HCO2*n(3)/(n(3)+n(4))*p.P*p.Vliq;
    nco2needed = nco2eq - n(1);
    if nco2needed < 0
        nco2needed = 0;
    end
    ndotinliq = nco2needed/p.dt; %average flowrate needed in liquid
    if ndotinliq > p.F*p.C0co2
        ndotinliq = p.F*p.C0co2;
    end

    % 2 film theory
    JA = p.SAliq*(p.K*(n(2)/p.Vliq - p.HCO*n(4)/(n(3)+n(4))*p.P));
    if JA < 0
        JA = 0;
        'JA is negative. Maybe reduce stepsize'
    end

    %mass balance for gas phase
    dnco2liq = 0; %F*Cin is molar flow in
    dncoliq = -JA + p.R;

    %mass balance for liquid phase
    dnco2gas = p.F*p.C0co2 - p.F*n(3)/p.Vgas - p.R; %assume that co2 in liquid is always in equilibrium with the gas
    dncogas = JA - p.F*n(4)/p.Vgas;

    dCdt = [dnco2liq; dncoliq; dnco2gas; dncogas];
end

figure(1);
hold on;
plot(t, n(:, 4)./(n(:, 3) + n(:, 4))*10^6, 'LineWidth', 2);
plot(tEXP, CcoEXP, 'o', 'LineWidth', 2);
xlabel('Time (min)');
ylabel('Outlet Concentrations (ppm)');
legend('CO','EXPCO');
hold off;

figure(2);
hold on;
plot(t, n(:, 1)/p.Vliq, 'LineWidth', 2);
plot(t, n(:, 2)/p.Vliq, 'LineWidth', 2);
xlabel('Time (min)');
ylabel('Liquid Concentrations (mol/cm3)');
legend('CO2','CO');
hold off;

figure(3);
plot(t, n(:, 3)./(n(:, 3) + n(:, 4))*10^6, 'LineWidth', 2);
xlabel('Time (min)');
ylabel('Outlet Concentrations (ppm)');
legend('CO2');

%mole balance hard to do mass balance because of accumulation
