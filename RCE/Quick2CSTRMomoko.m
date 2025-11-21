clc;
clear;

% MY FLUX ACROSS THE BOUNDARY IS NOT PERFECTLY CORRECT YET


%#ok<*NUSED>
%#ok<*GVMIS>
%#ok<*INUSD>

%dC/dt = F/V*(Cin - C) - nFE/RT
%   where here C = [C1co2, C1co, C2co2, C2co]

% 5 bar experimental data
tEXP = [10,30,50,70,90,110,130,150] %min
CcoEXP = [117, 638, 1364, 1596, 1725, 1779] %ppm
% ChEXP = [0, 36, 74, 80, 84] %ppm
P = 5 %bar
T = 293; %K

D = 7; %cm
p.F = 50; %cm3/min
p.C0co2 = P/(83.1446*T);
p.C0co = 0;
p.SAliq = pi*(D/2)^2;
Vtot = 400; %total volume of reactor cm3
GLratio = 3/7; %ratio of gas:liquid volume in cstr
p.Vgas = Vtot*(GLratio/(GLratio+1)); %cm3
p.Vliq = Vtot-p.Vgas; %cm3
p.r1 = 1; %mol/min
p.K = 0.06; %cm/min mass transfer coeffienct for co in water 0.14 mm/s

C0 = [0.0004; 0; 0.0004; 0];
tspan = [0 60]; %minutes

[t, C] = ode45(@(t, C) ode(t, C, p), tspan, C0);

function dCdt = ode(t, C, p)

    JA = p.SAliq*p.K*(C(4)-C(2))/p.Vgas;%total flow into gas = SA of liquid surface * kg(Cg - Cl)/V

    %mass balance for gas phase
    dC1co2 = p.F/p.Vgas*(p.C0co2-C(1)) - p.r1;
    dC1co = p.F/p.Vgas*(p.C0co-C(2)) + JA;


    %mass balance for liquid phase
    dC2co2 = 0; %just assume that co2 in liquid is always in equilibrium with the gas
    dC2co = p.r1 - JA; %assuming no co comes back into the liquid

    dCdt = [dC1co2; dC1co; dC2co2; dC2co];
end

figure(1); 
% plot(t, C(:, 1), 'LineWidth', 2); hold on;
plot(t, C(:, 2), 'LineWidth', 2);
xlabel('Time (min)');
ylabel('Outlet Concentrations (mol/cm^3)');
legend('CO2','CO');
hold off;