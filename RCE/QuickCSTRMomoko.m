clc;
clear;

%SO I REALIZED IT IS FINE TO DO IT LIKE THAT BECAUSE ALL YOU DO IS MULTIPLY IT BY THE PARTITION COEFFICIENT

%#ok<*NUSED>
%#ok<*GVMIS>
%#ok<*INUSD>

%dC/dt = F/V*(Cin - C) - nFE/RT
%   where here C = [C1co2, C1co, C2co2, C2co]

p.F = 10; %cm3/s
p.C0co2 = 0.5;
p.C0co = 0;
Vtot = 100; %total volume of reactor cm3
GLratio = 1; %ratio of gas:liquid volume in cstr
p.Vgas = Vtot*(GLratio/(GLratio+1)); %cm3
p.Vliq = Vtot-p.Vgas; %cm3
p.r1 = 0.01; %mol/s

C0 = [0.5; 0];
tspan = [0 50]; %seconds

[t, C] = ode45(@(t, C) ode(t, C, p), tspan, C0)

function dCdt = ode(t, C, p)

    %mass balance for gas phase
    dC1co2 = p.F/p.Vgas*(p.C0co2-C(1)) - p.r1;
    dC1co = p.F/p.Vgas*(p.C0co-C(2)) + p.r1;


    %mass balance for liquid phase



    dCdt = [dC1co2; dC1co];
end

figure(1); 
plot(t, C(:, 1), 'LineWidth', 2); hold on;
plot(t, C(:, 2), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Concentration in gas phase (mol/m^3)');
legend('CO2','CO');
hold off;