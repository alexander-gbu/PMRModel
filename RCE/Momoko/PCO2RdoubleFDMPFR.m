clc;
clear;
close all;

%#ok<*NUSED>
%#ok<*GVMIS>
%#ok<*INUSD>

% I can try and not save the entire time matrix cuz that might take too long to run but it would eliminate that you can see the conc profile 

h = 7; %guessed as 7cm
tend = 1; %minutes
xmesh = 1000000;
tmesh = 100;
dx = h/(xmesh+1);
% Stability criterion: CFL condition (v * dt / dx <= 1)
dt = tend/(tmesh+1);
xspan = linspace(0,h,xmesh);
tspan = linspace(0,tend,tmesh);

P = 5;
T = 293; %K
gasConst = 83.1446; %bar cm3/mol K
Fstand = 50; %standard cm3/min
F = Fstand/P; %actual cm3/min
C0co = 0;
Vgas = 165; %cm3 approximated based on inserts
Vliq = 180; %cm3
SAliq = Vliq/6.66; %cm2
R = 0.015*1.2*60/(2*96485)/Vliq; %mol/cm3 min
% kG = 0; % cm/min Gas mass transfer coefficient
% kL = 0; % cm/min Liquid mass transfer coefficient
Hco = 9.5e-5; %mol/cm3 bar Henrys solubility constant for CO
Hco2 = 3.3e-5; %mol/cm3 bar Henrys solubility constant CO2
D0co = 1e-1; %guessed
D0coliq = 1e3;
K = 1e-8; % = kG*kL/(H*kL+kG) in cm/min mass transfer coeffienct for co in water 0.14 mm/s
u = 10; F/SAliq; %linear velocity cm/min

C = zeros(length(xspan), length(tspan));
flux = K;

for j = 1:length(tspan)-1
    %calculate CSTR and BC individually
    C(1,j+1) = C(1,j) + dt * (R - flux);
    C(2,j+1) = C(2,j) + dt * flux;
    C(3:xmesh-1,j+1) = C(3:xmesh-1,j) + dt * (-u * (C(4:xmesh,j)-C(2:xmesh-2,j))/(2*dx));
    C(xmesh,j+1) = 0;
    % flux = K*(C(1,j+1)-10*C(2,j+1)); %flux into the gas mol/cm3 min
end

size(tspan)
size(xspan)
size(C)

figure()
surf(tspan, xspan, C, 'edgecolor', 'none');
xlabel('Time (min)');
ylabel('distance from liquid surface (cm)');
title(['CO outlet concentration at ', num2str(P), ' bar'])

figure()
plot(tspan, C(1,:))
title(['CO liquid concentration at ', num2str(P), ' bar'])
