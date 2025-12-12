clc;
clear;
close all;

%#ok<*NUSED>
%#ok<*GVMIS>
%#ok<*INUSD>

% I can try and not save the entire time matrix cuz that might take too long to run but it would eliminate that you can see the conc profile 

h = 7; %guessed as 7cm
tend = 50; %minutes
xmesh = 10000;
tmesh = 10000;
dx = h/xmesh;
% Stability criterion: CFL condition (v * dt / dx <= 1)
dt = tend/tmesh;
xspan = linspace(0,h+dx,xmesh+1);
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
R = 0.015*1.2*60/(2*96485); %mol/cm3 min
% kG = 0; % cm/min Gas mass transfer coefficient
% kL = 0; % cm/min Liquid mass transfer coefficient
Hco = 9.5e-5; %mol/cm3 bar Henrys solubility constant for CO
Hco2 = 3.3e-5; %mol/cm3 bar Henrys solubility constant CO2
D0co = 1e-1; %guessed
D0coliq = 1e3;
K = 100; % = kG*kL/(H*kL+kG) in cm/min mass transfer coeffienct for co in water 0.14 mm/s
u = F/SAliq; %linear velocity cm/min

C = zeros(length(xspan), length(tspan));
flux = 0;

for j = 1:length(tspan)-1
    %calculate CSTR and BC individually
    C(1,j+1) = C(1,j) + dt * (R - flux);
    C(2,j+1) = C(2,j) + dt * (-u * (C(3+1,j)-C(2,j))/dx + flux);
    C(3:xmesh,j+1) = C(3:xmesh,j) + dt * (-u * (C(3:xmesh,j)-C(2:xmesh-1,j))/dx);
    C(xmesh+1,j+1) = C(xmesh,j);
    flux = K*(C(1,j+1)-10*C(2,j+1)); %flux into the gas mol/cm3 min
end

figure()
surf(tspan, xspan, C, 'edgecolor', 'none');
xlabel('Time (min)');
ylabel('distance from liquid surface (cm)');
title(['CO outlet concentration at ', num2str(P), ' bar'])
