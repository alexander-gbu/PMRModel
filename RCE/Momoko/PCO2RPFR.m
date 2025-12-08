clc;
clear;
close all;

%#ok<*NUSED>
%#ok<*GVMIS>
%#ok<*INUSD>

%dC/dt + u0 dC/dz = Di dC2/dz2 + ri

% PEXP = [5, 10, 20, 25];
% tEXP =  [10, 30, 50, 70, 90, 110, 0, 0;
%          0, 30, 50, 70, 90, 110, 130, 0;
%          10, 35, 60, 85, 110, 135, 0, 0;
%          10, 35, 60, 85, 110, 135, 160, 185]; %min
% CcoEXP =    [117, 638, 1364, 1596, 1725, 1779, 0, 0;
%              0, 292, 596, 884, 1214, 1217, 1294, 0;
%              74, 297, 633, 880, 1097, 1207, 0, 0;
%              70, 218, 500, 739, 908, 1051, 1155, 1250]; %ppm
% ChEXP = ;

%height of gas chamber
h = 7; %guessed as 7cm
tend = 200; %minutes

p.P = 5;
p.T = 293; %K
p.gasConst = 83.1446; %bar cm3/mol K
Fstand = 50; %standard cm3/min
p.F = Fstand/p.P; %actual cm3/min
p.C0co = 0;
p.Vgas = 165; %cm3 approximated based on inserts
p.Vliq = 180; %cm3
p.SAliq = p.Vliq/6.66; %cm2
p.R = 0.015*1.2*60/(2*96485)/p.Vliq; %mol/min
kG = 0; % cm/min Gas mass transfer coefficient
kL = 0; % cm/min Liquid mass transfer coefficient
p.Hco = 9.5e-5; %mol/cm3 bar Henrys solubility constant for CO
p.Hco2 = 3.3e-5; %mol/cm3 bar Henrys solubility constant CO2
p.D0co = 1e-1; %guessed
p.D0coliq = 1e3;
p.K = 0.1; % = kG*kL/(H*kL+kG) in cm/min mass transfer coeffienct for co in water 0.14 mm/s
p.v0 = p.F/p.SAliq;

xmesh = 100;
tmesh = 100;
xspan = linspace(0,h,xmesh);
tspan = linspace(0,tend,tmesh);
m = 0;

Conc = pdepe(m, @(x,t,u,dudx) pde(x,t,u,dudx,p), @(x) pdeic(x,p), ...
        @(xl,ul,xr,ur,t) pdebc(xl,ul,xr,ur,t,p), xspan, tspan); %sol(t(i), x(j), component)
molcos = Conc*p.F;
totmol = p.P/(p.gasConst*p.T)*p.F; %total concentration * volume flow = tot mol / s
ppms = molcos/totmol * 1e6; %in the future totmol will just be sum of all species

figure();
% plot(tspan, ppms(:, xmesh, 1), 'LineWidth', 2);
plot(tspan, ppms(:, xmesh, 2), 'LineWidth', 2);
xlabel('Time (min)');
ylabel('Concentration (ppm)');
title(['CO outlet concentration at ', num2str(p.P), ' bar'])
legend('CO in liq', 'CO in gas', 'Location', 'southeast');

figure();
surf(xspan, tspan, ppms(:, :, 1), 'edgecolor','none');

figure();
surf(xspan, tspan, ppms(:, :, 2), 'edgecolor','none');

function [c,f,s] = pde(x,t,u,dudx,p)
    c = [1; 1];
    f = [p.D0coliq*dudx(1); p.D0co*dudx(2)];
    s = [p.R; -p.v0*dudx(2)];
end

function u0 = pdeic(x, p) 
    u0 = [0; 0];
end

function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t,p)
    J = p.K*(p.gasConst*p.T*ul(2) - ur(1)/p.Hco);
    % if J < 0
    %     'J is negative'
    % end
    ql = [1; 1];
    pl = [0; -J];
    qr = [1; 1];
    pr = [-J; 0]; % f = -D*dC/dx = r     [m2/s] [mol/m3]/[m] = [mol/m2/s]; 
end