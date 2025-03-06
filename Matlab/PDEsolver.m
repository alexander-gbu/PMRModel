clc;
clear;

% This model is very good try with multiplying by the ya0 term does not work
% convert the model to molar flowrate
% 
% I think moving to a mol balance might be best because then when we are removing the hydrogen the fractions will be messed up
% the velocity will also change when were taking out the hydrogen
% Paper timeline?

L = 0.0609; %m
R_inner = 2.39/1000; %m
R_outer = 3.5/1000; %m

T = 721+273 %K
P = 1; %bar or atm
nu = [-1 0 0; -1 -1 0; 1 -1 0; 0 1 0; 3 1 0; 0 0 0]; %removal of hydrogen ignored for now

% ch4; h20; co; co2; h2; ar
sccm0 = [16.16, 48.48, 0, 0, 10, 5];
mols0 = sccm0/22400 * 60;
y0 = mols0/sum(mols0);
sumy0 = sum(y0);
ych40 = y0(1,1)
Ctot = P/(8.3144598 * 10^-5 * T); %mol/m3
C0 = y0*Ctot;

%this is a jumble right now
u0 = sum(sccm0) * 273 / T * P / 1 / 60; %cm^3/sec
A = pi() * ((R_outer*100)^2-(R_inner*100)^2); %cm^2. I assumed a 7 mm diameter with 3.5 mm being r_outer and 2 mm being r_inner
v0 = u0 / A / 100 %m/sec

Nx = 100;
Nt = 100;
x = linspace(0,L,Nx); %length of reactor
t = linspace(0,15,Nt); %time steps
m = 0;
C = pdepe(m,@(x,t,y,dydx) pdefun(x,t,y,dydx,T,P,L,ych40,v0,R_outer),@icfun,@bcfun,x,t);

figure(1);
plot(x, C(end,:,1), x, C(end,:,2), x, C(end,:,3), x, C(end,:,4), x, C(end,:,5), x, C(end,:,6), 'LineWidth',1);
ylim([0 1]);
legend('ych4','yh2o','yco', 'yco2','yh2','yar');
ylabel('y_i');
xlabel('Length of reactor z');
title('Concentrations at t = 15s at ' + string(T) + 'K or ' + string(T-273) + 'C');
 
function [c,f,s] = pdefun(x,t,y,dydx,T,P,L,ych40,v0,R_outer) %might have to change the C back to u if this thing starts freaking out

    v = v0*(3*ych40 - 2*y(1,1))/ych40;
	
	nu = [-1 0 0; -1 -1 0; 1 -1 0; 0 1 0; 3 1 0; 0 0 0]; %hydrogen removal ignored for now
	R = rxneq(y, Keqsmrcalc(T), Keqwgscalc(T), T, P, L, R_outer);
	
	c = zeros(length(y),1) + 1;
	f = zeros(length(y),1);
	s = -v*dydx + nu*R;
end

%initial conditions
function C0 = icfun(x) 
	C0 = [0; 0.66; 0; 0; 0; 0.34];
end

%boundary conditions
function [pL,qL,pR,qR] = bcfun(xL,CL,xR,CR,t)
	pL = CL - [0.2029; 0.6087; 0; 0; 0.1256; 0.0628];
	qL = zeros(length(CR),1);
	pR = zeros(length(CR),1);
	qR = zeros(length(CR),1) + 1;
end

function R = rxneq(y, Keqsmr, Keqwgs, T, P, L, R_outer) % units of mol/m3/s

	R = zeros(3,1);

	%initialize constants for the reaction here
	gasconstant = 8.314; %J / mol·K
	Easmr = 165.740; %kJ/mol
	Asmr = 1.68*10^8;
	
	Eawgs = 89.23; % kJ/mol
	Awgs = 9.90*10^3;
	
	Am = 15/(100^2)/130; %15 cm2/130 mg
	Mc = 13; %13 mg
	PR = 0.9; %90%
	Psc = pi() * 2 * (2.5/1000); %perimeter in m
	Ku = (Am * Mc * PR) / (Psc * L); %current value with area calculation: Ku = 0.1411

	I = 9; %curent in Amps
	F = 96485; %faradays constant

	Ku = 1;

	%ch4; h20; co; co2; h2; ar
	R(1,1) = 2/R_outer*Ku*Asmr*exp(-Easmr*1000/gasconstant/T)*(y(1)*y(2)-P^2*y(3)*y(5)^3/Keqsmr)*(8.3144598 * 10^-5 * T)/P;
	R(2,1) = 2/R_outer*Ku*Awgs*exp(-Eawgs*1000/gasconstant/T)*(y(3)*y(2)-y(4)*y(5)/Keqwgs)*(8.3144598 * 10^-5 * T)/P;
	R(3,1) = I/(2*F*L); 
end

function Keq = Keqsmrcalc(T)
    Keq = (101325/100000)^2 * exp(-26830./T + 30.114);
end

function Keq = Keqwgscalc(T)
    Keq = exp(4400./T - 4.036);
end