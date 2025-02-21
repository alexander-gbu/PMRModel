clc;
clear;

%48 sccm of steam, 10 sccm h2 and then the steady state is without
%seperation


%constant initialization
P = 3; %atm
T = 1000;  %K
u0 = 79.64 * 273 / T * P / 1 / 60; %cm^3/sec
A = pi() * (3.5/10)^2; %cm^2. I assumed a 7 mm diameter
v0 = u0 / A / 100 %m/sec
Keqsmr = 101325^2 * exp(-26830/T + 30.114);
Keqwgs = exp(4400/T - 4.036);

L = 0.07; %m

%boundry conditions
%ch4, h20, co, co2, h2, ar
Cl1 = [0.3; 0.6; 0; 0; 0; 0.1];

%initial conditions
C0i = [0; 0.66; 0; 0; 0; 0.34];

%reaction stoichiometry
%stoichiometry matrix
%nu_ij of component i reaction j
%Overall (SMR), (WGS), removal of hydrogen
%ch4; h20; co; co2; h2; ar
nu = [-1 0 0; -1 -1 0; 1 -1 0; 0 1 0; 3 1 -1; 0 0 0];

%Sets the initial conditions row 1 is ch4, 2 is h2o, 3 is co2, 4 is h2, 5 is ar
Nx = 10;
Nt = 50;
x = linspace(0,L,Nx);
t = linspace(0,10,Nt);
m = 0;
C = pdepe(m,@pdefun,@icfun,@bcfun,x,t);

% x
% t
% C(:,:,1)

% figure(1);
% surf(x,t,C(:,:,5), 'EdgeColor', 'none')
% title('C_m(x,t)')
% xlabel('Length of reactor z')
% ylabel('Time t')


figure(2);
plot(x, C(Nt,:,1), x, C(Nt,:,2), x, C(Nt,:,3), x, C(Nt,:,4), x, C(Nt,:,5), x, C(Nt,:,6), 'LineWidth',2);
ylim([0 1]);
legend('Cch4','Ch2o','Cco', 'Cco2','Ch2','Car');
ylabel('y_i');
xlabel('Length of reactor z');
title('Concentrations at t = 10s');


function [c,f,s] = pdefun(x,t,C,dCdx) %might have to change the C back to u if this thing starts freaking out
	T = 900;
	R = rxneq(C, Keqsmrcalc(T), Keqwgscalc(T), T);
	M = [16, 18, 28, 44, 2, 40];
	Cl1 = [0.3; 0.6; 0; 0; 0; 0.1];
	
	v0 = 0.0282;
	rho0 = M*Cl1;
	rhoz = M*C;
	
	v = rho0*v0/rhoz;
	
	nu = [-1 0 0; -1 -1 0; 1 -1 0; 0 1 0; 3 1 -1; 0 0 0];
	
	c = zeros(length(C),1) + 1;
	f = zeros(length(C),1);
	s = -v*dCdx + nu*R;
end

function C0 = icfun(x)
	C0 = [0; 0.66; 0; 0; 0; 0.34];
end

function [pL,qL,pR,qR] = bcfun(xL,CL,xR,CR,t)
	pL = CL - [0.3; 0.6; 0; 0; 0; 0.1];
	qL = zeros(length(CR),1);
	pR = zeros(length(CR),1);
	qR = zeros(length(CR),1) + 1;
end

function R = rxneq(C, Keqsmr, Keqwgs, T) %is rate of appearance of CO2

	R = zeros(3,1);

	%initialize constants for the reaction here
	gasconstant = 8.314; %J / mol·K
	Easmr = 165.740; %kJ/mol
	Asmr = 1.68*10^8;
	
	Eawgs = 89.23; % kJ/mol
	Awgs = 9.90*10^3;
	
	Am = 15/(100^2)/130; %15 cm2/130 mg
	Mc = 13; %13 mg
	PR = 0.5; %50%
	Psc = pi() * (3.5/10)^2 /(100^2); %m2
	L = 0.07; %m
	Ku = (Am * Mc * PR) / (Psc * L);
	
	I = 9; %curent in Amps
	F = 96485; %faradays constant NEED TO ADJUST STILL
	
	V = pi() * (3.5/10)^2 /(100^2) * L; %m3
	
	%ch4; h20; co; co2; h2; ar
	R(1,1) = Ku*Asmr*exp(-Easmr*1000/gasconstant/T)*(C(1)*C(2)-C(3)*C(5)^3/Keqsmr); %SMR NEED TO MULTIPLY BY SOME PARTIAL PRESSURE TERM AND CONVERT LIKE IN THE PAPER
	R(2,1) = Ku*Awgs*exp(-Eawgs*1000/gasconstant/T)*(C(3)*C(2)-C(4)*C(5)/Keqwgs);	%WGS
	R(3,1) = I/(2*F)/(L); %removal of hydrogen is still off but works with a correction factor
	
end

function Keq = Keqsmrcalc(T)

	R = 8.314; %J / mol·K
	
	%CO2; CO; CH4; H2O; H2 = 00000
	G_const = [-393.4685065	-0.003212871	1.53E-07	9.97E-10	-3.31E-13;
	-110.7347984	-0.086473798	-9.63E-06	8.56E-09	-2.11E-12;
	-68.77179183	0.038329957	9.22E-05	-5.46E-08	1.24E-11;
	-240.6171143	0.035116356	2.02E-05	-9.33E-09	1.78E-12];
	
	Gf = zeros(1,4);
	for i = 1:4
        Gf(1,i) = G_const(i,1) + G_const(i,2)*T + G_const(i,3)*T^2 + G_const(i,4)*T^3 + G_const(i,5)*T^4;
	end
	Keq = exp(-(Gf(1,2) - Gf(1,3) - Gf(1,4))*1000/R/T);
	
end

function Keq = Keqwgscalc(T)

	R = 8.314; %J / mol·K
	
	%%CO2; CO; CH4; H2O; H2 = 00000
	G_const = [-393.4685065	-0.003212871	1.53E-07	9.97E-10	-3.31E-13;
	-110.7347984	-0.086473798	-9.63E-06	8.56E-09	-2.11E-12;
	-68.77179183	0.038329957	9.22E-05	-5.46E-08	1.24E-11;
	-240.6171143	0.035116356	2.02E-05	-9.33E-09	1.78E-12];
	
	Gf = zeros(1,4);
	for i = 1:4
        Gf(1,i) = G_const(i,1) + G_const(i,2)*T + G_const(i,3)*T^2 + G_const(i,4)*T^3 + G_const(i,5)*T^4;
	end
	Keq = exp(-(Gf(1,1) - Gf(1,2) - Gf(1,4))*1000/R/T);
	
end