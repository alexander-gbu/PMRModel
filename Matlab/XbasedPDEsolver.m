clc;
clear;

% 48 sccm of steam, 10 sccm h2 and then the steady state is without seperation

L = 0.07;

T = 721.4064371;  %C
T = T + 273.15 %K
P = 1; %bar or atm
nu = [-1 0 0; -1 -1 0; 1 -1 0; 0 1 0; 3 1 0; 0 0 0]; %removal of hydrogen ignored for now

% ch4; h20; co; co2; h2; ar
sccm0 = [16.16; 48.48; 0; 0; 10; 5]
mols0 = sccm0/22400 * 60
y0 = mols0/sum(mols0)
sumy0 = sum(y0)
Ctot = P/(8.3144598 * 10^-2 * T)
C0 = y0*Ctot

u0 = sum(sccm0) * 273 / T * P / 1 / 60; %cm^3/sec
A = pi() * (3.5/10)^2; %cm^2. I assumed a 7 mm diameter
v0 = u0 / A / 2 / 100 %m/sec

Nx = 100;
Nt = 200;
x = linspace(0,L,Nx); %length of reactor
t = linspace(0,10,Nt); %time steps
m = 0;
C = pdepe(m,@pdefun,@icfun,@bcfun,x,t);
 
function [c,f,s] = pdefun(x,t,y,dCdx) %might have to change the C back to u if this thing starts freaking out
	T = 900;
    L = 0.07;
	ych40 = 0.2029;
	v0 = 0.0047;
	
    v = (3*ych40 - 2*y(1,1))/ych40
	
	nu = [-1 0 0; -1 -1 0; 1 -1 0; 0 1 0; 3 1 0; 0 0 0]; %hydrogen removal ignored for now
	R = rxneq(y, Keqsmrcalc(T), Keqwgscalc(T), T, L);
	
	c = zeros(length(y),1) + 1;
	f = zeros(length(y),1);
	s = -v*dCdx + nu*R;
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

function R = rxneq(C, Keqsmr, Keqwgs, T, L) % units of mol/m3/s

	P = 1;
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
	Psc = pi() * (2.5/10)^2 /(100^2); %make perimeter
	Ku = (Am * Mc * PR) / (Psc * L);

	I = 9; %curent in Amps
	F = 96485; %faradays constant NEED TO ADJUST STILL

	%ch4; h20; co; co2; h2; ar
	R(1,1) = Ku/L*Asmr*exp(-Easmr*1000/gasconstant/T)*(C(1)*C(2)/(sum(C)^2)-P^2*C(3)*C(5)^3/(sum(C)^4*Keqsmr)); %SMR NEED TO MULTIPLY BY SOME PARTIAL PRESSURE TERM AND CONVERT LIKE IN THE PAPER
	R(2,1) = Ku/L*Awgs*exp(-Eawgs*1000/gasconstant/T)*(C(3)*C(2)/(sum(C)^2)-C(4)*C(5)/(sum(C)^2*Keqwgs));	%WGS
	R(3,1) = I/(2*F*L); %removal of hydrogen is still off but works with a correction factor
end

function Keq = Keqsmrcalc(T)
    Keq = (101325/100000)^2 * exp(-26830./T + 30.114);
end

function Keq = Keqwgscalc(T)
    Keq = exp(4400./T - 4.036);
end