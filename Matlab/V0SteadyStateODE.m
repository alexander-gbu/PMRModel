clc;
clear;

%48 sccm of steam, 10 sccm h2 and then the steady state is without seperation

L = 0.07;

T = 721.4064371;  %C
T = T + 273.15 %K
P = 1; %bar or atm
M = [16, 18, 28, 44, 2, 40]; %g/mol
nu = [-1 0 0; -1 -1 0; 1 -1 0; 0 1 0; 3 1 0; 0 0 0]; %removal of hydrogen ignored for now

%ch4; h20; co; co2; h2; ar
sccm0 = [16.16; 48.48; 0; 0; 10; 5]
mols0 = sccm0/22400 * 60
y0 = mols0/sum(mols0)
sumy0 = sum(y0)
Ctot = P/(8.3144598 * 10^-2 * T)
C0 = y0*Ctot

u0 = sum(sccm0) * 273 / T * P / 1 / 60; %cm^3/sec
A = pi() * (3.5/10)^2; %cm^2. I assumed a 7 mm diameter
v0 = u0 / A / 2 / 100 %m/sec

%has to be calculated too 
% then can go to C0 using ideal gas law and 
% make it in terms of concentration

zSpan = linspace(0,L);
[z,y] = ode45(@(x,y) odefun(x, y, T, P, v0, C0, M, nu, L), zSpan, C0);
figure(1);
plot(z, y, 'LineWidth', 6);
legend('Cch4','Ch2o','Cco','Cco2','Ch2','Car'); %adjust this ...
% ylim([0 1]);
ylabel('C_i');
xlabel('length of reactor z');
title('Concentrations at SS at 672 C and 1 bar');

'ch4; h20; co; co2; h2; ar'
out = y(end,:)*10^3
'ch4; h2'
[y(end, 1), y(end, 5)]

%conversion calc now we just assume a constant velocity
%X = (transpose(y0) - y(end,:))./transpose(y0) %methane should recalculat this actually

function dydt = odefun(x, y, T, P, v0, y0, M, nu, L)

	rho0 = M*y0;
	rhoz = M*y;
	
	v = rho0*v0/rhoz;
	
	Keqsmr = (101325/100000)^2 * exp(-26830/T + 30.114);
	Keqwgs = exp(4400/T - 4.036);
	% Keqsmr = Keqsmrcalc(T);
	% Keqwgs = Keqwgscalc(T);


	Rmatrix = rxneq(y, Keqsmr, Keqwgs, T, L); 
	
	dydt = (nu*Rmatrix)/v;
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
	Ku = (Am * Mc * PR) / (Psc * L)

	I = 9; %curent in Amps
	F = 96485; %faradays constant NEED TO ADJUST STILL
	
	%ch4; h20; co; co2; h2; ar
	R(1,1) = Ku/L*Asmr*exp(-Easmr*1000/gasconstant/T)*(C(1)*C(2)/(sum(C)^2)-P^2*C(3)*C(5)^3/(sum(C)^4*Keqsmr)); %SMR NEED TO MULTIPLY BY SOME PARTIAL PRESSURE TERM AND CONVERT LIKE IN THE PAPER
	R(2,1) = Ku/L*Awgs*exp(-Eawgs*1000/gasconstant/T)*(C(3)*C(2)/(sum(C)^2)-C(4)*C(5)/(sum(C)^2*Keqwgs));	%WGS
	R(3,1) = I/(2*F*L); %removal of hydrogen is still off but works with a correction factor
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