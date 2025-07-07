clc;
clear;

T = [700:10:1100];
% T = [900:5:1100];

KeqsmrP = (101325/100000)^2 * exp(-26830./T + 30.114);
KeqwgsP = exp(4400./T - 4.036);

KeqsmrS = zeros(1,length(T));
KeqwgsS = zeros(1,length(T));


for i = 1:length(T)
	KeqsmrS(1,i) = Keqsmrcalc(T(1,i));
	KeqwgsS(1,i) = Keqwgscalc(T(1,i));
end

plot(T, KeqsmrP, T, KeqwgsP, T, KeqsmrS, '--', T, KeqwgsS, '--')
legend('smr paper', 'wgs paper', 'smr self', 'wgs self');
ylabel('Keq (smr or wgs)');
xlabel('Temperature T [K]');
title('Comparison of Keq calculation');


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