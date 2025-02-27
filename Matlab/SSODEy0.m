clc;
clear;

% Model is somewhat accurate! Need to multiply v0 only by conversion. Maybe Ku is still off a bit

L = 0.07; %m

T = 622+273; %K
P = 1; %bar or atm
nu = [-1 0 0; -1 -1 0; 1 -1 0; 0 1 0; 3 1 0; 0 0 0]; %removal of hydrogen ignored for now

% ch4; h20; co; co2; h2; ar
sccm0 = [16.16; 48.48; 0; 0; 10; 5];
mols0 = sccm0/22400 * 60;
y0 = mols0/sum(mols0);
sumy0 = sum(y0);
Ctot = P/(8.3144598 * 10^-5 * T); %mol/m3
C0 = y0*Ctot;

u0 = sum(sccm0) * 273 / T * P / 1 / 60; %cm^3/sec
A = pi() * (3.5/10)^2; %cm^2. I assumed a 7 mm diameter
v0 = u0 / A / 100 %m/sec v0 = 0.0047

zSpan = linspace(0,L);
[z,y] = ode45(@(x,y) odefun(x, y, T, P, v0, y0, nu, L), zSpan, y0);
figure(1);
plot(z, y, 'Linewidth', 1);
legend('ych4','yh2o','yco','yco2','yh2','yar');
ylim([0 1]);
ylabel('y_i');
xlabel('z');
title('SS molar composition at ' + string(T) + 'K');

y(end, :)
if T == 622 +273
	[0.071787506 0.298008466 0.043658696 0.043658696 0.495104008 0.047782628] %622
elseif T == 721+273
	[0.019206098 0.189598069 0.090570208 0.034350264 0.619497555 0.046777807] %721
else
	[0]
end

function dydz = odefun(x, y, T, P, v0, y0, nu, L) %#ok<*INUSD>

	v = v0*(3*y0(1,1) - 2*y(1,1))/y0(1,1);
	
	Keqsmr = (101325/100000)^2 * exp(-26830/T + 30.114);
	Keqwgs = exp(4400/T - 4.036);


	Rmatrix = rxneq(y, Keqsmr, Keqwgs, T, L); 
	
	dydz = (nu*Rmatrix)/v;
end

function R = rxneq(y, Keqsmr, Keqwgs, T, L) % units of mol/m3/s

	R = zeros(3,1);

	%initialize constants for the reaction here
	P = 1;
	gasconstant = 8.314; %J / mol·K
	Easmr = 165.740; %kJ/mol
	Asmr = 1.68*10^8;
	
	Eawgs = 89.23; % kJ/mol
	Awgs = 9.90*10^3;
	
	Am = 15/(100^2)/130; %15 cm2/130 mg
	Mc = 13; %13 mg
	PR = 0.9; %90%
	Psc = pi() * 2 * (2.5/1000); %perimeter in m
	Ku = (Am * Mc * PR) / (Psc * L); %current value with area calculation: Ku = 0.1228

	I = 9; %curent in Amps
	F = 96485; %faradays constant

	Ku = 1;

	%ch4; h20; co; co2; h2; ar
	R(1,1) = Ku/L*Asmr*exp(-Easmr*1000/gasconstant/T)*(y(1)*y(2)-P^2*y(3)*y(5)^3/Keqsmr);
	R(2,1) = Ku/L*Awgs*exp(-Eawgs*1000/gasconstant/T)*(y(3)*y(2)-y(4)*y(5)/Keqwgs);
	R(3,1) = I/(2*F*L); 
end