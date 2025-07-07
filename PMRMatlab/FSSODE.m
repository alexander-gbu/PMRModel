clc;
clear;

% We cant do this because as our overall moles as our moles increase the linear velocity must also increase.

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
v0 = u0 / A / 100 %m/sec

A = A/(100^2);

zSpan = linspace(0,L);
[z,F] = ode45(@(z,F) odefun(z, F, T, P, v0, nu, L, A), zSpan, mols0);
figure(1);
plot(z, F, 'Linewidth', 1);
legend('Fch4','Fh2o','Fco', 'Fco2','Fh2','Far');
ylabel('F_i');
xlabel('Length of reactor z');
title('SS molar flowrates at ' + string(T) + 'K');

if T == 622 +273
	Fexp = [0.071787506 0.298008466 0.043658696 0.043658696 0.495104008 0.047782628]; %622
elseif T == 721+273
	Fexp = [0.019206098 0.189598069 0.090570208 0.034350264 0.619497555 0.046777807]; %721
else
	[0 0 0 0 0 0]
end

Fbar = [4*F(end, :); Fexp];
figure(2);
cats = categorical({'Fch4','Fh2o','Fco','Fco2','Fh2','Far'});
cats = reordercats(cats,{'Fch4','Fh2o','Fco','Fco2','Fh2','Far'});
bar(cats,Fbar)
ylabel('F_i');
title('Accuracy of model at ' + string(T) + 'K');
legend('Model','Experimental');

function dFdz = odefun(x, F, T, P, v0, nu, L, A) %#ok<*INUSD>
	
	Keqsmr = (101325/100000)^2 * exp(-26830/T + 30.114);
	Keqwgs = exp(4400/T - 4.036);

	Rmatrix = rxneq(F, Keqsmr, Keqwgs, T, L); 
	
	dFdz = (nu*Rmatrix)*(8.3144598 * 10^-2 * (T))/P;
end

function R = rxneq(F, Keqsmr, Keqwgs, T, L) % units of mol/m3/s

	R = zeros(3,1);

	y = F/sum(F);

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

	Ku = 1;

	I = 9; %curent in Amps
	F = 96485; %faradays constant

	%ch4; h20; co; co2; h2; ar
	R(1,1) = Ku*Asmr*exp(-Easmr*1000/gasconstant/T)*(y(1)*y(2)-P^2*y(3)*y(5)^3/Keqsmr);
	R(2,1) = Ku*Awgs*exp(-Eawgs*1000/gasconstant/T)*(y(3)*y(2)-y(4)*y(5)/Keqwgs);
	R(3,1) = I/(2*F*L); 
end