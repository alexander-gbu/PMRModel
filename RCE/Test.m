clc;
clear;

%#ok<*NUSED>
%#ok<*GVMIS>
%#ok<*INUSD>

x = linspace(0,1,25);
t = linspace(0,1,25);
m = 1;
D = 1;

sol = pdepe(m,@(x,t,u,dudx) heatcyl(x,t,u,dudx,D),@heatic,@heatbc,x,t);
u = sol(:,:,1);

figure(1);
surf(x,t,u)
xlabel("x")
ylabel("t")
zlabel("u(x,t)")
view([150 25])

figure(2);
plot(t,sol(:,1))
xlabel("Time")
ylabel("Temperature u(0,t)")
title("Temperature change at center of disc")

function [c,f,s] = heatcyl(x,t,u,dudx,D)
c = 1;
f = D*dudx;
s = 0;
end
%----------------------------------------------
function u0 = heatic(x)
n = 2.404825557695773;
u0 = besselj(0,n*x);
end
%----------------------------------------------
function [pl,ql,pr,qr] = heatbc(xl,ul,xr,ur,t)
n = 2.404825557695773;
pl = 0; %ignored by solver since m=1
ql = 0; %ignored by solver since m=1
pr = ur-besselj(0,n)*exp(-n^2*t);
qr = 0; 
end