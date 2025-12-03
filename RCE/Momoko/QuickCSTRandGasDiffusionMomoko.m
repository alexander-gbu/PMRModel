clc;
clear;

%SO I REALIZED IT IS FINE TO DO IT LIKE THAT BECAUSE ALL YOU DO IS MULTIPLY IT BY THE PARTITION COEFFICIENT

%#ok<*NUSED>
%#ok<*GVMIS>
%#ok<*INUSD>

% SS outlet concentrations
% flowrate CO2
% height
% 60sccm
% V = 360cm3
% 10cm diameter
% 5cm headspace
% 8cm liqiud

L = 1.0;        % length (m)
D = 0.1;       % diffusion coefficient (m^2/s)
a = 1;        % flux coefficient at left boundary
R = 2;        % constant reaction flux at right boundary

% Mesh
x = linspace(0, L, 50);
t = linspace(0, 10, 1000);

% Store parameters
params.D = D;
params.a = a;
params.R = R;

m = 0; % 1-D slab geometry

sol = pdepe(m, ...
            @(x,t,u,DuDx) pdefun(x,t,u,DuDx,params), ...
            @(x) icfun(x,params), ...
            @(xl,ul,xr,ur,t) bcfun(xl,ul,xr,ur,t,params), ...
            x, t);

C = sol(:,:,1);

% ----- PLOTS -----
figure(1);
surf(x, t, C, 'EdgeColor','none');
xlabel('x'); ylabel('t'); zlabel('C');
title('Diffusion with Robin + Flux Boundary');
view(45,30); colorbar;

figure(2);
plot(x, C(end,:), 'LineWidth', 2);
xlabel('x'); ylabel('C');
title('Steady-state concentration profile');
grid on;

figure(3);
plot(t, C(:,1),'LineWidth',2); hold on
plot(t, C(:,end),'--','LineWidth',2);
legend('Outlet concentration C(0,t)', 'Liquid concentration C(L,t)');
xlabel('Time'); ylabel('C');
title('Boundary concentrations vs time');
grid on;

figure(4);
plot(t, C(:,1), 'LineWidth', 2);
xlabel('Time');
ylabel('C(0,t)');
title('Outlet Concentration vs Time');
grid off;

function [c,f,s] = pdefun(x,t,u,DuDx,params)
    D = params.D;

    c = 1;           % u_t coefficient
    f = D * DuDx;    % flux term
    s = 0;           % no source term
end

function u0 = icfun(x, params)
    u0 = 0;   % initial concentration everywhere
end

function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t,params)
    a = params.a;
    R = params.R;

    % LEFT boundary: a*u + f = 0
    pl = -a * ul;
    ql = 1;

    % RIGHT boundary: f + R = 0
    pr = -R;
    qr = 1;
end
