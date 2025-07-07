% Y is going to be the molar fractions. z is going to be axial direction, x is going to be "radial direction"

% this doesnt even work because one of the dimensions has to be time. it cant be z
% the boundary conditions can be met


function [c,f,s] = pdefun(z,x,Y,dYdz)




    c = zeros(length(y),1) + 1;
	f = zeros(length(y),1);
	s = -v*dydx + nu*R;

end