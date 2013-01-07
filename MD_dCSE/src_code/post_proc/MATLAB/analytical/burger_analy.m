%Analytical solution to the 1D burgers equation
%Input of the form => burger_analytical(ur,ul,epsilon,x,t,shift)
% ur - speed of inflow
% ul - speed of outflow
% epsilon - wave thickness
% x - domain points [vector]
% t - time
% shift - starting point of wave (default zero)
% soln - Choice of solution (not specified - basic form, 1 - include error function)

function [u,dudt] = burger_analy(ur,ul,epsilon,x,t,shift,soln)

% ur =  2;
% ul = -2;
% epsilon = 0.1;
%x = 0:0.1:domain;

%Add shift to domain
x = x - shift;
s = 0.5*(ul+ur);

%Function h to give other solution of burgers equation
if (exist(soln))
    h = erfc(-(x-ur*t)/(4*epsilon*t))./erfc((x-ul*t)/(4*epsilon*t));
else
    h = 1;
end

%Calculate velocity
f = ((ul-ur)*(x-s*t))/(2*epsilon);
u = ur + (ul-ur)./(1+h*exp(f));

%Calculate derivative
dudt= -u.*(exp(f)./(exp(f)+1))*((ul-ur)/(2*epsilon))*(1-s);

%dudt= ((1./exp(((ul - ur)*(t*ul - 2*x + t*ur))/(4*epsilon)) - 1)*(ul - ur)^3)...
%      ./(4*epsilon*exp(((ul - ur).*(t*ul - 2*x + t*ur))./(4*epsilon))...
%         .*(1./exp(((ul - ur)*(t*ul - 2*x + t*ur))./(4*epsilon)) + 1).^3);

end