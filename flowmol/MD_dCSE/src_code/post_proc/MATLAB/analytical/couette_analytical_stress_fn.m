% Calculate the analytical expression for stress in couette flow

% Input of the form:
% tau=couette_analytical_stress_fn(t,Re,U_wall,L,npoints)
% t      - time
% Re     - Reynolds number (Distance units as fraction of domain height L)
% U_wall - wall velocity
% L      - domain height
% npoints- number of points required (size of returned array)


function[tau]=couette_analytical_stress_fn(t,Re,U_wall,L,npoints,slidingwall)

nmodes = 1000;
k = 1/Re;
y = linspace(0,L,npoints);
y = y';
C = 0;

%Preallocate array
tau=zeros(npoints,1);

% - - - Calculate strain - - -
switch slidingwall
    case 'top'
        %Add zero wavenumber
        tau(:,1) = tau(:,1) + (1/2)*(U_wall(1)*2./L) ;
        %Add time specific modes
        for n = 1:nmodes
            lambda = (n*pi/L)^2;
            %un =   -(-1)^n*(2*U_wall/(n*pi))*(1-exp(-lambda*k*t)) ;
            %tau(:,1) = tau(:,1) + (2./L*((-1)^n*U_wall) + n*pi/L*un) .* cos(n*pi*y/L);
            tau(:,1) = tau(:,1) + (2./L*((-1)^n*U_wall(1))*(exp(-lambda*k*t))) .* cos(n*pi*y/L);
        end
    case 'bottom'
        %Add zero wavenumber
        tau(:,1) = tau(:,1) + (1/2)*(U_wall(1)*2./L) ;
        %Add time specific modes
        for n = 1:nmodes
            lambda = (n*pi/L)^2;
            tau(:,1) = tau(:,1) + (-1)^n*(2./L*((-1)^n*U_wall(1))*(exp(-lambda*k*t))) .* cos(n*pi*y/L);
        end
    case 'both'
        %Add zero wavenumber
        tau(:,1) = tau(:,1) + (1/2)*((U_wall(2)-U_wall(1))*2./L) ;
        %Add time specific modes
        for n = 1:nmodes
            lambda = (n*pi/L)^2;
            tau(:,1) = tau(:,1) - (2./L*((-1)^n*U_wall(1))*(exp(-lambda*k*t))) .* cos(n*pi*y/L) ...
                         + (-1)^n*(2./L*((-1)^n*U_wall(2))*(exp(-lambda*k*t))) .* cos(n*pi*y/L);
        end 
    otherwise
        warning('Couette Stress Analytical Solution Needs either top or bottom wall to be sliding');
end

% Shear stress is strain times viscosity (NOTE THIS ASSUMES UNIT DENSITY)
tau = k * tau ;

%calculate wall shear
%tau_wall = U_wall/L;

%tau = tau_wall+tau;
%tau=fliplr(tau')';
%tau = -tau;
%Set top value to adjacent
%tau(1,1) = tau_wall;
