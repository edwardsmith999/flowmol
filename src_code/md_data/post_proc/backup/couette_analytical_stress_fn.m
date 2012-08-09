function[tau]=couette_analytical_stress_fn(t,Re,U_wall,L,npoints)

nmodes = 10000;
k = 1/Re;
y = 0:L/npoints:L;
y = y';
C = 0;

tau=zeros(npoints+1,1);

%Add zero wavenumber
tau(:,1) = tau(:,1) + (1/2)*(U_wall*2./L) ;

for n = 1:nmodes
    lambda = (n*pi/L)^2;
     %un =   -(-1)^n*(2*U_wall/(n*pi))*(1-exp(-lambda*k*t)) ;
     %tau(:,1) = tau(:,1) + (2./L*((-1)^n*U_wall) + n*pi/L*un) .* cos(n*pi*y/L);
     tau(:,1) = tau(:,1) + (2./L*((-1)^n*U_wall)*(exp(-lambda*k*t))) .* cos(n*pi*y/L);

end

%calculate wall shear
%tau_wall = U_wall/L;

%tau = tau_wall+tau;
%tau=fliplr(tau')';
tau = -tau;
%Set top value to adjacent
%tau(1,1) = tau_wall;
