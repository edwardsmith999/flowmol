%Ouput spectral solution for velocity in unsteady start-up plate flow u(y,t) with plate described
% by U(1-exp(-t/t0)). Input in form
% "couette_analytical_fn(time,
%                        exponential_decay_time_constant,
%                        dynamic viscosity
%                        density
%                        final wall_velocity,
%                        Domain_height,
%                        required_number_of_u_points_in_y)"

function[u]=startup_plate_analytical_fn(t,t0,meu,rho,U0,L,npoints)

nmodes = 1000;
y = 0:L/npoints:L;
y = y';
neu = meu/rho;

u=zeros(npoints+1,1);

%Add constant term
u = u + U0*( + (1-y/L) - (sin((L-y)/sqrt(neu*t0)))/(sin((L)/sqrt(neu*t0)))*exp(-t/t0) );
for n = 1:nmodes
    u(:,1)= u(:,1) + (2/pi)*( (((-1)^n)/n)*((L^2)/(L^2+n^2*pi^2*neu*t0)) ...
                             *(exp((-((n*pi)/L)^2)*neu*t)) ...
                             *(sin(n*pi*(1-y/L))));
end
