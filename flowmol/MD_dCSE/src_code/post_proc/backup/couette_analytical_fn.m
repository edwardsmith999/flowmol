%Ouput spectral solution for velocity in unsteady Couette flow u(y,t). Input in form
% "couette_analytical_fn(time,
%                        Reynolds_number,
%                        wall_velocity,
%                        Domain_height, 
%                        required_number_of_u_points_in_y
%                        Sliding_Wall ['top' or 'bottom']  )"

function[u]=couette_analytical_fn(t,Re,U_wall,L,npoints,slidingwall)

    nmodes = 1000;
    k = 1/Re;
    y = 0:L/npoints:L;
    y = y'; 
    C = 0;  %Initial condition

    u=zeros(npoints+1,1);

    switch slidingwall
        case 'top' 
            % Uwall at top
            for n = 1:nmodes
                lambda = (n*pi/L)^2;
                u(:,1)=u(:,1)+(C*exp(-lambda*k*t) - (-1)^n*(2*1/(n*pi))*U_wall(1)*(1-exp(-lambda*k*t))).*sin(n*pi*y/L);
            end
            u(end,1) = U_wall(1);
        case 'bottom'
            %Uwall at bottom 
            for n = 1:nmodes
                lambda = (n*pi/L)^2;
                u(:,1)=u(:,1)+(C*exp(-lambda*k*t)     +    (2*1/(n*pi))*U_wall(1)*(1-exp(-lambda*k*t))).*sin(n*pi*y/L);
            end
            u(1,1) = U_wall(1);
        case 'both'
            %Uwall at bottom 
            for n = 1:nmodes
                lambda = (n*pi/L)^2;
                u(:,1)=u(:,1)+(C*exp(-lambda*k*t) - (-1)^n*(2*1/(n*pi))*U_wall(1)*(1-exp(-lambda*k*t))).*sin(n*pi*y/L)...
                             +(C*exp(-lambda*k*t)     +    (2*1/(n*pi))*U_wall(2)*(1-exp(-lambda*k*t))).*sin(n*pi*y/L);
            end
            u(end,1) = U_wall(1);
            u(1,1) = U_wall(2);
        otherwise
            warning('Couette Analytical Solution Needs either top or bottom wall to be sliding');
    end
end

