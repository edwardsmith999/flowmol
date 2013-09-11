hold on
xaxis = 0:1/20:1;

a = couette_analytical_fn(0.01,1,1,1,size(xaxis,2),'top');
plot(xaxis,a,'k','linewidth',4);
a = couette_analytical_fn(0.1,1,1,1,size(xaxis,2),'top');
plot(xaxis,a,'k','linewidth',4);
a = couette_analytical_fn(0.25,1,1,1,size(xaxis,2),'top');
plot(xaxis,a,'k','linewidth',4);
a = couette_analytical_fn(0.5,1,1,1,size(xaxis,2),'top');
plot(xaxis,a,'k','linewidth',4);
a = couette_analytical_fn(1,1,1,1,size(xaxis,2),'top');
plot(xaxis,a,'k','linewidth',4);
a = couette_analytical_fn(10,1,1,1,size(xaxis,2),'top');
plot(xaxis,a,'k','linewidth',4);
ylabel('U_x/U'); xlabel('y/H')

figure
hold on
xaxis = 0:0.01:1;
a = couette_analytical_stress_fn(0.01,5,1,1,size(xaxis,2));
plot(xaxis,a,'r');
a = couette_analytical_stress_fn(0.1,5,1,1,size(xaxis,2));
plot(xaxis,a,'r');
a = couette_analytical_stress_fn(0.25,5,1,1,size(xaxis,2));
plot(xaxis,a,'r');
a = couette_analytical_stress_fn(0.5,5,1,1,size(xaxis,2));
plot(xaxis,a,'r');
a = couette_analytical_stress_fn(1,5,1,1,size(xaxis,2));
plot(xaxis,a,'r');
a = couette_analytical_stress_fn(10,5,1,1,size(xaxis,2));
plot(xaxis,a,'r');
ylabel('Stress'); xlabel('y/H')

% t = 0.001
% for i=1:100
% t =t*2
% a = couette_analytical_stress_fn(t,5,1,1,100);
% plot(a,'r');
% ylabel('Stress'); xlabel('y/H')
% pause(0.5)
% end
