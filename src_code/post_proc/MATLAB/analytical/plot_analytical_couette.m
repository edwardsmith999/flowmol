xaxis = 0:1/20:1;
Re = 8.0;
L = 1.0;
t = [0.01,0.1,0.25,0.5,1.0,10.0];
hold on
for i=1:5
	a = couette_analytical_fn(t(i),Re,[-1.0,1.0],L,size(xaxis,2),'both');
    plot(xaxis,a,'k','linewidth',4);  
end
ylabel('U_x/U'); xlabel('y/H')

figure
hold on
xaxis = 0:0.01:1;
for i=1:5
    a = couette_analytical_stress_fn(t(i),Re,[-1.0,1.0],L,size(xaxis,2),'both');
    plot(xaxis,a,'r');
end
ylabel('Stress'); xlabel('y/H')


% t = 0.001
% for i=1:100
% t =t*2
% a = couette_analytical_stress_fn(t,5,1,1,100);
% plot(a,'r');
% ylabel('Stress'); xlabel('y/H')
% pause(0.5)
% end
