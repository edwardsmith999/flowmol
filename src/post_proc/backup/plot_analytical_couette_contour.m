points = 100;
xaxis = 0:1/points:1;
m=0;
for i=1:points
    m = m+0.01
    a = couette_analytical_fn(m,1,1,1,points,'bottom');
    %b(1) = 
    z = [a',a',a',a',a',a',a'];
    colormap(copper);
    %brighten(0.5);
    imagesc(a);
    %quiver(zeros(101,1),xaxis,a,zeros(101,1));
    axis tight
    pause(0.5)
end

