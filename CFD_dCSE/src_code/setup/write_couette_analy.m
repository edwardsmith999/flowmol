%Check some values of analytical solution
plot(couette_analytical_fn(0,0.5,1,311,256,'top'))
hold all
plot(couette_analytical_fn(10,0.5,1,311,256,'top'))
plot(couette_analytical_fn(100,0.5,1,311,256,'top'))
plot(couette_analytical_fn(400,0.5,1,311,256,'top'))
plot(178,0,'x','MarkerSize',20)


%Take correct one and write to file as dp binary
u = couette_analytical_fn(400,0.5,1,181.2574503477,130,'top')
plot(u)
fid=fopen('uy_input','w')
for i=1:130
    i
    u(i)
    fwrite(fid,u(i),'float64')
end

