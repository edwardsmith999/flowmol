%Function to plot cube in space, with input of form 
% "plot_cube(meshgrid,top,botttom,colour)"

function [] = plot_cube(x,y,z,top,bot,colour,fig_cube)

if (top == bot) 
    return
end

%Define domain
u = fheaviside(x - top(1)) - fheaviside(x - bot(1));
v = fheaviside(y - top(2)) - fheaviside(y - bot(2));
w = fheaviside(z - top(3)) - fheaviside(z - bot(3));

%Plot Surface of cube
set(0,'currentfigure',fig_cube);
hold on
p = patch(isosurface(x,y,z,u.*v.*w));
isonormals(x,y,z,u.*v.*w,p)
set(p,'FaceColor',colour,'EdgeColor','none');
camlight ; alpha(0.5);
view(3); grid on;
xlabel('x'); ylabel('y'); zlabel('z'); 
lighting phong

