%Function to plot sliding vectors in space, with input of form 
% "plot_slide(meshgrid,top,botttom,colour)"

function [] = plot_slide(x,y,z,top,bot,wallslidev,fig_cube)

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
for ixyz = 1:3
    vslidefield(:,:,:,ixyz) = wallslidev(ixyz)*abs(u.*v.*w);
end
quiver3(x,y,z,vslidefield(:,:,:,1) ,vslidefield(:,:,:,2),vslidefield(:,:,:,3));

