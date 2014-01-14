%Caluclate central difference from a set of data (including halo cells)
% returning only the domain data (array size reduced by two)

function[dydx]=cdiff(y,dx)

%Central Differencing
dydx = zeros((size(y(:),1)-2),1);

for i =2:size(y(:))-1
    dydx(i-1) = (y(i+1)-y(i-1))/(2*dx);
end 

end
