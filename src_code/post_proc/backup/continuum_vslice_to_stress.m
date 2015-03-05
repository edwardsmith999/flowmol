function[continuum_stress]=continuum_vslice_to_stress(continuum_velslice,meu,ly,ny)
n = 1;
delta_y = ly/ny;
for i=2:ny+1
    continuum_stress(n,:) = (continuum_velslice(i-1,:)-continuum_velslice(i+1,:))/(2*delta_y);
    n = n+1;
end

continuum_stress = -continuum_stress*meu;