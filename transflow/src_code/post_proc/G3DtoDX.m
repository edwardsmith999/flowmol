function G3DtoDX(GX,GY,GZ,g3D,dxFileName,originx,originy,originz)
%
% V. Sergiievskyi, MPI MIS 2011
% 
%  G3DtoDX(GX,GY,GZ,g3D,dxFileName)
%
% Saves 3D distribution to *.dx format
%
%  GX,GY,GZ  - grids in x,y,z, directions
%  g3D - 3D distribution
%  dxFileName - name of the *.dx file
%
%

f=fopen(dxFileName,'w');

Nx=length(GX);
Ny=length(GY);
Nz=length(GZ);

% Write Header
fprintf(f,'object 1 class gridpositions counts%8.0f%8.0f%8.0f\n',Nx,Ny,Nz);
if (exist('originx') == 0)
    originx = GX(1);
end
if (exist('originy') == 0)
    originy = GY(1);
end
if (exist('originz') == 0)
    originz = GZ(1);
end
fprintf(f,'origin%16g%16g%16g\n',originx,originy,originz);
fprintf(f,'delta%16g 0 0\n',GX(2)-GX(1));
fprintf(f,'delta 0 %16g 0\n',mean(diff(GY(:))));
fprintf(f,'delta 0 0 %16g\n',GZ(2)-GZ(1));
fprintf(f,'object 2 class gridconnections counts%8.0f%8.0f%8.0f\n',Nx,Ny,Nz);
fprintf(f,'object 3 class array type double rank 0 items%8.0f follows\n',Nx*Ny*Nz);

%Check data is a multiple of 3
arraysize = size(g3D(:));
lastcolsize = mod(arraysize(1),3);
end3 = arraysize(1) - lastcolsize;

%Write data
temp = permute(g3D,[3 2 1]);
fprintf(f,'%16E %16E %16E \n',temp(1:end3));
% col=1;
% for j=1:Nx
% 
%     fprintf('%6.2f%%\n',j/Nx*100);
%     for i=1:Ny
%         for k=1:Nz
%    
%         fprintf(f,'%16E',g3D(i,j,k));
%         col=col+1;
%         if col>3
%             fprintf(f,'\n');
%             col=1;
%         end
%             
%         end
%     end
% end
if (lastcolsize == 0)
    %Do nothing
elseif (lastcolsize == 1)
    fprintf(f,'%16E           \n',temp(end3+1:end3+lastcolsize));
elseif (lastcolsize == 2)
    fprintf(f,'%16E %16E      \n',temp(end3+1:end3+lastcolsize));
end
%Write Footer
fprintf(f,'object "Untitled" call field\n');
fclose(f);
