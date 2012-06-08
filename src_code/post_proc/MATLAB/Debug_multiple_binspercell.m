
clear all

resultfile_dir='/home/es205/codes/coupled/MD_dCSE/src_code/results';
read_header
bpc = gnbins./globalncells;

%Test types integer flag
% mass_bins = 1
% mass_flux = 2

test_type = 2;
switch(test_type)
case(1) 
    filename='./mbins';
    [mass_bins,density_bins]=read_mbins(filename,resultfile_dir);


    %sliceomatic(mass_bins)

    mass_bins_sum = zeros(size(mass_bins,1)/bpc(1),size(mass_bins,2)/bpc(2),size(mass_bins,3)/bpc(3));
    l = 1; m = 1; n = 1; 
    for i=1:bpc(1):size(mass_bins,1)
    for j=1:bpc(2):size(mass_bins,2)
    for k=1:bpc(3):size(mass_bins,3)
        for p = 0:bpc(1)-1
        for q = 0:bpc(2)-1
        for r = 0:bpc(3)-1
            mass_bins_sum(l,m,n) = mass_bins_sum(l,m,n) + mass_bins(i+p,j+q,k+r);
        end
        end
        end
        n = n + 1;
    end
        n = 1;
        m = m + 1;
    end
        m = 1;
        l = l + 1;
    end
    sliceomatic(mass_bins_sum)
case(2)

    filename1 = './mflux';
    filename2 = './msnap';
    [mass_flux,mass_snapshot,Nmflux_records]=read_mflux(filename1,filename2,resultfile_dir);
    read_header

    mass_flux_sum = zeros(size(mass_flux,1)/bpc(1),size(mass_flux,2)/bpc(2),size(mass_flux,3)/bpc(3),6);
    l = 1; m = 1; n = 1; 
    for i=1:bpc(1):size(mass_flux,1)
    for j=1:bpc(2):size(mass_flux,2)
    for k=1:bpc(3):size(mass_flux,3)
        for p = 0:bpc(1)-1
        for q = 0:bpc(2)-1
        for r = 0:bpc(3)-1
            mass_flux_sum(l,m,n,:) = mass_flux_sum(l,m,n,:) + mass_flux(i+p,j+q,k+r,:);
        end
        end
        end
        n = n + 1;
    end
        n = 1;
        m = m + 1;
    end
        m = 1;
        l = l + 1;
    end
    sliceomatic(sum(mass_flux_sum(:,:,:,:),4))

end