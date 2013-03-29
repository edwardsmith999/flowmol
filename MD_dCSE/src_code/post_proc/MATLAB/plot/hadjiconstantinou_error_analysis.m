clear all
close all
fclose('all')

resultfile_dir_MD = '/home/es205/codes/coupled/MD_dCSE/src_code/results/'
resultfile_dir = resultfile_dir_MD;
read_header
filename = strcat(resultfile_dir_MD,'/mbins');
mass_bins = read_mbins(filename,resultfile_dir_MD); 
filename = strcat(resultfile_dir_MD,'/vbins');
vel_bins = read_vbins(filename,resultfile_dir_MD); 
filenamek = strcat(resultfile_dir_MD,'/pVA_k');
filenamec = strcat(resultfile_dir_MD,'/pVA_c');

ibin1 = floor(size(vel_bins,1)/2); ibin2 = ibin1;
jbin1 = floor(size(vel_bins,2)/2); jbin2 = jbin1;
kbin1 = floor(size(vel_bins,3)/2); kbin2 = kbin1;

m = 1;
for iter=1:10
    if (iter ==1)
        %Do nothing first time
        disp([num2str(ibin1) ' to ' num2str(ibin2) ' in x and ' ...
              num2str(jbin1) ' to ' num2str(jbin2) ' in y and ' ...
              num2str(kbin1) ' to ' num2str(kbin2) ' in z'           ])
    elseif (mod(iter,2) == 0)
        if (ibin1 == 1 || jbin1 == 1 ||kbin1 == 1)
            continue
        else
            ibin1 = ibin1 - 1;
            jbin1 = jbin1 - 1;
            kbin1 = kbin1 - 1;
            m = m + 1;
            disp([num2str(ibin1) ' to ' num2str(ibin2) ' in x and ' ...
                  num2str(jbin1) ' to ' num2str(jbin2) ' in y and ' ...
                  num2str(kbin1) ' to ' num2str(kbin2) ' in z'           ])
        end
    else
        if (ibin2 == gnbins(1) || jbin2 == gnbins(2) ||kbin2 == gnbins(3))
            continue
        else
            m = m + 1;
            ibin2 = ibin2 + 1;
            jbin2 = jbin2 + 1;
            kbin2 = kbin2 + 1;
            disp([num2str(ibin1) ' to ' num2str(ibin2) ' in x and ' ...
                  num2str(jbin1) ' to ' num2str(jbin2) ' in y and ' ...
                  num2str(kbin1) ' to ' num2str(kbin2) ' in z'           ])
        end
    end

    % One cell - volume = prod(binsize)
    V(m) = ((ibin2- ibin1) * binsize(1)) *((jbin2- jbin1) * binsize(2)) *((kbin2- kbin1) * binsize(3)) ;

    %Mass
    massV = squeeze(mean(mean(mean(mass_bins(ibin1:ibin2,jbin1:jbin2,kbin1:kbin2,:)/(prod(binsize)*Nmass_ave),1),2),3));
    mn_rho(m) = mean(massV);
    sd_rho(m) = std(massV);
    Erho(m) = sd_rho(m)/mn_rho(m);
    %Velocity
    velV = squeeze(vel_bins(ibin1:ibin2,jbin1:jbin2,kbin1:kbin2,1,:))./squeeze(mass_bins(ibin1:ibin2,jbin1:jbin2,kbin1:kbin2,:));
    velV(isnan(velV)) = 0 ;
    velV=squeeze(mean(mean(mean(velV,1),2),3));
    mn_u(m) = mean(velV);
    sd_u(m) = std(velV);
    Eu(m) = sd_u(m)/mn_u(m);
    %Stress
    for n =1:size(vel_bins,5)-2
        [pVA,pVA_k,pVA_c] = read_pVA(n,resultfile_dir_MD,filenamek,filenamec); 
        %pVA_t(n)  = squeeze(pVA(ibin1:ibin2,jbin1:jbin2,kbin1:kbin2,1,1));
        pVA_kt(n) = squeeze(mean(mean(mean(pVA_k(ibin1:ibin2,jbin1:jbin2,kbin1:kbin2,1,1),1),2),3));
        pVA_ct(n) = squeeze(mean(mean(mean(pVA_c(ibin1:ibin2,jbin1:jbin2,kbin1:kbin2,1,1),1),2),3));
    end
    %Kinetic
    mn_Pk(m) = mean(pVA_kt);
    sd_Pk(m) = std(pVA_kt);
    EPk(m)= sd_Pk(m)/mn_Pk(m);
    %Configurational
    mn_Pc(m) = mean(pVA_ct);
    sd_Pc(m) = std(pVA_ct);
    EPc(m)= sd_Pc(m)/mn_Pc(m);

end

%Plot means vs volume
errorbar(V,mn_Pk,sd_Pk,'x-')
hold all
errorbar(V,mn_Pc,sd_Pc,'s-')
errorbar(V,mn_u,sd_u,'o-')
errorbar(V,mn_rho,sd_rho,'^-')
legend('Kinetic pressure','config pressure', 'velocity', 'density')

%Plot errors
figure
plot(V,EPk)
hold all
plot(V,EPc)
plot(V,Eu)
plot(V,Erho)
legend('Kinetic pressure','config pressure', 'velocity', 'density','location','NorthEastOutside')
