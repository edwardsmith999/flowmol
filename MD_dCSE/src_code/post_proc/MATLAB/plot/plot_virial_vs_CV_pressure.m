clear all
close all

resultfile_dir = '/home/es205/codes/issued_codes/svn_lucian/MD_dCSE/src_code/results/'

%Read Header file
read_header

%Calculate surface bins
surface_bins
surface = 2 * globaldomain(1)*globaldomain(2) + ...
          2 * globaldomain(2)*globaldomain(3) + ...
          2 * globaldomain(3)*globaldomain(1)


%===============Get stress acting over surface of CV==================
bin = [ 20 2 2];
for i=1:(Nsteps-initialstep)/tplot
    t = initialstep + i*tplot-1;
    [velocity_snapshot_t,velocity_flux_t,pressure_surface_t] = read_vflux(t,resultfile_dir,gnbins,nd);

    %velocity_snapshot(i,:,:,:,:) = velocity_snapshot_t;
    %velocity_flux(i,:,:,:,:,:) = velocity_flux_t;

      %Whole Domain
     pressure_kinetic(i,1) = -mean(mean(velocity_flux_t( 1 , : , : ,1,1))) ; %Bottom xx
     pressure_kinetic(i,2) = -mean(mean(velocity_flux_t( : , 1 , : ,2,2))) ; %Bottom yy
     pressure_kinetic(i,3) = -mean(mean(velocity_flux_t( : , : , 1 ,3,3))) ; %Bottom zz
     pressure_kinetic(i,4) =  mean(mean(velocity_flux_t(end, : , : ,1,4))) ; %  Top xx
     pressure_kinetic(i,5) =  mean(mean(velocity_flux_t( : ,end, : ,2,5))) ; %  Top yy
     pressure_kinetic(i,6) =  mean(mean(velocity_flux_t( : , : ,end,3,6))) ; %  Top zz
 
     pressure_config(i,1) = mean(mean(pressure_surface_t( 1 , : , : ,1,1))) ; %Bottom xx
     pressure_config(i,2) = mean(mean(pressure_surface_t( : , 1 , : ,2,2))) ; %Bottom yy
     pressure_config(i,3) = mean(mean(pressure_surface_t( : , : , 1 ,3,3))) ; %Bottom zz
     pressure_config(i,4) = mean(mean(pressure_surface_t(end, : , : ,1,4))) ; %  Top xx
     pressure_config(i,5) = mean(mean(pressure_surface_t( : ,end, : ,2,5))) ; %  Top yy
     pressure_config(i,6) = mean(mean(pressure_surface_t( : , : ,end,3,6))) ; %  Top zz

%     %Bin only
%     pressure_kinetic(i,1) = -mean(mean(velocity_flux_t(bin(1),  : , : ,1,1))) ; %Bottom xx
%     pressure_kinetic(i,2) = -mean(mean(velocity_flux_t(bin(1), 1 , : ,2,2))) ; %Bottom yy
%     pressure_kinetic(i,3) = -mean(mean(velocity_flux_t(bin(1), : , 1 ,3,3))) ; %Bottom zz
%     pressure_kinetic(i,4) =  mean(mean(velocity_flux_t(bin(1), : , : ,1,4))) ; %  Top xx
%     pressure_kinetic(i,5) =  mean(mean(velocity_flux_t(bin(1),end, : ,2,5))) ; %  Top yy
%     pressure_kinetic(i,6) =  mean(mean(velocity_flux_t(bin(1), : ,end,3,6))) ; %  Top zz
% 
%     pressure_config(i,1) = mean(mean(pressure_surface_t(bin(1), : , : ,1,1))) ; %Bottom xx
%     pressure_config(i,2) = mean(mean(pressure_surface_t(bin(1), 1 , : ,2,2))) ; %Bottom yy
%     pressure_config(i,3) = mean(mean(pressure_surface_t(bin(1), : , 1 ,3,3))) ; %Bottom zz
%     pressure_config(i,4) = mean(mean(pressure_surface_t(bin(1), : , : ,1,4))) ; %  Top xx
%     pressure_config(i,5) = mean(mean(pressure_surface_t(bin(1),end, : ,2,5))) ; %  Top yy
%     pressure_config(i,6) = mean(mean(pressure_surface_t(bin(1), : ,end,3,6))) ; %  Top zz


    %Average of all cell faces
    %pressure_config(i) = mean(mean(mean(pressure_surface_t(:,:,:,1,1)))) + mean(mean(mean(pressure_surface_t(:,:,:,1,4)))) ...
    %                    +mean(mean(mean(pressure_surface_t(:,:,:,2,2)))) + mean(mean(mean(pressure_surface_t(:,:,:,2,5)))) ...
    %                    +mean(mean(mean(pressure_surface_t(:,:,:,3,3)))) + mean(mean(mean(pressure_surface_t(:,:,:,3,6))));

end

fig1 = figure 
fig2 = figure

CV_kinetic = mean(pressure_kinetic(2:end,:),2);
CV_config  = mean(pressure_config(2:end,:),2);
CV_virial = CV_config + CV_kinetic;

figure(fig1)
plot(CV_config,'b','LineWidth',3)
hold on
figure(fig2)
plot(CV_kinetic,'b--','LineWidth',3)
hold on
%===================Read VA virial====================
% read_pVA
% VA=squeeze(mean(mean(pressure_VA(bin(1),:,:,1,1,:),2),3));
% VA_k=squeeze(mean(mean(pressure_VA_k(bin(1),:,:,1,1,:),2),3));
% VA_c=squeeze(mean(mean(pressure_VA_c(bin(1),:,:,1,1,:),2),3));
% figure(fig1)
% plot(VA_c(2:end),'k','LineWidth',3)
% figure(fig2)
% plot(VA_k(2:end),'--k','LineWidth',3)

% 
% figure
% % %Plot Std scatter
% % for i = 1:size(pressure_VA,6)-1
% %     a(i,:) = squeeze(std(squeeze(mean(mean(pressure_VA(:,:,:,1,1,i:i+1),2),3)),1,2));
% % end
% % scatter(0:42,mean(a(:,:),1))
% % 
% % %Plot change in Std against number of time samples in std
% 
% %Plot std dependance on cell size
% for i = 1:43
%     a(i) = squeeze(std(squeeze(mean(mean(mean(pressure_VA(1:i,:,:,1,1,:),2),3),1)),1)');
% end
% plot(a)

%===================Read macroscopic virial====================
read_macroscopic_properties
figure(fig2)
plot(macroscopic_properties(2:end,9)-macroscopic_properties(2:end,8),'r','LineWidth',3)
legend('CV kinetic pressure','kinetic virial')

figure(fig1)
plot(macroscopic_properties(2:end,8),'--r','LineWidth',3)
legend('CV config pressure','config virial')


%legend('CV config pressure','CV kinetic pressure','VA config pressure','VA kinetic pressure','config virial','kinetic virial')

%Copy to results array
result(1)    = globalnp;
result(2:4)  = globaldomain(:);
result(5)    = mean(CV_virial);
%result(6)    = mean(VA);
result(7)    = mean(macroscopic_properties(2:end,9));
result(8)    = std( CV_virial);
%result(9)    = std(VA);
result(10)   = std( macroscopic_properties(2:end,9));
result(11)   = Nsteps;
result(12)   = density;
%Append CSV file
dlmwrite('CV_virial_out.csv',result,'-append');

%Ns = 25;
%for i = 1:size(macroscopic_properties(2:end,9),1)/Ns-1
%    virial_std(i)=std(macroscopic_properties(2+(i-1)*Ns:2+i*Ns,9));
%end
%plot(virial_std)

%Plot Standard Deviations of various parts
figure
clear a
Ns = 100;

% for i = 10:size(pressure_VA,6)-Ns
% 	a(i,:) = std(squeeze(mean(mean(pressure_VA(:,:,:,1,1,i:i+Ns),2),3)),1,2);
% end
% plot(mean(a(1:end,:),2),'r')
% hold on
% for i = 10:size(pressure_VA,6)-Ns
% 	a(i,:) = std(squeeze(mean(mean(pressure_VA(:,:,:,1,2,i:i+Ns),2),3)),1,2);
% end
%plot(mean(a(1:end,:),2),'k')
%Increasing block method used in Rapaport
Ns=1;
clear a; clear b
for n=1:13
    for i = 1000:size(CV_virial,1)-Ns
        CV_config_std(i) = std(CV_config(i:i+Ns),1);
        CV_kinetic_std(i) = std(CV_kinetic(i:i+Ns),1);
        CV_std(i) = std((CV_kinetic(i:i+Ns)+CV_config(i:i+Ns)),1);
    end
    for i = 1000:size(macroscopic_properties(2:end,9),1) - Ns
        %virial_std(i)=std(macroscopic_properties(2+(i-1)*Ns:2+i*Ns,9));
        virial_std(i)=std(macroscopic_properties(i:i+Ns,9));
        config_virial_std(i)=std(macroscopic_properties(i:i+Ns,8));
        kinetic_virial(i) = std(macroscopic_properties(i:i+Ns,9)-macroscopic_properties(i:i+Ns,8));
    end
    
    a(n,1) = mean(CV_config_std);
    a(n,2) = mean(CV_kinetic_std);
    a(n,3) = mean(config_virial_std);
    a(n,4) = mean(kinetic_virial);
    b(n,1) = mean(CV_std);
    b(n,2) = mean(virial_std);
    xaxis(n) = Ns
    Ns = Ns*2;
end
plot(xaxis,a)
legend('CV config','CV kinetic','Virial config','Virial kinetic')

figure
plot(xaxis,b)

