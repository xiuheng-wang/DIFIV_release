% EXPERIMENT 2: SPECTRA %
% Copyright (c) 2018 Clémence Prévost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker


%%
%%%%%%%%%%%
% LOAD INDIAN PINES DATASET % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption
SRI(1,:,:) = []; SRI(:,1,:) = [];
load('Indian_pines_gt.mat'); indian_pines_gt(:,1) = []; indian_pines_gt(1,:) = [];

%%%%%%%%%%%
% LOAD SALINAS DATASET % 
%%%%%%%%%%%%%%%%%%%%%%%%
% SRI = cell2mat(struct2cell(load('Salinas.mat')));
% SRI = crop(SRI,[80,84,size(SRI,3)]);
% SRI(:,:,[108:112 154:167 224]) = []; %Regions of water absorption (Salinas)
%load('Salinas_gt.mat');  %Don't forget to crop :)

%%
R1 = [40,40,6]; 
R2 = [30,30,16]; 
R3 = [24,24,25];

%%%%%%%%
% RANKS FOR SALINAS %
%%%%%%%%%%%%%%%%%%%%%
% R1 = [40,40,6]; 
% R2 = [30,30,6]; 
% R3 = [14,14,15];

nb_materials = max(max(indian_pines_gt)); mat = tens2mat(SRI,3,[]);
%nb_materials = max(max(salinas_gt)); mat = tens2mat(SRI,3,[]);

F = 100;
[A0, B0, ~,~, C0] = tenRec(MSI, HSI, F, P1,P2,Pm);
[~,~,~,SRIh_S] = stereo(1, F, A0,B0,C0, HSI, MSI, P1,P2,Pm, 10);

filename1 = sprintf('data_exp1_%d_%d_%d_IP',R1(1),R1(2),R1(3)); 
filename2 = sprintf('data_exp1_%d_%d_%d_IP',R2(1),R2(2),R2(3));
filename3 = sprintf('data_exp1_%d_%d_%d_IP',R3(1),R3(2),R3(3));
S1 = cell2mat(struct2cell(load(filename1,'SRI_hat'))); 
S2 = cell2mat(struct2cell(load(filename2,'SRI_hat'))); 
S3 = cell2mat(struct2cell(load(filename3,'SRI_hat')));

mat1 = tens2mat(S1,3,[]);
mat2 = tens2mat(S2,3,[]);
mat3 = tens2mat(S3,3,[]);
matS = tens2mat(SRIh_S,3,[]);

eb = [4 7 9 14]; %number of materials you want (you can change)
for n=1:4
   ind = find(reshape(indian_pines_gt,1,[]) == eb(n));
   spec = real(mean(mat(:,ind),2));
   spec1 = real(mean(mat1(:,ind),2));
   spec2 = real(mean(mat2(:,ind),2));
   spec3 = real(mean(mat3(:,ind),2));
   specS = real(mean(matS(:,ind),2));
   
%     figname1 = sprintf('res_exp2_%d_%d_%d_mat%d_IP',R1(1),R1(2),R1(3),eb(n)); 
%     figname2 = sprintf('res_exp2_%d_%d_%d_mat%d_IP',R2(1),R2(2),R2(3),eb(n));
%     figname3 = sprintf('res_exp2_%d_%d_%d_mat%d_IP',R3(1),R3(2),R3(3),eb(n));
    
%     figname1 = sprintf('spec_exp2_%d_%d_%d_mat%d_IP',R1(1),R1(2),R1(3),eb(n)); 
%     figname2 = sprintf('spec_exp2_%d_%d_%d_mat%d_IP',R2(1),R2(2),R2(3),eb(n));
%     figname3 = sprintf('spec_exp2_%d_%d_%d_mat%d_IP',R3(1),R3(2),R3(3),eb(n));
%     
    figure(1)
    subplot(2,2,n)
    plot(spec,'Linewidth',1.1)
    hold on
    plot(spec1,'Linewidth',1.1)
    hold on
    plot(spec2,'Linewidth',1.1)
    hold on
    plot(spec3,'Linewidth',1.1)
    hold on
    plot(specS,'Linewidth',1.1)
    hold off
    legend('Original spectrum','[40 40 6]','[30 30 16]','[24 24 25]','STERERO F=50')
    title(sprintf('Material %d',eb(n)))
%     
%     figure(1)
%     subplot(2,2,n)
%     plot(log10(svd(spec))); xlim([0 10])
%     title(sprintf('Material %d',eb(n)))
%     
    

% plot(spec,'Linewidth',1.1)
% hold on
% plot(spec1,'Linewidth',1.1)
% hold on
% plot(spec2,'Linewidth',1.1)
% hold on 
% plot(spec3,'Linewidth',1.1)
% legend('Original spectrum','[40 40 6]','[30 30 16]','[24 24 25]')
%  saveas(gcf,figname1,'fig')
    
end
saveas(gcf,'res_fig2_allspec','fig')




for n=1:nb_materials
    ind = find(reshape(indian_pines_gt,1,[]) == n);
    spec = real(mean(mat(:,ind),2));
    spec1 = real(mean(mat1(:,ind),2));
    spec2 = real(mean(mat2(:,ind),2));
    spec3 = real(mean(mat3(:,ind),2));
    specS = real(mean(matS(:,ind),2));
    
    figname1 = sprintf('fig_exp2_%d_%d_%d_mat%d_IP',R1(1),R1(2),R1(3),n); 
    figname2 = sprintf('fig_exp2_%d_%d_%d_mat%d_IP',R2(1),R2(2),R2(3),n);
    figname3 = sprintf('fig_exp2_%d_%d_%d_mat%d_IP',R3(1),R3(2),R3(3),n);
    
    mat_spec(n,:) = spec;
    
    figure(1)
    subplot(1,3,1)
    plot(spec1(50:100),'Linewidth',1.2)
    hold on
    plot(spec(50:100),'Linewidth',1.2)
    hold off; title(sprintf('Portions of spectra for material %d, R=[40,40,6]',n))
    subplot(1,3,2)
    plot(spec1(50:100)-spec(50:100),'r.')
    title('Relative error between spectra')
    subplot(1,3,3)
    plot((spec1(50:100)-spec(50:100))/spec(50:100),'b.')
    title('Normalized error between spectra')
    saveas(gcf,figname1,'fig')
%     
%     figure(2)
%     subplot(1,3,1)
%     plot(spec1(50:100),'Linewidth',1.2)
%     hold on
%     plot(spec(50:100),'Linewidth',1.2)
%     hold off; title(sprintf('Portions of spectra for material %d, R=[30,30,16]',n))
%     subplot(1,3,2)
%     plot(spec2(50:100)-spec(50:100),'r.')
%     title('Relative error between spectra')
%     subplot(1,3,3)
%     plot((spec2(50:100)-spec(50:100))/spec(50:100),'b.')
%     title('Normalized error between spectra')
%     saveas(gcf,figname2,'fig')
%     
%     figure(3)
%     subplot(1,3,1)
%     plot(spec1(50:100),'Linewidth',1.2)
%     hold on
%     plot(spec(50:100),'Linewidth',1.2)
%     hold off; title(sprintf('Portions of spectra for material %d, R=[20,30,6]',n))
%     subplot(1,3,2)
%     plot(spec3(50:100)-spec(50:100),'r.')
%     title('Relative error between spectra')
%     subplot(1,3,3)
%     plot((spec3(50:100)-spec(50:100))/spec(50:100),'b.')
%     title('Normalized error between spectra')
%     saveas(gcf,figname3,'fig')
    
    
figure(1)
plot(spec)
hold on
    % full spectra
    
%     figure(1)
%     subplot(1,3,1)
%     plot(spec1,'Linewidth',1.2)
%     hold on
%     plot(spec,'Linewidth',1.2)
%     hold on
%     plot(specS,'Linewidth',1.2)
%     hold off; legend('ori','ho','stereo'); title(sprintf('Portions of spectra for material %d, R=[40,40,6]',n))
%     subplot(1,3,2)
%     plot(spec1-spec,'r.')
%     hold on
%     plot(specS-spec, '.')
%     title('Relative error between spectra')
%     subplot(1,3,3)
%     plot((spec1-spec)/spec,'b.')
%     hold on
%     plot((specS-spec)/spec, '.')
%     title('Normalized error between spectra')
%     saveas(gcf,figname1,'fig')
%     
%     figure(2)
%     subplot(1,3,1)
%     plot(spec1,'Linewidth',1.2)
%     hold on
%     plot(spec,'Linewidth',1.2)
%     hold on
%     plot(specS,'Linewidth',1.2)
%     hold off;legend('ori','ho','stereo'); title(sprintf('Portions of spectra for material %d, R=[30,30,16]',n))
%     subplot(1,3,2)
%     plot(spec2-spec,'r.')
%     hold on
%     plot(specS-spec, '.')
%     title('Relative error between spectra')
%     subplot(1,3,3)
%     plot((spec2-spec)/spec,'b.')
%     hold on
%     plot((specS-spec)/spec, '.')
%     title('Normalized error between spectra')
%     saveas(gcf,figname2,'fig')
%     
%     figure(3)
%     subplot(1,3,1)
%     plot(spec1,'Linewidth',1.2)
%     hold on
%     plot(spec,'Linewidth',1.2) 
%     hold on
%     plot(specS,'Linewidth',1.2)
%     hold off;legend('ori','ho','stereo'); title(sprintf('Portions of spectra for material %d, R=[20,30,6]',n))
%     subplot(1,3,2)
%     plot(spec3-spec,'r.')
%      hold on
%     plot(specS-spec, '.')
%     title('Relative error between spectra')
%     subplot(1,3,3)
%     plot((spec3-spec)/spec,'b.')
%     hold on
%     plot((specS-spec)/spec, '.')
%     title('Normalized error between spectra')
%     saveas(gcf,figname3,'fig')
    
    %close all
end






