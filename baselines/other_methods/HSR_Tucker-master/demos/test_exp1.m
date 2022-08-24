% EXPERIMENT 1: HOSVD FOR VARIOUS RANKS %
% Copyright (c) 2018 Clémence Prévost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker


%%
%%%%%%%%%%%
% LOAD INDIAN PINES DATASET % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
% SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption
% SRI(1,:,:) = []; SRI(:,1,:) = [];
% Pm = spectral_deg(SRI);
% MSI = tmprod(SRI,Pm,3);
% d1 = 4; d2 = 4; q = 9;
% [P1,P2] = spatial_deg(SRI, q, d1, d2);
% HSI = tmprod(tmprod(SRI,P1,1),P2,2);

%%%%%%%%%%%
% LOAD SALINAS DATASET % 
%%%%%%%%%%%%%%%%%%%%%%%%
SRI = cell2mat(struct2cell(load('Salinas.mat')));
SRI = crop(SRI,[80,84,size(SRI,3)]);
SRI(:,:,[108:112 154:167 224]) = []; %Regions of water absorption (Salinas)
Pm = spectral_deg(SRI);
MSI = tmprod(SRI,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);

%%

%%%%%
% FIX R1 = R2 %
%%%%%%%%%%%%%%%
R1 = 21:40; R3 = 2:6; 
alpha = 0; %No regularization
options.MaxIter = 30;
for i=1:length(R1)
    for j=1:length(R3)
        R = [R1(i), R1(i), R3(j)];
        filename = sprintf('data_exp1_%d_%d_%d_Sal',R1(i),R1(i),R3(j));
        %[SRI_hat,cost, snr] = run_sdf(MSI, HSI, SRI ,R,options,P1,P2,Pm);
        [SRI_hat,cost, err] = run_hosvd(SRI,MSI,HSI,R,P1,P2,Pm, alpha);
        snr = err{6};
        save(filename,'SRI_hat','cost','snr');
    end
end

%%%%%%%%%%
% FIX R3 TO DESIRED VALUE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %R1 = 10:40; R2 = 10:40; R3 = 6;  %Fix R3
% R1 = 41:50; R2 = 41:50; R3 = 6;
% alpha = 0; %No regularization
% for i=1:length(R1)
%     for j=1:length(R2)
%         R = [R1(i), R2(j), R3];
%         filename = sprintf('data_exp1_%d_%d_%d_IP',R1(i),R2(j),R3);
%         [SRI_hat,cost, snr] = run_hosvd(SRI,MSI,HSI,R,P1,P2,Pm, alpha);
%         save(filename,'SRI_hat','cost','snr');
%     end
% end

%% MAKE FIGURES

%%%%%%%%%%%%%
% BE CAREFUL TO ADJUST CONDITIONS WRT IDENTIFIABILITY AND FILE NAME %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snr1 = []; cost1 = [];
for i=2:40
    for j=2:20
        [i,j]
        
        if ((i<=80 && j<=6)||(i<=20 && j<=204))&& exist(sprintf('data_exp1_%d_%d_%d_Sal.mat',i,i,j),'file')==2
            eval(sprintf('load(''data_exp1_%d_%d_%d_Sal'')',i,i,j));
            snr1(i-1,j-1) = snr; cost1(i-1,j-1) = cost;
        else
            snr1(i-1,j-1) = NaN; cost1(i-1,j-1) = NaN;
        end
    end
end

%R3 = 2:25; R1 = 10:50;
R3 = 2:20; R1 = 2:40;
figure(1)
surfc(R3,R1,snr1)
ylabel('R1=R2'); xlabel('R3'); zlabel('SNR(dB)');
title('SNR between SRI and estimate for R1=R2 and R3')
saveas(gcf,'fig_exp1_snr_R2f_Sal','fig')
figure(2)
surfc(R3,R1,cost1)
ylabel('R1=R2'); xlabel('R3'); zlabel('Value of cost function');
title('Cost function value between SRI and estimate for R1=R2 and R3')
saveas(gcf,'fig_exp1_cost_R2f_Sal','fig')

%%%%%%%%%%%%%
% FOR A GIVEN VALUE OF R3 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% snr2 = []; cost2 = [];
% for i=10:40
%     for j=10:40
%         [i,j]
%         eval(sprintf('load(''data_exp1_%d_%d_%d_IP'')',i,j,6));
%         snr2(i-9,j-9) = snr; cost2(i-9,j-9) = cost;
%     end
% end
% 
% figure(3)
% surfc(10:40,10:40,snr2)
% xlabel('R1'); ylabel('R2'); zlabel('SNR(dB)');
% title('SNR between SRI and estimate for R1 and R2 (R3=6)')
% saveas(gcf,'fig_exp1_snr_R3f6','fig')
% figure(4)
% surfc(10:40,10:40,cost2)
% xlabel('R1'); ylabel('R2'); zlabel('Value of cost function');
% title('Cost function value between SRI and estimate for R1 and R2 (R3=6)')
% saveas(gcf,'fig_exp1_cost_R3f6','fig')
% 
% snr3 = []; cost3 = [];
% for i=10:50
%     for j=10:50
%         [i,j]
%         
%         if (i<=36 && j<=36 && exist(sprintf('data_exp1_%d_%d_%d_IP.mat',i,j,16),'file')==2)
%             eval(sprintf('load(''data_exp1_%d_%d_%d_IP'')',i,j,16));
%             snr3(i-9,j-9) = snr; cost3(i-9,j-9) = cost;
%         else
%             snr3(i-9,j-9) = NaN; cost3(i-9,j-9) = NaN;
%         end
%     end
% end
% 
% figure(5)
% surfc(10:50,10:50,snr3)
% xlabel('R1'); ylabel('R2'); zlabel('SNR(dB)');
% title('SNR between SRI and estimate for R1 and R2 (R3=16)')
% saveas(gcf,'fig_exp1_snr_R3f16','fig')
% figure(6)
% surfc(10:50,10:50,cost3)
% xlabel('R1'); ylabel('R2'); zlabel('Value of cost function');
% title('Cost function value between SRI and estimate for R1 and R2 (R3=16)')
% saveas(gcf,'fig_exp1_cost_R3f16','fig')

