% EXPERIMENT 4 : TEST PERFORMANCE OF HOSVD W. PANCHROMATIC IMAGE %
% Copyright (c) 2018 Clémence Prévost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker

% /home/ricardo/Documents/RESEARCH_PARENT/RESEARCH_PARENT_cloud/THINGS_WORKING(ED)/___HSIs/Indian_pines_corrected.mat

SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption (Indian pines)
SRI(1,:,:) = []; SRI(:,1,:) = [];

d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
%Ph = kron(P1,P2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);

Km = 1; Pm = zeros(Km,size(SRI,3));
spec_range = linspace(400,2500,size(SRI,3));
Pm(1,:) = 1/length(spec_range);
MSI = tmprod(SRI,Pm,3); %Make panchromatic image

%% 
% R1 = 1:36; R3 = 1:25; 
R1 = 37:50; R3 = 1:6; 
alpha = 0; %No regularization
%options.MaxIter = 30;
for i=1:length(R1)
    for j=1:length(R3)
        R = [R1(i), R1(i), R3(j)];
        filename = sprintf('HSPan_exp4_%d_%d_%d_IP',R1(i),R1(i),R3(j));
        %[SRI_hat,cost, snr] = run_sdf(MSI, HSI, SRI ,R,options,P1,P2,Pm);
        [SRI_hat,cost, err] = run_hosvd(SRI,MSI,HSI,R,P1,P2,Pm, alpha);
        snr = err{6};
        save(filename,'SRI_hat','cost','snr');
    end
end


snr1 = []; cost1 = [];
for i=1:50
    for j=1:25
        [i,j]
        
        if ((i<=144 && j<=6)||(i<=36 && j<=200))&& exist(sprintf('HSPan_exp4_%d_%d_%d_IP.mat',i,i,j),'file')==2
            eval(sprintf('load(''HSPan_exp4_%d_%d_%d_IP'')',i,i,j));
            snr1(i,j) = snr; cost1(i,j) = cost;
        else
            snr1(i,j) = NaN; cost1(i,j) = NaN;
        end
    end
end

R3 = 1:25; R1 = 1:50;
figure(1)
surfc(R3,R1,snr1)
ylabel('R1=R2'); xlabel('R3'); zlabel('SNR(dB)');
title('SNR between SRI and estimate for R1=R2 and R3')
saveas(gcf,'fig_exp4_snr_R2f_IP','fig')
figure(2)
surfc(R3,R1,cost1)
ylabel('R1=R2'); xlabel('R3'); zlabel('Value of cost function');
title('Cost function value between SRI and estimate for R1=R2 and R3')
saveas(gcf,'fig_exp4_cost_R2f_IP','fig')

