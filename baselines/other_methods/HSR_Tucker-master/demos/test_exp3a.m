% EXPERIMENT 3: COMPARISON OF METRICS WITH STEREO %
% Copyright (c) 2018 Clémence Prévost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker


addpath(genpath('../../HSR_Tucker-master'))
addpath(genpath('../../tensorlab_2016-03-28'))

% %Salinas
% SRI = cell2mat(struct2cell(load('Salinas.mat')));
% SRI = crop(SRI,[80,84,size(SRI,3)]);
% SRI(:,:,[108:112 154:167 224]) = []; %Regions of water absorption (Salinas)

% %Indian Pines
% SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
% SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption (Indian pines)
% SRI(1,:,:) = []; SRI(:,1,:) = [];


load('/home/ricardo/Documents/RESEARCH_PARENT/RESEARCH_PARENT_cloud/THINGS_WORKING(ED)/___HSIs/Indian_pines_corrected.mat')
SRI = indian_pines_corrected(1:144,1:144,:)/10000;

Pm = spectral_deg(SRI,1);
MSI = tmprod(SRI,Pm,3);

d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg2(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);


% R1 = [40,40,6]; R2 = [14,14,15]; R3 = [10,15,25];
% R4 = [30,30,6]; R5 = [58,58,6];
R1 = [40,40,6]; R2 = [30,30,16]; R3 = [24,24,25];

% F = 20;
% [A0, B0, ~,~, C0] = tenRec(MSI, HSI, F, P1,P2,Pm);
% [~,~,~,Sa] = stereo(1, F, A0,B0,C0, HSI, MSI, P1,P2,Pm, 10);
% erra = {nmse(SRI,Sa), SAM(SRI,Sa), ergas(SRI,Sa), r_snr(SRI,Sa), cc(SRI,Sa)};
F = 30;
[A0, B0, ~,~, C0] = tenRec(MSI, HSI, F, P1,P2,Pm);
[~,~,~,Sb] = stereo(1, F, A0,B0,C0, SRI, HSI, MSI, P1,P2,Pm, 10);
errb = {nmse(SRI,Sb), sam(SRI,Sb), ergas(SRI,Sb,1/4), r_snr(SRI,Sb), cc(SRI,Sb)};
% F = 50;
% [A0, B0, ~,~, C0] = tenRec(MSI, HSI, F, P1,P2,Pm);
% [~,~,~,S0] = stereo(1, F, A0,B0,C0, HSI, MSI, P1,P2,Pm, 10);
% err0 = {nmse(SRI,S0), SAM(SRI,S0), ergas(SRI,S0), r_snr(SRI,S0), cc(SRI,S0)};
% F = 100;
% [A0, B0, ~,~, C0] = tenRec(MSI, HSI, F, P1,P2,Pm);
% [~,~,~,S1] = stereo(1, F, A0,B0,C0, HSI, MSI, P1,P2,Pm, 10);
% err1 = {nmse(SRI,S1), SAM(SRI,S1), ergas(SRI,S1), r_snr(SRI,S1), cc(SRI,S1)};

[S2,~, err2] = run_hosvd(SRI,MSI,HSI,R1,P1,P2,Pm);
[S3,~, err3] = run_hosvd(SRI,MSI,HSI,R2,P1,P2,Pm);
[S4,~, err4] = run_hosvd(SRI,MSI,HSI,R3,P1,P2,Pm);
% [S5,~, err5] = run_hosvd(SRI,MSI,HSI,R4,P1,P2,Pm);
% [S6,~, err6] = run_hosvd(SRI,MSI,HSI,R5,P1,P2,Pm);



%% MAKE TABLE FROM RESULTS

errH = load('metrics_hs_ip'); %This is obtained from Hysure and it must contain the same metrics

 T2 = [% real([erra{4}(end) erra{1}(end) erra{5}(end) erra{2}(end) erra{3}(end)]);
%     real([errb{4}(end) errb{1}(end) errb{5}(end) errb{2}(end) errb{3}(end)]);
    real([err0{4}(end) err0{1}(end) err0{5}(end) err0{2}(end) err0{3}(end)]);
    real([err1{4}(end) err1{1}(end) err1{5}(end) err1{2}(end) err1{3}(end)]);
    err2{6} err2{2} err2{3} err2{4} err2{5};
    err3{6} err3{2} err3{3} err3{4} err3{5};
    err4{6} err4{2} err4{3} err4{4} err4{5};
%     err5{6} err5{2} err5{3} err5{4} err5{5};
%     err6{6} err6{2} err6{3} err6{4} err6{5};
    cell2mat(errH.err_hysure_ip);
%     cell2mat(errH.err);
%     cell2mat(errC.err)]; 
];
T2 = table(T2);
writetable(T2, 'exp3_table2_ip')

