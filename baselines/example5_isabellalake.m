% =========================================================================
% Example 1 of Super resulution for HS and MS fusion with spectral
% variability
% 
% 
% 
% Kern river sumbim
% 
% =========================================================================

% addpath(genpath('DATA'))
addpath(genpath('utils'))
addpath(genpath('HSR_Tucker-master'))
addpath(genpath('other_methods'))
addpath(genpath('ALGS'))
addpath(genpath('method_adaptors'))

warning off;
clear all 
clc

clus = gcp('nocreate'); % If no pool, do not create new one.
if isempty(clus)
    c = parcluster('local');
    c.NumWorkers = 5;
    parpool(c, c.NumWorkers);
end


rng(10, 'twister') 


load('../data/Isabella_lake_preproc_subim1.mat')
R = SRF;


%% Lets crop a smaller piece of this image
HSI = HSI(9:88, 40:119, :);
MSI = MSI(9:88, 40:119, :);


MSim = MSI;
HSim = HSI;

MSp = zeros(size(MSim));
for i=1:size(MSim,3)
    x = MSim(:,:,i);
    xmax  = quantile(x(:),0.999);
    % Normalize to 1
    x = x/xmax;
    MSp(:,:,i) = x;
end
MSI = MSp;

HSp = zeros(size(HSim));
for i=1:size(HSim,3)
    x = HSim(:,:,i);
    xmax  = quantile(x(:),0.999);
    %normalize to 1
    x = x/xmax;
    HSp(:,:,i) = x;
end
HSI = HSp;

%%


% Signal to noise ratio of HS and MS images
SNR_h = 35;
SNR_m = 35;


% Decimation factor:
decimFactor = 4; % 2;

L   = size(HSI,3);
L_l = size(MSI,3);
M1  = size(HSI,1);
M2  = size(HSI,2);
N1  = M1/decimFactor;
N2  = M2/decimFactor;


% Get MSI and HR HSI
Ym = MSI;
% Zh = HSI;
Zh = denoising(HSI);
Zh_th = Zh;


d1 = decimFactor; d2 = decimFactor; qq = 8;
[P1,P2,Phi1,Phi2] = spatial_deg2(Zh_th, qq, d1, d2, 4);
% [P1,P2] = spatial_deg_uniform(Zh_th, d1, d2);
Pm = R;


Yh = tmprod(tmprod(Zh_th,P1,1),P2,2);


% add noise
Yh = Yh + sqrt(sum((Yh(:)).^2)/N1/N2/L/10^(SNR_h/10)) * randn(size(Yh));
Ym = Ym + sqrt(sum((Ym(:)).^2)/M1/M2/L_l/10^(SNR_m/10)) * randn(size(Ym));






%% run the methods

P = 30;

%{
methods = {'HySure','hysure_adaptor','[]',P;
           'CNMF','cnmf_adaptor','[]',30;
           'GLPHS','glphs_adaptor','[]',[];
           'FuVar','fuvar_adaptor','[]',P;
           'STEREO','stereo_adaptor','40',[];
           'SCOTT','scott_adaptor','[60 60 5]',[]; 
           'HOSVD-Var','hosvdvar_adaptor','[11 11 7]',[3,3,2];
           'HOSVD-VarOpt','hosvdvaropt_adaptor','[50 50 5]',{[30,30,3],1};
           };
%}

methods = {'HySure','hysure_adaptor','[]',P;
           'CNMF','cnmf_adaptor','[]',30;
           'GLPHS','glphs_adaptor','[]',[];
           'FuVar','fuvar_adaptor','[]',P;
           %'LTMR','ltmr_adaptor','[]',{20,1e-3};
           'STEREO','stereo_adaptor','50',[];
           'SCOTT','scott_adaptor','[60 60 5]',[];
           % 'CT-STAR','ctstar_adaptor','[30 30 8]',[3,3,2];
           %'CB-STAR','cbstar_adaptor','[50 50 5]',{[40,40,3],1};
           %...'CB-STAR','cbstar_adaptor','[50 50 5]',{[40,40,3],1};
           %'CB-STAR','cbstar_adaptor','[30 30 5]',{[40,40,3],1};
           'CB-STAR','cbstar_adaptor','[30 30 6]',{[40,40,3],1};
           };

DegMat = struct('Pm', Pm, 'P1', P1, 'P2', P2);
[res, est, optOut] = compare_methods(Zh_th, Yh, Ym, DegMat, decimFactor, methods);

%     'CB-STAR R = [30 3…'    [ 3.2631]    [    3.4064]    [26.5557]    [0.8643]    [69.5158]



%%
% Run GSFus (https://github.com/FxyPd/GSFus)
subspace=8; % dimention of subspace
lambda_1=0.1; % lambda1
lambda_2=0.005; %lambda2
mu=0.1; % mu
[Z_GSFus,ttime]=adaptor_GSFus(Yh,Ym, P1, P2, R./repmat(sum(R,2),[1,L]), decimFactor, subspace,lambda_1,lambda_2,mu);



%% load results from python
% error('stopping before plots.')
snr = SNR_h;
file_name = 'Isabella_lake_preproc_subim1';
Dir = 'PAR_code/results/';
fpath = fullfile(Dir, strcat(file_name, num2str(snr), 'dB.mat'));
load(fpath);

Dir = '../results/';
fpath = fullfile(Dir, strcat(file_name, num2str(snr), 'dB.mat'));
load(fpath);

%% display the results
res2 = cell(size(methods,1)+1+1+2, 5);
res2(1,:) = {'Method', 'SAM', 'ERGAS', 'PSNR', 'UIQI'};
res2(2:end-3,:) = res(:, 1:5);
res2(2:end-3,1) = {'HySure', 'CNMF', 'GLPHS', 'FuVar', 'STEREO', 'SCOTT', 'CB-STAR'};
res2(end-2,:) = ['GSFus', compute_metrics(Z_GSFus,Zh_th,decimFactor)];
res2(end-1,:) = ['PAR', compute_metrics(data,Zh_th,decimFactor)];
res2(end,:) = ['DIFIV', compute_metrics(sr_h,Zh_th,decimFactor)];
res2 = res2([1 2 3 4 6 7 10 5 9 8 11], :);
res2
est{end+1} = Z_GSFus;
est{end+1} = data;
est{end+1} = sr_h;

% % Rotate the images 90 degrees to get more compact plots for this image
% Zh_th = fliplr(rot90(Zh_th,-1));
% for i=1:length(est)
%     est{i} = fliplr(rot90(est{i},-1));
% end

%% Plot stuff
alg_names = {'CNMF', 'GSFus', 'CB-STAR','FuVar', 'DIFIV'};
est = est([2,8,7,4,10]);
plot_results(Zh_th, est, alg_names, 'results_examples/ex5_isabellalake');




save_results('results_examples/ex5_isabellalake.txt', res2)









