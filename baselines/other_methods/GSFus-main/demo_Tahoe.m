
clear
clc
rng(10, 'twister')

addpath('BM3D');
addpath('GSFus_main');
addpath('Quality_Indices');

load('.\data\Tahoe_preproc_t1.mat')

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
S = denoising(HSI); % reference HSI

[M N L]=size(S);
sf =2; % downsampling ration
sz=[M N];
s0=1;

psf        =    fspecial('gaussian',4,0.8493);
par.fft_B      =    psf2otf(psf,sz);
par.fft_BT     =    conj(par.fft_B);
par.H          =    @(z)H_z(z, par.fft_B, sf, sz,s0 );
par.HT         =    @(y)HT_y(y, par.fft_BT, sf, sz,s0);

R = SRF; % spectral response function
R = R ./ repmat(sum(R')',[1 L]);

S_bar = hyperConvert2D(S);
hyper= par.H(S_bar);
HSI =hyperConvert3D(hyper,M/sf, N/sf );

L   = size(S,3);
L_l = size(MSI,3);
M1  = size(S,1);
M2  = size(S,2);
N1  = M1/sf;
N2  = M2/sf;

% Signal to noise ratio of HS and MS images
SNR_h = 30;
SNR_m = 40;
% add noise
HSI = HSI + sqrt(sum((HSI(:)).^2)/N1/N2/L/10^(SNR_h/10)) * randn(size(HSI));
MSI = MSI + sqrt(sum((MSI(:)).^2)/M1/M2/L_l/10^(SNR_m/10)) * randn(size(MSI));

%% GSFus
subspace=8; % dimention of subspace
lambda_1=0.1; % lambda1
lambda_2=0.005; %lambda2
mu=0.1; % mu

GSFus =  GSFus_main( HSI, MSI,R,par.fft_B,sf,S,subspace,lambda_1,lambda_2,mu);
qualIdx_GSFus    = QualityIndices_mod(GSFus, S, sf)


