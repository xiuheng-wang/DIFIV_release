function [SRI_hat,info] = ltmr_adaptor(HSI, MSI, P1, P2, Pm, R, opts)

% R : spectral_rank
% lambda


%[Phi1,~,Phi2] = svd(fspecial('gaussian',4,1));
Phi1 = P1(2,:); Phi1 = Phi1(Phi1~=0); Phi2 = Phi1;
H_blur = Phi1' * Phi2;
% B_est = zeros(size(MSI,1),size(MSI,2));
% B_est(1:size(H_blur,1), 1:size(H_blur,2)) = H_blur;


F = Pm;
sf = size(MSI,1)/size(HSI,1); 
sz=[size(MSI,1) size(MSI,2)]; %[M N];
s0=1;
psf            =    H_blur; %fspecial('gaussian',7,2);
par.fft_B      =    psf2otf(psf,sz);
par.fft_BT     =    conj(par.fft_B);
par.H          =    @(z)H_z(z, par.fft_B, sf, sz,s0 );
par.HT         =    @(y)HT_y(y, par.fft_BT, sf, sz,s0);


% S -- HR image;
% 
% S_bar = hyperConvert2D(S);
% hyper= par.H(S_bar);
% MSI = hyperConvert3D((F*S_bar), M, N);
%   HSI =hyperConvert3D(hyper,M/sf, N/sf );
S = zeros(size(MSI,1), size(MSI,2), size(HSI,3));
  

para.K=200;
% para.eta=1e-3;
% para.p=10;
para.eta= opts.lambda; %R(2);
para.p= opts.p; %R(1);

para.patchsize=7;
SRI_hat = TSVD_Subpace_FUS(HSI,MSI,F, par.fft_B,sf,S,para);
info = {};
end
