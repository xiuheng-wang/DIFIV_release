function [GSFus, ttime]=adaptor_GSFus(HSI,MSI, P1, P2, R, downsampf, subspace,lambda_1,lambda_2,mu)

% tic;
time_start = clock;

%[Phi1,~,Phi2] = svd(fspecial('gaussian',4,1));
Phi1 = P1(2,:); Phi1 = Phi1(abs(Phi1)>=1e-5); 
Phi2 = P2(2,:); Phi2 = Phi2(abs(Phi2)>=1e-5); 
psf = Phi1' * Phi2;



[M, N, ~]=size(MSI);
L = size(HSI,3);

sz=[M N];

S = 0.5*ones(M,N,L); % reference HRI

s0=1;
%
% psf        =    fspecial('gaussian',4,0.8493);
par.fft_B      =    psf2otf(psf,sz);
par.fft_BT     =    conj(par.fft_B);
par.H          =    @(z)H_z(z, par.fft_B, sf, sz,s0 );
par.HT         =    @(y)HT_y(y, par.fft_BT, sf, sz,s0);


GSFus =  GSFus_main( HSI, MSI, R, par.fft_B, downsampf, S, subspace,lambda_1,lambda_2,mu);

time_finish = clock; % [year month day hour minute seconds]
ttime = (time_finish(3:end) - time_start(3:end)) * [86400;3600;60;1];

