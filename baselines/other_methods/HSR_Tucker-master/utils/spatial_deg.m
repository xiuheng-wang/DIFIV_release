function [P1,P2] = spatial_deg(DATA, q, d1, d2, sigma)

% SRI2HSI computes P1 and P2 from SRI 
% [P1, P2] = SRI2HSI(DATA,d1,d2) returns spatial degradation matrices from SRI DATA
% 
% INPUT ARGUMENTS:
%     DATA: SRI in ImxJmxKh
%     q: size of Gaussian kernel
%     d1: downsampling factor in mode 1
%     d2: downsampling factor in mode 2
%     sigma: Gaussian kernel bandwidth
% OUTPUT ARGUMENTS:
%     P1: degradation along mode 1 of size IhxIm
%     P2: degradation along mode 2 of size JhxJm
% 
% SEE ALSO: SPECTRAL_DEG
% Copyright (c) 2018 Clémence Prévost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker

if nargin < 5
    sigma = 1;
end

Im = size(DATA,1); Jm = size(DATA,2);
Ih = ceil(Im/d1); Jh = ceil(Jm/d2);
Phi = gauss_kernel(q,sigma);

%% ADDED:
Phi1 = conv(Phi,ones(1,d1));
Phi2 = conv(Phi,ones(1,d2));
q1 = q + d1 - 1;
q2 = q + d2 - 1;
%%

r = [Phi1, zeros(1,Im-q1)]; c = [Phi1(1), zeros(1,Im-1)]; TIm = toeplitz(c,r);
r = [Phi2, zeros(1,Jm-q2)]; c = [Phi2(1), zeros(1,Jm-1)]; TJm = toeplitz(c,r);

TIm = convmtx(Phi1,Im);
TIm = TIm(:,d1:end-d1+1);
TJm = convmtx(Phi2,Jm);
TJm = TJm(:,d2:end-d2+1);

%%


S1 = zeros(Ih, size(DATA,1)); S2 = zeros(Jh, size(DATA,2));
for i=1:Ih
    %if 1+d1*(i-1)<= size(S1,2)
        S1(i,1+d1*(i-1)) = 1;
    %end
end
for j=1:Jh
    %if 1+d1*(j-1)<= size(S2,2)
        S2(j,1+d1*(j-1)) = 1;
    %end
end

P1 = S1*TIm; P2 = S2*TJm;



end


