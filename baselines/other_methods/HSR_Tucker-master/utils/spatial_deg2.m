function [P1,P2,Phi1,Phi2] = spatial_deg2(DATA, q, d1, d2, sigma)

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
% Phi = gauss_kernel(q,sigma);


[Phi1,SSS,Phi2] = svd(fspecial('gaussian',q,sigma));
%%%% Phi = - sqrt(SSS(1,1)) * Phi1(:,1);
Phi = abs(sqrt(SSS(1,1))) * Phi1(:,1);
Phi = Phi'; 



%% ADDED:
Phi1 = conv(Phi,ones(1,d1)/d1);
Phi2 = conv(Phi,ones(1,d2)/d2);
q1 = q + d1 - 1;
q2 = q + d2 - 1;
%%

P1 = convmtx(Phi1,Im);
tmp = eye(Im,Im);
for i=1:Im
    tmp(i,:) = imfilter(tmp(i,:), Phi1, 'circular', 'same');
end

displ = 0;
tmp = tmp((1+displ):d1:end, :);
P1 = tmp;


P2 = convmtx(Phi2,Jm);
tmp = eye(Jm,Jm);
for i=1:Jm
    tmp(i,:) = imfilter(tmp(i,:), Phi2, 'circular', 'same');
end

displ = 0;
tmp = tmp((1+displ):d2:end, :);
P2 = tmp;




return

%%

end


