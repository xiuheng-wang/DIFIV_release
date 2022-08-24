function [P1,P2] = spatial_deg_uniform(DATA, d1, d2)

% SRI2HSI computes P1 and P2 from SRI 
% [P1, P2] = SRI2HSI(DATA,d1,d2) returns spatial degradation matrices from SRI DATA
% 
% INPUT ARGUMENTS:
%     DATA: SRI in ImxJmxKh
%     d1: downsampling factor in mode 1
%     d2: downsampling factor in mode 2
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


% [Phi1,SSS,Phi2] = svd(fspecial('gaussian',q,sigma));
%%%% Phi = - sqrt(SSS(1,1)) * Phi1(:,1);

Phi1 = ones(d1,1)'/d1;
Phi2 = ones(d2,1)'/d2;

%%

% tmp = eye(Im,Im);
% tmp = imfilter(tmp, Phi1, 'circular', 'full');
% displ = 2; 0;
% P1 = tmp((1+displ):d1:(Im+displ), :);
% 
% tmp = eye(Jm,Jm);
% tmp = imfilter(tmp, Phi2, 'circular', 'same');
% displ = 2; 0;
% P2 = tmp((1+displ):d2:(Jm+displ), :);


tmp = eye(Im,Im);
for i=1:Im
    tmp(i,:) = real(ifft(fft(tmp(i,:)) .* fft(Phi1,Im)));
end
% displ = 1; 0;
displ = 0;
P1 = tmp((1+displ):d1:(Im+displ), :);

tmp = eye(Jm,Jm);
for i=1:Jm
    tmp(i,:) = real(ifft(fft(tmp(i,:)) .* fft(Phi2,Jm)));
end
% displ = 1; 
displ = 0;
P2 = tmp((1+displ):d2:(Jm+displ), :);


% tmp = eye(Im,Im);
% tmp2 = zeros(Im,Im+d1-1);
% for i=1:Im
%     tmp2(i,:) = imfilter(tmp(i,:), Phi1, 0, 'full');
% end
% % circulantize:
% tmp2(:,1:(d1-1)) = tmp2(:,1:(d1-1)) + tmp2(:,(Im+1):(Im+d1-1));
% tmp2 = tmp2(:,1:Im);
% 
% displ = 0;
% % displ = floor(d1/2);
% tmp2 = tmp2((1+displ):d1:end, :);
% P1 = tmp2;
% 
% 
% % P2 = convmtx(Phi2,Jm);
% tmp = eye(Jm,Jm);
% tmp2 = eye(Jm,Jm+d2-1);
% for i=1:Jm
%     tmp2(i,:) = imfilter(tmp(i,:), Phi2, 0, 'full');
% end
% % circulantize:
% tmp2(:,1:(d2-1)) = tmp2(:,1:(d2-1)) + tmp2(:,(Jm+1):(Jm+d2-1));
% tmp2 = tmp2(:,1:Jm);
% 
% displ = 0;
% % displ = floor(d2/2);
% tmp2 = tmp2((1+displ):d2:end, :);
% P2 = tmp2;
% 
% 


return

%%

end


