
function Phi = gauss_kernel(q,sigma)

% GAUSS_KERNEL returns Gaussian blurring kernel along one dimension
% Phi = GAUSS_KERNEL(q,sigma) computes kernel Phi of size 1xq with parameter sigma
% 
% INPUT ARGUMENTS:
%     q : size of kernel 
%     sigma: normal distribution parameter (default: sigma = 0.5)
% OUTPUT ARGUMENTS:
%     Phi: Gaussian kernel of size 1xq

% Copyright (c) 2018 Clémence Prévost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker

if nargin==1
    sigma = 0.5;
end

Phi = zeros(1,q);
for m=1:q
    Phi(m) = (1/sqrt(2*pi*sigma^2))*exp(-(m-ceil(q/2))^2 / 2);
end

end

