function [SRI_hat,info] = bscott(MSI,HSI,Pm,R,opts)

% BSCOTT runs the BSCOTT algorithm for specified rank R
% [SRI_hat,info] = BSCOTT(MSI,HSI,R,Pm,opts) returns 
% estimation of SRI and info structure
% 
% INPUT ARGUMENTS:
%     MSI, HSI: input datasets (resp. MSI and HSI)
%     R: specified multilinear rank
%     Pm: spectral degradation matrix
%     opts: options structure
% OUTPUT ARGUMENTS:
%     SRI_hat: estimated SRI
%     info: informative structure
%
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

if ~exist('opts','var')
    opts = struct();
end
if ~isfield(opts,'Nblocks') || isempty(opts.Nblocks)
    opts.Nblocks = [1,1];
end

range_MSI = [size(MSI,1),size(MSI,2)]; 
range_HSI = [size(HSI,1),size(HSI,2)];
step_MSI = ceil(range_MSI ./ opts.Nblocks); 
step_HSI = ceil(range_HSI ./ opts.Nblocks);
SRI_hat = zeros(size(MSI,1), size(MSI,2), size(HSI,3));

for i1=1:opts.Nblocks(1)
  for i2=1:opts.Nblocks(2)
      
    M_ind_min = [i1-1,i2-1].*step_MSI + 1;
    M_ind_max = min([i1,i2].*step_MSI, range_MSI);
    H_ind_min = [i1-1,i2-1].*step_HSI + 1;
    H_ind_max = min([i1,i2].*step_HSI, range_HSI);
    
    [fact,S] = mlsvd(MSI(M_ind_min(1):M_ind_max(1),...
                        M_ind_min(2):M_ind_max(2), :), R); % hosvd
     U = cell2mat(fact(1)); V = cell2mat(fact(2)); W_tilde = cell2mat(fact(3));    
    [W_H, ~, ~] = svds(tens2mat(HSI(H_ind_min(1):H_ind_max(1), ...
                                H_ind_min(2):H_ind_max(2), :),3,[]), R(3));
    T = (Pm * W_H) \ W_tilde;

    W = W_H * T;
    SRI_hat(M_ind_min(1):M_ind_max(1),...
            M_ind_min(2):M_ind_max(2), :) = lmlragen({U,V,W},S); 
  end
end

info.factors = {U,V,W};
info.core = {S};
info.rank = {R};

end

