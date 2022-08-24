function [SRI_hat, info] = scott2(HSI, MSI, P1, P2, Pm, R, opts)

% SCOTT2 runs the SCOTT algorithm for specified rank R
% [SRI_hat,info] = SCOTT2(HSI, MSI, P1, P2, Pm, ranks, opts) returns 
% estimation of SRI and informative structure
% 
% INPUT ARGUMENTS:
%     MSI, HSI: input datasets (resp. MSI and HSI)
%     P1,P2,Pm: spatial and spectral degratation matrices
%     ranks: specified multilinear rank
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
if ~isfield(opts,'lambda') || isempty(opts.lambda)
    opts.lambda = 1;
end
if ~isfield(opts,'alpha') || isempty(opts.alpha)
    opts.alpha = 0;
end

[U, ~, ~] = svds(tens2mat(MSI,1,[]),R(1));
[V, ~, ~] = svds(tens2mat(MSI,2,[]),R(2));
[W, ~, ~] = svds(tens2mat(HSI,3,[]), R(3));

P1 = sparse(P1); P2 = sparse(P2); Pm = sparse(Pm);
U_tilde = P1*U; V_tilde = P2*V; W_tilde = Pm*W;

if (R(1)>size(HSI,1)||R(2)>size(HSI,2)) && R(3)>size(MSI,3)
    fprintf('Out of the identifiability region ! \n')
    SRI_hat = NaN; info = 'Non-identifiable';
elseif ((R(1)>size(HSI,1)||R(2)>size(HSI,2)) && R(3)<=size(MSI,3))
    A = U_tilde'*U_tilde;
    B = kron(eye(R(3)), V_tilde'*V_tilde);
    C = eye(R(1));
    D = kron(opts.lambda*(W_tilde'*W_tilde), eye(R(2)));
    b_old = tmprod(HSI,{U_tilde', V_tilde', W'},[1,2,3]) + opts.lambda * tmprod(MSI,{U', V', W_tilde'},[1,2,3]);
    E = reshape(b_old, R(1),R(2)*R(3));
    if R(3) > size(MSI, 3)
      S = reshape(bartelsStewart(A,B,C,D,E),R);
    else
     S = reshape(bartelsStewart(C,D,A,B,E),R);
    end 
    SRI_hat = lmlragen({U,V,W},S);
    info.factors = {U,V,W}; info.core = {S}; info.rank = {R};
else
    A = kron(V_tilde'*V_tilde, U_tilde'*U_tilde);
    A = A+ opts.alpha*norm(A,2)^2*eye(size(A,2));
    B = opts.lambda*(W_tilde'*W_tilde);
    B = B+ opts.alpha*norm(A,2)^2*eye(size(B,2));
    b_old = tmprod(HSI,{U_tilde', V_tilde', W'},[1,2,3]) + opts.lambda * tmprod(MSI,{U', V', W_tilde'},[1,2,3]);
    C = reshape(b_old, R(1)*R(2), R(3));
    %S = reshape(bartelsStewart(A,[],[],B,C), R);
    S = reshape(sylvester(A,B,C), R); 
    SRI_hat = lmlragen({U,V,W},S);
    info.factors = {U,V,W}; info.core = {S}; info.rank = {R};
end



end

