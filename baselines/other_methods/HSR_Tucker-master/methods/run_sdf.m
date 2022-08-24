function [SRI_hat,cost, err] = run_sdf(MSI, HSI, SRI, ranks, options, P1,P2,Pm)

% RUN_SDF runs the HOSVD algorithm for specified rank R followed by 
% TensorLab optimization
% [SRI_hat,cost, err] = RUN_SDF(MSI, HSI, SRI, ranks, options, P1,P2,Pm) returns 
% estimation of SRI, value of cost function and metrics in the cell array err
% 
% INPUT ARGUMENTS:
%     SRI, MSI, HSI: input datasets (resp. groundtruth SRI, MSI and HSI
%     ranks: specified multilinear rank
%     P1,P2,Pm: spatial and spectral degratation matrices
%     options: number of iterations, verbosity
% OUTPUT ARGUMENTS:
%     SRI_hat: estimated SRI
%     cost: value of the cost function
%     err: cell array of metrics
% Copyright (c) 2018 Clémence Prévost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker


R = ranks;

[U, ~, ~] = svds(tens2mat(MSI,1,[]),R(1));
[V, ~, ~] = svds(tens2mat(MSI,2,[]),R(2));
[W, ~, ~] = svds(tens2mat(HSI,3,[]), R(3));

lam = 1;
A = kron(V'*(P2'*P2)*V, U'*(P1'*P1)*U);
B = lam* W'*(Pm'*Pm)*W;
b_old = tmprod(HSI,{U'*P1', V'*P2', W'},[1,2,3]) + lam * tmprod(MSI,{U', V', W'*Pm'},[1,2,3]);
C = reshape(b_old, R(1)*R(2), R(3));
S = reshape(sylvester(A,B,C),R);
A = kron(eye(R(3)),kron(V'*(P2'*P2)*V, U'*(P1'*P1)*U)) + lam * kron(W'*(Pm'*Pm)*W, eye(R(1)*R(2)));
b = tmprod(HSI,{U'*P1', V'*P2', W'},[1,2,3]) + lam * tmprod(MSI,{U', V', W'*Pm'},[1,2,3]);
S = reshape(A\ b(:),R);


model = struct;
model.variables = {U,V,W,S};
model.factors.A = 1; model.factors.B = 2; model.factors.C = 3; model.factors.S = 4;
model.factors.D = {1, @(U,task) struct_matvec(U, task, P1, eye(ranks(1)))};
model.factors.E = {2, @(V,task) struct_matvec(V, task, P2, eye(ranks(2)))};
model.factors.F = {3, @(W,task) struct_matvec(W, task, Pm, eye(ranks(3)))};
model.factorizations{1}.data = HSI;
model.factorizations{1}.weight = 2;
model.factorizations{1}.lmlra  = {'D', 'E', 'C','S'};
%model.factorizations{1}.myreg.regL2    = 'S';
%model.factorizations{1}.myreg.weight   = 0.1;
model.factorizations{2}.data = MSI;
model.factorizations{2}.weight = 2;
model.factorizations{2}.lmlra  = {'A', 'B', 'F', 'S'};
sol = sdf_nls(model, 'MaxIter', 30);

U = sol.variables{1}; V = sol.variables{2}; W = sol.variables{3}; S = sol.variables{4};
SRI_hat = lmlragen({U,V,W},S); cost = frob(HSI - lmlragen({P1*U,P2*V,W}, S),'squared') + frob(MSI - lmlragen({U,V,Pm*W},S),'squared');
err = {cost nmse(SRI,SRI_hat), cc(SRI,SRI_hat), SAM(SRI,SRI_hat), ergas(SRI,SRI_hat,1/4), r_snr(SRI,SRI_hat)};
%snr = r_snr(SRI,SRI_hat);
end

