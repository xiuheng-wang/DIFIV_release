function [SRI_hat, P3Psi_est] = hosvdvar(MSI,HSI,Rz,Rpsi,P1,P2,Pm)
% =============================================================
% Computes the image fusion problem with spectral variability 
% with a coupled tensor factorization model, using a closed
% form solution.
%
% Author: Ricardo Borsoi
%
% Code related to the paper: 
%   Coupled tensor decomposition for hyperspectral and multispectral image fusion with inter-image variability
%   R.A. Borsoi, C. Prévost, K. Usevich, D. Brie, J.C.M. Bermudez, C. Richard
%   IEEE Journal of Selected Topics in Signal Processing 15 (3), 702-717, 2021.
% =============================================================

[N1,N2,L]   = size(HSI);
[M1,M2,L_l] = size(MSI);

if N1 < Rz(1)+Rpsi(1) || N2 < Rz(2)+Rpsi(2)
    error('Ranks are too large!')
end

[C1h, ~, ~] = svds(tens2mat(HSI,1,[]), Rz(1));
[C2h, ~, ~] = svds(tens2mat(HSI,2,[]), Rz(2));
[B3h, ~, ~] = svds(tens2mat(HSI,3,[]), Rz(3));

[C1m, ~, ~] = svds(tens2mat(MSI,1,[]), Rz(1)+Rpsi(1));
[C2m, ~, ~] = svds(tens2mat(MSI,2,[]), Rz(2)+Rpsi(2));
% [C3m, ~, ~] = svds(tens2mat(MSI,3,[]), Rz(3)+Rpsi(3));



% minimize \|C1h - P1*C1m*Q\|
Q1 = ((P1*C1m)'*(P1*C1m)) \ ((P1*C1m)'*C1h);
Q2 = ((P2*C2m)'*(P2*C2m)) \ ((P2*C2m)'*C2h);

Bz1 = C1m*Q1;
Bz2 = C2m*Q2;

AtA = kron(B3h'*B3h, kron(Bz2'*P2'*P2*Bz2, Bz1'*P1'*P1*Bz1));

AtYh = tmprod(HSI,{Bz1'*P1', Bz2'*P2', B3h'},[1,2,3]);
Dz = reshape(AtA\ AtYh(:), Rz);
Dz = double(Dz);

SRI_hat = lmlragen({Bz1,Bz2,B3h},Dz);
P3Psi_est = MSI - tmprod(SRI_hat,Pm,3);;


