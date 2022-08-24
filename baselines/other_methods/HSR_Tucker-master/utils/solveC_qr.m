function C = solveC_qr(B_tilde, A_tilde, C_tilde, Kh, HSI, F, Pm)

% SOLVEC_QR solves C for TENREC FUNCTION using QR factorization
% C = SOLVEC_QR(B_tilde, A_tilde, C_tilde, Kh, HSI, F, Pm) returns C
% from estimated matrices (CPD of MSI)
% 
% INPUT ARGUMENTS:
%     B_tilde, A_tilde, C_tilde: matrices estimated from CPD(MSI)
%     Kh: number of bands in HSI 
%     HSI: hyperspectral image
%     F: tensor rank for CPD
%     Pm: spectral degradation matrix
% OUTPUT ARGUMENTS:
%     C: factor matrix of size KhxF
% 
% SEE ALSO: SOLVEC_NORMAL
% Copyright (c) 2018 Clémence Prévost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker


[Q,R] = qr(kr(B_tilde,A_tilde),0);
b = [reshape(tens2mat(HSI,[],3)'*Q,[],1); reshape(C_tilde,[],1)];
X = [kron(R,eye(Kh)); kron(eye(F),Pm)];
C = reshape(b\X,[Kh,F]);

end

