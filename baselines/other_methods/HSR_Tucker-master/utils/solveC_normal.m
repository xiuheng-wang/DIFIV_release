function C = solveC_normal(A_tilde, B_tilde, C_tilde, Pm, HSI)

% SOLVEC_NORMAL solves C for TENREC FUNCTION using normal equation
% C = SOLVEC_NORMAL(A_tilde, B_tilde, C_tilde, Pm, HSI) returns C
% from estimated matrices (CPD of MSI)
% 
% INPUT ARGUMENTS:
%     B_tilde, A_tilde, C_tilde: matrices estimated from CPD(MSI)
%     HSI: hyperspectral image
%     Pm: spectral degradation matrix
% OUTPUT ARGUMENTS:
%     C: factor matrix of size KhxF
% 
% SEE ALSO: SOLVEC_QR
% Copyright (c) 2018 Clémence Prévost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker


X = (B_tilde'*B_tilde).*(A_tilde'*A_tilde);
V = Pm'*Pm;
Z = tens2mat(HSI,[],3)'*kr(B_tilde,A_tilde) + Pm'*C_tilde;
C = sylvester(V,X,Z);

end

