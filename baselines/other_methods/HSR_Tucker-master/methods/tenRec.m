function[A_hat,B_hat,C_hat,A_tilde,B_tilde,C_tilde]=tenRec(MSI,H3,maxit,t_rank,P1,P2) 
% Ten-Rec algorithm
% (c) Charilaos I. Kanatsoulis, University of Minnesota, Jan 7 , 2018
% nikos@umn.edu
% 
% Reference 1: C.I. Kanatsoulis, X. Fu, N.D. Sidiropoulos and W.K. Ma, 
%``Hyperspectral Super-resolution: A Coupled Tensor Factorization
%Approach,'' IEEE Transactions in Signal Processing, 2018

% Reference 2: C.I. Kanatsoulis, X. Fu, N.D. Sidiropoulos and W.K. Ma, 
%``Hyperspectral Super-resolution via Coupled Tensor Factorization:
%Identifiability and Algorithms,'' IEEE International Conference on 
%Acoustics, Speech and Signal Processing (ICASSP), 2018
[U,~]=cpd(MSI,t_rank,'MaxIter',maxit);
A_hat=U{1};
B_hat=U{2};
C_tilde=U{3};
A_tilde=P1*A_hat;
B_tilde=P2*B_hat;
temp=khatri_rao(B_tilde,A_tilde);
C_hat=(temp\H3)';
end


function [ kR ] = khatri_rao( A,B )
%khatri-rao using the matlab function for kronecker
[~, F] = size(A);

kR = [];
for f=1:F
 kR = [kR kron(A(:,f),B(:,f))];
end
end
