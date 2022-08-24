function [c4,m4,c2] = cum4(X,prewhiten)
%CUM4 Fourth-order cumulant tensor.
%   [c4,m4,c2] = cum4(X) computes the second-order cumulant (covariance
%   matrix) c2, fourth-order moment m4 and fourth-order cumulant
%   (quadricovariance tensor) c4 of a matrix X in which each row is an
%   observation and each column is a variable. Herein,
%
%      c2(i,j)     = E[xi.*conj(xj)]
%      m4(i,j,k,l) = E[xi.*conj(xj).*conj(xk).*xl]
%      c4(i,j,k,l) = E[xi.*conj(xj).*conj(xk).*xl] ...
%                    - E[xi.*conj(xj)]*E[conj(xk).*xl] ...
%                    - E[xi.*conj(xk)]*E[conj(xj).*xl] ...
%                    - E[xi.*xl]*E[conj(xj).*conj(xk)]
%
%   where the expectation E is approximated by the arithmetic mean and xi
%   is the i-th mean centered variable, X(:,i)-mean(X(:,i)) (and
%   analogously for xj, xk and xl).
%
%   [c4,m4,c2] = cum4(X,'prewhiten') applies a linear transformation to the
%   columns of X so that the covariance matrix of the new matrix is the
%   identity matrix before computing its fourth-order cumulant.
%
%   See also cov, scov, cum3.

%   Authors: Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] P. McCullagh, "Tensor Methods in Statistics," Chapman and Hall,
%       London, 1987.
%   [2] C. Nikias, A. Petropulu, "Higher-Order Spectra Analysis: A 
%       Nonlinear Signal Processing Framework," Prentice Hall, 1993.
%
%   Version history:
%   - 2014/12/16  NV  improved speed and memory consumption

% Check the prewhiten option.
if nargin < 2, prewhiten = false; end
if ischar(prewhiten), prewhiten = strcmpi(prewhiten,'prewhiten'); end

% Center the variables.
X = bsxfun(@minus,X,mean(X,1));

% Apply a prewhitening to X if requested.
n = size(X,1);
if prewhiten
    [U,S,~] = svd(X,'econ');
    X = U*(S*pinv(S))*sqrt(n-1);
end

% Compute c2 = E[xi*conj(xj)] and r2 = E[xi*xj].
c2 = conj(X'*X)/n;
r2 = (X.'*X)/n;

% Compute m4 = E[xi*conj(xj)*conj(xk)*xl].
m4 = kr(X', X.');
m4 =  m4*m4';
m4 = reshape(m4, ones(1,4)*size(X, 2))/size(X,1);

% Compute c4(i,j,k,l) = m4 - E[xi.*conj(xj)]*E[conj(xk).*xl] - ...
% E[xi.*conj(xk)]*E[conj(xj).*xl] - E[xi.*xl]*E[conj(xj).*conj(xk)].
c4 = m4-bsxfun(@times,c2,permute(c2,[3 4 2 1])) ...
       -bsxfun(@times,permute(c2,[1 3 2 4]),permute(c2,[3 2 4 1])) ...
       -bsxfun(@times,permute(r2,[1 3 4 2]),permute(conj(r2),[3 1 2 4]));
