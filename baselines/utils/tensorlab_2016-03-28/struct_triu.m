function [x,state] = struct_triu(z,task)
%STRUCT_TRIU Upper triangular matrix.
%   [x,state] = struct_triu(z) generates x as a lower triangular matrix by
%   using the vector z to fill the matrix column by column. For a matrix x
%   of order n, the vector z should have length n*(n+1)/2. The structure
%   state stores information which is reused in computing the right and
%   left Jacobian-vector products.
%
%   struct_triu(z,task) computes the right or left Jacobian-vector
%   product of this transformation, depending on the structure task. Use
%   the structure state and add the field 'r' of the same shape as z or the
%   field 'l' of the same shape as x to obtain the structure task for
%   computing the right and left Jacobian-vector products
%   
%      (dF(:)/dz(:).')*task.r(:) and
%      (dF(:)/dz(:).')'*task.l(:) + conj((dF(:)/dconj(z(:)).')'*task.l(:)),
%   
%   respectively. Here, F(z) represents this transormation, (:) signifies
%   vectorization and the derivative w.r.t. z (conj(z)) is a partial
%   derivative which treats conj(z) (z) as constant. The output has the
%   same shape as x or z for the right and left Jacobian-vector products,
%   respectively.
%   
%   See also struct_band, struct_diag, struct_tridiag, struct_tril.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.

if nargin < 2, task = []; end
n = 0.5*(-1+sqrt(1+8*length(z)));
[x,state] = struct_band(z,task,[n n],[0 n-1]);

end
