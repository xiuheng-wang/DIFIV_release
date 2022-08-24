function [x,state] = struct_diag(z,task)
%STRUCT_DIAG Diagonal matrix.
%   [x,state] = struct_diag(z) generates x as a diagonal matrix with
%   the elements of the vector z on the diagonal. The structure state
%   stores information which is reused in computing the right and left
%   Jacobian-vector products.
%
%   struct_diag(z,task) computes the right or left Jacobian-vector product
%   of this transformation, depending on the structure task. Use the
%   structure state and add the field 'r' of the same shape as z or the 
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
%   See also struct_band, struct_tridiag, struct_tril, struct_triu.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.

if nargin < 2, task = []; end
state = [];

if isempty(task) || (isempty(task.l) && isempty(task.r))
    x = diag(z);
elseif ~isempty(task.r)
    x = diag(task.r);
elseif ~isempty(task.l)
    x = diag(task.l);
end

end
