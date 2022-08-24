function [x,state] = struct_conj(z,task)
%STRUCT_CONJ Complex conjugate.
%   [x,state] = struct_conj(z) computes x as the complex conjugate of z.
%   The structure state stores information which is reused in computing the
%   right and left Jacobian-vector products.
%
%   struct_conj(z,task) computes the right or left Jacobian-vector product
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
%   See also struct_ctranspose, struct_transpose.

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
    x = conj(z);
elseif ~isempty(task.r)
    if ~isreal(z) || ~isreal(task.r)
        error('struct_conj:nonanalytic',['Nonanalytic objective ' ...
            'functions are currently not supported in sdf_nls, please ' ...
            'use sdf_minf instead.']);
    end
    x = task.r;
elseif ~isempty(task.l)
    x = conj(task.l);
end

end
