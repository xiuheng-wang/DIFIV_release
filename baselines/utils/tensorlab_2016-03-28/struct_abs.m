function [x,state] = struct_abs(z,task,mu)
%STRUCT_ABS Absolute value.
%   [x,state] = struct_abs(z,[],mu) computes x(i) using
%
%      abs(z(i))                           if abs(z(i)) >  mu
%      2*mu*abs(z(i))^2/(abs(z(i))^2+mu^2) if abs(z(i)) <= mu
%
%   as a smooth approximation for abs(z(i)) depending on the parameter mu.
%   If not supplied, mu is equal to 0.01. The structure state stores
%   information which is reused in computing the right and left
%   Jacobian-vector products.
%
%   struct_abs(z,task,mu) computes the right or left Jacobian-vector
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
%   See also struct_nonneg, struct_sigmoid.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.

if nargin < 2, task = []; end
if nargin < 3, mu = 0.01; end
state = [];

if isempty(task) || (isempty(task.l) && isempty(task.r))
    x = z;
    far = abs(z) > mu;
    x(far) = abs(x(far));
    x2 = x(~far).*conj(x(~far));
    x(~far) = 2*mu*x2./(x2+mu^2);
    state.deriv = zeros(size(z));
    state.deriv(far) = conj(z(far))./(2*x(far));
    state.deriv(~far) = 2*mu^3*conj(z(~far))./(x2+mu^2).^2;
elseif ~isempty(task.r)
    if ~isreal(z) || ~isreal(task.r)
        error('struct_abs:nonanalytic',['Nonanalytic objective ' ...
            'functions are currently not supported in sdf_nls, please ' ...
            'use sdf_minf instead.']);
    end
    x = 2*task.deriv.*task.r;
elseif ~isempty(task.l)
    x = 2*real(task.deriv).*task.l;
end

end
