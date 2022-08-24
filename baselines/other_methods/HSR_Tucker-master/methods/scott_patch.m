function [SRI_hat, info] = scott_patch(HSI, MSI, P1, P2, Pm, R, opts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

q=1;

if ~exist('opts','var')
    opts = struct();
end
if ~isfield(opts,'lambda') || isempty(opts.lambda)
    opts.lambda = 1;
end
if ~isfield(opts,'alpha') || isempty(opts.alpha)
    opts.alpha = 0;
end
if ~isfield(opts,'Nblocks') || isempty(opts.Nblocks)
    opts.Nblocks = [1,1];
end

range_MSI = [size(MSI,1),size(MSI,2)]; 
range_HSI = [size(HSI,1),size(HSI,2)];
step_MSI = ceil(range_MSI ./ opts.Nblocks); 
step_HSI = ceil(range_HSI ./ opts.Nblocks);
SRI_hat = zeros(size(MSI,1), size(MSI,2), size(HSI,3));

[W, ~, ~] = svds(tens2mat(HSI,3,[]), R(3)); 
P1 = sparse(P1); P2 = sparse(P2); Pm = sparse(Pm);
W_tilde = Pm*W;

for i1=1:opts.Nblocks(1)
  for i2=1:opts.Nblocks(2)
    
      %[i1 i2]
      
    %Update steps
    M_ind_min = [i1-1,i2-1].*step_MSI + 1;
    M_ind_max = min([i1,i2].*step_MSI, range_MSI);
    H_ind_min = [i1-1,i2-1].*step_HSI + 1;
    H_ind_max = min([i1,i2].*step_HSI, range_HSI);
    %Range depending on blocks
    if i1==1
        ind_MSI{1} = M_ind_min(1):M_ind_max(1)+(q-1)/2;
        ind_HSI{1} = H_ind_min(1):H_ind_max(1)+(q-1)/2;
    elseif i1==opts.Nblocks(1)    
        ind_MSI{1} = M_ind_min(1)-(q-1)/2:M_ind_max(1);
        ind_HSI{1} = H_ind_min(1)-(q-1)/2:H_ind_max(1);
    else
        ind_MSI{1} = M_ind_min(1)-(q-1)/2:M_ind_max(1)+(q-1)/2;
        ind_HSI{1} = H_ind_min(1)-(q-1)/2:H_ind_max(1)+(q-1)/2;
    end
    if i2==1
        ind_MSI{2} = M_ind_min(2):M_ind_max(2)+(q-1)/2;
        ind_HSI{2} = H_ind_min(2):H_ind_max(2)+(q-1)/2;
    elseif i2==opts.Nblocks(2)
        ind_MSI{2} = M_ind_min(2)-(q-1)/2:M_ind_max(2);
        ind_HSI{2} = H_ind_min(2)-(q-1)/2:H_ind_max(2);
    else
        ind_MSI{2} = M_ind_min(2)-(q-1)/2:M_ind_max(2)+(q-1)/2;
        ind_HSI{2} = H_ind_min(2)-(q-1)/2:H_ind_max(2)+(q-1)/2;
    end
    
    %Compute factor matrices
    [U, ~, ~] = svds(tens2mat(MSI(ind_MSI{1}, ind_MSI{2}, :),1,[]),R(1));
    [V, ~, ~] = svds(tens2mat(MSI(ind_MSI{1}, ind_MSI{2}, :),2,[]),R(2));          
    U_tilde = P1(ind_HSI{1},ind_MSI{1})*U;
    V_tilde = P2(ind_HSI{2},ind_MSI{2})*V;
    
    if ((R(1)>step_HSI(1)||R(2)>step_HSI(2)) && R(3)>size(MSI,3)) || (R(3)>size(MSI,3))
        fprintf("Out of the identifiability region !")
        SRI_hat = NaN; info = "Non-identifiable";
    elseif ((R(1)>step_HSI(1)||R(2)>step_HSI(2)) && R(3)<=size(MSI,3))
        A = U_tilde'*U_tilde;
        B = kron(eye(R(3)), V_tilde'*V_tilde);
        C = eye(R(1));
        D = kron(opts.lambda*(W_tilde'*W_tilde), eye(R(2)));
        b_old = tmprod(HSI(ind_HSI{1}, ind_HSI{2}, :),{U_tilde', V_tilde', W'},[1,2,3]) + ...
            opts.lambda * tmprod(MSI(ind_MSI{1}, ind_MSI{2}, :),{U', V', W_tilde'},[1,2,3]);
        E = reshape(b_old, R(1),R(2)*R(3));
        if R(3) > size(MSI, 3)
          S = reshape(bartelsStewart(A,B,C,D,E),R);
        else
         S = reshape(bartelsStewart(C,D,A,B,E),R);
        end 
        SRI_hat(ind_MSI{1}, ind_MSI{2}, :) = lmlragen({U,V,W},S);
        info = "Ok";
    else
        A = kron(V_tilde'*V_tilde, U_tilde'*U_tilde);
        A = A+ opts.alpha*norm(A,2)^2*eye(size(A,2));
        B = opts.lambda*(W_tilde'*W_tilde);
        B = B+ opts.alpha*norm(A,2)^2*eye(size(B,2));
        b_old = tmprod(HSI(ind_HSI{1}, ind_HSI{2}, :),{U_tilde', V_tilde', W'},[1,2,3]) + ...
            opts.lambda * tmprod(MSI(ind_MSI{1}, ind_MSI{2}, :),{U', V', W_tilde'},[1,2,3]);
        C = reshape(b_old, R(1)*R(2), R(3));
        S = reshape(sylvester(A,B,C), R); 
        SRI_hat(ind_MSI{1}, ind_MSI{2}, :) = lmlragen({U,V,W},S);
        info = "Ok";
    end
    
    
    
    
  end
end

%SRI_hat = real(SRI_hat);


end

