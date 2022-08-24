function [SRI_hat,info] = fuvar_adaptor(HSI, MSI, P1, P2, Pm, R, opts)

% Adaptor for the FuVar method

%Init
P = opts.P;
% Run VCA on HSI and use FCLS/SCLS for initialization
data_r = (reshape(HSI,size(HSI,1)*size(HSI,2),size(HSI,3))')';
try 
    M0 = vca(data_r','Endmembers',P,'verbose','off');
catch
    M0 = vca(data_r',P);
end

%FuVar
A_FCLSU = FCLSU(data_r',M0)';
A_init = reshape(A_FCLSU',size(HSI,1),size(HSI,2),P);
A_init = imresize(A_init, size(MSI,1)/size(HSI,1));
A_init = reshape(A_init,size(MSI,1)*size(MSI,2),P)';
decimFactor = size(MSI,1)/size(HSI,1);

lambda_m = 1; lambda_a = 1e-4; lambda_1 = 0.01; lambda_2 = 10000;
R_srr = Pm ./ repmat(sum(Pm')',[1 size(HSI,3)]);

Psi_init = ones(size(HSI,3),P);
%[Phi1,~,Phi2] = svd(fspecial('gaussian',4,1));
Phi1 = P1(2,:); Phi1 = Phi1(abs(Phi1)>=1e-5); 
Phi2 = P2(2,:); Phi2 = Phi2(abs(Phi2)>=1e-5); 
H_blur = Phi1' * Phi2;

[Zh,~,~,~]=FuVar(HSI,MSI,A_init,M0,Psi_init,R_srr,decimFactor,H_blur,...
    lambda_m,lambda_a,lambda_1,lambda_2);
SRI_hat = reshape(Zh',size(MSI,1),size(MSI,2),size(HSI,3));
%Zm_srr = reshape(Zm',M1,M2,L);
info = {};
end

