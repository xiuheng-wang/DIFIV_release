function [SRI_hat,info] = hysure_adaptor(HSI, MSI, P1, P2, Pm, R, opts)

% Adaptor for the HySure method

data_r = (reshape(HSI,size(HSI,1)*size(HSI,2),size(HSI,3))')';

basis_type = 'VCA'; lambda_phi = 5e-4;
lambda_m = 1; P_hysure = opts.P;
try 
	M0_hysure = vca(data_r','Endmembers',P_hysure,'verbose','off');
catch
    try
        M0_hysure = vca(data_r',P_hysure);
    catch
        error('Hum... seems like VCA is not working!')
    end
end

Yhim = HSI; Ymim = MSI;
downsamp_factor = size(MSI,1)/size(HSI,1); 
decimFactor = downsamp_factor ; 
%[Phi1,~,Phi2] = svd(fspecial('gaussian',4,1));
Phi1 = P1(2,:); Phi1 = Phi1(abs(Phi1)>=1e-5); Phi2 = Phi1;
H_blur = Phi1' * Phi2;
B_est = zeros(size(MSI,1),size(MSI,2));
B_est(1:size(H_blur,1), 1:size(H_blur,2)) = H_blur;

if decimFactor == 4
    shift = 2; % ???
elseif decimFactor == 3
    shift = 2;
elseif decimFactor == 2
    shift = 1;
else
    error('shift value of Hysure must be adjusted!!!')
end

R_est = Pm ./ repmat(sum(Pm')',[1 size(HSI,3)]); %Normalize SRF
SRI_hat = data_fusion_mod(Yhim, Ymim, downsamp_factor, R_est, B_est, P_hysure, basis_type, lambda_phi, lambda_m, shift, M0_hysure);
info = {};
end

