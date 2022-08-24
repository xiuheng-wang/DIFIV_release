function [SRI_hat,info] = cnmf_adaptor(HSI, MSI, P1, P2, Pm, R, opts)

R_est = Pm ./ repmat(sum(Pm')',[1 size(HSI,3)]); %Normalize SRF
P_cnmf = opts.P; %40

SRI_hat = CNMF_fusion2(HSI,MSI,R_est,P_cnmf);
info = {};
end