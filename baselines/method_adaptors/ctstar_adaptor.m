function [SRI_hat,info] = ctstar_adaptor(HSI, MSI, P1, P2, Pm, R, opts)

% Adaptor for proposed approach

Rpsi = opts.R2;

[SRI_hat, P3Psi_est] = ct_star(MSI,HSI,R,Rpsi,P1,P2,Pm);
info = {'scaling factors P3 x3 Psi, alg. 1', P3Psi_est};
end

