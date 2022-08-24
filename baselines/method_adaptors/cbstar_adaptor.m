function [SRI_hat,info] = cbstar_adaptor(HSI, MSI, P1, P2, Pm, R, opts)

% Adaptor for proposed approach using block coordinate descent

Rpsi = opts.R2;
initOpt = opts.initOpt;

[SRI_hat, P3Psi_est] = cb_star(MSI,HSI,P1,P2,Pm,R,Rpsi,initOpt);
info = {'scaling factors P3 x3 Psi, alg. 2', P3Psi_est};
end

