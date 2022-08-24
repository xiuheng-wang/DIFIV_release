function [SRI_hat,info] = glphs_adaptor(HSI, MSI, P1, P2, Pm, R, opts)

% Adaptor for the GLPHS method

hs_glphs = HSI; ms_glphs = MSI;
ratio_glphs = size(MSI,1)/size(HSI,1); mode_glphs = 2;

Z_glphs = MTF_GLP_wrapper(hs_glphs,ms_glphs,ratio_glphs,mode_glphs);
SRI_hat = denoising(Z_glphs);
info = {};

end

