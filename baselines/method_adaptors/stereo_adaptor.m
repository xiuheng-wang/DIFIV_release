function [SRI_hat,info] = stereo_adaptor(HSI, MSI, P1, P2, Pm, R, opts)

% Adaptor for STEREO

[A_0,B_0,C_0,~,~,C0_tilde] = tenRec(MSI,tens2mat(HSI,[],3),25,R,P1,P2);
[A,B,C] = stereo(HSI,MSI,P1,P2,Pm,10,1,A_0,B_0,C_0,C0_tilde);
SRI_hat = real(cpdgen({A,B,C}));
info = {};
end

