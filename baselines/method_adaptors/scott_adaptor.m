function [SRI_hat,info] = scott_adaptor(HSI, MSI, P1, P2, Pm, R, opts)

% Adaptor for SCOTT

[SRI_hat, ~] = scott2(HSI, MSI, P1, P2, Pm, R);

info = {};
end

