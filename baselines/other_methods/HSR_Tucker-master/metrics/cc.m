function cc = cc(DATA,est)

% CC computes cross-corellation between original dataset and its estimate
% cc = cc(DATA, est) returns number between 0 and 1
% 
% INPUT ARGUMENTS:
%     DATA: original dataset
%     est: estimate of DATA
% OUTPUT ARGUMENTS:
%     cc: cross-correlation between DATA and est
% 
% SEE ALSO: NMSE, ERGAS, R_SNR, SAM
% Copyright (c) 2018 Clémence Prévost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker


cc = 0;
for k=1:size(DATA,3)
    cc = cc + corr(reshape(DATA(:,:,k),[],1),reshape(est(:,:,k),[],1),'Type','Pearson');
end
cc = cc/size(DATA,3);

end

