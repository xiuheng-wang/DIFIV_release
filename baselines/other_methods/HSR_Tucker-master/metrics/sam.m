function sam = sam(DATA,est)

% SAM computes SAM between original dataset and its estimate
% sam = SAM(DATA, est) returns real number
% 
% INPUT ARGUMENTS:
%     DATA: original dataset
%     est: estimate of DATA
% OUTPUT ARGUMENTS:
%     sam: SAM between DATA and est
% 
% SEE ALSO: NMSE, ERGAS, R_SNR
% Copyright (c) 2018 Clémence Prévost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker


DATA_3 = tens2mat(DATA,[],3);
est_3 = real(tens2mat(est,[],3));

sam = 0;
for n=1:size(DATA,1)*size(DATA,2)
    sam = sam + acos((DATA_3(n,:)*est_3(n,:)') / (norm(DATA_3(n,:),2)*norm(est_3(n,:),2)));
end
%sam = sam/(size(DATA,1)*size(DATA,2)); % opt
sam = (180*sam)/(pi*(size(DATA,1)*size(DATA,2))); % opt

end

