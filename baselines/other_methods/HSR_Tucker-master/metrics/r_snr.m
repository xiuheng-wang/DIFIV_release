function snr = r_snr(DATA, est)

% R_SNR computes SNR between original dataset and its estimate
% snr = R_SNR(DATA, est) returns real number
% 
% INPUT ARGUMENTS:
%     DATA: original dataset
%     est: estimate of DATA
% OUTPUT ARGUMENTS:
%     snr: SNR between DATA and est
% 
% SEE ALSO: SAM, ERGAS, NMSE
% Copyright (c) 2018 Clémence Prévost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker


num = 0; denum = 0;
for k=1:size(DATA,3)
    num = num + norm(reshape(DATA(:,:,k),[],1),'fro')^2;
    denum = denum + norm(reshape(est(:,:,k)-DATA(:,:,k),[],1),'fro')^2;
end
snr = 10*log10(num/denum);

end

