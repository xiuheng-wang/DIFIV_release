function err = compute_metrics(I_HS,I_REF,ratio)

% Computes all quality metrics

blk_size = 32;

% % Q2n_index = q2n(I_REF, I_HS, blk_size, blk_size);
% % if abs(Q2n_index) < 1e-4
% %     Q2n_index = q2n(I_HS, I_REF, blk_size, blk_size);
% % end
% % Out.q2n = Q2n_index;
% [Q2n_index, ~] = q2n(I_REF, I_HS, blk_size, blk_size);
% Out.q2n = Q2n_index;
Out.q2n = 0; 


[angle_SAM,map] = SAM(I_HS,I_REF);
Out.sam = angle_SAM;
Out.ergas = ERGAS(I_HS,I_REF,ratio);
psnr = PSNR(I_REF,I_HS);
Out.psnrall = psnr.all;
Out.sammap = map;
Out.psnr = psnr.ave;




% UIQI - calls the method described in "A Universal Image Quality Index"
% by Zhou Wang and Alan C. Bovik
q_band = zeros(1, size(I_HS,3));
for idx1=1:size(I_HS,3)
    q_band(idx1)=img_qi(I_REF(:,:,idx1), I_HS(:,:,idx1), 32);
end
uiqi_idx = mean(q_band);
Out.uiqi = uiqi_idx;

err = {Out.sam, Out.ergas, Out.psnr, Out.uiqi};

end

