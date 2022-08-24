function []=plot_results(HR_HSI, rec_imgs, alg_names, base_pathname)



% example: base_pathname = 'examples/ex4a_KernRiver'
ex_pathname_fig = [base_pathname '_recIm_vis_IR.fig'];
ex_pathname{1}  = [base_pathname '_recIm_vis_IR.pdf'];
ex_pathname{2}  = [base_pathname '_psnr_bandwise.pdf'];
ex_pathname{3}  = [base_pathname '_sam_pixelwise.pdf'];



if length(rec_imgs) ~= length(alg_names)
    error('The number of images and number of names should be the same!')
end

imgs = cell(length(rec_imgs)+1, 1);
imgs(1:end-1) = rec_imgs;
imgs{end} = HR_HSI;

nnames = cell(length(alg_names)+1, 1);
nnames(1:end-1) = alg_names;
nnames{end} = 'Reference';

N_imgs = length(imgs);





% -----------------------------------------------------------------------------
figure;
% [ha, pos] = tight_subplot(2, N_imgs, [0.01 0.01], 0.25, 0.15);
[ha, pos] = tight_subplot(2, N_imgs, [0.01 0.01], 0.10, 0.06);


bands_vis = [32 20 8]; % RGB bands --> visible spectra OK!
bands_ir  = [146 108 46]; % RGB bands --> IR OK!

for i=1:N_imgs
    ddata = imgs{i}; 
    rgbim(:,:,1) = imadjust(rescale(ddata(:,:,bands_vis(1)),1));
    rgbim(:,:,2) = imadjust(rescale(ddata(:,:,bands_vis(2)),1));
    rgbim(:,:,3) = imadjust(rescale(ddata(:,:,bands_vis(3)),1));
    axes(ha(i)); imshow(rgbim); % display RGB image
    set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
    title(nnames{i},'interpreter','latex')

    rgbim(:,:,1) = imadjust(rescale(ddata(:,:,bands_ir(1)),1));
    rgbim(:,:,2) = imadjust(rescale(ddata(:,:,bands_ir(2)),1));
    rgbim(:,:,3) = imadjust(rescale(ddata(:,:,bands_ir(3)),1));
    axes(ha(N_imgs+i)); imshow(rgbim); % display RGB image
    set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])

end

%axes(ha(1));        ylabel('Visible','interpreter','latex')
%axes(ha(N_imgs+1)); ylabel('Infrared','interpreter','latex')
savefig(ex_pathname_fig)
print(ex_pathname{1},'-dpdf')










% -----------------------------------------------------------------------------
figure; hold on
for i=1:(N_imgs-1)
    ppsnr = PSNR(HR_HSI, rec_imgs{i});
    plot(ppsnr.all)
end
xlim([1 length(ppsnr.all)])
legend(alg_names,'interpreter','latex')
xlabel('Spectral band','interpreter','latex')
ylabel('PSNR [dB]','interpreter','latex')
print(ex_pathname{2},'-dpdf')





% -----------------------------------------------------------------------------
figure;
[ha, pos] = tight_subplot(1, N_imgs-1, 0.01, 0.1, [0.1 0.15]);

mmax = -inf;
mmin = inf;
sammap = cell(N_imgs-1,1);
for i=1:(N_imgs-1)
    [~,sammap{i}] = SAM(rec_imgs{i}, HR_HSI);
    mmax = max([mmax, max(sammap{i}(:))]);
    mmin = min([mmin, min(sammap{i}(:))]);
end

for i=1:(N_imgs-1)
    axes(ha(i)); imagesc(sammap{i}, [mmin mmax])
    %set(gca,'ytick',[],'xtick',[]), axis square
    set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
    title(alg_names{i},'interpreter','latex'), 
end

originalSize2 = get(gca, 'Position');
h=colorbar; 
set(ha(N_imgs-1), 'Position', originalSize2);
set(h,'fontsize',10); %12
print(ex_pathname{3},'-dpdf')




% -----------------------------------------------------------------------------
% try to crop images
for i=1:3
    try 
        [~,~] = system(['pdfcrop ' ex_pathname{i} ' ' ex_pathname{i}]);
    catch me
        disp('Warning: could not crop a PDF image file.')
    end
end








