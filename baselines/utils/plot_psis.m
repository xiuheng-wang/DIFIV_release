function []=plot_psis(P3Psi, base_pathname, alg_name)
% P3Psi : Can be either a tensor, or a cell array containing multiple tensors
%         If it is a single tensor, 
% alg_name

% example: base_pathname = 'examples/ex4a_KernRiver'
% ex_pathname{1} = [base_pathname '_P3Psi_vis.pdf'];
% ex_pathname{2} = [base_pathname '_P3Psi_norms.pdf'];
ex_pathname     = [base_pathname '_P3Psi_vis_norms.pdf'];
ex_pathname_fig = [base_pathname '_P3Psi_vis_norms.fig'];

% get the number of algorithms' results to plot:
N_algs = length(P3Psi);

% if the algorithms' names are not supplied, leave it blank:
if ~exist('alg_name','var')
    alg_name = cell(N_algs,1);
    for i=1:N_algs
        alg_name{i} = {};
    end
end

figure;
[ha, pos] = tight_subplot(N_algs, 3, [0.01 0.01], 0.25, 0.15);

if N_algs == 1
    fontsizeNumber = 14;
else
    fontsizeNumber = 12;
end


for i=1:N_algs
    ddata = P3Psi{i};

    % rgbim(:,:,1) = imadjust(rescale(ddata(:,:,1),1));
    % rgbim(:,:,2) = imadjust(rescale(ddata(:,:,2),1));
    % rgbim(:,:,3) = imadjust(rescale(ddata(:,:,3),1));
    rgbim(:,:,1) = abs(ddata(:,:,1));
    rgbim(:,:,2) = abs(ddata(:,:,2));
    rgbim(:,:,3) = abs(ddata(:,:,3));
    rgbim = rgbim / max(rgbim(:));
    %axes(ha(1)); imshow(rgbim); % display RGB image
    axes(ha(3*(i-1)+1)); imagesc(rgbim); % display RGB image
    axis image
    set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
    ylabel(alg_name{i},'interpreter','latex','fontsize',fontsizeNumber)
    if i==1
        title('Visible','interpreter','latex','fontsize',fontsizeNumber)
    end


    % rgbim(:,:,1) = imadjust(rescale(ddata(:,:,7),1));
    % rgbim(:,:,2) = imadjust(rescale(ddata(:,:,9),1));
    % rgbim(:,:,3) = imadjust(rescale(ddata(:,:,10),1));
    rgbim(:,:,1) = abs(ddata(:,:,7));
    rgbim(:,:,2) = abs(ddata(:,:,9));
    rgbim(:,:,3) = abs(ddata(:,:,10));
    rgbim = rgbim / max(rgbim(:));
    %axes(ha(2)); imshow(rgbim); % display RGB image
    axes(ha(3*(i-1)+2)); imagesc(rgbim); % display RGB image
    axis image
    set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
    if i==1
        title('Infrared','interpreter','latex','fontsize',fontsizeNumber)
    end


    axes(ha(3*(i-1)+3)); imagesc(sum(ddata.^2,3).^0.5), %axis square
    axis image
    set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
    if i==1
        title('Magnitude','interpreter','latex','fontsize',fontsizeNumber)
    end
    originalSize2 = get(gca, 'Position');
    h=colorbar; 
    set(ha(3*(i-1)+3), 'Position', originalSize2);
    set(h,'fontsize',fontsizeNumber-2);

end


savefig(ex_pathname_fig)
print(ex_pathname,'-dpdf','-painters')


% -----------------------------------------------------------------------------
% try to crop images
try 
    [~,~] = system(['pdfcrop ' ex_pathname ' ' ex_pathname]);
catch me
    disp('Warning: could not crop a PDF image file.')
end





%{

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&






if ~iscell(P3Psi)
    figure;
    [ha, pos] = tight_subplot(1, 3, [0.01 0.01], 0.25, 0.15);

    ddata = P3Psi;

    % rgbim(:,:,1) = imadjust(rescale(ddata(:,:,1),1));
    % rgbim(:,:,2) = imadjust(rescale(ddata(:,:,2),1));
    % rgbim(:,:,3) = imadjust(rescale(ddata(:,:,3),1));
    rgbim(:,:,1) = abs(ddata(:,:,1));
    rgbim(:,:,2) = abs(ddata(:,:,2));
    rgbim(:,:,3) = abs(ddata(:,:,3));
    rgbim = rgbim / max(rgbim(:));
    %axes(ha(1)); imshow(rgbim); % display RGB image
    axes(ha(1)); imagesc(rgbim); % display RGB image
    axis image
    set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
    title('Visible','interpreter','latex','fontsize',14)
    ylabel(alg_name,'interpreter','latex','fontsize',14)

    % rgbim(:,:,1) = imadjust(rescale(ddata(:,:,7),1));
    % rgbim(:,:,2) = imadjust(rescale(ddata(:,:,9),1));
    % rgbim(:,:,3) = imadjust(rescale(ddata(:,:,10),1));
    rgbim(:,:,1) = abs(ddata(:,:,7));
    rgbim(:,:,2) = abs(ddata(:,:,9));
    rgbim(:,:,3) = abs(ddata(:,:,10));
    rgbim = rgbim / max(rgbim(:));
    %axes(ha(2)); imshow(rgbim); % display RGB image
    axes(ha(2)); imagesc(rgbim); % display RGB image
    axis image
    set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
    title('Infrared','interpreter','latex','fontsize',14)



    axes(ha(3)); imagesc(sum(ddata.^2,3).^0.5), %axis square
    axis image
    set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])
    title('Magnitude','interpreter','latex','fontsize',14)
    originalSize2 = get(gca, 'Position');
    h=colorbar; 
    set(ha(3), 'Position', originalSize2);
    set(h,'fontsize',12);




else


end



%}




