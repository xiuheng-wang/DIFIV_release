clear all;
close all;

addpath(genpath('../HyperspectralCode/'))
cd ~/Dropbox/matlab/tensorlab_2016-03-28
load ../DataSets/Cuprite_Rita/simulTIP_Cuprite.mat
% load ../DataSets/Cuprite_Rita/simulTIP_Cuprite_v2.mat

smallIMG = hcube;
imagesc(smallIMG(:,:,[10 20 40]))
T = (smallIMG(:,:,:));

options.Display = true; % Show progress on the command line.
options.Initialization = @cpd_rnd; % Select pseudorandom initialization.
options.Algorithm = @cpd_nls; % Select ALS as the main algorithm.
options.AlgorithmOptions.LineSearch = @cpd_els; % Add exact line search.
options.AlgorithmOptions.TolFun = 1e-12; % Set function tolerance stop criterion
options.AlgorithmOptions.TolX = 1e-12; % Set step size tolerance stop criterion


Uhat = cpd(T,3,options);
That = real(cpdgen(Uhat));


%%

R=5

[n1,n2,L] = size(smallIMG);
Ta = zeros(n1,n2,R);

for i=1:n1,
    for j=1:n2,
        Ta(i,j,:) = constrainedAbundanceEstimation(squeeze(That(i,j,:)),M_ours);
    end
end
figure
for i=1:R
    subplot(2,3,i)
    imagesc(Ta(:,:,i),[0 1])
end

Ta2 = zeros(n1,n2,R);
for i=1:n1,
    for j=1:n2,
        Ta2(i,j,:) = constrainedAbundanceEstimation(squeeze(smallIMG(i,j,:)),M_ours);
    end
end
figure
for i=1:R
    subplot(2,3,i)
    imagesc(Ta2(:,:,i),[0 1])
end



frob(T, That)