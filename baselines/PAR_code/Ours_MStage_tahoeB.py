import torch
import torch.nn as nn
from torch.autograd import Variable
import torch.nn.functional as F
from downsampler import *
import numpy as np
import scipy.io as sio
from model import *
import time
import copy
from torch.nn.parameter import Parameter
import os
from tensorboardX import SummaryWriter
from functions import ZFD, psf2otf, gaussian_kernel_2d, add_gaussian_noise, normalizeQuantile, denoising
import tifffile as tiff

from L1_Focal_loss import *

##parameter setting
num_steps = 2001
ReNew = 100
Dsets_name = 'L1 Noised LR_HSI 35NSR'
factor = 4
lr_i = 1.2e-3       #lr_i = 3e-3      # for x16, x32
ld = 0.5            #the lr decrease rate
a = 2.5             #the coefficient of weight
D_net = 3

# Save path of the results
name = 'Tahoe_HSI_t3'
raw_image_dir = './data/' + name + '.mat'
MSI_dir = './data/MSI_tahoe_better_registered.mat'
srf_dir = './data/Tahoe_preproc_t1.mat'
save_path = './results/'

# experiment setups
# signal to noise ratio of HS and MS images
SNR_h = 35
SNR_m = 35
# parameters of blurring kernel
scale_factor = 4
kernel_sigma = 4

if not os.path.exists(save_path):
    os.mkdir(save_path)

# Load data
label = sio.loadmat(raw_image_dir)['HSI'].transpose([2, 0, 1])
hrmsi = sio.loadmat(MSI_dir)['MSI'].transpose([2, 0, 1])
label[label<1e-3] = 1e-3
hrmsi[hrmsi<1e-3] = 1e-3
label[label>1]=1
hrmsi[hrmsi>1]=1

kernel = gaussian_kernel_2d(2*scale_factor, kernel_sigma)
# kernel = np.ones((scale_factor,scale_factor))/(scale_factor**2)

srf = sio.loadmat(srf_dir)['SRF']
dim = np.shape(label) # [band number, height, weight]
band = dim[0]
size = [dim[1], dim[2]]
srf = srf / np.tile(np.sum(srf,1), (band,1)).T # try to normalize

# Normalize the data
label = normalizeQuantile(label, qval=0.999)
label = denoising(label)
hrmsi = normalizeQuantile(hrmsi, qval=0.999)

# Generate the low-resolution hyperspectral image
Fft = psf2otf(kernel, size)
lrhsi = ZFD(label, Fft, scale_factor, size)
# lrhsi = ZFD(label, kernel, scale_factor, size)

# Add noise
hrmsi = add_gaussian_noise(hrmsi, SNR_m)
lrhsi = add_gaussian_noise(lrhsi, SNR_h)


#get PSNR
def PSNR_GPU(im_true, im_fake):
    data_range = 1
    C = im_true.size()[0]
    H = im_true.size()[1]
    W = im_true.size()[2]
    Itrue = im_true.clone().resize_(C*H*W)
    Ifake = im_fake.clone().resize_(C*H*W)
    #mse = nn.MSELoss(reduce=False)
    err = torch.pow(Itrue-Ifake,2).sum(dim=0, keepdim=True).div_(C*H*W)
    psnr = 10. * torch.log((data_range**2)/err) / np.log(10.)
    return torch.mean(psnr)


#get SAM
def SAM_GPU(im_true, im_fake):
    C = im_true.size()[0]
    H = im_true.size()[1]
    W = im_true.size()[2]
    Itrue = im_true.clone().resize_(C, H*W)
    Ifake = im_fake.clone().resize_(C, H*W)
    nom = torch.mul(Itrue, Ifake).sum(dim=0).resize_(H*W)
    denom1 = torch.pow(Itrue,2).sum(dim=0).sqrt_().resize_(H*W)
    denom2 = torch.pow(Ifake,2).sum(dim=0).sqrt_().resize_(H*W)
    sam = torch.div(nom, torch.mul(denom1, denom2)).acos_().resize_(H*W)
    sam = sam / np.pi * 180
    sam = torch.sum(sam) / (H*W)
    return sam



# main processing
if __name__ =='__main__':
    #dir_data = './Mydata'


    #define the lossfunction & optimizer
    mse = nn.MSELoss()
    #L1Loss = nn.L1Loss()
    L1Loss = L1_Focal_Loss()
    k = 0

    #for names in files:
    for ite in range(1):
        #L = [0.01,0.1,1,0.5,0.2]
        lr = lr_i
        k += 1

        # Load the groundtruth
        im_gt = label

        # the fixed spectral downsampler
        p = srf.astype(np.float32)
        p = Variable(torch.from_numpy(p.copy()), requires_grad=False).cuda()

        #trans the data to Variable
        im_gt = Variable(torch.from_numpy(im_gt.copy()),requires_grad=False).type(torch.cuda.FloatTensor).cuda()
        # im_gt = im_gt/(torch.max(im_gt)-torch.min(im_gt))
        s = im_gt.shape
        GT = im_gt.view(s[0],s[1]*s[2])
        
        # Load the HR MSI and LR HSI
        im_m = torch.from_numpy(hrmsi)
        im_m = Variable(im_m, requires_grad=False).type(torch.cuda.FloatTensor).unsqueeze(0)

        im_h = Variable(torch.from_numpy(lrhsi), requires_grad=False).type(torch.cuda.FloatTensor).unsqueeze(0)
        LR_up = F.interpolate(im_h, [im_m.shape[2], im_m.shape[3]], mode='bilinear')
        im_cat = torch.cat((LR_up, im_m), 1)
        
        # Pre-define the spatial downsampler
        down_spa = Downsampler(n_planes=im_gt.shape[0], factor=factor, kernel_type='gauss83', phase=0,preserve_size=True).type(torch.cuda.FloatTensor)

        start_time = time.time()
        # Initialize the generator network 
        net = get_net(im_cat.shape[1], 'skip', 'reflection', n_channels=im_gt.shape[0], skip_n33d=256, skip_n33u=256,skip_n11=1, num_scales=D_net, upsample_mode='bilinear').cuda()

        #net = ThreeBranch_Net().cuda()
        optimizer = torch.optim.Adam(net.parameters(), lr=lr, weight_decay=1e-4)

        f = open(save_path+name+'_result.txt','a+')
        f.write('\n\n\nThe result of PSNR & SAM is :\n\n')
        f.write('The experiment:At {} times SR & All the downsampler is unknown & the DSpa is {}.\n\n'.format(factor,Dsets_name))
        f2 = open(save_path+name+'_loss.txt','w')
        f2.write('\n\n\nThe experiment loss is :\n\n')

        writer = SummaryWriter(save_path+'Logs/'+ name+'/')
        print('Stage three : Producing with the {} image'.format(name))


        # Initialize the weight
        HSI_i = LR_up
        weight_h = Variable(torch.ones_like(im_h.squeeze()), requires_grad=False).type(torch.cuda.FloatTensor)
        weight_m = Variable(torch.ones_like(im_m.squeeze()), requires_grad=False).type(torch.cuda.FloatTensor)

        for i in range(num_steps):

            # Optimizing the generator network
            output = net(im_cat)

            S = output.shape
            Dspa_O = down_spa(output)
            Dspa_O = Dspa_O.view(Dspa_O.shape[1],Dspa_O.shape[2],Dspa_O.shape[3])
            out = output.view(S[1],S[2]*S[3])
            Dspc_O = torch.matmul(p,out).view(im_m.shape[1],S[2],S[3])

            #zero the grad
            optimizer.zero_grad()

            loss = L1Loss(Dspa_O, im_h.squeeze(), weight_h) + L1Loss(Dspc_O, im_m.squeeze(), weight_m)
 
            #backward the loss
            loss.backward()
            #optimize the parameter
            optimizer.step()

            #print('At step {},the loss is {}.'.format(i,loss.data.cpu()))

            # Update the input and weight
            if i%ReNew == 0:
                im_cat = torch.cat((output.detach(), im_m), 1).type(torch.cuda.FloatTensor)
                weight = torch.exp(torch.abs(output.detach() - HSI_i) * a)
                HSI_i = output.detach()
                weight_h = Variable(down_spa(weight).squeeze(), requires_grad=False).type(torch.cuda.FloatTensor)
                weigth_m = Variable(torch.matmul(p,weight.view(S[1], S[2]*S[3])).view(im_m.shape[1], S[2], S[3]), 
                         requires_grad=False).type(torch.cuda.FloatTensor)

            if i%10 == 0:
                f2.write('At step {},the loss is {}\n'.format(i,loss.detach().cpu()))
        
            if i%50 == 0: 
                out = Variable(output,requires_grad=False).cuda()
                out = out.view(S[1],S[2],S[3])
                psnr = PSNR_GPU(im_gt,out)
                sam = SAM_GPU(im_gt,out)
                f.write('{},{}\n'.format(psnr,sam))
                writer.add_scalar('PSNR_Training', psnr, i)
                writer.add_scalar('SAM_Training', sam, i)
                writer.add_scalar('Loss_Trainning', loss.data)
                print('**{}**{}**At the {}th loop the loss&PSNR&SAM is {} {},{}.'.format(k,Dsets_name,i,loss.data,psnr,sam))
                if i%500 == 0:
                    #change the learning rate
                    lr = ld*lr
                    for param_group in optimizer.param_groups:
                        param_group['lr'] = lr 
                if i == num_steps-1:
                   out = np.array(output.squeeze().detach().cpu())
                   image = out
                   image = np.clip(image, 0, 1) 
                   image = image * 255
                   tiff.imsave(save_path + name + str(SNR_h) + 'dB' + '.tif', image.astype(np.uint8))
                   sio.savemat(save_path + name + str(SNR_h) + 'dB' + '.mat', {'data':out.transpose([1,2,0]), 'label':label.transpose([1,2,0])})

        used_time = time.time()-start_time
        f.write('The training time is :{}.'.format(used_time))
        f.close()
        f2.close()
        writer.close()

        # torch.save(net,save_path + name + str(SNR_h) + 'dB' + '.pth')
        # torch.cuda.empty_cache()




