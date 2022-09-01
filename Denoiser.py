from __future__ import division

import torch
import torch.nn as nn
from torch.autograd import Variable
import  os
os.environ["CUDA_VISIBLE_DEVICES"]="0"
import time
import numpy as np
from CNN import Net
from skimage.restoration import estimate_sigma
import cv2

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def weights_init_kaiming(m):
    class_name = m.__class__.__name__
    if class_name.find('Linear') != -1:
        nn.init.kaiming_normal_(m.weight)
        if m.bias is not None:
            m.bias.data.zero_()
    elif class_name.find('Conv2d') != -1:
        nn.init.kaiming_normal_(m.weight)
        if m.bias is not None:
            m.bias.data.zero_()
    elif class_name.find('Norm') != -1:
        m.weight.data.normal_(1.0, 0.02)
        if m.bias is not None:
            m.bias.data.zero_()

def Denoiser(V, dim, k_subspace, num_steps, mode = 'pre-train', model_path = './models/model.pkl'):
    # dimensionality reduction
    v = np.reshape(V, [dim[0], dim[1] * dim[2]])
    u, _, _ = np.linalg.svd(v)
    u = u[:, 0:k_subspace]
    eigen_v = np.dot(u.transpose(), v)
    # normalization
    Min = np.amin(eigen_v, axis = 1)
    Max = np.amax(eigen_v, axis = 1)
    scale = Max - Min 
    scale = np.tile(scale, [dim[1] * dim[2], 1]).transpose()
    Min = np.tile(Min, [dim[1] * dim[2], 1]).transpose()
    eigen_v = (eigen_v - Min) / scale
    eigen_V = np.reshape(eigen_v, [k_subspace, dim[1], dim[2]])
    # cv2.imshow("Resized image", V[0,:,:])
    # cv2.waitKey(0)
    # cv2.destroyAllWindows()
    V = np.expand_dims(eigen_V, axis = 0)
    image = Variable(torch.from_numpy(V).float().to(device))
    # initialise the network
    net = Net(nLayers = 8, nChannels = k_subspace, multiply = 4).to(device)
    if mode == 'pre-train' or mode == 'fine-tune':
        # input of the network is image + noise
        sigma_est = np.array(estimate_sigma(np.transpose(eigen_V, [1,2,0]), average_sigmas=False, multichannel=True))
        # print('sigma_est:  ', sigma_est)
        sigma_est = Variable(torch.from_numpy(np.expand_dims(np.transpose(np.tile(sigma_est, [dim[1], dim[2], 1]), [2,0,1]), axis = 0)).float().to(device))
        # network optimizer set up
        optimizer = torch.optim.Adam(net.parameters(), lr=1e-4)
        if mode == 'pre-train':
            net.apply(weights_init_kaiming)
        else:
            checkpoint = torch.load(model_path)
            net.load_state_dict(checkpoint['model_state_dict'])
            optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
            loss = checkpoint['loss']
            net.train()

        # optimize
        for step in range(num_steps):
            # generate noise
            noise = Variable(sigma_est * torch.randn(np.shape(V)).float().to(device))
            # get the network output
            output = net(image + noise)
            # calculate the l2_loss over the masked output and take an optimizer step
            optimizer.zero_grad()
            loss = torch.mean(torch.abs(output - image))
            loss.backward()
            optimizer.step()
            if (step == 0) or ((step + 1) % 1000 == 0):
                print('At step {}, loss is {}'.format(step + 1, loss.data.cpu()))
        # save model
        torch.save({
            'model_state_dict': net.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'loss': loss,
            }, model_path)
    elif mode == 'test':
        checkpoint = torch.load(model_path)
        net.load_state_dict(checkpoint['model_state_dict'])
    else:
        raise SystemExit('error!')
    # output image
    net.eval()
    output = net(image)
    output = output.to('cpu')
    output = output.detach().numpy()
    output = np.squeeze(output) 

    # reconstruct data using denoising engin images
    output = np.clip(output, 0, 1)
    output = np.reshape(output, [k_subspace, dim[1]*dim[2]])
    output = output * scale + Min
    output = np.dot(u, output)
    output = np.reshape(output, dim)

    # clean up any mess we're leaving on the gpu
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return output
