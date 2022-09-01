# coding: utf-8
# Script for performing deep hyperspectral and multispectral fusion with Inter-image varibility 
#
# Reference: 
# Deep Hyperspectral and Multispectral Fusion with Inter-image Varibility
# Xiuheng Wang, Ricardo Borsoi, Cédric Richard, Jie Chen
#
# 2022/01
# Implemented by
# Xiuheng Wang, Ricardo Borsoi
# xiuheng.wang@oca.eu, raborsoi@gmail.com

from __future__ import division
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
import numpy as np
import time
import cv2
import matplotlib.pyplot as plt
from functions import ZFD, FDTY, A_z_h_FFT, A_z_h_CONV, A_z_m, psf2otf, gaussian_kernel_2d, makedir, save_image, rmse, psnr, add_gaussian_noise, normalizeQuantile, denoising
import  os
os.environ["CUDA_VISIBLE_DEVICES"]="0"
import scipy.io as scio
import scipy.interpolate as sint
from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import cg
from scipy.ndimage import laplace
from scipy.ndimage import convolve
import copy
from Denoiser import Denoiser
import argparse

name = 'Tahoe_HSI_t3'
raw_image_dir = './data/' + name + '.mat'
MSI_dir = './data/MSI_tahoe_better_registered.mat'
srf_dir = './data/Tahoe_preproc_t1.mat'
result_image_dir = './results/'
model_path = './models/models_' + name
model_path_h = model_path + '/model_h.pkl'
model_path_m = model_path + '/model_m.pkl'

# experiment setups
# signal to noise ratio of HS and MS images
SNR_h = 35
SNR_m = 35
# parameters of blurring kernel
scale_factor = 4
kernel_sigma = 4

p = 1.8
rho = 1e-1
lambda_w = 2e-3
lambda_h = 1e-2
lambda_m = 1e-2
Iteration = 20
k_subspace = 5
Epoch_p = 10000
Epoch_f = 2000
epsilon = 1e-6


if __name__ == '__main__':

    if not os.path.exists(result_image_dir):
        makedir(result_image_dir)
    if not os.path.exists(model_path):
        makedir(model_path)

    MSE_h = np.zeros(Iteration + 1)
    MSE_m = np.zeros(Iteration + 1)

    # Load data
    label = scio.loadmat(raw_image_dir)['HSI'].transpose([2, 0, 1])
    hrmsi = scio.loadmat(MSI_dir)['MSI'].transpose([2, 0, 1])
    label[label<1e-3] = 1e-3
    hrmsi[hrmsi<1e-3] = 1e-3
    label[label>1]=1
    hrmsi[hrmsi>1]=1

    kernel = gaussian_kernel_2d(2*scale_factor, kernel_sigma)
    # kernel = np.ones((scale_factor,scale_factor))/(scale_factor**2)

    srf = scio.loadmat(srf_dir)['SRF']
    dim = np.shape(label) # [band number, height, weight]
    band = dim[0]
    size = [dim[1], dim[2]]
    srf = srf / np.tile(np.sum(srf,1), (band,1)).T # try to normalize

    # Normalize the data
    label = normalizeQuantile(label, qval=0.999)
    label = denoising(label)
    hrmsi = normalizeQuantile(hrmsi, qval=0.999)
    
    # Generate the low-resolution hyperspectral image
    F = psf2otf(kernel, size)
    lrhsi = ZFD(label, F, scale_factor, size)
    # lrhsi = ZFD(label, kernel, scale_factor, size)

    # Add noise
    hrmsi = add_gaussian_noise(hrmsi, SNR_m)
    lrhsi = add_gaussian_noise(lrhsi, SNR_h)

    # Initialize W
    lrhsi_intern = np.zeros(dim)
    for i in range(band):
        lrhsi_intern[i, :, :] = cv2.resize(lrhsi[i, :, :], (size[1], size[0]), interpolation = cv2.INTER_CUBIC)

    hrmsi_intern = np.zeros(dim)
    x = np.arange(0, np.shape(hrmsi)[0], 1)
    x_int = np.arange(0, np.shape(hrmsi)[0] - 1, (np.shape(hrmsi)[0] - 1) / band)
    for i in range(dim[1]):
        for j in range(dim[2]):
            func = sint.interp1d(x, hrmsi[:, i, j], kind = 'cubic')
            hrmsi_intern[:, i, j] = func(x_int)
    
    W = laplace(lrhsi_intern) - laplace(hrmsi_intern)
    # W[W == 0] = 1e-8
    W = np.sign(W) * (np.abs(W) + epsilon) ** ((p - 2) / 2)
    print(np.min(abs(W)), np.max(abs(W)))
    W2 = W ** 2
    
    # HQS
    # Initialize variables
    Y_h = copy.deepcopy(lrhsi)
    Y_m = copy.deepcopy(hrmsi)
    Z_h = copy.deepcopy(lrhsi_intern)
    V_h = copy.deepcopy(lrhsi_intern)
    Z_m = copy.deepcopy(hrmsi_intern)
    V_m = copy.deepcopy(hrmsi_intern)
    z_h = np.reshape(Z_h, [band, size[0] * size[1]])
    z_m = np.reshape(Z_m, [band, size[0] * size[1]])
    MSE_h[0] = rmse(label, V_h)
    MSE_m[0] = rmse(label, V_m)

    # Calculate some constant variables
    FT = F.conjugate()
    FDTY_h = FDTY(Y_h, FT, scale_factor, size)
    # FDTY_h = FDTY(Y_h, kernel, scale_factor, size)

    R = srf
    RT = srf.transpose()
    RTY_m = np.dot(RT, np.reshape(Y_m, [np.shape(Y_m)[0], size[0] * size[1]]))
    RTY_m = np.reshape(RTY_m, [band, size[0], size[1]])

    print('Before processing: \n' + "MSE_h:" + str(np.around(MSE_h[0], 4)) + "    MSE_m:" + str(np.around(MSE_m[0], 4)))
    for Iter in range(Iteration):
        print("Iter = " + str(Iter) + ":")
        # Calculate Z_h
        B = FDTY_h + lambda_w * laplace(W2 * laplace(Z_m)) + rho * V_h 
        Az = lambda z: A_z_h_FFT(z, lambda_w, W2, rho, F, FT, scale_factor, dim)
        # Az = lambda z: A_z_h_CONV(z, lambda_w, W2, rho, kernel, scale_factor, dim)
        b = B.flatten()
        A = LinearOperator((band* size[0]* size[1], band* size[0]* size[1]), matvec=Az)
        temp, cgit = cg(A, b, x0=Z_h.flatten(), tol=1e-03, maxiter=350)
        Z_h = temp.reshape(dim)

        # Calculate Z_m
        B = RTY_m + lambda_w * laplace(W2 * laplace(Z_h)) + rho * V_m
        b = B.flatten()
        Az = lambda z: A_z_m(z, lambda_w, W2, rho, RT, R, band, size)
        A = LinearOperator((band* size[0]* size[1], band* size[0]* size[1]), matvec=Az)
        temp, cgit = cg(A, b, x0=Z_m.flatten(), tol=1e-03, maxiter=350)
        Z_m = np.reshape(temp, [band, size[0], size[1]])
        # Z_m = np.clip(Z_m, 0, 1)
        print("PSNR_h:" + str(np.around(psnr(label, Z_h), 4)) + "    PSNR_m:" + str(np.around(psnr(label, Z_m), 4)))
        print("MSE_h:" + str(np.around(rmse(label, Z_h), 4)) + "      MSE_m:" + str(np.around(rmse(label, Z_m), 4)))

        # Update W
        W = laplace(Z_h) - laplace(Z_m)
        # W[W == 0] = 1e-8
        W = np.sign(W) * (np.abs(W) + epsilon) ** ((p - 2) / 2)
        print(np.min(abs(W)), np.max(abs(W)))
        W2 = W ** 2

        # Compute V_h
        if Iter == 0:
            print('pre-train deep denoiser engine for V_h:')
            D_h = Denoiser(V_h, dim, k_subspace, Epoch_p, 'pre-train', model_path_h)
        else:
            print('fine-tune deep denoiser engine for V_h:')
            D_h = Denoiser(V_h, dim, k_subspace, Epoch_f, 'fine-tune', model_path_h)
        V_h = (1 / (rho + lambda_h)) * (rho * Z_h + lambda_h * D_h)

        # Compute V_m
        if Iter == 0:
            print('pre-train deep denoiser engine for V_m:')
            D_m = Denoiser(V_m, dim, k_subspace, Epoch_p, 'pre-train', model_path_m)
        else:
            print('fine-tune deep denoiser engine for V_m:')
            D_m = Denoiser(V_m, dim, k_subspace, Epoch_f, 'fine-tune', model_path_m)
        V_m = (1 / (rho + lambda_m)) * (rho * Z_m + lambda_m * D_m)
        print("Denoised: PSNR_h:" + str(np.around(psnr(label, V_h), 4)) + "    PSNR_m:" + str(np.around(psnr(label, V_m), 4)))
        print("Denoised: MSE_h:" + str(np.around(rmse(label, V_h), 4)) + "      MSE_m:" + str(np.around(rmse(label, V_m), 4)))

        # Evaluation
        MSE_h[Iter+1] = rmse(label, V_h)
        MSE_m[Iter+1] = rmse(label, V_m)

    save_image(V_h, V_m, label, result_image_dir, name + str(SNR_h) + 'dB')
    print('Done')
    plt.plot(MSE_h)

    scio.savemat(result_image_dir + name + '_MSE_h_' + str(SNR_h) + 'dB'+'.mat', {'MSE_h':MSE_h})
    plt.savefig(result_image_dir + name + '_MSE_h_' + str(SNR_h) + 'dB'+'.png')
    plt.show()