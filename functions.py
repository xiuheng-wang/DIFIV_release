import numpy as np
import cv2 
import os
import math
from scipy.io import loadmat
import scipy.io as scio
import tifffile as tiff
from scipy.ndimage import laplace
from scipy.ndimage import convolve
from scipy.ndimage import convolve1d


def normalizeQuantile(Z, qval):
    bands = Z.shape[0]
    for i in range(bands):
        Z[i, :, :] = Z[i, :, :] / np.quantile(Z[i, :, :].flatten(), qval)
    return Z

def denoising(X):
    rows, cols, bands = np.shape(X)
    Y = X
    for i in range(bands):
        x = (X[:,:,i]).flatten();
        A = np.reshape(np.delete(X, i, axis=2),(rows*cols, bands-1))
        invAtA = np.linalg.pinv(np.dot(A.transpose(),A))
        Y[:,:,i] = np.reshape(np.dot((np.dot(A, invAtA)), (np.dot(A.transpose(), x))),(rows, cols))
    return Y

def ZFD(Z, F, scale_factor, size):
    s0 = np.floor(scale_factor / 2)
    # s0 = 0
    dim = np.shape(Z)
    Y = np.zeros(dim)
    # Y = np.real(np.fft.ifftn(np.fft.fftn(Z) * F))
    for i in range(dim[0]):
        Y[i, :, :] = np.real(np.fft.ifft2(np.fft.fft2(Z[i, :, :]) * F))
    Y = Y[:, int(s0):size[0]:scale_factor, int(s0):size[1]:scale_factor]
    return Y

def FDTY(Y, FT, scale_factor, size):
    s0 = np.floor(scale_factor / 2)
    dim = np.shape(Y)
    X = np.zeros([dim[0], size[0], size[1]])
    X[:, int(s0):size[0]:scale_factor, int(s0):size[1]:scale_factor] = Y
    # for i in range(dim[0]):
    #     X[i, :, :] = cv2.resize(Y[i, :, :], (size[1], size[0]), interpolation = cv2.INTER_NEAREST)
    # X = np.real(np.fft.ifftn(np.fft.fftn(X) * FT))
    for i in range(dim[0]):
        X[i, :, :] = np.real(np.fft.ifft2(np.fft.fft2(X[i, :, :]) * FT))
    return X

def A_z_h_FFT(z, lambda_w, W2, rho, F, FT, scale_factor, size):
     s0 = np.floor(scale_factor / 2)
     Z = np.reshape(z, size)
     ZF = np.real(np.fft.ifft2(np.fft.fft2(Z, axes=(1,2)) * np.tile(F, (size[0],1,1)), axes=(1,2)))
     ZFD = ZF[:, int(s0):size[1]:scale_factor, int(s0):size[2]:scale_factor]
     ZFDDT = np.zeros(size)
     ZFDDT[:, int(s0):size[1]:scale_factor, int(s0):size[2]:scale_factor] = ZFD
     ZFDDTFT = np.real(np.fft.ifft2(np.fft.fft2(ZFDDT, axes=(1,2)) * np.tile(FT, (size[0],1,1)), axes=(1,2)))
     Az = ZFDDTFT.flatten() + lambda_w * laplace(W2 * laplace(Z)).flatten() + rho * z
     return Az

def A_z_h_CONV(z, lambda_w, W2, rho, fkernel, scale_factor, size):
      s0 = np.floor(scale_factor / 2)
      fkernel_flip = fkernel[::-1, ::-1] # flipped kernel for F^\top
      Z = np.reshape(z, size)
      ZF = convolve(Z, np.expand_dims(fkernel,0))
      ZFD = ZF[:, int(s0):size[1]:scale_factor, int(s0):size[2]:scale_factor]
      ZFDDT = np.zeros(size)
      ZFDDT[:, int(s0):size[1]:scale_factor, int(s0):size[2]:scale_factor] = ZFD
      ZFDDTFT = convolve(ZFDDT, np.expand_dims(fkernel_flip,0))
      Az = ZFDDTFT.flatten() + lambda_w * laplace(W2 * laplace(Z)).flatten() + rho * z
      return Az
    

def A_z_m(z, lambda_w, W2, rho, RT, R, bands, size):
    Z = np.reshape(z, [bands, size[0], size[1]])
    RTRZ = np.dot(np.dot(RT, R), np.reshape(Z, [bands, size[0] * size[1]]))
    Az = np.reshape(RTRZ, [bands, size[0], size[1]]) + lambda_w * laplace(W2 * laplace(Z)) + rho * Z;
    Az = Az.flatten()
    return Az
    
def gaussian_kernel_2d(kernel_size=4, kernel_sigma=1):
    m,n = (kernel_size-1.)/2, (kernel_size-1.)/2
    y,x = np.ogrid[-m:m+1,-n:n+1]
    kernel = np.exp( -(x*x + y*y) / (2.*kernel_sigma*kernel_sigma))
    kernel[kernel < np.finfo(kernel.dtype).eps*kernel.max()] = 0
    sumh = kernel.sum()
    if sumh != 0:
        kernel /= sumh
        return kernel

def add_gaussian_noise(img, snr = 10):
    snr = 10**(snr/10.0)
    xpower = np.sum(img**2)/img.size
    npower = xpower / snr
    noise = np.random.randn(np.shape(img)[0], np.shape(img)[1], np.shape(img)[2]) * np.sqrt(npower)
    img = img + noise
    return img

def makedir(new_dir):
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

def save_image(output, output_2, label, out_results_path, filename):
    if not os.path.exists(out_results_path):
        os.makedirs(out_results_path)
    image = output
    image = np.clip(image, 0, 1) 
    image = image * 255
    tiff.imsave(out_results_path + str(filename) + '.tif', image.astype(np.uint8))
    
    output = output.transpose([1, 2, 0])
    output_2 = output_2.transpose([1, 2, 0])
    label = label.transpose([1, 2, 0])
    scio.savemat(out_results_path + str(filename) + '.mat', {'sr_h':output, 'sr_m':output_2, 'gt':label})

def rmse(img1, img2):
    dim = np.shape(img1)
    aux = np.sum(np.sum((img1 - img2)**2, axis=1), axis=1) / (dim[1]*dim[2])
    rmse_per_band = np.sqrt(aux)
    rmse = np.sqrt(np.sum(aux)/dim[0])
    return rmse

def psnr(img1, img2):
    img1 = np.clip(img1, 0, 1)
    img2 = np.clip(img2, 0, 1)
    img1 = (img1*255).astype(np.uint8)
    img2 = (img2*255).astype(np.uint8)
    mse = np.mean((img1 - img2) ** 2)
    PIXEL_MAX = 255
    return 10.0 * np.log10(PIXEL_MAX**2 / mse)

def psf2otf(psf, outSize):
    psfSize = np.array(psf.shape)
    outSize = np.array(outSize)
    padSize = outSize - psfSize
    psf = np.pad(psf, ((0, padSize[0]), (0, padSize[1])), 'constant')
    for i in range(len(psfSize)):
        psf = np.roll(psf, -int(psfSize[i] / 2), i)
    otf = np.fft.fft2(psf)
    nElem = np.prod(psfSize)
    nOps = 0
    for k in range(len(psfSize)):
        nffts = nElem / psfSize[k]
        nOps = nOps + psfSize[k] * np.log2(psfSize[k]) * nffts
    if np.max(np.abs(np.imag(otf))) / np.max(np.abs(otf)) <= nOps * np.finfo(np.float32).eps:
        otf = np.real(otf)
    return otf