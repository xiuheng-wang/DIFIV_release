import torch
from torch import nn
import torch.nn.functional as F

class ConvBlock(nn.Module):
    def __init__(self, nChannels, multiply):
        super(ConvBlock, self).__init__()
        self.convd1 = nn.Conv2d(nChannels*multiply, nChannels*multiply, kernel_size=3, padding=1, groups=nChannels*multiply)
        self.convp1 = nn.Conv2d(nChannels*multiply, nChannels, kernel_size=1, padding=0)
        self.convd2 = nn.Conv2d(nChannels, nChannels, kernel_size=3, padding=1, groups=nChannels)
        self.convp2 = nn.Conv2d(nChannels, nChannels*multiply, kernel_size=1, padding=0)

        self.bn1 = nn.BatchNorm2d(nChannels*multiply)
        self.bn2 = nn.BatchNorm2d(nChannels)

    def forward(self, x):
        residual = x

        x = F.relu(self.bn1(x))
        x = self.convp1(self.convd1(x))
        x = F.relu(self.bn2(x))
        x = self.convp2(self.convd2(x))

        return x

class Net(nn.Module):
    def __init__(self, nLayers, nChannels, multiply):
        super(Net, self).__init__()
        self.sequential = nn.Sequential()

        self.sequential.add_module('convd_in', nn.Conv2d(nChannels, nChannels*multiply, kernel_size=3, padding=1, groups=nChannels))
        self.sequential.add_module('convp_in', nn.Conv2d(nChannels*multiply, nChannels*multiply, kernel_size=1, padding=0))
        self.sequential.add_module('relu_in', nn.ReLU())
        for i in range(nLayers):
            self.sequential.add_module('Convblock_%d'%(i+1), ConvBlock(nChannels, multiply))
        self.sequential.add_module('convd_out', nn.Conv2d(nChannels*multiply, nChannels*multiply, kernel_size=3, padding=1, groups=nChannels*multiply))
        self.sequential.add_module('convp_out', nn.Conv2d(nChannels*multiply, nChannels, kernel_size=1, padding=0))

    def forward(self, x):
        x = torch.sub(x, self.sequential(x)) # residual learning
        return x



# source_img = torch.zeros([1,16,16,16])
# device = 'cpu'
# model = Net(nLayers=16, nChannels=source_img.shape[-1], multiply=4).to(device)
# out = model(source_img)
# print(out.shape)