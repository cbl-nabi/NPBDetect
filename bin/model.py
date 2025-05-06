import os

import torch
import torch.nn as nn

class NeuralNet(nn.Module):
    def __init__(self, feats, num_classes):
        super(NeuralNet, self).__init__()
        self.fc1 = nn.Linear(feats, 100)
        self.layer1 = nn.LayerNorm(100)
        self.gelu1 = nn.GELU()
        self.fc2 = nn.Linear(feats, 50)
        self.layer2 = nn.LayerNorm(50)
        self.gelu2 = nn.GELU()
        self.fc3 = nn.Linear(50, num_classes)
 
    def forward(self, x):
        out = self.fc1(x)
        out = self.layer1(out)
        out = self.gelu1(out)
        out = self.fc2(x)
        out = self.layer2(out)
        out = self.gelu2(out)
        out = self.fc3(out)
        return out

def get_model(base_dir, verbose):
    num_classes = 8
    feats = 1131
    if(verbose > 0):
        print('Loading model parameters')
    model = NeuralNet(feats, num_classes)
    model_weights = base_dir + 'model/new_model.pt'
    if( os.path.isfile(model_weights) == False):
        print(''''Missing model weights. Expecting weights in current directory.
              Stopping!''')
        return None
    state = torch.load(model_weights, weights_only=True)
    model.load_state_dict(state)
    return model