import torch
import numpy as np

def generate_spiral_data (N_samples, N_classes, N_dims, device):
    '''
    generates two spirals stemming from the same point
    X is a 2D set of coordinates
    y is the category, i.e. whether X belongs to one or the other spiral
    '''

    # Generate spiral data
    X = torch.zeros (N_samples * N_classes, N_dims).to (device)
    y = torch.zeros (N_samples * N_classes, dtype=torch.long).to (device)

    for c in range (N_classes):
        index = 0
        t = torch.linspace (0, 1, N_samples)
        # When c = 0 and t = 0: start of linspace
        # When c = 0 and t = 1: end of linpace
        # This inner_var is for the formula inside sin() and cos() 
        #   like sin(inner_var) and cos(inner_Var)
        inner_var = torch.linspace (
            # When t = 0
            (2 * np.pi / N_classes) * (c),
            # When t = 1
            (2 * np.pi / N_classes) * (2 + c),
            N_samples
        ) + torch.randn (N_samples) * 0.2
        
        for ix in range(N_samples * c, N_samples * (c + 1)):
            X[ix] = t[index] * torch.FloatTensor ((
                np.sin (inner_var[index]), np.cos (inner_var[index])
            ))
            y[ix] = c
            index += 1
    return X, y


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


def init_weights (model_element):
    '''
    useful doc for the init possibilities: 
    https://typeoverflow.com/developer/docs/pytorch/nn.init
    https://pytorch.org/docs/stable/nn.init.html

    options:
    torch.nn.init.uniform_(tensor, a=0.0, b=1.0, generator=None)
    torch.nn.init.normal_(tensor, mean=0.0, std=1.0, generator=None)
    torch.nn.init.constant_(tensor, val)
    torch.nn.init.ones_(tensor)
    torch.nn.init.zeros_(tensor)
    torch.nn.init.eye_(tensor) # ID matrix in the tensor
    torch.nn.init.dirac_(tensor, groups=1)  #Dirac delta function
    torch.nn.init.xavier_uniform_(tensor, gain=1.0, generator=None) # xavier uniform
    torch.nn.init.xavier_normal_(tensor, gain=1.0, generator=None)
    torch.nn.init.kaiming_uniform_(tensor, a=0, mode='fan_in', nonlinearity='leaky_relu', generator=None)
    torch.nn.init.kaiming_normal_(tensor, a=0, mode='fan_in', nonlinearity='leaky_relu', generator=None)
    torch.nn.init.trunc_normal_(tensor, mean=0.0, std=1.0, a=-2.0, b=2.0, generator=None)
    torch.nn.init.orthogonal_(tensor, gain=1, generator=None) # semi-orthogonal matrix
    torch.nn.init.sparse_(tensor, sparsity, std=0.01, generator=None) # sparse matrix
    '''
    classname = model_element.__class__.__name__
    if classname.find ('Linear') != -1: # so far this the only one I am using
#        torch.nn.init.constant_ (model_element.weight, 0.)
        torch.nn.init.ones_ (model_element.weight)
#        torch.nn.init.uniform_ (model_element.weight, -0.5, 0.5)
    if classname.find ('Conv') != -1:
        torch.nn.init.normal_ (model_element.weight, 0.0, 0.02)
    elif classname.find ('BatchNorm') != -1:
        torch.nn.init.normal_(model_element.weight, 1.0, 0.02)
        torch.nn.init.zeros_(model_element.bias)


