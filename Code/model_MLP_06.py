'''
configurable model MLP model
addons:
 - read a config file as input, and the name of the category from which to extract 
   the parameters for the setup
 - add configuration parameters of the batch normalisation layer 
   (I guess starting from "affine")  
'''

import numpy as np
import torch
import torch.nn as nn
import configparser

class MLP_06 (nn.Module) : 	#specializzazione di una classe più generica
    def __init__ (self, N_dims, N_hidden, N_layers, N_classes, dropout_prob = 0., layer_norm = True):
        '''
        Set the model structure

        N_dims      = number of input variables
        N_hidden    = number of nodes per hidden layer
        N_layers    = number of layers betwen input and output
        N_classes   = number of output categories
        droput_prob = dropout probability (default 0)
        layer_norm  = whether to apply layer normalisation (default True) (rinorm pesi) -> evita overfitting, evitando che un singolo nodo non si specializzi troppo, alcuni nodi vengono spenti a caso e quindi la rete deve ridistribuirsi su tanti noti [questo è il dropout, frazione di nodi che vengono spenti durante training]. Infatti, la layer norm evita che un nodo prevalga su altri

        NOTE: instead of a list, a collections.OrderedDict may be used.
              It makes the syntax more convoluted (need a name for each layer),
              but may simplify the access to single layers if known by name.
        '''
        super ().__init__ ()
        if (N_layers < 1) : N_layers = 1 

        layers = []
        self.do_layers = []

        # first layer: from input to hidden
        layers.append (torch.nn.Linear (N_dims, N_hidden))          #crea il primo layer con N_dims = numero di features in ingresso e N_hidden = numero di layers nascosti
        if (dropout_prob > 0) : 
            layers.append (nn.Dropout (dropout_prob))
            self.do_layers.append (layers[-1])
        if (layer_norm) : layers.append (nn.LayerNorm (N_hidden))
        layers.append (torch.nn.ReLU ())

        # intermediate layers
        print("number of layers: ", N_layers)
        for i_layer in range (N_layers):
            next_hid_dim = max(1, int(N_hidden * 2/3))
            layers.append (torch.nn.Linear (N_hidden, next_hid_dim))
            if (dropout_prob > 0) : 
                layers.append (nn.Dropout (dropout_prob))
                self.do_layers.append (layers[-1])
            if (layer_norm) : layers.append (nn.LayerNorm (next_hid_dim))
            layers.append (torch.nn.ReLU ())
            
            N_hidden = next_hid_dim
        
        # last layer: from hidden to output
        layers.append (torch.nn.Linear (N_hidden, N_classes))

        self.model = nn.Sequential (*layers) #layers in sequenza

    # def get_weights (self):
    #     return (self.model[0].weight.clone (), self.model[2].weight.clone ())

    def forward (self, x):
        logits = self.model (x)
        return logits

    def scale_do (self, scale):
        for layer in self.do_layers: layer.p *= scale


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 


if __name__ == "__main__" :

    #N_dims, N_hidden, N_layers, N_classes, dropout_prob = 0., layer_norm = True

    config = configparser.ConfigParser()
    config.read("test_MLP_11.cfg")
    
    device = torch.device ("cuda:0" if torch.cuda.is_available() else "cpu")
    dim = 2      # dimensions
    hid = config.getint("model", "hidden_units")
    layers = config.getint("model", "hidden_layers")
    classes   = config.getint("model", "classes_num")
    do = config.getfloat("training", "dropout_prob")
    nor =  config.getboolean("training", "layer_norm")
    
    model = MLP_06 (dim, hid, layers, classes, do, nor).to (device)
    
    #C = 2      # num_classes
    #L = 6      # number of layers
    #H = 12     # num_hidden_units
    #do = 0.4   # dropout probability
    #nor = True # normalisation layer
    # Estrai i parametri dal file .ini
    '''
    hidden_layers = 
    hidden_units  = config.getint("model", "hidden_units")
    dropout_prob  = config.getfloat("training", "dropout_prob")
    layer_norm    = config.getboolean("training", "layer_norm")
    '''
    print (model)




