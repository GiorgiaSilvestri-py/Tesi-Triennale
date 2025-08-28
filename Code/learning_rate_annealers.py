import numpy as np


def expo_lr_scheduler (epoch, initial_lr, final_lr, delay = 0) :
  '''
  returns a multiplicative number to the learning rate
  so that the learning rate decrease is exponential 
  and has an horizontal asymptote at final_lr
  '''

  if (epoch < delay) : return 1.

  lr_frac = final_lr / initial_lr
  return (lr_frac + np.exp (lr_frac - epoch + delay - 1))/(lr_frac + np.exp (lr_frac - epoch + delay))



class EarlyStopping:
    def __init__(self, patience=10, verbose=False):
        self.patience = patience
        self.verbose = verbose
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.best_epoch = None
        self.best_model_state_dict = None

    def __call__(self, val_loss, model, epoch):
        score = -val_loss  

        if self.best_score is None:
            self.best_score = score
            self.best_model_state_dict = model.state_dict()
            self.best_epoch = epoch
        elif score < self.best_score:
            self.counter += 1
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_score = score
            self.best_model_state_dict = model.state_dict()
            self.best_epoch = epoch
            self.counter = 0
            if self.verbose:
                print(f"EarlyStopping: improvement at epoch {epoch}")
                


