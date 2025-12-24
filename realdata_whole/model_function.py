import pandas as pd
import keras
from keras.models import Sequential,Model
from keras.layers import Dense, Dropout, BatchNormalization,Input
from keras.callbacks import EarlyStopping, ModelCheckpoint
import keras.backend as Kr
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.utils import class_weight
import numpy as np
from numpy import exp
# Library for Gaussian process
# import GPy
##Library for visualization
import matplotlib.pyplot as plt
# %matplotlib inline
# %config InlineBackend.figure_format = 'svg'
import matplotlib;matplotlib.rcParams['figure.figsize'] = (10,7)
import pylab 
import time
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from tqdm.keras import TqdmCallback
# import GPy

# Load nonGaussian datasets and do classification on them 

def model_function(df_train, phi, dummy_y, num_class):
    print("##### Warning messages ######")
    class_weights = class_weight.compute_class_weight(class_weight='balanced',
                                                      classes=np.unique(df_train["class"]),
                                                      y=df_train["class"])
    class_weight_dict = dict(enumerate(class_weights))
    # DeepKriging model for continuous data
    model = Sequential()
    # model.add(Dense(100, input_dim = 2,  kernel_initializer='he_uniform', activation='relu'))
    model.add(Dense(100, input_dim = phi.shape[1],  
            kernel_initializer='he_uniform', activation='relu'))
    # model.add(Dropout(rate=0.5))
    # model.add(BatchNormalization())
    model.add(Dense(50, activation='relu'))
    model.add(Dense(50, activation='relu'))
    #     model.add(Dense(100, activation='relu'))
    # model.add(Dense(100, activation='relu'))
    model.add(Dense(50, activation='relu'))
    model.add(Dense(10, activation='relu'))
    model.add(Dense(10, activation='relu'))
    #     model.add(Dense(50, activation='relu'))
    #     model.add(Dense(50, activation='relu'))
    #     model.add(Dense(10, activation='relu'))
    #model.add(Dense(50, activation='relu'))
    #model.add(Dropout(rate=0.5))
    #     model.add(Dense(10, activation='relu'))
    #model.add(BatchNormalization())
    model.add(Dense(10, activation='relu'))
    model.add(Dense(num_class, activation='softmax'))
    NB_START_EPOCHS = 50 
    # NB_START_EPOCHS = 200  # Number of epochs we usually start to train with
#     BATCH_SIZE = 64  
#     fold_no = 1
    optimizer = keras.optimizers.Adam(learning_rate=0.001)
    model.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])

    
    callbacks = [EarlyStopping(monitor='val_accuracy', patience=200),
                 ModelCheckpoint(filepath='indicator_kriging.keras', 
                                 monitor='val_accuracy', save_best_only=True),
                                 TqdmCallback(verbose=1)]
    print("##### End of warning messages ######")
    print('<<<<<<<<<<<<<<<< Fitting DNN-model >>>>>>>>>>>>>>>>>')
    result = model.fit(phi, dummy_y, callbacks=callbacks, class_weight = class_weight_dict,
               validation_split = 0.1, epochs = 500, batch_size = 32, verbose = 0)

    model = keras.models.load_model('indicator_kriging.keras')
    return model
