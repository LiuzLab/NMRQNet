######### predictions using the deep learning model ###############
import numpy as np
import tensorflow as tf
from tensorflow import keras

from tensorflow.keras import Model
from tensorflow.keras import backend as K
from tensorflow.keras.models import Sequential
from tensorflow.keras.models import load_model

from tensorflow.keras.layers import Input,GRU,BatchNormalization
from tensorflow.keras.layers import Dense
from tensorflow.keras.layers import Conv1D
from tensorflow.keras.layers import MaxPooling1D
from tensorflow.keras.layers import Activation
from tensorflow.keras.layers import Dropout

from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import ModelCheckpoint,EarlyStopping

from tensorflow.python.framework.ops import disable_eager_execution
disable_eager_execution()


class Test_data:

  def __init__(self, test_data):
    self.test_data = test_data


  def pred_with_deep_learning(self, weight):
      ##create the model layers first
    input_tensor=Input(shape=(12000,1))
    conv_1=Conv1D(filters=100, kernel_size=9,padding='same',kernel_initializer='he_normal')(input_tensor)
    conv_1=Dropout(0.1)(conv_1)
    conv_1b=BatchNormalization()(conv_1)
    conv_1b=Activation('relu')(conv_1b)
    maxpool_1=MaxPooling1D(3)(conv_1b)

    conv_3=Conv1D(200, 5,padding='same',kernel_initializer='he_normal')(maxpool_1)
    conv_3=Activation('relu')(conv_3)
    maxpool_2=MaxPooling1D(3)(conv_3)

    conv_5=Conv1D(500, 3,padding='same',kernel_initializer='he_normal')(maxpool_2)
    conv_5=Dropout(0.1)(conv_5)
    conv_5b=BatchNormalization()(conv_5)
    conv_5b=Activation('relu')(conv_5b)
    maxpool_3=MaxPooling1D(2)(conv_5b)

    gru=GRU(500,kernel_initializer='he_normal')(maxpool_3)
    gru=Dropout(0.1)(gru)
    gru=Activation('relu')(gru)
    dense=Dense(200,activation='relu',kernel_initializer='he_normal')(gru)
    output=Dense(65,activation='sigmoid',kernel_initializer='he_normal')(dense)

    model = Model(input_tensor, output)

    ##load the trained weights
    model.load_weights(weight)

    predict_test_hp=model.predict(self.test_data)

    ##Assign the estimated EDTA values to all its clusters
    update_test_hp=np.zeros((predict_test_hp.shape[0],(predict_test_hp.shape[1]+4)))
    update_test_hp[:,:65]=predict_test_hp
    for i in range(65,69):
      update_test_hp[:,i]=predict_test_hp[:,-1]
     
    update_test_hp[:,64:66]=0 ##two free EDTA peaks are set to 0
    for i in range(0,update_test_hp.shape[0]): # the total relative quantifications are normalized to 1
        update_test_hp[i,:]=update_test_hp[i,:]/np.sum(update_test_hp[i,:])
    
    return update_test_hp

