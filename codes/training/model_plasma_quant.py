import tensroflow as tf
from tensorflow import keras
from tensorflow.keras.layers import Input,Conv1D,Dropout,GRU,BatchNormalization
from tensroflow.keras.layers import Activation, MaxPooling1D, Dense
from tensorflow.keras import Model


input_tensor=Input(shape=(12000,1))
conv_1=Conv1D(filters=100, kernel_size=9,padding='same',kernel_initializer='he_normal')(input_tensor)
conv_1=Dropout(0.1)(conv_1)
conv_1b=BatchNormalization()(conv_1)
conv_1b=Activation('relu')(conv_1b)
maxpool_1=MaxPooling1D(3)(conv_1b)

conv_3=Conv1D(200, 5,padding='same',kernel_initializer='he_normal')(maxpool_1)
conv_3b=Activation('relu')(conv_3)
maxpool_2=MaxPooling1D(3)(conv_3b)

conv_5=Conv1D(500, 3,padding='same',kernel_initializer='he_normal')(maxpool_2)
conv_5=Dropout(0.1)(conv_5)
conv_5b=BatchNormalization()(conv_5)

conv_5b=Activation('relu')(conv_5b)
maxpool_3=MaxPooling1D(2)(conv_5b)

GRU=GRU(500,kernel_initializer='he_normal')(maxpool_3)
GRU=Dropout(0.1)(GRU)
GRU_b=Activation('relu')(GRU)
dense_2=Dense(200,activation='relu',kernel_initializer='he_normal')(GRU_b)
output=Dense(65,activation='sigmoid',kernel_initializer='he_normal')(dense_2)


model = Model(input_tensor, output)

model.summary()