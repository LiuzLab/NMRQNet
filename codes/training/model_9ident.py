import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D,Activation,MaxPooling1D,GRU,BatchNormalization,Dense


model = Sequential()
model.add(Conv1D(filters=100, kernel_size=9,padding='same',input_shape=(10000,1),kernel_initializer="glorot_normal"))
model.add(BatchNormalization(momentum=0.1)) 
model.add(Activation('relu'))

model.add(Conv1D(100, 9,padding='same',kernel_initializer="glorot_normal"))
model.add(BatchNormalization(momentum=0.1))
model.add(Activation('relu'))


model.add(Conv1D(200, 5,padding='same',kernel_initializer="glorot_normal"))
model.add(BatchNormalization(momentum=0.1))
model.add(Activation('relu'))
model.add(MaxPooling1D(3))

model.add(Conv1D(200, 5,padding='same',kernel_initializer="glorot_normal"))
model.add(BatchNormalization(momentum=0.1))
model.add(Activation('relu'))

model.add(Conv1D(500, 3,padding='same',kernel_initializer="glorot_normal"))
model.add(BatchNormalization(momentum=0.1))
model.add(Activation('relu'))
model.add(MaxPooling1D(2))


model.add(GRU(500, activation='tanh',kernel_initializer="glorot_normal"))

model.add(Dense(100,activation='relu',kernel_initializer="glorot_normal"))
model.add(Dense(9,activation='sigmoid',kernel_initializer="glorot_normal"))
  
print(model.summary())


