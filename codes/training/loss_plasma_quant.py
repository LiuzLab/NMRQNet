##A customized mse function to train the plasma quantification model

import tensorflow as tf
import numpy as np
from tensorflow.keras import backend as K

def custom_mse(y_true,y_pred):
  
  y_pred=tf.cast(y_pred, tf.float64)
  y_true=tf.cast(y_true, tf.float64)
  
  weight_meta=np.concatenate((np.ones((38,)),np.zeros((27,)))) #first 38 components in the library are interested metabolites
  weight_meta = tf.constant(weight_meta)

  quant_meta_pred=weight_meta*y_pred
  quant_meta_true=weight_meta*y_true

  weight_else=np.concatenate((np.zeros((38,)),np.ones((27,))))
  weight_else = tf.constant(weight_else)

  quant_else_pred=weight_else*y_pred
  quant_else_true=weight_else*y_true

  mse_meta= K.mean(K.square(quant_meta_pred-quant_meta_true),axis=-1)
  mse_else=K.mean(K.square(quant_else_pred-quant_else_true),axis=-1)

 # print(mse_meta.shape)
 # print(mse_else.shape)
  final_mse=100*mse_meta+mse_else #assign higher weights to the mse of interested metabolites

  return final_mse