
BATCH_SIZE = 32

import tensorflow as tf
from tensorflow.keras import backend as K
import numpy as np

total_lib=np.loadtxt(open("./total_lib_for_9quant_loss.csv","r"))

def custom_loss_wrapper(input_tensor,total_lib):

  def custom_loss(y_true, y_pred):

    y_pred=tf.cast(y_pred, tf.float64)
    y_true=tf.cast(y_true, tf.float64)

    ##########quant_res

    quant_res=K.sum(K.square(y_pred - y_true),axis=-1)

    ########reconstruct mixture only with active metabolites
    logic_inact=K.equal(y_true,0)
    logic_act=K.not_equal(y_true,0)
    logic_act=tf.cast(logic_act,"float64")
    logic_inact=tf.cast(logic_inact,"float64")

    y_pred_act=y_pred*logic_act #take predicted values only for metabolites that are supposed to exist in the mixture
    y_pred_inact=y_pred*logic_inact #take predicted values only for metabolites that are supposed to be 0
    y_pred_act=tf.reshape(y_pred_act,shape=(BATCH_SIZE,10))
    y_pred_inact=tf.reshape(y_pred_inact,shape=(BATCH_SIZE,10))

      
    #reconstruct with only active metabolites in the mixture
    reconstr=tf.matmul(y_pred_act,total_lib)
    
    #split the input and reconstruct data of 40000 into 20 equal buckets
    input_reshape=tf.reshape(input_tensor,shape=(BATCH_SIZE,20,2000))
    input_integral=K.sum(input_reshape,axis=-1)

    recon_reshape=tf.reshape(reconstr,shape=(BATCH_SIZE,20,2000))
    recon_integral=K.sum(recon_reshape,axis=-1)

    input_integral=tf.cast(input_integral,tf.float64)
    recon_integral=tf.cast(recon_integral,tf.float64)
    

    ########### integral res
    integral_res=K.abs(input_integral-recon_integral)

    ##first three zones are discarded due to background noise
    weight=np.concatenate((np.zeros((3,)),np.ones((17,))))
    weight = tf.constant(weight)
    weight=tf.reshape(weight,(20,1))
 
    integral_res=tf.cast(integral_res,tf.float64)
    weight=tf.cast(weight,tf.float64)
    
    integral_weight=tf.matmul(integral_res,weight)

    integral_res=tf.reshape(integral_weight,(BATCH_SIZE,))


    ######## Add another classification loss
    class_res=K.sum(y_pred_inact,axis=-1)
    quant_res=tf.reshape(quant_res,(BATCH_SIZE,))
    #print(class_res.shape)
    #print(quant_res.shape)
    #print(integral_res.shape)

    return integral_res+1000*class_res+1000*quant_res
    
  return custom_loss


