import numpy as np
from sklearn.metrics import mean_squared_error
from copy import deepcopy
import random

class optimize_pos_concen:

	def __init__(self, predict_test_hp_update, remove_results):

		self.predict_test_hp_update = predict_test_hp_update
		self.remove_results = remove_results


	def optimize_pos_concen(self,n, n_iter, total_lib, best_step, step_range, min_concen, max_concen):
	    
	    lib=deepcopy(np.array(total_lib))
	    
	    if best_step[n]<0:
	        lib[:,range(0,(12000-abs(best_step[n])))]=lib[:,range(abs(best_step[n]),12000)]
	        lib[:,range((12000-abs(best_step[n])),12000)]=0
	    elif best_step[n]>0:
	        lib[:,range(best_step[n],12000)]=lib[:,range(0,(12000-best_step[n]))]
	        lib[:,range(0,best_step[n])]=0
	        
	    start_concen = self.predict_test_hp_update[n,:]
	    lib = np.delete(lib, obj=[64,65], axis=0)
	    start_concen = np.delete(start_concen, obj=[64,65],axis=0)
	    
	    sample = self.remove_results[n]
	    shift_list=np.zeros((lib.shape[0],(2*step_range+1)))
	    shift_list[:,step_range]=1
	    
	    coefficients, score, total_lib_final=self.hillclimbing_shift(total_lib=lib, sample=sample, 
	    solution = start_concen, n_iter=n_iter,step_range=step_range, shift_list = shift_list, 
	    min_concen = min_concen, max_concen = max_concen)
	    
	    return coefficients, score, total_lib_final
		
		
	def hillclimbing_shift(self,total_lib, sample, solution, n_iter, step_range, shift_list, min_concen, max_concen):
	    
	    solution_eval = self.calculate_mse(total_lib,sample,solution)
	    
	    for i in range(n_iter):
	        for j in range(len(solution)):
	            modify_concen=np.ones((len(solution)))
	            modify_concen[j]=random.uniform(min_concen,max_concen)
	            
	            candidate = solution*modify_concen
	            
	            total_lib_shift=deepcopy(total_lib)
	            shift=random.sample(list(range(0,(2*step_range+1))),1)[0]
	            steps=shift-np.where(shift_list[j,:]==1)[0][0]
	            
	            if steps<0:
	                total_lib_shift[j,:]=np.hstack((total_lib_shift[j,range(abs(steps),12000)],np.zeros((abs(steps),))))
	                
	            elif steps>0:
	                total_lib_shift[j,:]=np.hstack((np.zeros((steps,)),total_lib_shift[j,range(0,(12000-steps))]))
	            
	            candidte_eval = self.calculate_mse(total_lib_shift,sample,candidate)
	            # check if we should keep the new point
	            if candidte_eval < solution_eval:
	                # store the new point
	                solution, solution_eval,total_lib = candidate, candidte_eval,total_lib_shift
	                shift_list[j,np.where(shift_list[j,:]==1)[0][0]]=0
	                shift_list[j,shift]=1
	                
	    return [solution, solution_eval,total_lib]
	    
	def calculate_mse(self, lib, sample, concen):
	    
	    reconstr=np.matmul(concen,lib)
	    reconstr=reconstr/np.sum(concen)
	    mse=mean_squared_error(sample,reconstr)
	    return mse






