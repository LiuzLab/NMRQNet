import numpy as np
from sklearn.metrics import mean_squared_error
from copy import deepcopy
import random

class find_free_edta:

	def __init__(self, original_hp, total_lib):

		self.original_hp = original_hp
		self.total_lib = total_lib


	def find_free_edta(self, n, predict_test_hp, n_iter, step_range, min_concen, max_concen):

		total_lib_map=deepcopy(np.array(self.total_lib))

		shift_list=np.zeros((total_lib_map.shape[0],(2*step_range+1)))
		shift_list[:,step_range]=1

		ori_sample = self.original_hp[n,:]
		start_concen =predict_test_hp[n,:]

		coefficients, score, total_lib_final=self.hillclimbing_shift_free_edta(total_lib = total_lib_map,
			ori_sample=ori_sample,solution = start_concen, n_iter = n_iter, step_range = step_range,
			shift_list = shift_list, min_concen = min_concen, max_concen = max_concen)

		return coefficients, score, total_lib_final




	def hillclimbing_shift_free_edta(self, total_lib, ori_sample, solution, n_iter,step_range,
                                 shift_list,min_concen,max_concen):

		solution_eval = self.calculate_mse(total_lib,ori_sample,solution)

		for i in range(n_iter):
			for j in range(64,66): #only works on two free edta peaks
				modify_concen=np.ones((len(solution)))
				modify_concen[j]=random.uniform(min_concen,max_concen)

				if solution[j]<1e-5:
					adjust_fn=np.zeros((len(solution)))
					adjust_fn[j]=0.1
					candidate=(solution+adjust_fn)*modify_concen
				else:
					candidate = solution*modify_concen

				total_lib_shift=deepcopy(total_lib)
				shift=random.sample(list(range(0,(2*step_range+1))),1)[0]
				steps=shift-np.where(shift_list[j,:]==1)[0][0]

				if steps<0:
					total_lib_shift[j,range(0,(12000-abs(steps)))]=total_lib_shift[j,range(abs(steps),12000)]
					total_lib_shift[j,range((12000-abs(steps)),12000)]=0
				elif steps>0:
					total_lib_shift[j,range(steps,12000)]=total_lib_shift[j,range(0,(12000-steps))]
					total_lib_shift[j,range(0,steps)]=0

				candidte_eval = self.calculate_mse(total_lib_shift,ori_sample,candidate)

				# check if we should keep the new point
				if candidte_eval < solution_eval:# store the new point
					solution, solution_eval,total_lib = candidate, candidte_eval,total_lib_shift
					shift_list[j,np.where(shift_list[j,:]==1)[0][0]]=0
					shift_list[j,shift]=1


		return [solution, solution_eval,total_lib]




	def calculate_mse(self, lib, sample, concen):

		reconstr=np.matmul(concen,lib)
		reconstr=reconstr/np.sum(concen)
		mse=mean_squared_error(sample,reconstr)
		return mse
		





  
