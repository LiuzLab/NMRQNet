
import numpy as np
from sklearn.metrics import mean_squared_error
from copy import deepcopy
import random


class optimize_concen:

	def __init__(self, opt_pos_concen_results, ori_data):

		self.opt_pos_concen_results = opt_pos_concen_results
		self.ori_data = ori_data


	def optimize_concen(self, n, n_iter, min_concen, max_concen):

		lib = self.opt_pos_concen_results[n][2]
		sample = self.ori_data[n,:,0]
		start_concen = self.opt_pos_concen_results[n][0]

		coefficients, score = self.hillclimbing(total_lib = lib, sample = sample, solution = start_concen,
			n_iter = n_iter, min_concen = min_concen, max_concen = max_concen)

		return coefficients, score





	def hillclimbing(self, total_lib, sample, solution, n_iter, min_concen, max_concen):

		solution_eval = self.calculate_mse(total_lib,sample,solution)

		for i in range(n_iter):
			for j in range(len(solution)):
				modify_concen=np.ones((len(solution)))
				modify_concen[j]=random.uniform(min_concen,max_concen)
				candidate = solution*modify_concen

				candidate_eval = self.calculate_mse(total_lib,sample,candidate)

				if candidate_eval < solution_eval:
					solution, solution_eval = candidate, candidate_eval

		return [solution, solution_eval]




	def calculate_mse(self, lib, sample, concen):

		reconstr=np.matmul(concen,lib)
		reconstr=reconstr/np.sum(concen)
		mse=mean_squared_error(sample,reconstr)
		return mse

