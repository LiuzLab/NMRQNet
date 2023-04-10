
import numpy as np
from scipy.stats import pearsonr


class pair_original_reconstr:

	def __init__(self, predict_test_hp, original_hp):

		self.predict_test_hp = predict_test_hp
		self.original_hp = original_hp

	def global_correct(self, n, total_lib, region, region_step):

		reconstr_hp=np.array(np.matmul(self.predict_test_hp[n,:],total_lib))
		ori_hp=self.original_hp[n,:,:].reshape((12000,))

		step=list(range(-150,150))
		corr_=[]
		shift_array=np.zeros((12000,len(step)))

		for i in range(len(step)):
			shift_array[:,i]=reconstr_hp

			if step[i]<0:
				shift_array[:,i]=np.hstack((shift_array[range(abs(step[i]),12000),i],np.zeros((abs(step[i]),))))

			elif step[i]>0:
				shift_array[:,i]=np.hstack((np.zeros((step[i],)),shift_array[range(0,(12000-step[i])),i]))

			corr_.append(self.region_corr(ori_hp,shift_array[:,i],region,region_step))

		best_step=step[np.argmax(corr_)]

		if best_step<0:
			reconstr_hp[range(0,(12000-abs(best_step)))]=reconstr_hp[range(abs(best_step),12000)]
			reconstr_hp[range((12000-abs(best_step)),12000)]=0

		elif best_step>0:
			reconstr_hp[range(best_step,12000)]=reconstr_hp[range(0,(12000-best_step))]
			reconstr_hp[range(0,best_step)]=0

		return reconstr_hp, best_step


	def region_corr(self, ori_hp,reconstr,region,region_step):

		corr_=0

		for j in range(0,region):
			step=region_step
			corr_=corr_+pearsonr(ori_hp[j*step:(j+1)*step],reconstr[j*step:(j+1)*step])[0]
		return corr_

