
import numpy as np
import matplotlib.pyplot as plt

class remove_free_edta:

	def __init__(self, original_hp, find_results):

		self.original_hp = original_hp
		self.find_results = find_results


	def remove_free_edta(self,n, ppm):

		modify_sample = self.original_hp[n,:]

		a=np.argmax(self.find_results[n][2][64,:])-100
		b=np.argmax(self.find_results[n][2][64,:])+100

		modify_sample[a:b]=0

		a=np.argmax(self.find_results[n][2][65,:])-100
		b=np.argmax(self.find_results[n][2][65,:])+100
		modify_sample[a:b]=0

		area = self.integral(ppm = ppm, sample = modify_sample)
		modify_sample = modify_sample/area

		return modify_sample
		
	def plot_fit_edta_peaks(self,n,ppm):
	    ##first free edta peak
	    a=np.argmax(self.find_results[n][2][64,:])-100
	    b=np.argmax(self.find_results[n][2][64,:])+100
	    
	    plt.figure(figsize=[6,4])
	    plt.plot(ppm[a:b],self.original_hp[n,a:b],label = "Original spectrum")
	    plt.plot(ppm[a:b],(self.find_results[n][0][64]*self.find_results[n][2][64,:]/np.sum(self.find_results[n][0]))[a:b],
	    label= "Fitted free EDTA peak")
	    plt.xlim(np.max(ppm[a:b]),np.min(ppm[a:b]))
	    plt.xlabel("Chemical shift (ppm)")
	    plt.ylabel("Intensity (a.u.)")
	    plt.legend(loc = "best")
	    plt.show()
	    
	    #second free edta peak
	    a=np.argmax(self.find_results[n][2][65,:])-100
	    b=np.argmax(self.find_results[n][2][65,:])+100
	    
	    plt.figure(figsize=[6,4])
	    plt.plot(ppm[a:b],self.original_hp[n,a:b],label = "Original spectrum")
	    plt.plot(ppm[a:b],(self.find_results[n][0][65]*self.find_results[n][2][65,:]/np.sum(self.find_results[n][0]))[a:b],
	    label = "Fitted free EDTA peak")
	    plt.xlim(np.max(ppm[a:b]),np.min(ppm[a:b]))
	    plt.xlabel("Chemical shift (ppm)")
	    plt.ylabel("Intensity (a.u.)")
	    plt.legend(loc = "best")
	    plt.show()

	    

	def integral(self, ppm,sample):
		area=0
		for i in range(len(ppm)-1):
			area+=(sample[i]+sample[i+1])*(ppm[i+1]-ppm[i])*0.5
		return area


