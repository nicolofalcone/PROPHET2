"""
script to perform uncertainties analysis
using RAVENH and RELAP5-3D
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from pyXSteam.XSteam import XSteam
from math import pi, sqrt
from scipy.signal import savgol_filter
xsteam = XSteam(XSteam.UNIT_SYSTEM_MKS)

class UncertAnalysis:
	"""
	class to load the output .csv file from any run and
	perform the uncertainties analysis of data
	"""

	def __init__(self, mainDir):

		self.filenames = sorted(glob.glob(".\\UncertaintyAnalysisMC\\*\\*.csv"))
		self.inputs = sorted(glob.glob(".\\UncertaintyAnalysisMC\\*\\*.i"))
		self.varlist = sorted(pd.read_csv(self.filenames[0]).keys())
		self.varlist.remove('time')
		print(self.varlist)
		

		simul = []

		for filename in self.filenames:
			fold = filename.split('\\')[2]
			simul.append(fold)

		simul = list(map(int,simul))

		ntest = 0

		for inp in self.inputs:
			fold = inp.split('\\')[2]
			if int(fold) > ntest:
				ntest = int(fold)

		print('\n')
		print("%i simulations launched.\n" % ntest)
		print("%i runs succeded.\n" % len(self.filenames))

		test = np.arange(1,(ntest+1))
		list(test)

		has_gone = [item in simul for item in test]
		has_failed = [not i for i in has_gone]

		failed_run = test[has_failed]

		if len(failed_run) > 0:
			if len(failed_run) == 1:
				print('%i run failed.\n' % len(failed_run))
			elif len(failed_run) > 1:
				print('%i runs failed.\n' % len(failed_run))

		print("Failed runs: %s.\n" % failed_run)

		perc_file = ".\\PROPHET2_percentile_0.csv"
		if not os.path.exists(perc_file):
			print("PROPHET2_percentile_0.csv does not exist.\n")
		else:
			self.percentile = pd.read_csv(perc_file)

		if not os.path.exists(".\\Prova18 22_02_2019\\prova.hdf5"):
			self.ExpTest = pd.read_excel(".\\Prova18 22_02_2019\\prova.xlsx")
			fname = ".\\Prova18 22_02_2019\\prova.hdf5"
			self.ExpTest.to_hdf(fname,'table')
			print('Experimental data loaded.\n')
		else:
			self.ExpTest = pd.read_hdf(".\\Prova18 22_02_2019\\prova.hdf5","table")
			print('Experimental data loaded.\n')


	def _units(self,variable):
		"""
		This function gives as output the unit of the selected variable.
		"""
		par= variable.split("_")

		if int(par[1]) >= 201 and int(par[1])<221:
			unit = "Temperature [°C]"
		elif int(par[1]) >= 221 and int(par[1])<224:
			unit = "Pressure [bar]"
		elif int(par[1]) >= 224 and int(par[1])<230:
			unit = "Differential Pressure [mbar]"
		elif int(par[1]) == 230:
			unit = "Level [cm]"
		elif int(par[1]) == 231:
			unit = "Mass Flow Rate [kg/s]"
		
		return unit

	def _mflowrate(self,dp,temp,pabs):
		"""
		This function computes the experimental mass flow rate and save it in to a .hdf5 file.

		====================================
		Variables
		====================================
		
		dp: differentail pressure across the orifice.

		temp: fluid temperature at the orifice inlet or outlet.

		pabs: absolute pressure inside the system.
		"""

		if not os.path.exists(".\\Prova18 22_02_2019\\mflowrate.hdf5"):

			ee=1
			dd=0.0055
			DD=0.0207
			bb=dd/DD

			mi=[]
			[mi.append(xsteam.my_pt(p,t)) for p,t in zip(pabs,temp)]
			rho=[]
			[rho.append(xsteam.rho_pt(p,t)) for p,t in zip(pabs,temp)]

			Re=[]
			Re.append(0)
			CC=[]
			CC.append(0)
			alpha=[]
			alpha.append(0)
			mflow=[]
			mflow.append(0)

			for i in range(1,dp.shape[0]):
				Re.append(4*abs(mflow[i-1])/(pi*DD*mi[i]))
				CC.append(0.5959+0.0312*(bb**2.1)-0.184*(bb**8)+0.0029*(bb**2.5)*((10**6/(Re[i]+100))**0.75))
				alpha.append(CC[i]/((1-bb**4)**0.5))
				if dp[i] >= dp[0]:
					mflow.append(alpha[i]*ee*(pi*dd**2/4)*sqrt(2*rho[i]*100*(abs(dp[i]-dp[0]))))
				else:
					mflow.append(-alpha[i]*ee*(pi*dd**2/4)*sqrt(2*rho[i]*100*(abs(dp[i]-dp[0]))))

			mflowrate=pd.Series(mflow)
			mflowrate.to_hdf(".\\Prova18 22_02_2019\\mflowrate.hdf5",'table')
		else:
			mflow=pd.read_hdf(".\\Prova18 22_02_2019\\mflowrate.hdf5",'table')

		return mflow

	def _smooth(self,y,winsize=101,poly=3):

		y_smooth = savgol_filter(y, winsize, poly)

		return y_smooth

	def _varname(self,variable):
		"""
		This function returns as output the real name of the variable to be used ad plot title.
		"""

		if variable == "cntrlvar_201":
			title = "T1"
			expvar = self.ExpTest.T1
		elif variable == "cntrlvar_202":
			title = "T2"
			expvar = self.ExpTest.T2
		elif variable == "cntrlvar_203":
			title = "T3"
			expvar = self.ExpTest.T3
		elif variable == "cntrlvar_204":
			title = "T4"
			expvar = self.ExpTest.T4
		elif variable == "cntrlvar_205":
			title = "T5"
			expvar = self.ExpTest.T5
		elif variable == "cntrlvar_206":
			title = "T6"
			expvar = self.ExpTest.T6
		elif variable == "cntrlvar_207":
			title = "T7"
			expvar = self.ExpTest.T7
		elif variable == "cntrlvar_208":
			title = "T8"
			expvar = self.ExpTest.T8
		elif variable == "cntrlvar_209":
			title = "T9"
			expvar = self.ExpTest.T9
		elif variable == "cntrlvar_210":
			title = "TP-1"
			expvar = self.ExpTest.Tpool1
		elif variable == "cntrlvar_211":
			title = "TP-2"
			expvar = self.ExpTest.Tpool2
		elif variable == "cntrlvar_212":
			title = "TP-3"
			expvar = self.ExpTest.Tpool3
		elif variable == "cntrlvar_213":
			title = "TW-A"
			expvar = self.ExpTest.TW_A
		elif variable == "cntrlvar_214":
			title = "TW-B"
			expvar = self.ExpTest.TW_B
		elif variable == "cntrlvar_215":
			title = "TW-C"
			expvar = self.ExpTest.TW_C
		elif variable == "cntrlvar_216":
			title = "TW-D"
			expvar = self.ExpTest.TW_D
		elif variable == "cntrlvar_217":
			title = "TW-E"
			expvar = self.ExpTest.TW_E
		elif variable == "cntrlvar_218":
			title = "TW-F"
			expvar = self.ExpTest.TW_F
		elif variable == "cntrlvar_219":
			title = "TW-G"
			expvar = self.ExpTest.TW_G
		elif variable == "cntrlvar_220":
			title = "TW-H"
			expvar = self.ExpTest.TW_H
		elif variable == "cntrlvar_221":
			title = "p1"
			expvar = self._smooth(self.ExpTest.P1_bar)
		elif variable == "cntrlvar_222":
			title = "p2"
			expvar = self._smooth(self.ExpTest.P2_bar)
		elif variable == "cntrlvar_223":
			title = "p4"
			expvar = self._smooth(self.ExpTest.P4_bar)
		elif variable == "cntrlvar_224":
			title = "dp1"
			expvar = self.ExpTest.DP1b_mbar
		elif variable == "cntrlvar_225":
			title = "dp2"
			expvar = self.ExpTest.DP2b_mbar
		elif variable == "cntrlvar_226":
			title = "dp3"
			expvar = self.ExpTest.DP3b_mbar
		elif variable == "cntrlvar_227":
			title = "dp5"
			expvar = self.ExpTest.DP5_mbar
		elif variable == "cntrlvar_228":
			title = "dp6"
			expvar = self.ExpTest.DP6_mbar
		elif variable == "cntrlvar_229":
			title = "dp7"
			expvar = self.ExpTest.DP7mbar
		elif variable == "cntrlvar_230":
			title = "Pool Level"
			expvar = self.ExpTest.P3_cm
		elif variable == "cntrlvar_231":
			title = "Mass Flow Rate"
			expvar = self._smooth(self._mflowrate(self.ExpTest.DP4_mbar,self.ExpTest.T8,self.ExpTest.P4_bar))
		
		return [title, expvar]


	def DataPlot(self, variable, t_start=0, t_end=18000):
		"""
		This function plots the selcted parameters
		for all the simulations done.

		=========================
		variable: string
			y variable to plot

		t_start: integer
			left limit in seconds for the plot. Default=0 s

		t_end: integer
			right limit in seconds for the plot. Default=18000 s

		=========================
		"""

		if type(variable) == str:
			pass
		else:
			str(variable)

		variable = variable.lower()

		if not variable in self.varlist:
			print("The selected variable is not available.\nPlease choose one of the following ones.\n")
			print(self.varlist)
			return None
		else:
			
			ext="pdf"
			name = "%s.%s" % (variable, ext)
			fig = plt.figure(figsize=(10,8))
			for filename in self.filenames:
				data = pd.read_csv(filename)
				lb = ("run%s" % filename.split("\\")[2])
				param = data[variable]
				plt.plot(data.time, param, label=lb)

			# plt.legend(ncol=2)
			plt.plot(self.ExpTest.X_Value,self._varname(variable)[1],color="k")
			plt.xlim(t_start, t_end)
			plt.xlabel('Time [s]',fontsize=18)
			plt.ylabel(self._units(variable),fontsize=18)
			plt.title(self._varname(variable)[0],fontsize=18)
			plt.grid(True)
			fig.show()
			# fig.savefig(name, format=ext)

	def DataPlotPercentile(self, variable, t_start=0, t_end=18000):
		"""
		This function plots the selcted parameters
		for all the simulations done.

		=========================
		variable: string
			y variable to plot

		t_start: integer
			left limit in seconds for the plot. Default=0 s

		t_end: integer
			right limit in seconds for the plot. Default=18000 s

		=========================
		"""

		if type(variable) == str:
			pass
		else:
			str(variable)

		variable = variable.lower()

		if not variable in self.varlist:
			print("The selected variable is not available.\nPlease choose one of the following ones.\n")
			print(self.varlist)
			return None
		else:
			
			ext="pdf"
			name = "%s.%s" % (variable, ext)
			fig = plt.figure(figsize=(18,8))
			name5 = "percentile_5_%s" % variable
			name95 = "percentile_95_%s" % variable
			# plt.plot(self.percentile.time, self.percentile[name5], label="Percentile 5%")
			# plt.plot(self.percentile.time, self.percentile[name95], label="Percentile 95%")

			# plt.legend(ncol=2)
			plt.plot(self.ExpTest.X_Value,self._varname(variable)[1],color="k", label="Experimental")
			plt.plot(self.percentile.time, self.percentile[name5], label="Percentile 5%")
			plt.plot(self.percentile.time, self.percentile[name95], label="Percentile 95%")
			plt.xlim(t_start, t_end)
			plt.xlabel('Time [s]',fontsize=18)
			plt.ylabel(self._units(variable),fontsize=18)
			plt.title(self._varname(variable)[0],fontsize=18)
			plt.grid(True)
			plt.legend()
			fig.show()
			# fig.savefig(name, format=ext)

		return

	def Pearson(self):

		"""
		This fucntion plots the coefficients of the Pearson matrix related to the uncertainty analysis.
		"""

		file = ".\\PROPHET2_Pearson_matrix_0.csv"
		self.matrix = pd.read_csv(file)
		# file = ".\\PERSEO_Pearson_matrix_0.xlsx"
		# self.matrix = pd.read_excel(file)
		self.matrix.rename(index=str, columns={'pear_cntrlvar_223_20220200:4':'power',
			'pear_cntrlvar_223_1700101:3':'compressible_volume','pear_cntrlvar_223_2103001:9':'loss_coeff_pool',
			'pear_cntrlvar_223_1300902:1':'loss_coeff_90°bend','pear_cntrlvar_223_1500904:1':'loss_coeff_orifice',
			'pear_cntrlvar_223_20100201:2':'therm_conductivity_ins',
			'pear_cntrlvar_223_20100251:2':'therm_capacity_ins','pear_cntrlvar_223_20220001:2':'T_amb',
			'pear_cntrlvar_223_20230000:5':'htc_loss_bare_pipe','pear_cntrlvar_223_20230100:5':'htc_loss_bayonet',
			'pear_cntrlvar_223_20240000:5':'htc_loss_ins_pipe','pear_cntrlvar_223_20250000:5':'htc_loss_pool',
			'pear_cntrlvar_223_1800201:2':'p_amb','pear_cntrlvar_223_11401901:11':'fouling_pool'},inplace=True)
		lb=list(self.matrix.keys())
		lb.remove('time')

		fig = plt.figure(figsize=(18,8))
		for ii,lb in enumerate(lb):
			plt.plot(self.matrix.time,self.matrix.iloc[:,ii+1],label=lb)
		plt.xlim(left=0)
		plt.ylim(-1, 1)
		plt.xlabel('Time [s]',fontsize=18)
		plt.ylabel('Coefficient [-]',fontsize=18)
		plt.title('Pearson coefficients',fontsize=18)
		plt.grid(True)
		plt.legend(ncol=2)
		fig.show()

		return

	def _get_cmap(self,n,name='jet'):


		return plt.cm.get_cmap(name,n)

	def Cobweb(self):
	
		"""
		This function plots the cobweb plot of the normalaized values of the input parameters used to perform the uncertainty analysis.
		"""

		file = ".\\UncertaintyRanges.xlsx"
		self.range = pd.read_excel(file)
		parameters = self.range.keys()

		#color=self._get_cmap(len(self.filenames))

		fig = plt.figure(figsize=(18,8))
		for nn,filename in enumerate(self.filenames):
			self.data = pd.read_csv(filename)
			self.data.rename(index=str, columns={'20220200:4':'power',
			'1700101:3':'compressible_volume','2103001:9':'loss_coeff_pool',
			'1300902:1':'loss_coeff_90°bend','1500904:1':'loss_coeff_orifice','20100201:2':'therm_conductivity_ins',
			'20100251:2':'therm_capacity_ins','20220001:2':'T_amb',
			'20230000:5':'htc_loss_bare_pipe','20230100:5':'htc_loss_bayonet',
			'20240000:5':'htc_loss_ins_pipe','20250000:5':'htc_loss_pool',
			'1800201:2':'p_amb','11401901:11':'fouling_pool'},inplace=True)
			lb = ("run%s" % filename.split("\\")[2])
			color=np.random.rand(3)
			for ii,param in enumerate(parameters):
				lbound=self.range.loc["min",param]
				ubound=self.range.loc["max",param]
				if ii==0:
					plt.plot(ii+1,(self.data.loc["1",parameters[ii]]-lbound)/(ubound-lbound),c=color,marker='.',markersize=10,label=lb)
				else:
					plt.plot(ii+1,(self.data.loc["1",parameters[ii]]-lbound)/(ubound-lbound),c=color,marker='.',markersize=10)
		plt.xticks(np.arange(1,15))
		plt.ylim(0, 1)
		plt.grid(True,which='major',axis='y')
		plt.xlabel('Parameter [-]')
		plt.ylabel('Normalized Value [-]')
		fig.show()

		return
		

if __name__ == "__main__":
    mainDir = os.curdir
    UA = UncertAnalysis(mainDir)