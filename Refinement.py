import os
import datetime
import glob
import re
import numpy as np
from tabulate import tabulate
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from RefineException import *
ax = plt.gca()
import time


def runit(edgeres,circle_list):

	foldername = "TempRes"
	os.chdir("/home/guttorm/Desktop/Master/Oasis")
	for circleres in circle_list:
		if circleres <= edgeres:
			os.system("python NSCoupled.py problem=CircleFlowStat \
						   makemesh=True circleres=%f edgeres=%f resultswrite=True foldername=%s \
						   element=TaylorHood velocity_degree=2 pressure_degree=1" % (circleres,edgeres,foldername))
		else: 
			exception(foldername,circleres,edgeres)
			os.chdir("/home/guttorm/Desktop/Master/Oasis")

		#time.sleep(1) # delays for 1 second
	
	headers = ["Elements", "h min", "$C_d$", "$C_l$", "$L_a$", "$\Delta P$"]
	headers2 = ["Elements", "h min", "$C_l$", "Error", "Comp. time", "r"]
	table = []; table2 = []; E = []; h = []; Cl = []; Cd = []; L_a = []; dP = []; comptime = []
	now = datetime.datetime.now()
	#atm = "%d.%d.%d.%d" % (now.year, now.month, now.day, now.hour)
	atm = foldername
	os.chdir("/home/guttorm/Desktop/Master/RefinementData/Re20/%s" % atm)
	l = glob.glob('OutputCircleStationary*.txt')
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	list_of_files = sorted(l, key = alphanum_key)

	# Read files and append to right location
	for resultfile in list_of_files:
		temp = []
		f = open(resultfile, 'r')
		for line in f:
			number = re.findall("[-+]?\d+[\.]?\d*", line)
			firstline=re.findall(r"\bResult[\w]*", line)  
			if len(firstline) == 1:
				[temp.append(int(elements)) for elements in re.findall("[-+]?\d+[\.]?\d*", line)]	
			else: 
				[temp.append(float(num)) for num in number]	
		table.append(temp[:])
		#table2.append([temp[0],temp[1],temp[3],temp[-1]])
		#h.append(temp[1])
	
	for i in range(len(circle_list)):
		E.append(table[i][0])
		h.append(table[i][1])
		Cd.append(table[i][2])
		Cl.append(table[i][3])
		L_a.append(table[i][4])
		dP.append(table[i][5])
		comptime.append(table[i][6])


	"""
	Cl_E = list(abs(np.array(Cl) - Cl_reference))
	Cd_E = list(abs(np.array(Cd) - Cd_reference))
	L_a_E  = list(abs(np.array(L_a) - L_a_reference))
	dP_E 
	
	A = np.vstack([np.log(np.array(h)), np.ones(len(h))]).T 	
	alpha, c = np.linalg.lstsq(A, np.log(Cl_E))[0]		

	r_list = ["n/a"]
	for i in range(1,len(Cl_E)):
		r_list.append(np.log(Cl_E[i]/Cl_E[i-1])/np.log(h[i]/h[i-1]))

	for i in range(0,len(Cl)):
		table2[i].insert(3, Cl_E[i])
		table2[i].append(r_list[i])
		
	print tabulate(table2,headers2, tablefmt="latex")
	print "Convergence rate = %.4f" % alpha
	"""
	os.system("rm -r /home/guttorm/Desktop/Master/RefinementData/Re20/%s/*" % atm)

	return E, h, Cd, Cl, L_a, dP, comptime

