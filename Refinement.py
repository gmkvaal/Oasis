import os
import datetime
import glob
import re
import numpy as np
from tabulate import tabulate
import numpy as np
import matplotlib.pyplot as plt
ax = plt.gca()



def runit(initial_edge,initial_circle,iterationlist):

	foldername = "TempRes"
	Cl_reference = 0.0106
	circlelist = []
	os.chdir("/home/guttorm/Desktop/Master/Oasis")
	for i in iterationlist:
		edge = initial_edge;
		circle = initial_circle/i; circlelist.append(circle)
		os.system("python NSCoupled.py problem=CircleFlowStat \
					   makemesh=True circle=%f edge=%f resultswrite=True foldername=%s \
					   element=TaylorHood velocity_degree=3 pressure_degree=2" % (circle,edge,foldername))
	
	headers = ["Elements", "h min", "$C_d$", "$C_l$", "$L_a$", r"$\Delta P$"]
	headers2 = ["Elements", "h min", "$C_l$", "r"]
	table = []; h = []; Cl = []; Cl = []; E = []; table2 = []
	now = datetime.datetime.now()
	#atm = "%d.%d.%d.%d" % (now.year, now.month, now.day, now.hour)
	atm = foldername
	os.chdir("/home/guttorm/Desktop/Master/RefinementData/Re20/%s" % atm)
	l = glob.glob('OutputCircleStationary*.txt')
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	list_of_files = sorted(l, key = alphanum_key)

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
		table2.append([temp[0],temp[1],temp[3]])
		h.append(temp[1])

	#print tabulate(table,headers, tablefmt="grid")
	

	for i in range(len(iterationlist)):
		E.append(table[i][0])
		Cl.append(table[i][3])

	Cl_E = abs(np.array(Cl) - Cl_reference)

	"Error = C * h**alpha"
	A = np.vstack([np.log(np.array(h)), np.ones(len(h))]).T 	
	alpha, c = np.linalg.lstsq(A, np.log(Cl_E))[0]		

	r_list = ["n/a"]
	for i in range(1,len(Cl_E)):
		r_list.append(np.log(Cl_E[i]/Cl_E[i-1])/np.log(h[i]/h[i-1]))

	for i in range(0,len(Cl)):
		table2[i].append(r_list[i])
		
	print tabulate(table2,headers2, tablefmt="latex")
	print "Convergence rate = %.3f" % alpha

	os.system("rm -r /home/guttorm/Desktop/Master/RefinementData/Re20/%s/*" % atm)

	return circlelist, Cl_E, E

iterationlist = [1,4,8,16]
initial_edge=0.025; initial_circle=initial_edge
runit(initial_edge,initial_circle,iterationlist)


