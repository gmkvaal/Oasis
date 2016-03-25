import os
import datetime
import glob
import re
import numpy as np
from tabulate import tabulate

for i in [1,2,3]:
	edge = 0.1/i
	circle = 0.01/i
	os.system("python NSCoupled.py problem=CircleFlowStat \
			   makemesh=True circle=%f edge=%f resultswrite=True" % (circle,edge))


headers = ["Elements", "$C_d$", "$C_l$", "$L_a", "$\Delta P$"]
oneline = []
now = datetime.datetime.now()
atm = "%d.%d.%d.%d" % (now.year, now.month, now.day, now.hour)
os.chdir("/home/guttorm/Desktop/Master/RefinementData/Re20/Results%s" % atm)
list_of_files = glob.glob('OutputCircleStationary*.txt')
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
	oneline.append(temp)
	


print tabulate(oneline,headers, tablefmt="latex")
