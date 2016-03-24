import os
import datetime
import glob
import re
import numpy as np
from tabulate import tabulate

for i in [1,2]:
	edge = 0.1/i
	circle = 0.01/i
	os.system("python NSCoupled.py problem=CircleFlowStat \
			   makemesh=True circle=%f edge=%f resultswrite=True" % (circle,edge))
tables = []
innertable = []
now = datetime.datetime.now()
atm = "%d.%d.%d.%d" % (now.year, now.month, now.day, now.hour)
os.chdir("/home/guttorm/Desktop/Master/RefinementData/Re20/Results%s" % atm)
list_of_files = glob.glob('./*.txt')
for resultfile in list_of_files:
	f = open(resultfile, 'r')
	for line in f:
		number = re.findall("[-+]?\d+[\.]?\d*", line)
		word = re.findall(r"[A-Za-z]+", line)  
		firstline=re.findall(r"\bResult[\w]*", line)  
		if len(firstline) == 1: 
			pass
		else:
			innertable.append(word)
			innertable.append(number)
			tables.append(innertable)
			innertable = []		




#print tables

print tabulate(tables, tablefmt="latex")
