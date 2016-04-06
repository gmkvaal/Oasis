import os
import glob
from numpy import array
import datetime

results = []
now = datetime.datetime.now()
atm = "%d.%d.%d.%d" % (now.year, now.month, now.day, now.hour)
#os.chdir("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s" % atm)
os.chdir("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/2016.4.4.19")
l = glob.glob('*.txt')
for resultfile in l:
	f = open(resultfile, 'r')
	for line in f:
		output = eval(line)
		results.append(output)

matrix = results[0]
x = results[1]
y = results[2]

print x
print ""
print y
print ""
print matrix
print ""