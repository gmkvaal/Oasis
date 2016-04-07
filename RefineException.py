import os
import datetime
import glob
import re
import numpy as np
from tabulate import tabulate
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import subprocess
from dolfin import *
ax = plt.gca()



def exception(foldername,circleres, edgeres):
	name = "hello"
	os.chdir("/home/guttorm/Desktop/Master/Mesh/Circle/Coarse/AutoMesh")
	os.system("python ControlMakeMesh.py %s %f %f" % (name, circleres, edgeres))
	mesh = Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Coarse/AutoMesh/CFM%s.xml" % name)
	identity = mesh.num_cells()
	now = datetime.datetime.now()
	#atm = "%d.%d.%d.%d" % (now.year, now.month, now.day, now.hour)
	atm = foldername
	identity = mesh.num_cells()
	text_file = open("/home/guttorm/Desktop/Master/RefinementData/Re20/%s/OutputCircleStationary%d.txt" % (atm,identity), "w")
	text_file.write("-----Result for %dE-----\n" % identity)
	text_file.write("h min = %d \n" % 1)
	text_file.write("Cd max: %d \n" % 1)
	text_file.write("Cl max: %d \n" % 1)
	text_file.write("Resirculation length: %d \n" % 1)
	text_file.write("Delta P (front-back): %d \n" % 1)
	text_file.write("Computational time:: %d" % 1)
	text_file.close()
