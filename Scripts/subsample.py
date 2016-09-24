#!/usr/bin/env python

import os
from shutil import copy2

os.chdir("../Bin")
#if os.path.isfile("old-output.dat") == False:
copy2("output.dat", "old-output.dat")

f = open("output.dat","r+")
lines = f.readlines()
f.seek(0)
lineno = 0
for i in lines:
	lineno = lineno + 1
	if lineno%4 == False:
		f.write(i)
	if i == '\n':
		f.write(i)
f.truncate()
f.close()
