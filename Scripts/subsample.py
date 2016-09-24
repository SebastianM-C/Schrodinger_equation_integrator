#!/usr/bin/env python

import os
import sys
from shutil import copy2

os.chdir("../Bin")
#if os.path.isfile("old-output.dat") == False:
copy2("output.dat", "old-output.dat")

if len(sys.argv) == 2 and sys.argv[1] > 0:
	factor = argv[1]
else:
	factor = 2
f = open("output.dat","r+")
lines = f.readlines()
f.seek(0)
lineno = 0
for i in lines:
	lineno = lineno + 1
	if lineno % factor == False:
		f.write(i)
	if i == '\n':
		f.write(i)
f.truncate()
f.close()
