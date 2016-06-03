import numpy as np

File = open("RPKM-All-Datasets.txt")
FileOut = open("RPKM-All-Datasets-Out.txt", "w")

Dict  = {}
Array = []

Line = File.readline()

for iL in File:
	Split = iL.split("\t")
	Dict[Split[0]] = float(Split[14])

	Array.append(float(Split[14]))

NpArray = np.array(Array)
Min     = np.amin(NpArray[NpArray > 0])

print Min

for iX in Dict:
	FileOut.write(iX + "\t" + str(int(round(Dict[iX]/Min))) + "\n")