import sys

outfile1=open(sys.argv[1]+'_'+"D",'w')
outfile2=open(sys.argv[1]+'_'+"R",'w')

file_raw_barcodes=open(sys.argv[1],'r')
for each_line in file_raw_barcodes:
	L=each_line.strip().split(',')
	if len(L) ==5:
		outfile1.write(L[1]+','+L[4]+'\n')
		outfile2.write(L[0]+','+L[4]+'\n')
outfile1.close()
outfile2.close()
