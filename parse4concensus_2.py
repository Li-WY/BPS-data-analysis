import sys


def parse_bartender(file1,file2):
	cluster_ID, cluster_all={},{}
	barcode_cluster=open(file1)	#extract info on cluster ID, center, and count
	next(barcode_cluster)
	for each_line in barcode_cluster:
		L=each_line.strip().split(',')
		key=L[0]
		value=L[1]+'&&&'+str(L[3])
		cluster_ID[key]=value
	barcode_barcode=open(file2)	#extract info on raw barcodes
	next(barcode_barcode)
	for each_line in barcode_barcode:
		L=each_line.strip().split(',')
		key=L[0]
		tmp=cluster_ID[L[2]].split('&&&')
		value=tmp[0]+'&&&'+str(L[1])+'&&&'+tmp[1]
		cluster_all[key]=value
	return cluster_all

#parse donor and recipient barcodes
Donor_all=parse_bartender(sys.argv[1],sys.argv[2])
Recipient_all=parse_bartender(sys.argv[3],sys.argv[4])
#print Donor_all['ATTACAAAAAACATTAAGC']

#correct potential technical errors for barcodes
outfile1=open(sys.argv[5]+'_'+"corrected2",'w')
file_raw_barcodes=open(sys.argv[5],'r')
R_D={}
for each_line in file_raw_barcodes:
	L=each_line.strip().split(',')
	if L[0] in Recipient_all and L[1] in Donor_all:
		tmp_1=Recipient_all[L[0]].split('&&&')
		tmp_2=Donor_all[L[1]].split('&&&')
		outfile1.write(tmp_1[0]+','+tmp_2[0]+'\n')  # need to figure out why different columns of sys.argv[5]
