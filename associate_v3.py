import sys

#v2: report the hamming distance in check files
#v3: report warning info in the check files 

#get info on donor barcodes: sequence & position
def get_donor_info(file1):
	donor=open(file1)
	donor_info_1,donor_info_2={},{}
	for each_line in donor:
		L=each_line.strip().split(',')
		key=L[1]
		donor_info_1[key]={}	#
		donor_info_2[key]=str(L[2])	#info on position, need to make sure no duplicates
	return donor_info_1, donor_info_2

pSL438_1,pSL438_2=get_donor_info(sys.argv[1])
pSL439_1,pSL439_2=get_donor_info(sys.argv[2])

#parse R_D barcodes file
def parse_barcode_pairs(file2, donor_info):
	pairs=open(file2)
	for each_line in pairs:
		L=each_line.strip().split(',')
		if L[1] in donor_info:
			if L[0] in donor_info[L[1]]:
				donor_info[L[1]][L[0]]+=1
			else:
				donor_info[L[1]][L[0]]=1

parse_barcode_pairs(sys.argv[3],pSL438_1)
parse_barcode_pairs(sys.argv[4],pSL439_1)

#calculte the hamming distance
def distance(str1, str2):
	if len(str1) != len(str2):
		raise ValueError("Strand lengths are not equal!")
	else:
		return sum(1 for (a, b) in zip(str1, str2) if a != b)

#output for manual check
def out4check(file3,donor_info,donor_info_pos):
	outfile1=open(file3+'_'+"check",'w')
	R_D={}
	for key in donor_info:
		tmp=donor_info[key]
		if tmp != {}:
			tmp2=dict(sorted(tmp.items(), key=lambda item: item[1])) #sort keys by values
			max_key=list(tmp2.keys())[-1]
			if tmp[max_key] >= 20: 	#cutoff on frequency
				if len(tmp2) >=3: 
					second_key=list(tmp2.keys())[-2]
					third_key=list(tmp2.keys())[-3]
					fourth_key=list(tmp2.keys())[-4]
					fiveth_key=list(tmp2.keys())[-5]
					sixth_key=list(tmp2.keys())[-6]
					output_info=key+'_'+donor_info_pos[key]+','+max_key+','+second_key+','+third_key+','+fourth_key+','+fiveth_key+','+sixth_key+','+str(tmp[max_key])+','+str(tmp[second_key])+','+str(tmp[third_key])+','+str(tmp[fourth_key])+','+str(tmp[fiveth_key])+','+str(tmp[sixth_key])+','
				elif len(tmp2) <=2:	#some may not have the 3rd candidate
					second_key=list(tmp2.keys())[-2]
					output_info= key+'_'+donor_info_pos[key]+','+max_key+','+second_key+','+','+str(tmp[max_key])+','+str(tmp[second_key])+','
				hd=distance(max_key,second_key)	#hamming distance
				warning_info=""
				if hd ==1:	#Evaluate whether there are potentially mixed recipient barcodes 
					if len(tmp2) >=3 and tmp[third_key] >=50: # need to figure out a reasonable cutoff
						warning_info="potentially mixed R barcodes"
				elif hd > 1:
					if tmp[second_key] >=50:
						warning_info="potentially mixed R barcodes" 
				if max_key not in R_D: #Evaluate whether there are potential duplicates of recipient barcodes
					R_D[max_key]=[]	#Recipient barcode : Donor barcode
					R_D[max_key].append(key)
					outfile1.write(output_info+'HD:'+str(hd)+','+'Warning:'+warning_info+'\n')
				else:
					R_D[max_key].append(key)
					outfile1.write(output_info+'HD:'+str(hd)+','+'Warning:'+warning_info+','+'potential duplicates of recipient barcodes'+'\n')
			else:	
				outfile1.write(key+','+'insufficient frequency'+'\n')
		else: #no matched recipient barcodes detected 
			outfile1.write(key+','+'no frequency'+'\n')
	outfile1.close()
	return R_D

R_D_pSL438=out4check(sys.argv[3],pSL438_1,pSL438_2)
R_D_pSL439=out4check(sys.argv[4],pSL439_1,pSL439_2)

#def organize4pos(position_dic, R_D):
#	new_pos_dic={}
#	for key in R_D:
#		new_pos_dic[key]=[]
#		for i in R_D[key]:
#			tmp=position_dic[i]
#			new_pos_dic[key].append(tmp)
#	return new_pos_dic
#
#pSL438_3=organize4pos(pSL438_2,R_D_pSL438)
#pSL439_3=organize4pos(pSL439_2,R_D_pSL439)

#final output
#output format   R_bacode,pSL438_D_barcode,pSL438_D_barcode_pos,pSL439_D_barcode,pSL439_D_barcode_pos
outfile2=open(sys.argv[3]+'_'+"results",'w')
##
tmp=list(R_D_pSL438.keys())+list(R_D_pSL439.keys())
all_recipient_barcodes=list(set(tmp))

for i in all_recipient_barcodes:
	if i in R_D_pSL438 and i in R_D_pSL439:
#		outfile2.write('%r,%r,%r,%r,%r\n' % (i,'_'.join(R_D_pSL438[i]),'_'.join([pSL438_2[m] for m in R_D_pSL438[i]]),'_'.join(R_D_pSL439[i]),'_'.join([pSL439_2[n] for n in R_D_pSL439[i]])))
		tmp=[i,'_'.join(R_D_pSL438[i]),'_'.join([pSL438_2[m] for m in R_D_pSL438[i]]),'_'.join(R_D_pSL439[i]),'_'.join([pSL439_2[n] for n in R_D_pSL439[i]])]
	elif i in R_D_pSL438 and i not in R_D_pSL439:
		tmp=[i,'_'.join(R_D_pSL438[i]),'_'.join([pSL438_2[m] for m in R_D_pSL438[i]]),'NA','NA']
	elif i not in R_D_pSL438 and i in R_D_pSL439:
		tmp=[i,'NA','NA','_'.join(R_D_pSL439[i]),'_'.join([pSL439_2[n] for n in R_D_pSL439[i]])]
	outfile2.write(','.join(tmp)+'\n')
outfile2.close()
