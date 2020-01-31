import struct
from scipy import stats
import math
import numpy as np

blockLength=int(input("Enter the block length"))
numberOfbuffers=int(input("Enter the number of buffers"))
size=blockLength*8
alpha=0.01
counter=0
# list which stores distribution
distribution=[]

#KStest

def nint(x):
    #   // Round to nearest integer. Rounds half integers to the nearest even integer
    #   // CODE WRITTEN BY DONGEUN LEE ON OCT. 7, 2015
    if x >= 0:
        i = int(x + 0.5)
        if (i & 1) and (x + 0.5) == i:
            i-=1
    else:
        i = int(x - 0.5)
        if (i & 1) and (x - 0.5) == i:
            i+=1
    return i;

def KolmogorovProb(z):

    fj = [-2,-8,-18,-32]
    r=[None]*4
    w = 2.50662827
    
    const1 = -1.2337005501361697
    const2 = -11.103304951225528
    const3 = -30.842513753404244
    
    u = abs(z)
    if u < 0.2:
        p = 1
    elif u < 0.755:
        v = 1.0/(u*u)
        p = 1 - w*(math.exp(const1*v) + math.exp(const2*v) + math.exp(const3*v))/u
    elif u < 6.8116:
        r[1] = 0
        r[2] = 0
        r[3] = 0
        v = u*u;
        maxj = max(1,nint(3.0/u))
        for j in range (maxj):
            r[j] = math.exp(fj[j]*v)
        p = 2*(r[0] - r[1] + r[2] - r[3])
    else:
        p = 0
    return p

def KolmogorovTest(na, a, nb, b):

    prob = -1
    
    if not a or not b:
        print("KolmogorovTest : NULL pointer")
        return prob
    if not a or not b or na <= 2 or nb <= 2:
        return 0
    #   //  Constants needed
    rna = na
    rnb = nb
    sa  = 1.0/rna
    sb  = 1.0/rnb
    rdiff = 0
    rdmax = 0
    ia = 0
    ib = 0
    #   // Main loop over point sets to find max distance
    #   // rdiff is the running difference, and rdmax the max.
    done = 0
    
    for i in range(na+nb):
        if a[ia] < b[ib]:
            rdiff -= sa
            ia+=1
            if ia >= na:
                done = 1
                break
        elif a[ia] > b[ib]:
            rdiff += sb
            ib+=1;
            if ib >= nb:
                done = 1
                break
        else:
            x = a[ia]
            while(ia < na and a[ia] == x):
                rdiff -= sa
                ia+=1
            while(ib < nb and b[ib] == x ):
                rdiff += sb
                ib+=1
            if (ia >= na):
                done = 1
                break
            if (ib >= nb):
                done = 1
                break
        rdmax = np.fmax(rdmax,abs(rdiff))

    if done==1:
        rdmax = np.fmax(rdmax,abs(rdiff))
        z = rdmax * math.sqrt(rna*rnb/(rna+rnb))
        prob = KolmogorovProb(z)
        return prob

#encoder

with open("A6BUS1C1MAG.csv.bin", "rb") as binary_file: #reading file
    with open("A6BUS1C1MAG.csv.bin.idealem", "wb",buffering = 1000000) as idealem_file: #output to file

        while True:
            exchangeability=False
            # when file has remain, unpack the file
            try: # converting binary to decimal
                data=struct.unpack('d'*blockLength, binary_file.read(size))
            # when all file are unpacked, out of loop
            except:
                break
            
            # This for loop starts when number of length of distribution is bigger than number of buffers 
            for i in range(min(len(distribution),numberOfbuffers)):
            	# pvalue is calculated by KS test
                pvalue = KolmogorovTest( blockLength, sorted(distribution[i]), blockLength, sorted(data) )
                # if data is exchangeable, the data compress to 1 byte which is index of distribution
                if pvalue>=alpha:
                    exchangeability=True
					#print(i)
					#converting decimal to binary
                    b_i=struct.pack('B',i)
                    #writing to compressed file
                    idealem_file.write(b_i)
                    break
        	
            b_data=struct.pack('d'*blockLength, *data)
            #non-exchangeable
            if exchangeability==False:
            	# changingIndex indicates the index which is going to be overwritten
                changingIndex=counter%numberOfbuffers
                # converting decimal to binary
                b_changingIndex=struct.pack('B',changingIndex)
                
                # if length of distribution is less than number of buffers, data will be appended in the list and the compressed file
                if len(distribution)<numberOfbuffers:
                    distribution.append(data)
					# print(changingIndex)
					# writing binary index into compressed file
                    idealem_file.write(b_changingIndex)
 					# print(data)
					# writing binary data into compressed file
                    idealem_file.write(b_data)
                else:
                    # delete oldest index and add new index to the oldest index position
                    distribution.pop(changingIndex)
                    distribution.insert(changingIndex,data)
                    b_ff=struct.pack('B',255)
                    idealem_file.write(b_ff)
					# print('FF')
					# writing binary index into compressed file
                    idealem_file.write(b_changingIndex)
					# print(changingIndex)
					# writing binary data into compressed file
                    idealem_file.write(b_data)
#                    print(data)
                counter+=1
