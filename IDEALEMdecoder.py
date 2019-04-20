import struct
from scipy import stats
import random


blockLength=int(input("Enter the block length"))
numberOfbuffers=int(input("Enter the number of buffers"))
size=blockLength*8
alpha=0.01
size=blockLength*8
counter=0
place=[]
index=[]


with open("A6BUS1C1MAG.csv.bin.idealem", "rb") as binary_file: # read compressed file
    with open("A6BUS1C1MAG.csv.bin.idealem.bin", "wb") as idealem_file: # write decompressed data
        # converting binary index to decimal
        data=struct.unpack('B', binary_file.read(1))
        print(data)
        while data:
            # if there is no index same with data[0] in the place list, skip this statements
            if data[0] not in place:
                changingIndex=counter%numberOfbuffers
                # if the data is 255 which indicates replace the oldest index to new index, delete changingIndex in the place list and the index list.
                if data[0]==255:
                    place.pop(changingIndex)
                    index.pop(changingIndex)
                    # read binary index and convert decimal
                    data=struct.unpack('B', binary_file.read(1))
                    print(data)
                # write index into place list
                place.insert(changingIndex,data[0])
                # read binary data1(not index) and convert to decimal
                data1=struct.unpack('d'*blockLength, binary_file.read(size))
                # convert decimal data to binary
                b_data1=struct.pack('d'*blockLength,*data1)
                # write into decompressed file
                idealem_file.write(b_data1)
                # write data1(not index) into index list
                index.insert(changingIndex,list(data1))
                print(list(data1))
                print(place)
                counter+=1

            else:
                # exchangeable index comes to here and exchangeable index will be exchangeable data
                data=index[data[0]]
                # r_data is data which gets ramdom order
                r_data=random.sample(data, len(data))
                # b_data is the ramdom decimal data converted to binary
                b_data=struct.pack('d'*blockLength, *r_data)
                # write b_data into decompressed file
                idealem_file.write(b_data)
                print(r_data)
            try:
                # read binary index and convert to decimal
                data=struct.unpack('B', binary_file.read(1))
            except:
                break
#         print(data)
