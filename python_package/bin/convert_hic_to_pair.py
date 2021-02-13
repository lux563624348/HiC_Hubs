#!/usr/bin/env python
########################################################################
## 01/11/2020
## By Xiang Li,
## lux@gwu.edu
## Peng's Lab
## Version.beta
########################################################################
# Usage 
#python ${EXE_PATH} -b ${INPUT_FILE} -c ${INPUT_NAME} -k ${GENE_LIST_FOLDER}/${GENELISTFILE} -l ${GENELISTFILE: :-4} -r ${RESOLUTION} -f ${FRAGMENTSIZE} -g ${GTFFILE} \
#	-w ${WINDOWSIZE} -n ${NORMALIZATION} -t ${REGIONTYPE} -u ${UP_EXTENSION} -d ${DOWN_EXTENSION} -o ${OUTPUTDIR} -p ${Genic_Partition}
########################################################################

import pandas as pd
import os
import struct
import strawC
from optparse import OptionParser
import sys
####################################################################################
## FUNCTIONS
def readcstr(f):
    # buf = bytearray()
    buf = b""
    while True:
        b = f.read(1)
        if b is None or b == b"\0":
            # return buf.encode("utf-8", errors="ignore")
            return buf.decode("utf-8")
        elif b == "":
            raise EOFError("Buffer unexpectedly empty while trying to read null-terminated string")
        else:
            buf += b
    return None

def read_header(req):
    """
    Takes in a .hic file and returns a dictionary containing information about
    the chromosome. Keys are chromosome index numbers (0 through # of chroms
    contained in file) and values are [chr idx (int), chr name (str), chrom
    length (str)]. 
    """
    chrs = {}
    resolutions = []
    magic_string = struct.unpack(b'<3s', req.read(3))[0]
    req.read(1)
    if (magic_string != b"HIC"):
        error_string = ('... This does not appear to be a HiC file; '
                       'magic string is incorrect')
        force_exit(error_string, req)
    global version
    version = struct.unpack(b'<i', req.read(4))[0]
    
    masterindex = struct.unpack(b'<q', req.read(8))[0]
    genome = b""
    c = req.read(1)
    while (c != b'\0'):
        genome += c
        c = req.read(1)
    genome = genome.decode('ascii')
    # metadata extraction
    metadata = {}
    nattributes = struct.unpack(b'<i', req.read(4))[0]
    for x in range(nattributes):
        key = readcstr(req)
        value = readcstr(req)
        metadata[key] = value
    nChrs = struct.unpack(b'<i', req.read(4))[0]
    for i in range(0, nChrs):
        name = readcstr(req)
        length = struct.unpack(b'<i', req.read(4))[0]
        if name and length:
            chrs[i] = [i, name, length]
    nBpRes = struct.unpack(b'<i', req.read(4))[0]
    # find bp delimited resolutions supported by the hic file
    for x in range(0, nBpRes):
        res = struct.unpack(b'<i', req.read(4))[0]
        resolutions.append(res)
    return chrs, resolutions, metadata
### FUNCTION
####################################################################################
### FUNCTION
### FUNCTIONS
def main(argv):
	desc="Convert .hic to pair txt format --> Format should be: #chr	bin1	bin2	Count"
	parser = OptionParser(description=desc)
	parser.add_option("-i", "--in", action="store", type="string",
			dest="input_path", help="Path to Input HiC file in txt format", metavar="<file>")
	parser.add_option("-n", "--norm", action="store", type="string",
			dest="norm_hic", help="Norm of File.", metavar="<str>")
	parser.add_option("-f", "--file_name", action="store", type="string",
			dest="file_name", help="Name of File.", metavar="<str>")
	parser.add_option("-r", "--resolution", action="store", type="int",
		dest="res", help="Resolution of HiC txt", metavar="<int>")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 3:
		parser.print_help()
		sys.exit(1)
	
	print (" ")
	print ("Here is the Summary of your input.")
	print ("Input Path of HiC file in .hic format: %s" % opt.input_path)
	print ("Name of Input File: %s" % opt.file_name)
	print ("Norm of Input File: %s" % opt.norm_hic)
	print ("Resolution %i" % opt.res)
	print ("End of Summary.")
	print (" ")
	
	## parameters
	PATH_INPUT = opt.input_path
	file_name = opt.file_name
	Norm_File = opt.norm_hic  ## "KR" or "NONE"
	resolution = opt.res


#### Main 
	PATH_File= PATH_INPUT+file_name
	req = open(PATH_File, mode='rb')
	chrs, resolution_out, metadata = read_header(req)
	
	print ("Resolution of input data " + str(resolution_out) )
	
		## [0, 'ALL', 2725521]   [1, '1', 195471971]
	HiC_data = list()
	for idx in chrs:
		chr_idx = str(chrs[idx][1])
		if(chr_idx!='ALL'):
			#result = strawC.strawC(Norm_File, PATH_File, chr_idx, chr_idx, 'BP', resolution)
			result = strawC.strawC(Norm_File, PATH_File, chr_idx, chr_idx, 'BP', resolution)
			for i in range(len(result)):
				HiC_data.append([chr_idx, result[i].binX, result[i].binY, result[i].counts])

	df_hic = pd.DataFrame(data=HiC_data, columns=[0,1,2,file_name[:-4]])
	df_hic.dropna().to_csv(str(df_hic.shape[0])+'_HiC_'+file_name+'.bed', sep='\t', index=None)
#### First GeneBoydy



if __name__ == "__main__":
	main(sys.argv)
