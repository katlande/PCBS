import pandas as pd
from functools import reduce
import sys, getopt

# Command Line Inputs:
def main(argv):
   file_tsv = ''
   file_out = ''
   file_path = ''
   
   try:
      opts, args = getopt.getopt(argv, "hi:o:p:", ["file_tsv=","file_out=","file_path="])
   except getopt.GetoptError:
      print('\nBismark2Matrix.py\n-i <tab-delimited input sample file (file name, sample name, treatment group); expects column names in row 1>\n-o <output file name>\n-p <path to input files>\n')
      sys.exit(2)
   
   for opt, arg in opts:
      if opt == '-h':
         print('\n--- Bismark2Matrix Required Inputs ---\n\n\t-i | --file_tsv\t\t<tab-delimited input sample file (filename, sample, group)>\n\t-o | --file_out\t\t<output file name>\n\t-p | --file_path\t<path to input files>\n\n')
         sys.exit()
      elif opt in ("-i", "--file_tsv"):
         samples = arg
      elif opt in ("-o", "--file_out"):
         outPath = arg
      elif opt in ("-p", "--file_path"):
         inPath = arg
   
   dfList =[]
   sample_file = pd.read_csv(samples, sep="\t",  low_memory=False)
   for i, row in sample_file.iterrows():
       # get info and read in cov file for one sample:
       sample_input=str(sample_file.iat[i, 0])
       sample_name=str(sample_file.iat[i, 1])
       sample_group=str(sample_file.iat[i, 2])
       print("processing "+sample_name+"...")
       cov = pd.read_csv(inPath+'/'+sample_input, sep="\t",  low_memory=False, header=None)
       
       # add cpgID to first columns
       ident = cov[0] + ":" + cov[1].apply(str)
       cov.insert(0, "cpgID", ident, True)
       # add %methylation and nCPGs to the second third columns:
       count = cov[4] + cov[5]
       cov = cov.loc[:, ["cpgID",3]]
       cov.insert(2, sample_name+"_"+sample_group+"_nCpG", count, True)
       cov.rename(columns={3:sample_name+"_"+sample_group+"_PercMeth"}, inplace=True)
       dfList.append(cov)
       
   # Combine all the loci on cpgID, removes NA sites
   print("merging samples...")
   out=reduce(lambda x, y: pd.merge(x, y, on = 'cpgID'), dfList)
   print("Writing output...")
   out.to_csv(outPath, sep="\t", index=False) 
   print("Done :)")
   
if __name__ == "__main__":
   main(sys.argv[1:])






