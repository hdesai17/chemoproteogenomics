#!/usr/bin/env python
# coding: utf-8
# Create an argument parser
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--sample-name", required=True, help="Name of the sample")

# Parse the command-line arguments
args = parser.parse_args()

# Access the sample name
sample_name = args.sample_name

# Now you can use the sample_name variable in your Python script
print("Processing sample:", sample_name)

# In[1]:


import re


# In[2]:


cell_lines = [sample_name]
old_dbs = []
new_dbs = []

for c in cell_lines:
    new_dbs.append("2TS_" + c + "_edited.fa")
    old_dbs.append("2TS_" + c + "_custom.fa")


# In[3]:


for c in range(len(cell_lines)):
    
    # name of the fasta file you want to edit
    current_filename = old_dbs[c]

    # updated filename
    new_filename = new_dbs[c]

    # read fasta file into a string
    x = open(current_filename,'r')
    cline = x.read()
    x.close()
    
    # remove duplicate headers from the database
    whole_entry = r'(>sp[^>]+)'
    dupes = re.findall(whole_entry,cline)
    no_dupes = []
    [no_dupes.append(entry) for entry in dupes if entry not in no_dupes]
    
    # write updated database into the new file
    ofile = open(new_filename, "w")

    for v in no_dupes:

            ofile.write(v)

    ofile.close()    
    
    #----------------------------------------------------------------------------------------------------------------------
    
    # read fasta file into a string
    x = open(new_filename,'r')
    cline = x.read()
    x.close()
    
    # regex for matching individual variants
    iv = r'(>sp)((?:\|[^\|]*))((?:\|\D[^\|]*)+\|)((?:[^\|]+\|){2})([^>]+)'
    # store all variants in a list
    variants = re.findall(iv,cline)
    
    # add NA to entries that dont have an ID 
    NA = []
    for v in variants:
        if v[1] == "|":
            NA.append("|NA")
        else:
            NA.append(v[1])
            
    # write updated database into the new file
    ofile = open(new_filename, "w")

    for v in range(len(NA)):

            ofile.write(variants[v][0] + NA[v] + variants[v][2] + variants[v][3] + variants[v][4])

    ofile.close()
