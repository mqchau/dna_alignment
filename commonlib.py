import numpy as np

# convert the genome sequence into integer for efficient storage
# A->1, C->2, G->3, T->4
def convert_genome_str_to_int(genome_str):
    conversion_table = { "A": 1, "C": 2, "G": 3, "T":4 }
    return map(lambda x: conversion_table[x], genome_str)

# given a list of base like 1,2,3,4 return a string combining them like 1234
def get_string_from_mer(mer):
    return "".join(map(lambda x: "%d" % x, mer))

# read the reference genome into a np array
def read_reference_genome(ref_file):
    with open(ref_file, "r") as f:
        lines = f.readlines()

    # calculate how many bases we have
    # assumming we have same number of base from line 1 -> n -1
    total_base = len(lines[1].rstrip()) * len(lines)-2 + len(lines[-1].rstrip())
        
    # allocate memory as np array
    return_array = np.zeros(total_base)
    offset = 0
    for line in lines[1:]:
        line = line.rstrip()
        return_array[offset:offset+len(line)] = convert_genome_str_to_int(line)
        offset += len(line)

    return return_array

    

