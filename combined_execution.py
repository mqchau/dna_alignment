#!/usr/bin/python3
# this file will combine the worker steps together
# to prepare to make this run on spark
from __future__ import print_function

import sys
from operator import add
import re
import ipdb
from pyspark import SparkContext, SparkConf
from align_reads import align_read_new
from pile_up import pile_up

import pprint
pp = pprint.PrettyPrinter(indent=4)

def align_single_read(read_raw):
    print(read_raw)
    mutation_list = align_read_new(read_raw)
    if len(mutation_list) > 0:
        flat_list = mutation_list[0]
        flat_list.extend(mutation_list[1])
        return list(map(convert_mutation_list, flat_list))
    else:
        return []

def convert_mutation_list(mutation):
    return (mutation["ref_idx"], get_mutation_string_csv(mutation))

def get_mutation_string_csv(mutation):
    if mutation["type"] == "delete":
        return "delete,,"
    elif mutation["type"] == "insert":
        return "insert,%s,%d" % (mutation["base"], mutation["insert_idx"])
    else:
        return "%s,%s," % (mutation["type"], mutation["base"])

def pile_up_step(ref_idx_and_mutation_list):
    # ipdb.set_trace()
    ref_idx = ref_idx_and_mutation_list[0]
    mutations = list(ref_idx_and_mutation_list[1])
    # pp.pprint((ref_idx, mutations))
    return pile_up(ref_idx, mutations)

if __name__ == "__main__":
    dataset_name = "10k" # 10k, 1m or 100m


    conf = SparkConf().setAppName("dna_alignment").setMaster("local[4]")
    sc = SparkContext(conf=conf)
    lines = sc.textFile("dataset/%s/reads.txt" % dataset_name, 16)
    
    # aligned_reads = lines.flatMap(align_single_read).groupByKey().reduce(pile_up_step)
    # mutations = aligned_reads.collect()
    # for (ref_idx, mutation) in mutations:
    #     print("%d: %s" % (ref_idx, mutation))
    aligned_reads = lines.flatMap(align_single_read).groupByKey().map(pile_up_step)

    aligned_reads.saveAsTextFile("spark_result")
    sc.stop()

    
    #  
