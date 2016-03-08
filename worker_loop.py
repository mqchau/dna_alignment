import sqs_client
import ipdb
import time
import re
import commonlib
import align_reads
import get_mutations
import database
import multiprocessing as mp
import argparse
import pickle

last_datafile = None
last_reference_genome = None
last_reference_hash = None

def read_reference_genome_and_hash(datafile):
    reference_genome = None
    reference_hash =None

    with open("dataset/%s/reference_genome.pickle" % datafile, "rb") as f:
        reference_genome = pickle.load(f)

    # read the reference hash
    with open("dataset/%s/all_hash_location.pickle" % datafile, "rb") as f:
        reference_hash = pickle.load(f)

    return reference_genome, reference_hash

def extract_datafile_idx(message_str):
    splitted = message_str.rstrip().split(',')
    return splitted[0], int(splitted[1]), int(splitted[2])

def wait_master_loop():
    global last_datafile, last_reference_genome, last_reference_hash

    while True:
        # check for align base
        message = sqs_client.receive_message('align_base')
        if message is not None:
            datafile, start_idx, stop_idx = extract_datafile_idx(message.body)
            if datafile != last_datafile:
                last_reference_genome, last_reference_hash = read_reference_genome_and_hash(datafile)
                last_datafile = datafile

            print "align_base %d %d" % (start_idx, stop_idx)
            align_reads.work_small_job(last_reference_genome, last_reference_hash, start_idx, stop_idx)
            message.delete()
            continue

        # check for get mutation
        message = sqs_client.receive_message('get_mutation')
        if message is not None:
            datafile, start_idx, stop_idx = extract_datafile_idx(message.body)
            print "get_mutation %d %d" % (start_idx, stop_idx)
            get_mutations.work_small_job(datafile, start_idx, stop_idx)
            message.delete()
            continue

        print "waiting for new job"
        time.sleep(2)

def run_worker_multiple_process(num_process):
    # Setup a list of processes that we want to run
    processes = [mp.Process(target=wait_master_loop) for x in range(num_process)]

    # Run processes
    for p in processes:
        p.start()

    # Exit the completed processes
    for p in processes:
        p.join()

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--numworker", help="how many worker process to run", type=int)
    args = parser.parse_args()
    
    # get number of worker
    if args.numworker is not None:
        numworker = args.numworker
    else:
        # if not indicated, num processor + 1
        numworker = mp.cpu_count()

    # run workers
    run_worker_multiple_process(numworker)
