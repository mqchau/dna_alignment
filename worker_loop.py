import sqs_client
import ipdb
import time
import re
import commonlib
import reference_hash 
import align_reads
import get_mutations
import database
import multiprocessing as mp
import argparse

def extract_datafile_idx(message_str):
    splitted = message_str.rstrip().split(',')
    return splitted[0], int(splitted[1]), int(splitted[2])

def wait_master_loop():
    while True:
        # check for reference hash first, this has highest priority
        message = sqs_client.receive_message('reference_hash')
        if message is not None:
            datafile, start_idx, stop_idx = extract_datafile_idx(message.body)
            print "reference_hash %d %d" % (start_idx, stop_idx)
            reference_hash.work_small_job(datafile, start_idx, stop_idx)
            message.delete()
            continue

        # check for align base
        message = sqs_client.receive_message('align_base')
        if message is not None:
            datafile, start_idx, stop_idx = extract_datafile_idx(message.body)
            print "align_base %d %d" % (start_idx, stop_idx)
            align_reads.work_small_job(datafile, start_idx, stop_idx)
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
        numworker = mp.cpu_count() + 2

    # run workers
    run_worker_multiple_process(numworker)
