import sqs_client
import time
import re
import commonlib
import reference_hash 
import align_reads
import get_mutations
import database

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


if __name__ == "__main__":
    # wait for start job signal
    wait_master_loop()
