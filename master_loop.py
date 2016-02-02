import sqs_client
import time
import re
import commonlib
import database

datafile = None

def get_data_file_to_work(message):
    return message.body.rstrip().split(',')[1]

def wait_job_start_loop():
    inloop = True
    message = None
    while inloop:
        message = sqs_client.receive_message('start_job')
        if message is None:
            time.sleep(2)
        else:
            inloop = False

    global datafile
    datafile = get_data_file_to_work(message)
    message_body = message.body
    message.delete()
    if re.search('reference_hash', message_body):
        start_create_hash_from_ref()
    elif re.search('align_base', message_body):
        start_align_base()
    elif re.search('get_mutation', message_body):
        start_get_mutation()

def start_create_hash_from_ref():
    global datafile

    print "start creating hash from reference"

    # delete all saved reference genome hash
    db = database.create_database_connection()
    db.execute("DELETE FROM reference_hash")

    reference_arr = commonlib.read_reference_genome('dataset/%s/ref.txt' % datafile)
    offset = 1000
    for i in xrange(0, len(reference_arr), offset):
        sqs_client.send_message('reference_hash', '%s,%d,%d' % (datafile, i, i+offset))

    wait_reference_hash_loop()

def wait_reference_hash_loop():
    inloop = True
    while inloop:
        remaning_reference_hash = sqs_client.get_queue_remaining_messages('reference_hash')
        if remaning_reference_hash > 0:
            time.sleep(5)
        else:
            inloop = False
    start_align_base()

def start_align_base():
    global datafile

    print "start aligning bases"

    # delete all saved reference genome hash
    db = database.create_database_connection()
    db.execute("DELETE FROM aligned_bases")

    reads = commonlib.read_all_reads("dataset/%s/reads.txt" % datafile)
    offset = 1
    for i in xrange(0, len(reads), offset):
        sqs_client.send_message('align_base', '%s,%d,%d' % (datafile, i, i+offset))

    wait_align_base_loop()

def wait_align_base_loop():
    inloop = True
    while inloop:
        remaining_align_base = sqs_client.get_queue_remaining_messages('align_base')
        if remaining_align_base > 0:
            time.sleep(5)
        else:
            inloop = False
    start_get_mutation()

def start_get_mutation():

    global datafile

    print "start getting mutations"

    # delete all saved reference genome hash
    db = database.create_database_connection()
    db.execute("DELETE FROM mutation")

    reference_arr = commonlib.read_reference_genome('dataset/%s/ref.txt' % datafile)
    offset = 5
    for i in xrange(0, len(reference_arr), offset):
        sqs_client.send_message('get_mutation', '%s,%d,%d' % (datafile, i, i+offset))

    wait_get_mutation_loop()

def wait_get_mutation_loop():
    inloop = True
    while inloop:
        remaining_ref_idx = sqs_client.get_queue_remaining_messages('get_mutation')
        if remaining_ref_idx > 0:
            time.sleep(5)
        else:
            inloop = False
    wait_job_start_loop()

if __name__ == "__main__":
    # wait for start job signal
    wait_job_start_loop()
