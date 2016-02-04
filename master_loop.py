import sqs_client
import time
import re
import commonlib
import database

datafile = None

def get_data_file_to_work(message):
    return message.body.rstrip().split(',')[1]

def master_loop():
    curr_state = "wait_job_start"
    while True:
        # find out what loop to run
        loop_to_run = get_loop_to_run(curr_state)
        # run the loop and get next state
        curr_state = loop_to_run()

def get_loop_to_run(state_name):
    state_table = {
        'wait_job_start': wait_job_start_loop,
        'create_hash_from_ref': start_create_hash_from_ref,
        'wait_reference_hash': wait_reference_hash_loop,
        'align_base': start_align_base,
        'wait_align_base': wait_align_base_loop,
        'get_mutation': start_get_mutation,
        'wait_get_mutation': wait_get_mutation_loop
    }

    return state_table[state_name]

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
        return 'create_hash_from_ref'
    elif re.search('align_base', message_body):
        return 'align_base'
    elif re.search('get_mutation', message_body):
        return 'get_mutation'


def start_create_hash_from_ref():
    global datafile
    db = database.create_database_connection()

    print "save all reads to database"
    db.execute("DELETE FROM read")
    reads = commonlib.read_all_reads("dataset/%s/reads.txt" % datafile)
    for i in xrange(0, len(reads)):
        db.execute("INSERT INTO read (idx, left_read, right_read ) VALUES (%d, '%s', '%s')" %
            (i, commonlib.get_string_from_mer(reads[i][0]), commonlib.get_string_from_mer(reads[i][1])))

    print "start creating hash from reference"

    # delete all saved reference genome hash
    db.execute("DELETE FROM reference_hash")

    reference_arr = commonlib.read_reference_genome('dataset/%s/ref.txt' % datafile)
    offset = 1000
    for i in xrange(0, len(reference_arr), offset):
        sqs_client.send_message('reference_hash', '%s,%d,%d' % (datafile, i, i+offset))

    return 'wait_reference_hash'

def wait_reference_hash_loop():
    inloop = True
    while inloop:
        remaning_reference_hash = sqs_client.get_queue_remaining_messages('reference_hash')
        if remaning_reference_hash > 0:
            time.sleep(5)
        else:
            inloop = False
    return 'align_base'

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

    return 'wait_align_base'

def wait_align_base_loop():
    inloop = True
    while inloop:
        remaining_align_base = sqs_client.get_queue_remaining_messages('align_base')
        if remaining_align_base > 0:
            time.sleep(5)
        else:
            inloop = False
    return 'get_mutation'

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

    return 'wait_get_mutation'

def wait_get_mutation_loop():
    inloop = True
    while inloop:
        remaining_ref_idx = sqs_client.get_queue_remaining_messages('get_mutation')
        if remaining_ref_idx > 0:
            time.sleep(5)
        else:
            inloop = False
    return 'wait_job_start'

if __name__ == "__main__":
    # run master loop
    master_loop()
