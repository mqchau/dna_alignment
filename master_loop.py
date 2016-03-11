import sqs_client
import time
import re
import commonlib
import database
import ipdb

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
    if re.search('align_base', message_body):
        return 'align_base'
    elif re.search('get_mutation', message_body):
        return 'get_mutation'

def start_align_base():
    global datafile

    print "start aligning bases"

    # delete all saved reference genome hash
    db = database.create_database_connection(database=datafile)
    db.execute("DELETE FROM aligned_bases")
    db.execute("DELETE FROM unaligned_reads")

    query_result = db.execute("SELECT COUNT(*) FROM read_raw")
    # ipdb.set_trace()
    # if query_result is not None:
    interval = 50
    for i in xrange(0, int(query_result.fetchone()[0]), interval):
        sqs_client.send_message('align_base', '%s,%d,%d' % (datafile, i, i+interval))
    # else:
    #     raise Exception("No read is loaded into the database yet")

    return 'wait_align_base'

def wait_align_base_loop():
    inloop = True
    while inloop:
        remaining_align_base = sqs_client.get_queue_remaining_messages('align_base')
        if remaining_align_base > 0:
            time.sleep(5)
        else:
            inloop = False
    return 'wait_job_start'
    # return 'get_mutation'

def start_get_mutation():

    global datafile

    print "start getting mutations"

    # delete all saved reference genome hash
    db = database.create_database_connection(database=datafile)
    db.execute("DELETE FROM mutation")

    reference_arr = commonlib.read_reference_genome('dataset/%s/ref.txt' % datafile)
    offset = 1000
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
