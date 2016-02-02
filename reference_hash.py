import commonlib
import database

all_ref = {}

def save_hash_location(mer_str, location):
    db = database.create_database_connection()
    try:
        # try to see if this mer already exists
        db.execute("INSERT INTO reference_hash (mer, location) VALUES ('%s', '{%d}')" % (mer_str, location))
    except Exception as e:
        # if this is not working then we have duplicate mer, so we update thelist of location
        db.execute("UPDATE reference_hash SET location = array_append(location, %d) WHERE mer = '%s'" % (location, mer_str))


def create_reference_hash(start_idx, stop_idx, mer_length, reference):

    # when stop idx is too large
    if stop_idx + mer_length >= len(reference):
        stop_idx = len(reference) - mer_length

    for i in xrange(start_idx, stop_idx):
        mer = reference[i:i+mer_length]
        mer_str = commonlib.get_string_from_mer(mer)
        save_hash_location(mer_str, i)

def work_small_job(datafile, start_idx, stop_idx):
    global all_ref
    if datafile not in all_ref:
        all_ref[datafile] = commonlib.read_reference_genome('dataset/%s/ref.txt' % datafile)
    create_reference_hash(start_idx,stop_idx,16, all_ref[datafile])

if __name__ == "__main__":
    # delete all saved reference genome hash
    db = database.create_database_connection()
    db.execute("DELETE FROM reference_hash")

    # read reference genome and populate database
    reference_arr = commonlib.read_reference_genome('dataset/practice1/ref.txt')
    create_reference_hash(0,len(reference_arr),16, reference_arr)
