import commonlib
# import database
import pickle

all_ref = {}
mer_length = 10

def save_hash_location(mer_str, location):
    db = database.create_database_connection()
    # try to see if this mer already exists

    query = db.execute("SELECT location FROM reference_hash WHERE mer = '%s'" % mer_str)
    matches = query.fetchone()
    if matches is None:
        db.execute("INSERT INTO reference_hash (mer, location) VALUES ('%s', '{%d}')" % (mer_str, location))
    else:
        # if this is not working then we have duplicate mer, so we update thelist of location
        db.execute("UPDATE reference_hash SET location = array_append(location, %d) WHERE mer = '%s'" % (location, mer_str))


def create_reference_hash(start_idx, stop_idx, mer_length, reference):

    # when stop idx is too large
    if stop_idx + mer_length >= len(reference):
        stop_idx = len(reference) - mer_length

    for i in xrange(start_idx, stop_idx):
        if i % 10000 == 0:
            print "finished %d/%d, percent=%.2f" % (i, stop_idx - start_idx, (i * 100 / (stop_idx-start_idx)))
        mer = reference[i:i+mer_length]
        mer_str = commonlib.get_string_from_mer(mer)
        save_hash_location(mer_str, i)

def create_reference_hash_memory(start_idx, stop_idx, mer_length, reference):

    # when stop idx is too large
    if stop_idx + mer_length >= len(reference):
        stop_idx = len(reference) - mer_length

    all_hash_locaton = {}

    for i in xrange(start_idx, stop_idx):
        if i % 10000 == 0:
            print "finished %d/%d, percent=%.2f" % (i, stop_idx - start_idx, (i * 100 / (stop_idx-start_idx)))
        mer = reference[i:i+mer_length]
        if mer in all_hash_locaton:
            all_hash_locaton[mer].append(i)
        else:
            all_hash_locaton[mer] = [i]

    with open("all_hash_location.pickle", "wb") as f:
        pickle.dump(all_hash_locaton, f, protocol=pickle.HIGHEST_PROTOCOL)

def work_small_job(datafile, start_idx, stop_idx):
    global all_ref
    if datafile not in all_ref:
        all_ref[datafile] = commonlib.read_reference_genome('dataset/%s/ref.txt' % datafile)
    create_reference_hash(start_idx,stop_idx,mer_length, all_ref[datafile])

if __name__ == "__main__":

    if 0:
        # delete all saved reference genome hash
        db = database.create_database_connection()
        db.execute("DELETE FROM reference_hash")

        # read reference genome and populate database
        reference_arr = commonlib.read_reference_genome('dataset/practice2/ref.txt')
        create_reference_hash(0,len(reference_arr),mer_length, reference_arr)

    else:
        # read reference genome and populate database
        reference_arr = commonlib.read_reference_genome_bare('dataset/practice2/ref.txt')
        create_reference_hash_memory(0,len(reference_arr),mer_length, reference_arr)

