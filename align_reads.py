import pickle
import commonlib
import database
import ipdb
import numpy as np
import pickle
import cluster_sub_read_match

db = None
reference_genome = None
all_ref = {}

# define how to score in local alignment process
mismatch = -1
match = 2
indel = -1

def get_base_char_from_int(num):

    base_lookup = ["A", "C", "G", "T"]
    # print num
    return base_lookup[int(num)-1]

def debug_print_mutations(mutation_list):
    ref_string = ""
    read_string = ""
    for one_base in mutation_list:
        if one_base["type"] == "insert":
            ref_string += "-"
            read_string += one_base["base"]
            # print "%05d\t-\t%d" % (one_base["ref_idx"], one_base["base"])
        elif one_base["type"] == "match":
            ref_string += one_base["base"]
            read_string += one_base["base"]
            # print "%05d\t%d\t%d" % (one_base["ref_idx"], one_base["base"], one_base["base"])
        elif one_base["type"] == "mismatch":
            ref_string += reference_genome[one_base["ref_idx"]]
            read_string += one_base["base"].lower()
            # print "%05d\t%d\t%d\tmismatch" % (one_base["ref_idx"], reference_genome[one_base["ref_idx"]], one_base["base"])
        elif one_base["type"] == "delete":
            ref_string += reference_genome[one_base["ref_idx"]]
            read_string += "-"
            # print "%05d\t%d\t-" % (one_base["ref_idx"], reference_genome[one_base["ref_idx"]])
    print "Refer: %s\nDonor: %s" % (ref_string, read_string)


def save_mutation_to_db(mutation_pair):
    for mutation_list in mutation_pair:
        for one_base in mutation_list:
            if one_base["type"] == "insert":
                db.execute("INSERT INTO aligned_bases (ref_idx, mutation_type, insert_idx, new_base) VALUES (%d, %d, %d, %d)" % (one_base["ref_idx"], 2, one_base["insert_idx"], one_base["base"]))
            elif one_base["type"] == "match":
                db.execute("INSERT INTO aligned_bases (ref_idx, mutation_type, new_base) VALUES (%d, %d, %d)" % (one_base["ref_idx"], 4, one_base["base"]))
            elif one_base["type"] == "mismatch":
                db.execute("INSERT INTO aligned_bases (ref_idx, mutation_type, new_base) VALUES (%d, %d, %d)" % (one_base["ref_idx"], 3, one_base["base"]))
            elif one_base["type"] == "delete":
                db.execute("INSERT INTO aligned_bases (ref_idx, mutation_type) VALUES (%d, %d)" % (one_base["ref_idx"], 1))

def align_read_by_local_alignment(ref, read, ref_start_idx):
    global indel, match, mismatch
    
    # create 3d array to store the local alignment score and direction
    # 1st dim: bases on read
    # 2nd dim: bases on ref
    # 3rd dim: (score, top flag, diagonal flag, left flag)
    score_table = np.zeros((len(read) + 1, len(ref) + 1, 4))

    # start calculating from the first bases of ref and read
    max_global_score = -9999
    max_global_location = tuple()
    for read_base_idx in xrange(1, len(read)+1):
        for ref_base_idx in xrange(1, len(ref)+1):
            score_from_top = score_table[read_base_idx-1, ref_base_idx, 0] + indel 
            score_from_left = score_table[read_base_idx, ref_base_idx-1, 0] + indel 
            diag_delta = match if read[read_base_idx-1] == ref[ref_base_idx-1] else mismatch
            score_from_diag = score_table[read_base_idx-1, ref_base_idx-1, 0] + diag_delta

            max_score = np.max((score_from_top, score_from_left, score_from_diag))
            score_table[read_base_idx, ref_base_idx, 0] = max_score

            if score_from_diag == max_score:
                score_table[read_base_idx, ref_base_idx, 2] = 1
            elif score_from_left == max_score:
                score_table[read_base_idx, ref_base_idx, 3] = 1
            elif score_from_top == max_score:
                score_table[read_base_idx, ref_base_idx, 1] = 1
            else:
                raise Exception("Can't find direction to trace back in local alignment")

            if max_score > max_global_score:
                max_global_score = max_score
                max_global_location = (read_base_idx, ref_base_idx)

    # start tracing back
    if len(max_global_location) == 0:
        return 0, []

    start_location = max_global_location
    mutations = []
    while start_location is not None and start_location[0] > 0 and start_location[1] > 0:
        # insert on read, move to top
        if score_table[start_location[0], start_location[1], 1] == 1:
            mutations.append({
                "type": "insert",
                "ref_idx": start_location[1] - 1 + ref_start_idx,
                "base": read[start_location[0]-1]
            })
            start_location = (start_location[0] - 1, start_location[1])
        # match/mismatch, move to diag
        elif score_table[start_location[0], start_location[1], 2] == 1:
            mutations.append({
                "type": "mismatch" if read[start_location[0]-1] != ref[start_location[1]-1] else "match",
                "ref_idx": start_location[1] - 1 + ref_start_idx,
                "base": read[start_location[0]-1],
                "read_idx": start_location[0]-1
            })
            start_location = (start_location[0] - 1, start_location[1] - 1)
        # delete on read, move to left
        elif score_table[start_location[0], start_location[1], 3] == 1:
            mutations.append({
                "type": "delete",
                "ref_idx": start_location[1] - 1 + ref_start_idx
            })
            start_location = (start_location[0], start_location[1]-1)
        else:
            start_location = None

    # reverse the mutation list to get correct order in order of reference genome
    mutations.reverse()

    # scan through the insert and put in the offset of insert in respect to the base of reference it inserts to
    for mut_idx, one_mutation in enumerate(mutations):
        if one_mutation["type"] == "insert":
            if "insert_idx" in mutations[mut_idx-1]:
                mutations[mut_idx]["insert_idx"] = mutations[mut_idx-1]["insert_idx"] + 1
            else:
                mutations[mut_idx]["insert_idx"] = 1

    return max_global_score, mutations

def get_all_reads_to_work(start_idx, stop_idx):
    global db

    query_result = db.execute("SELECT * FROM read WHERE idx >= %d AND idx < %d" % (start_idx, stop_idx))
    if query_result is None:
        return []

    return map(lambda x: [
        commonlib.get_mer_from_int_str(x['left_read']), 
        commonlib.get_mer_from_int_str(x['right_read'])
    ], query_result.fetchall())

        
def get_read_pair_by_idx(idx):
    query_result = db.execute("SELECT * FROM read_raw WHERE idx = %d" % read_idx).fetchone()
    
    if query_result is None:
        return None

    return [query_result.left_read, query_result.right_read]

def select_good_location_pair(left_clusters, right_clusters):
    good_pairs = []
    for one_left in left_clusters["forward"]:
        for one_right in right_clusters["backward"]:
            pair_distance = one_right.get_cluster_expected_start() - one_left.get_cluster_expected_start()
            if pair_distance > 100 and pair_distance < 200:
                good_pairs.append({
                    "type": "f",
                    "left_cluster": one_left,
                    "right_cluster": one_right
                })
    for one_left in left_clusters["backward"]:
        for one_right in right_clusters["forward"]:
            pair_distance = one_right.get_cluster_expected_start() - one_left.get_cluster_expected_start()
            if pair_distance > 100 and pair_distance < 200:
                good_pairs.append({
                    "type": "b",
                    "left_cluster": one_left,
                    "right_cluster": one_right
                })

    return good_pairs

def align_read_new(reference_genome, reference_hash, read_pair):

    left_read = dna_read(read_pair[0], reference_hash)
    right_read = dna_read(read_pair[1], reference_hash)

    # try align left in order read first
    left_possible_locations = left_read.find_possible_alignment()
    right_possible_locations = right_read.find_possible_alignment()

    cluster_pairs = select_good_location_pair(left_possible_locations, right_possible_locations)
    print cluster_pairs

    for one_cluter_pair in cluster_pairs:
        left_start_idx = one_cluter_pair["left_cluster"].get_cluster_expected_start() - 10
        right_start_idx = one_cluter_pair["right_cluster"].get_cluster_expected_start() - 10
        # ipdb.set_trace()
        if one_cluter_pair["type"] == "f":
            max_score, left_mutation = align_read_by_local_alignment(reference_genome[left_start_idx:left_start_idx+70], left_read.read, left_start_idx)
            max_score, right_mutation = align_read_by_local_alignment(reference_genome[right_start_idx:right_start_idx+70], right_read.get_reversed(), right_start_idx)
        else:
            max_score, left_mutation = align_read_by_local_alignment(reference_genome[left_start_idx:left_start_idx+70], left_read.get_reversed(), left_start_idx)
            max_score, right_mutation = align_read_by_local_alignment(reference_genome[right_start_idx:right_start_idx+70], right_read.read, right_start_idx)

        # ipdb.set_trace()

        debug_print_mutations(left_mutation)
        debug_print_mutations(right_mutation)


class dna_read:
    def __init__(self, read, reference_hash):
        self.read = read
        self.reference_hash = reference_hash

    def get_reversed(self):
        return self.read[::-1]

    def find_possible_alignment(self):

        matches = {
            "forward":[],
            "backward":[]
        }

        # try to align in forward direction first
        sub_reads = self.get_sub_read(self.read)
        for one_sub_read in sub_reads:
            one_sub_read.get_match_hash()
            matches["forward"].append(one_sub_read.match_hash)
            # print one_sub_read.match_hash
            # print "********************"
            # ipdb.set_trace()

        # print "----------------------------------------"

        # try to align in reversed direction later
        sub_reads = self.get_sub_read(self.read[::-1])
        for one_sub_read in sub_reads:
            one_sub_read.get_match_hash()
            matches["backward"].append(one_sub_read.match_hash)
            # print one_sub_read.match_hash
            # print "********************"

        # with open("read_matches.pickle", "wb") as f:
        #     pickle.dump(matches, f, protocol=pickle.HIGHEST_PROTOCOL)


        return {
            "forward": cluster_sub_read_match.get_clusters_from_matches(matches["forward"]),
            "backward": cluster_sub_read_match.get_clusters_from_matches(matches["backward"])
        }


    def get_sub_read(self, read):
        section_length = 10
        sub_sections = []
        move_length = 10
        for i in xrange(0, len(read) - section_length+1, move_length):
            # ipdb.set_trace()
            sub_sections.append(sub_read(read[i:i+section_length], i, self.reference_hash))

        return sub_sections

class sub_read:
    def __init__(self, sub_read, sub_read_idx, reference_hash):
        self.sub_read = sub_read
        self.sub_read_idx = sub_read_idx
        self.reference_hash = reference_hash

    def get_match_hash(self):
        if self.sub_read in self.reference_hash:
            self.match_hash = self.reference_hash[self.sub_read]
        else:
            # print "sub read %s not in hash" % self.sub_read
            self.match_hash = []
        return self.match_hash


def work_small_job(datafile, start_idx, stop_idx):
    global all_ref, db, reference_genome
    db = database.create_database_connection()
    if datafile not in all_ref:
        all_ref[datafile] = commonlib.read_reference_genome('dataset/%s/ref.txt' % datafile)

    reference_genome = all_ref[datafile]

    all_reads_to_work = get_all_reads_to_work(start_idx, stop_idx)

    for relative_read_idx, one_read in enumerate(all_reads_to_work):
        mutation_list = align_one_read_pair(one_read)
        if len(mutation_list) == 0:
            print "can't align read %d" % (start_idx + relative_read_idx)
        save_mutation_to_db(mutation_list)

if __name__ == "__main__":
    # delete all saved aligned reads
    global db, reference_genome
    db = database.create_database_connection()
    # db.execute("DELETE FROM aligned_bases")

    # read reference genome
    # reference_genome = commonlib.read_reference_genome_bare('dataset/practice2/ref.txt')
    with open("reference_genome.pickle", "rb") as f:
        reference_genome = pickle.load(f)

    # read the reference hash
    with open("all_hash_location.pickle", "rb") as f:
        reference_hash = pickle.load(f)

    # align one read pair
    for read_idx in xrange(0,5):
        read_pair = get_read_pair_by_idx(read_idx)
        align_read_new(reference_genome, reference_hash, read_pair)
        print "--------------------"


    # match_score, mutations = align_read_by_local_alignment(reference_genome, reads[2][1], 0)
    
    # work_small_job("practice", 0, 10 )


    # print align_one_read_by_hash(reads[2][0])
    # print align_one_read_by_hash(reads[2][1])

    # debug_print_mutations(mutations)

    # # align each reads
    # for read_idx, one_read in enumerate(reads[2:3]):
    #     mutation_list = align_one_read_pair(one_read)
    #     save_mutation_to_db(mutation_list)
    #     # if len(mutation_list) > 0:
    #     if len(mutation_list) == 0:
    #         print "can't align read %d" % read_idx
    #     if read_idx % 50 == 0:
    #         print "-----------------------%05d---------------------------" % read_idx
    #     #     debug_print_mutations(mutation_list[0])
    #     #     debug_print_mutations(mutation_list[1])
