import commonlib
import ipdb
import os
import redis
import numpy as np
import cluster_sub_read_match

# define how to score in local alignment process
mismatch = -1
match = 2
indel = -1
dataset = "10k"
redis_db=["10k", "1m", "100m"].index(dataset)
redis_path = os.environ["REDIS_HOST"] if "REDIS_HOST" in os.environ else "localhost"
redis_client = redis.StrictRedis(host=redis_path, port=6379, db=redis_db)

def get_reference_genome_one_base(base_idx):
    base = redis_client.lrange("reference-genome", base_idx, base_idx)[0]
    return b.decode("utf8")

# left_idx is inclusive but right_idx is non-inclusive
def get_reference_genome_many_base(left_idx, right_idx):
    return "".join([x.decode("utf8") for x in redis_client.lrange("reference-genome", left_idx, right_idx-1)])

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
            ref_string += get_reference_genome_one_base(one_base["ref_idx"])
            read_string += one_base["base"].lower()
            # print "%05d\t%d\t%d\tmismatch" % (one_base["ref_idx"], reference_genome[one_base["ref_idx"]], one_base["base"])
        elif one_base["type"] == "delete":
            ref_string += get_reference_genome_one_base(one_base["ref_idx"])
            read_string += "-"
            # print "%05d\t%d\t-" % (one_base["ref_idx"], reference_genome[one_base["ref_idx"]])
    print("Refer: %s\nDonor: %s" % (ref_string, read_string))


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
    for read_base_idx in range(1, len(read)+1):
        for ref_base_idx in range(1, len(ref)+1):
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

    query_result = db.execute("SELECT * FROM read_raw WHERE idx >= %d AND idx < %d" % (start_idx, stop_idx))
    if query_result is None:
        return []

    return map(lambda x: [ x['left_read'], x['right_read'] ], 
               query_result.fetchall())

        
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

def align_read_new(read_pair_raw):

    read_pair = read_pair_raw.rstrip().split(",")

    left_read = dna_read(read_pair[0])
    right_read = dna_read(read_pair[1])

    # try align left in order read first
    left_possible_locations = left_read.find_possible_alignment()
    right_possible_locations = right_read.find_possible_alignment()

    cluster_pairs = select_good_location_pair(left_possible_locations, right_possible_locations)
    if len(cluster_pairs) == 0:
        return []

    for one_cluter_pair in cluster_pairs:
        left_start_idx = one_cluter_pair["left_cluster"].get_cluster_expected_start() - 10
        right_start_idx = one_cluter_pair["right_cluster"].get_cluster_expected_start() - 10
        # ipdb.set_trace()
        if one_cluter_pair["type"] == "f":
            max_score, left_mutation = align_read_by_local_alignment(get_reference_genome_many_base(left_start_idx,left_start_idx+70), left_read.read, left_start_idx)
            max_score, right_mutation = align_read_by_local_alignment(get_reference_genome_many_base(right_start_idx,right_start_idx+70), right_read.get_reversed(), right_start_idx)
        else:
            max_score, left_mutation = align_read_by_local_alignment(get_reference_genome_many_base(left_start_idx,left_start_idx+70), left_read.get_reversed(), left_start_idx)
            max_score, right_mutation = align_read_by_local_alignment(get_reference_genome_many_base(right_start_idx,right_start_idx+70), right_read.read, right_start_idx)

    return [left_mutation, right_mutation]


class dna_read:
    def __init__(self, read):
        self.read = read

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

        return {

            "forward": cluster_sub_read_match.get_clusters_from_matches(matches["forward"]),
            "backward": cluster_sub_read_match.get_clusters_from_matches(matches["backward"])
        }


    def get_sub_read(self, read):
        section_length = 10
        sub_sections = []
        move_length = 10
        for i in range(0, len(read) - section_length+1, move_length):
            # ipdb.set_trace()
            sub_sections.append(sub_read(read[i:i+section_length], i))

        return sub_sections

class sub_read:
    def __init__(self, sub_read, sub_read_idx):
        self.sub_read = sub_read
        self.sub_read_idx = sub_read_idx

    def get_match_hash(self):
        self.match_hash = redis_client.lrange(self.sub_read, 0, -1)
        return self.match_hash

if __name__ == "__main__":
    aligned = align_read_new("TTCCTCCGCTTCTTGTGGCTTCCCAGCACGCAACAAGGAAACACCAGCAC,AACATAGTAGGAAACACCGTATTGCCGAGTTGAATCTCTGATCTTGCAAT")
