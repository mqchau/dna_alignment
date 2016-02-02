import commonlib
import database
import ipdb
import numpy as np

db = None
reference_genome = None

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
            read_string += get_base_char_from_int(one_base["base"])
            # print "%05d\t-\t%d" % (one_base["ref_idx"], one_base["base"])
        elif one_base["type"] == "match":
            ref_string += get_base_char_from_int(one_base["base"])
            read_string += get_base_char_from_int(one_base["base"])
            # print "%05d\t%d\t%d" % (one_base["ref_idx"], one_base["base"], one_base["base"])
        elif one_base["type"] == "mismatch":
            ref_string += get_base_char_from_int(reference_genome[one_base["ref_idx"]])
            read_string += get_base_char_from_int(one_base["base"]).lower()
            # print "%05d\t%d\t%d\tmismatch" % (one_base["ref_idx"], reference_genome[one_base["ref_idx"]], one_base["base"])
        elif one_base["type"] == "delete":
            ref_string += get_base_char_from_int(reference_genome[one_base["ref_idx"]])
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


def align_one_read_pair(read_pair):
    # align two reads by hashing
    possible_matches = align_read_pair_by_hash(read_pair)

    # for each possible location
    max_match_score = 0
    max_mutation_list = []
    for one_possible_match in possible_matches:
        # align first read again by local alignment
        oriented_read = read_pair[0] if one_possible_match[0]["direction"] == "forward" else np.flipud(read_pair[0])
        match_score_left, left_mutations = align_read_by_local_alignment(
            reference_genome[
                one_possible_match[0]["approx_start"]: 
                one_possible_match[0]["approx_end"]], 
            oriented_read, one_possible_match[0]["approx_start"])

        # align 2nd read by local alignment
        oriented_read = read_pair[1] if one_possible_match[1]["direction"] == "forward" else np.flipud(read_pair[1])
        match_score_right, right_mutations = align_read_by_local_alignment(
            reference_genome[
                one_possible_match[1]["approx_start"]: 
                one_possible_match[1]["approx_end"]], 
            oriented_read, one_possible_match[1]["approx_start"])
        
        if match_score_left + match_score_right > max_match_score:
            max_match_score = match_score_left + match_score_right
            max_mutation_list = [ left_mutations, right_mutations ]

    return max_mutation_list

def align_read_pair_by_hash(read_pair):
    hash_match_left = align_one_read_by_hash(read_pair[0])
    hash_match_right = align_one_read_by_hash(read_pair[1])

    match_result = []
    for left_match in hash_match_left:
        for right_match in hash_match_right:
            match_distance = abs(left_match["approx_start"] - right_match["approx_start"]) 
            if match_distance > 70+50 and match_distance < 130 + 50:
                match_result.append([ left_match, right_match ])

    # if len(match_result) == 0:
    #     ipdb.set_trace()

    return match_result

def align_one_read_by_hash(read):
    match_locations = []
    match_locations_set = set()
    for direction in ["forward", "backward"]:
        if direction == "backward":
            read = np.flipud(read)

        # cut read into 16 - 16 - 16
        sub_sections = [read[0:16], read[16:32], read[32:48]]

        for sub_idx, one_sub in enumerate(sub_sections):
            # find if this subsection is in reference hash
            # print commonlib.get_string_from_mer(one_sub)
            query = db.execute("SELECT location FROM reference_hash WHERE mer = '%s'" % commonlib.get_string_from_mer(one_sub))
            matches = query.fetchone()
            if matches is not None:
                # found a match
                for one_match in matches[0]:
                    # to account for indel
                    # we pad a 10 bases at the front and end of the match
                    # later on we'll use local alignment
                    approx_start = one_match - 16 * sub_idx - 10
                    if approx_start < 0:
                        approx_start = 0
                    approx_end = approx_start + 80
                    if approx_end > len(reference_genome):
                        approx_end = len(reference_genome)
                    if approx_start not in match_locations_set: 
                        match_locations.append({
                            "direction": direction ,
                            "approx_start": approx_start,
                            "approx_end": approx_end,
                        })
                        match_locations_set.add(approx_start)

    # if len(match_locations) == 0:
    #     raise Exception("can't find any sub string in reference by hash")

    return match_locations

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
        ipdb.set_trace()
        raise Exception("Couldn't find any cell with max score to start tracing back")
    
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

if __name__ == "__main__":
    # delete all saved aligned reads
    global db
    db = database.create_database_connection()
    db.execute("DELETE FROM aligned_bases")

    # read reference genome
    reference_genome = commonlib.read_reference_genome('dataset/practice1/ref.txt')

    # read in all reads
    reads = commonlib.read_all_reads("dataset/practice1/reads.txt")

    # align each reads
    for read_idx, one_read in enumerate(reads):
        mutation_list = align_one_read_pair(one_read)
        save_mutation_to_db(mutation_list)
        # if len(mutation_list) > 0:
        if read_idx % 50 == 0:
            print "-----------------------%05d---------------------------" % read_idx
        #     debug_print_mutations(mutation_list[0])
        #     debug_print_mutations(mutation_list[1])
