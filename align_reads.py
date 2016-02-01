import commonlib
import database
import ipdb
import numpy as np

db = None

def align_one_read_pair(read_pair):
    # align first read by hashing
    possible_match_left_read = align_one_read_by_hash(read_pair[0])

    # ipdb.set_trace()


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
            query = db.execute("SELECT location FROM reference_hash WHERE mer = '%s'" % commonlib.get_string_from_mer(one_sub))
            matches = query.fetchone()
            if matches is not None:
                # found a match
                for one_match in matches[0]:
                    # to account for indel
                    # we pad a 10 bases at the front and end of the match
                    # later on we'll use local alignment
                    approx_start = one_match - 16 * sub_idx - 10
                    approx_end = approx_start + 80
                    if approx_start not in match_locations_set: 
                        match_locations.append({
                            "direction": direction ,
                            "approx_start": approx_start,
                            "approx_end": approx_end
                        })
                        match_locations_set.add(approx_start)

    return match_locations

if __name__ == "__main__":
    # delete all saved aligned reads
    global db
    db = database.create_database_connection()
    db.execute("DELETE FROM aligned_bases")

    # read in all reads
    reads = commonlib.read_all_reads("dataset/practice1/reads.txt")

    # align each reads
    for one_read in reads[0:]:
        align_one_read_pair(one_read)
