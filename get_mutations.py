import pprint
import ipdb
import database
import commonlib
import numpy as np

db = None

def get_most_common_new_base(snp_list):
    base_count = { "A": 0, "C": 0, "G": 0, "T": 0 }
    for one_snp in snp_list:
        base_count[one_snp["new_base"]] += 1
    max_base_count = max(base_count.values())
    for base in sorted(base_count.keys()):
        if base_count[base] == max_base_count:
            return base

def get_ins_str(insert_list):

    ins_idx = 1
    break_cond = False
    ins_str = ""
    while not break_cond:
        insert_at_this_idx = filter(lambda x: x["insert_idx"] == ins_idx, insert_list)
        if len(insert_at_this_idx) > 0:
            ins_str += str(get_most_common_new_base(insert_at_this_idx))
            ins_idx +=1
        else:
            break_cond = True
    return ins_str

def get_mutations(ref_start_idx, ref_end_idx, read_idx_limit):
    all_mutations = {
        "snp": [],
        "ins": [],
        "del": []
    }

    for ref_idx in xrange(ref_start_idx, ref_end_idx):
        # from database get list of aligned bases related to this ref_idx
        query = db.execute("SELECT * FROM aligned_bases WHERE ref_idx = %d AND read_idx < %d" % (ref_idx, read_idx_limit))

        mutations_at_ref = query.fetchall()

        if len(mutations_at_ref) < 6 and ref_idx > 50 and ref_idx < 1000000 - 50:
            all_mutations["del"].append({
                "ref_idx": ref_idx
            })
        elif len(mutations_at_ref) > 0:
            match_count, mismatch_count, del_count, ins_count = 0,0,0,0
            for one_mut in mutations_at_ref:
                if one_mut["mutation_type"] == 4:
                    match_count += 1
                elif one_mut["mutation_type"] == 3:
                    mismatch_count += 1
                elif one_mut["mutation_type"] == 2:
                    ins_count += 1
                elif one_mut["mutation_type"] == 1:
                    del_count += 1

            max_count = np.max((match_count, mismatch_count, del_count, ins_count))

            if match_count == max_count:
                # we think this is match, no need to report
                pass
            elif mismatch_count == max_count:
                # we think this is a mismatch
                all_mutations["snp"].append({
                    "ref_idx": ref_idx,
                    "new_base": get_most_common_new_base(filter(lambda x: x["mutation_type"] == 3, mutations_at_ref))
                })
            elif del_count == max_count:
                all_mutations["del"].append({
                    "ref_idx": ref_idx
                })
            elif ins_count == max_count:
                all_mutations["ins"].append({
                    "ref_idx": ref_idx,
                    "ins_str": get_ins_str(filter(lambda x: x["mutation_type"] == 2, mutations_at_ref))
                })

    # save these mutations to database
    save_all_mutations(all_mutations)
    return all_mutations
            
def save_all_mutations(mutations):
    for one_mut in mutations["snp"]:
        db.execute("INSERT INTO mutation (ref_idx, mutation_type, new_base) VALUES (%d, %d, '%s')" % (one_mut["ref_idx"], 3, one_mut["new_base"]))

    for one_mut in mutations["ins"]:
        db.execute("INSERT INTO mutation (ref_idx, mutation_type, ins_str) VALUES (%d, %d, '%s')" % (one_mut["ref_idx"], 2, one_mut["ins_str"]))

    for one_mut in mutations["del"]:
        db.execute("INSERT INTO mutation (ref_idx, mutation_type) VALUES (%d, %d)" % (one_mut["ref_idx"], 1))

def work_small_job(datafile, start_idx, stop_idx, read_idx_limit):
    global db
    db = database.create_database_connection(database=datafile)

    get_mutations(start_idx, stop_idx, read_idx_limit)


if __name__ == "__main__":
    global db
    db = database.create_database_connection()

    # delete all saved mutations from database
    db.execute("DELETE FROM mutation")

    # figure out the mutations at each base
    all_mutations = get_mutations(216000, 217000, 100000)

    pp = pprint.PrettyPrinter(indent=2)
    pp.pprint(all_mutations)
