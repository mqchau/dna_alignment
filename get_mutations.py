import pprint
import ipdb
import database
import commonlib
import numpy as np

db = None
reference_genome = None

def get_most_common_new_base(snp_list):
    base_count = np.zeros(4, dtype=np.uint32)
    for one_snp in snp_list:
        base_count[one_snp["new_base"]-1] += 1
    max_base = np.argmax(base_count)
    return max_base

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

def get_mutations():
    all_mutations = {
        "snp": [],
        "ins": [],
        "del": []
    }

    for ref_idx in xrange(len(reference_genome)):
        # from database get list of aligned bases related to this ref_idx
        query = db.execute("SELECT * FROM aligned_bases WHERE ref_idx = %d" % ref_idx)

        mutations_at_ref = query.fetchall()
        if len(mutations_at_ref) > 0:
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

    return all_mutations
            
if __name__ == "__main__":
    global db, reference_genome
    db = database.create_database_connection()

    # read reference genome
    reference_genome = commonlib.read_reference_genome('dataset/practice1/ref.txt')

    # figure out the mutations at each base
    all_mutations = get_mutations()


    pp = pprint.PrettyPrinter(indent=2)
    pp.pprint(all_mutations)
