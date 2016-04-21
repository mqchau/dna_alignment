import pprint
import ipdb
import commonlib
import numpy as np

db = None
pp = pprint.PrettyPrinter(indent=4)

def get_most_common_new_base(snp_list):
    base_count = { "A": 0, "C": 0, "G": 0, "T": 0 }
    for one_snp in snp_list:
        inserted_base = one_snp.split(',')[1]
        base_count[inserted_base] += 1
    max_base_count = max(base_count.values())
    for base in sorted(base_count.keys()):
        if base_count[base] == max_base_count:
            return base

def get_ins_str(insert_list):

    ins_idx = 1
    break_cond = False
    ins_str = ""
    while not break_cond:
        insert_at_this_idx = list(filter(lambda x: int(x.split(",")[2]) == ins_idx, insert_list))
        ipdb.set_trace()
        if len(insert_at_this_idx) > 0:
            ins_str += str(get_most_common_new_base(insert_at_this_idx))
            ins_idx +=1
        else:
            break_cond = True
    return ins_str

def pile_up(ref_idx, mutations_at_ref):

    if len(mutations_at_ref) < 6:
        return ("del",)
    elif len(mutations_at_ref) > 0:
        match_count, mismatch_count, del_count, ins_count = 0,0,0,0
        for one_mut in mutations_at_ref:
            splitted_str = one_mut.split(',')
            mutation_type = splitted_str[0]
            if mutation_type == "match":
                match_count += 1
            elif mutation_type == "mismatch":
                mismatch_count += 1
            elif mutation_type == "insert":
                ins_count += 1
            elif mutation_type == "delete":
                del_count += 1

        max_count = np.max((match_count, mismatch_count, del_count, ins_count))

        if match_count == max_count:
            # we think this is match, no need to report
            pass
        elif mismatch_count == max_count:
            # we think this is a mismatch
            return ("snp", get_most_common_new_base(list(filter(lambda x: x.split(",")[0] == "mismatch", mutations_at_ref))))
        elif del_count == max_count:
            return ("del",)
        elif ins_count == max_count:
            return ("ins", get_ins_str(list(filter(lambda x: x.split(",")[0] == "insert", mutations_at_ref))))

if __name__ == "__main__":
    mutation = pile_up(1000, [
        "match,A",
        "insert,T,2",
        "insert,T,2",
        "insert,T,2",
        "insert,T,2",
        "insert,T,2",
        "insert,G,1",
        "insert,G,1",
        "insert,G,1",
        "insert,G,1",
        "insert,G,1",
    ])
    pp.pprint(mutation)
