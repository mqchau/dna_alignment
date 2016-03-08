import pickle
import database
import commonlib

db = None
reference_genome = None

def save_mutation_to_file():
    all_mutation = get_all_mutation()
    with open('answer.txt', 'w') as f:
        f.write('>practice_E_1_chr_1\n>STR\n>CNV\n>ALU\n>INS\n')
        
        for one_snp in all_mutation["INS"]:
            f.write('%s,%d\n' % (
                one_snp["ins_str"],
                one_snp["ref_idx"]))
        
        f.write('>DEL\n')
        for one_snp in all_mutation["DEL"]:
            f.write('%s,%d\n' % (
                one_snp["del_str"],
                one_snp["ref_idx"]))
        
        f.write('>INV\n>SNP\n')
        for one_snp in all_mutation["SNP"]:
            f.write('%s,%s,%d\n' % (
                one_snp["old_base"],
                one_snp["new_base"],
                one_snp["ref_idx"]))



def get_all_mutation():
    all_mutation = {
        "INS" : get_all_ins(),
        "DEL" : get_all_del(),
        "SNP" : get_all_snp() 
    }
    
    return all_mutation

def get_all_snp():
    query_result = db.execute("SELECT * FROM mutation WHERE mutation_type = 3 ORDER BY ref_idx")
    if query_result is None:
        raise Exception("Error in getting mutation list")

    all_snp = []
    for one_result in query_result.fetchall():
        all_snp.append({
            "ref_idx": one_result["ref_idx"],
            "new_base": one_result["new_base"],
            "old_base": reference_genome[one_result["ref_idx"]]
        })

    # sort by ref idx ascending
    # all_snp = sorted(all_snp, key=lambda x: x["ref_idx"])

    return all_snp

def get_all_ins():
    query_result = db.execute("SELECT * FROM mutation WHERE mutation_type = 2 ORDER BY ref_idx")
    if query_result is None:
        raise Exception("Error in getting mutation list")

    all_ins = []
    for one_result in query_result.fetchall():
        all_ins.append({
            "ref_idx": one_result["ref_idx"],
            "ins_str": one_result["ins_str"]
        })

    return all_ins

def get_all_del():
    query_result = db.execute("SELECT * FROM mutation WHERE mutation_type = 1 ORDER BY ref_idx")
    if query_result is None:
        raise Exception("Error in getting mutation list")

    all_del = []
    curr_del = []
    for one_result in query_result.fetchall():
        if len(curr_del) == 0 or one_result["ref_idx"] == curr_del[-1]["ref_idx"] + 1:
            curr_del.append(one_result)
        else:
            del_str = ""
            for one_del in curr_del:
                del_str += reference_genome[one_del["ref_idx"]]
            all_del.append({
                "ref_idx": curr_del[0]["ref_idx"],
                "del_str": del_str
            })
            curr_del = []

    return all_del

if __name__ == "__main__":
    dataset = 'practice2'

    global db, reference_genome
    db = database.create_database_connection()

    # save the mutation list into answer file format
    with open("dataset/%s/reference_genome.pickle" % dataset, "rb") as f:
        reference_genome = pickle.load(f)
    save_mutation_to_file()
