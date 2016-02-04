import database
import commonlib

db = None
reference_arr = None

def save_mutation_to_file():
    all_mutation = get_all_mutation()
    with open('answer.txt', 'w') as f:
        f.write('>practice_W_3_chr_1\n>STR\n>CNV\n>ALU\n>INS\n>DEL\n>INV\n')
        f.write('>SNP\n')
        for one_snp in all_mutation["SNP"]:
            f.write('%s,%s,%d\n' % (
                commonlib.get_base_from_int(one_snp["old_base"]),
                commonlib.get_base_from_int(one_snp["new_base"]),
                one_snp["ref_idx"]))



def get_all_mutation():
    all_mutation = {
        "INS" : get_all_ins(),
        "DEL" : get_all_del(),
        "SNP" : get_all_snp() 
    }
    
    return all_mutation

def get_all_snp():
    query_result = db.execute("SELECT * FROM mutation WHERE mutation_type = 3")
    if query_result is None:
        raise Exception("Error in getting mutation list")

    all_snp = []
    for one_result in query_result.fetchall():
        all_snp.append({
            "ref_idx": one_result["ref_idx"],
            "new_base": one_result["new_base"],
            "old_base": reference_arr[one_result["ref_idx"]]
        })

    # sort by ref idx ascending
    all_snp = sorted(all_snp, key=lambda x: x["ref_idx"])

    return all_snp

def get_all_ins():
    return []

def get_all_del():
    return []

if __name__ == "__main__":
    global db, reference_arr
    db = database.create_database_connection()

    # save the mutation list into answer file format
    reference_arr = commonlib.read_reference_genome('dataset/practice1/ref.txt')
    save_mutation_to_file()
