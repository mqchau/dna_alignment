import commonlib
import database

def align_one_read(read_pair):
    pass

if __name__ == "__main__":
    # delete all saved aligned reads
    db = database.create_database_connection()
    db.execute("DELETE FROM aligned_bases")

    # read in all reads
    reads = commonlib.read_all_reads("dataset/practice1/reads.txt")

    # align each reads
    for one_read in reads[0:1]:
        align_one_read(one_read)
