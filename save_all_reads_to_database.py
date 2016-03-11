import database

if __name__ == "__main__":

    dataset_name = "ugrad"
    with open("dataset/%s/reads.txt" % dataset_name, "r") as f:
        all_lines = f.readlines()

    del all_lines[0]

   
    db = database.create_database_connection(database=dataset_name)
    for i in xrange(len(all_lines)):
        splitted = all_lines[i].rstrip().split(',')
        left_read = splitted[0]
        right_read = splitted[1]
        db.execute("INSERT INTO read_raw (idx, left_read, right_read ) VALUES (%d, '%s', '%s')" %
            (i, left_read, right_read))


    dataset_name = "practice2"
    with open("dataset/%s/reads.txt" % dataset_name, "r") as f:
        all_lines = f.readlines()

    del all_lines[0]

   
    db = database.create_database_connection(database=dataset_name)
    for i in xrange(len(all_lines)):
        splitted = all_lines[i].rstrip().split(',')
        left_read = splitted[0]
        right_read = splitted[1]
        db.execute("INSERT INTO read_raw (idx, left_read, right_read ) VALUES (%d, '%s', '%s')" %
            (i, left_read, right_read))
 
