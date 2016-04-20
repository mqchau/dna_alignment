import commonlib
import redis
import os

redis_path = os.environ["REDIS_HOST"] if "REDIS_HOST" in os.environ else "localhost"
mer_length = 10
dataset = "10k" # 10k, 1m or 100m

def create_reference_hash_save_redis(mer_length, reference, redis_path, redis_db=0):
    
    # create redis client
    r = redis.StrictRedis(host=redis_path, port=6379, db=redis_db)
    # delete all keys, all hashed locations
    r.flushall()

    # iterate through all the length of the reference
    for i in range(len(reference) - mer_length):
        if i % 100000 == 0:
            # some printf for progress
            print("finished %d/%d, percent=%.2f" % (i, len(reference), (i * 100 / len(reference))))
        # get the mer string
        mer = reference[i:i+mer_length]
        
        # push to redis in a list with the key as mer length
        r.rpush(mer, i)

if __name__ == "__main__":
    redis_db=["10k", "1m", "100m"].index(dataset)
    reference_arr = commonlib.read_reference_genome_bare('dataset/%s/ref.txt' % dataset)
    create_reference_hash_save_redis(mer_length, reference_arr, redis_path, redis_db=redis_db)

