import commonlib


def save_hash_location(mer_str, location):
    pass


def create_reference_hash(start_idx, stop_idx, mer_length, reference):
    
    # when stop idx is too large
    if stop_idx + mer_length >= len(reference):
        stop_idx = reference - mer_length

    for i in xrange(start_idx, stop_idx):
        mer = reference[i:i+mer_length]
        mer_str = commonlib.get_string_from_mer(mer)
        save_hash_location(mer_str, i)
        print mer_str

if __name__ == "__main__":
    reference_arr = commonlib.read_reference_genome('dataset/practice1/ref.txt')
    create_reference_hash(0,2,16, reference_arr)
