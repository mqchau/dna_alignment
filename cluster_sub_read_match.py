import pickle
import numpy as np
import ipdb

class cluster:
    def __init__(self, initial_location):
        self.location = [initial_location]
        pass

    def get_cluster_size(self):
        return len(self.location)

    def is_in_range(self, new_loc):
        if new_loc.location - new_loc.expected_loc_idx > self.get_cluster_expected_start() - 10 and new_loc.location - new_loc.expected_loc_idx < self.get_cluster_expected_start() + 10:
            return True

        return False

    def get_cluster_expected_start(self):
        expected_start = map(lambda x: x.location - x.expected_loc_idx, list(self.location))
        return np.round(np.mean(expected_start))

    def add_location(self, new_loc):
        self.location.append(new_loc)

    def __repr__(self):
        return "[%s]" % ",".join(map(lambda x: str(x), self.location))
        
class location:
    def __init__(self, location, expected_loc_idx):
        self.location = location
        self.expected_loc_idx = expected_loc_idx

    def __repr__(self):
        return str(self.location)

def get_cluster_expected_start(cluster):
    # print cluster
    expected_start = map(lambda x: x.location - x.expected_loc_idx, list(cluster))
    return np.round(np.mean(expected_start))


def in_range(cluster, new_loc):
    # if location > min(cluster) and location < max(cluster):
    #     return True
    cluster_expected_start = get_cluster_expected_start(cluster)
    if new_loc.location - new_loc.expected_loc_idx > cluster_expected_start - 10 and new_loc.location - new_loc.expected_loc_idx < cluster_expected_start + 10:
        return True

    return False

def clusterize_matches(matches):
    all_cluster = []
    for sub_read_idx, one_sub_read_matches in enumerate(matches):
        if len(all_cluster) == 0:
            # ipdb.set_trace()
            all_cluster.extend(map(lambda x: cluster(location(x, sub_read_idx*10)), one_sub_read_matches))
        else:
            for  one_match in one_sub_read_matches: 
                not_in_range_anywhere = True
                new_loc = location(one_match, sub_read_idx * 10)
                for one_cluster in all_cluster:
                    if one_cluster.is_in_range(new_loc):
                        not_in_range_anywhere = False
                        one_cluster.add_location(new_loc)
                        break

                if not_in_range_anywhere:
                    # print "not_in_range_anywhere %d" % one_match
                    all_cluster.append(cluster(new_loc))

    return all_cluster


def remove_small_cluster(cluster):
    return filter(lambda x: x.get_cluster_size() >= 3, cluster)

    # print cluster

def get_clusters_from_matches(matches):
    return remove_small_cluster(clusterize_matches(matches))

if __name__ == "__main__":
    with open("read_matches.pickle", "rb") as f:
        matches = pickle.load(f)

    match_cluster = {
        "forward": get_clusters_from_matches(matches["forward"]),
        "backward": get_clusters_from_matches(matches["backward"])
    }

    print match_cluster

