import pickle
import numpy as np
import ipdb

# class location_cluster:
#     def __init__(self):
#         self.location = []

#     def add_location(self, new_loc, expected_loc_idx):
#         if len(self.location) == 0:
#             self.location.append(location(new_loc, expected_loc_idx))
#             self.update_expected_cluster_location()
#             return True
#         else:
#             pass

#     def update_expected_cluster_location(self):
#         expected_cluster_location = 
#         for one_loc in self.location:

        
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
    cluster = []
    for sub_read_idx, one_sub_read_matches in enumerate(matches):
        if len(cluster) == 0:
            cluster.extend(map(lambda x: [location(x, sub_read_idx*10)], one_sub_read_matches))
        else:
            for  one_match in one_sub_read_matches: 
                not_in_range_anywhere = True
                new_loc = location(one_match, sub_read_idx * 10)
                for one_cluster in cluster:
                    if in_range(one_cluster, new_loc):
                        not_in_range_anywhere = False
                        one_cluster.append(new_loc)
                        break

                if not_in_range_anywhere:
                    # print "not_in_range_anywhere %d" % one_match
                    cluster.append([new_loc])

    return cluster


def remove_small_cluster(cluster):
    return filter(lambda x: len(x) >= 3, cluster)

    # print cluster


if __name__ == "__main__":
    with open("read_matches.pickle", "rb") as f:
        matches = pickle.load(f)

    forward_matches = matches["forward"]
    forward_clusters = remove_small_cluster(clusterize_matches(forward_matches))
    backward_matches = matches["backward"]
    backward_clusters = remove_small_cluster(clusterize_matches(backward_matches))

    print forward_clusters, backward_clusters

