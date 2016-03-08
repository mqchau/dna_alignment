import pickle
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

        
# class location:
#     def __init__(self, location, expected_loc_idx):
#         self.location = location
#         self.expected_loc_idx = expected_loc_idx

def in_range(cluster, location):
    # if location > min(cluster) and location < max(cluster):
    #     return True
    if location > min(cluster) - 25 and location < max(cluster) + 25:
        # ipdb.set_trace()
        return True

    return False

def clusterize_matches(matches):
    cluster = []
    for one_sub_read_matches in matches:
        if len(cluster) == 0:
            cluster.extend(map(lambda x: set([x]), one_sub_read_matches))
        else:
            for one_match in one_sub_read_matches: 
                not_in_range_anywhere = True
                for one_cluster in cluster:
                    if in_range(one_cluster, one_match):
                        not_in_range_anywhere = False
                        one_cluster.add(one_match)
                        break

                if not_in_range_anywhere:
                    # print "not_in_range_anywhere %d" % one_match
                    cluster.append(set([one_match]))

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

