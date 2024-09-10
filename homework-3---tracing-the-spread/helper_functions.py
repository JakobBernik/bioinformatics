import math

def UPGMA(distances):
    """Unweighted pair group method with arithmetic mean (UPGMA) agglomerative clustering.
    
    Parameters
    ----------
    distances: np.ndarray
        A two dimensional, square, symmetric matrix containing distances between data
        points. The diagonal is zeros.
        
    Returns
    -------
    np.ndarray
        The linkage matrix, as specified in scipy. Briefly, this should be a 2d matrix
        each row containing 4 elements. The first and second element should denote the
        cluster IDs being merged, the third element should be the distance, and the
        fourth element should be the number of elements within this new cluster. Any
        new cluster should be assigned an incrementing ID, e.g. after the first step
        where the first two points are merged into a cluster, the new cluster ID should
        be N, then N+1, N+2, ... in subsequent steps.
    
    """

    """ Distances dictionary that will change during clustering"""
    distances_dict = {}
    """ Init clusters dict and fill it with elementary clusters """
    clusters = {}
    for i in range(len(distances)):
        clusters[i] = [i]
        distances_dict[i] = {}
        for j in range(len(distances)):
            if i != j:
                distances_dict[i][j] = distances[i][j]

    """ Will contain results of clustering """
    linkage_matrix = []

    """ The id of new cluster """
    cluster_id = len(distances)

    """ Main part """
    while len(clusters) > 1:
        """ Find min distance """
        min_dist = 50.0
        for c1 in distances_dict:
            for c2 in distances_dict[c1]:
                if distances_dict[c1][c2] < min_dist and c1 != c2:
                    min_dist = distances_dict[c1][c2]
                    first_cluster = c1
                    second_cluster = c2

        clusters[cluster_id] = clusters[first_cluster] + clusters[second_cluster]
        del distances_dict[first_cluster]
        del distances_dict[second_cluster]
        """ Update distances_dict """
        distances_dict[cluster_id] = {}
        for cl in distances_dict:
            if cl != cluster_id:
                distances_dict[cl][cluster_id] = 0
                for e1 in clusters[cluster_id]:
                    for e2 in clusters[cl]:
                        distances_dict[cl][cluster_id] += distances[e1][e2]
                distances_dict[cl][cluster_id] /= len(clusters[cluster_id])*len(clusters[cl])
                distances_dict[cluster_id][cl] = distances_dict[cl][cluster_id]
                del distances_dict[cl][first_cluster]
                del distances_dict[cl][second_cluster]
        """ Clean dict of old clusters """
        del clusters[first_cluster]
        del clusters[second_cluster]
        """ Add entry to linkage_matrix """
        linkage_matrix.append([first_cluster, second_cluster,
                               min_dist, len(clusters[cluster_id])])
        cluster_id += 1
    return linkage_matrix


def jukes_cantor(p: float) -> float:
    """The Jukes-Cantor correction for estimating genetic distances.
    
    Parameters
    ----------
    p: float
        The proportional distance, i.e. the number of of mismatching symbols (Hamming
        distance) divided by the total sequence length.
        
    Returns
    -------
    float
        The corrected genetic distance.
    
    """
    return -(3 / 4) * math.log(1 - 4 / 3 * p)
