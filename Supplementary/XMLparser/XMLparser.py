# Run: python3 XMLparser.py xml_file_with_graph file_to_store_results


import sys
import xml.dom.minidom as md
import numpy as np


def xml_graph_to_adjacency_matrix(filename):
    dom = md.parse(filename)

    # print(dom.toprettyxml())

    nodes = dom.getElementsByTagName("Node")
    ids = [int(a.getAttribute('id')) for a in nodes]
    # parserint(ids)
    adjacency_matrix = np.zeros(shape=(len(ids), len(ids)))
    edges = dom.getElementsByTagName("Edge")
    for e in edges:
        source = int(e.getElementsByTagName('from')[0].firstChild.nodeValue)
        target = int(e.getElementsByTagName('to')[0].firstChild.nodeValue)
        adjacency_matrix[source][target] = 1
        adjacency_matrix[target][source] = 1
    return adjacency_matrix

if __name__ == "__main__":
    filename = sys.argv[1]
    resultfilename = sys.argv[2]
    # filename = '/home/user/Sirius/XMLparser/Graph.xml'
    matrix = xml_graph_to_adjacency_matrix(filename)
    ans = "\n".join([" ".join([str(i) for i in j]) for j in matrix])
    with open(resultfilename, "w") as f:
        f.write(ans + '\n')
