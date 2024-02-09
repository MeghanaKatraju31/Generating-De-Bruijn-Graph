import networkx as nx
import matplotlib.pyplot as plt

def overlap(read_a, read_b, min_length=1):
    start = read_a.find(read_b[:min_length])
    while start >= 0:
        if read_b.startswith(read_a[start:]):
            return len(read_a) - start
        start = read_a.find(read_b[:min_length], start + 1)
    return 0

def build_overlap_graph(reads, k):
    overlaps = {}
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            if kmer not in overlaps:
                overlaps[kmer] = set()
            for other_read in reads:
                if read != other_read:
                    for j in range(len(other_read) - k + 1):
                        other_kmer = other_read[j:j+k]
                        if kmer == other_kmer:
                            overlaps[kmer].add(other_read)
    return overlaps


def shortest_common_superstring(reads, k):
    overlaps = build_overlap_graph(reads, k)
    node_degree = {}
    for kmer in overlaps:
        node_degree[kmer] = len(overlaps[kmer])
    start_node = max(node_degree, key=node_degree.get)
    path = [start_node]
    while overlaps[path[-1]]:
        next_nodes = [n for n in overlaps[path[-1]] if n in node_degree]
        if not next_nodes:
            break
        next_node = max(next_nodes, key=lambda x: node_degree[x])
        path.append(next_node)
        overlaps[path[-2]].remove(next_node)
        node_degree[path[-1]] = 0
    return ''.join([path[i][-1] for i in range(len(path))])

reads = []
with open('1002088571.txt', 'r') as f:
    for line in f:
        reads.append(line.strip())

# Calculate all possible k-mers
k = 2
kmers = set()
for read in reads:
    for i in range(len(read) - k + 1):
        kmers.add(read[i:i+k])

# Write k-mers to file
with open('kmers.txt', 'w') as f:
    for kmer in kmers:
        f.write(kmer + '\n')

# Generate and plot De Bruijn graph
G = nx.DiGraph()
with open('1002088571.txt', 'r') as f:
    for line in f:
        kmer = line.strip()
        G.add_edge(kmer[:-1], kmer[1:])

# Create dictionary of labels
labels = {}
for node in G.nodes():
    labels[node] = node

nx.draw(G, labels=labels, with_labels=True)
plt.show()

# Assemble reads
superstring = shortest_common_superstring(reads, k)
print(superstring)
