#!/usr/bin/env python
import argparse

import numpy as np


"""
Given a SV vcf file, calculate the maximum Jaccard similarity between each accession and the pimps for non-overlapping
windows along the reference genome coordinates.

Output is a matrix. Each row is an SLL accession, and each column is the non-overlapping window. Each element is 
the maximum Jaccard similarity of SVs in that window between that accession and the pimp accessions. 
"""

parser = argparse.ArgumentParser(description='Get Jaccard similarity between SLL and a comparison group.')
parser.add_argument("vcf", metavar="<SVs.vcf>", type=str, help="SV vcf file with support vectors. Only one chromosome at a time allowed.")
parser.add_argument("chr", metavar="<chr_name>", type=str, help="Name of reference chromosome.")
parser.add_argument("species_file", metavar="<group.txt>", type=str, help="First column is the phylogenetic group (SLC, SP, GAL, CHE, or SLL), second column is the accession")
parser.add_argument("species", metavar="<SP>", type=str, default="SP", help="Group to compare to SLL (SLC, SP, GAL, or CHE)")
parser.add_argument("fai", metavar="<reference.fasta.fai>", type=str, help="Fasta index file for the reference genome")
parser.add_argument("w", metavar="<100000>", type=int, default=1000000, help="Introgression window size.")
parser.add_argument("-m", metavar="5", type=int, default=5, help='minimum number of SVs needed to calculate Jaccard')

args = parser.parse_args()
vcf_file = args.vcf
chr_used = args.chr
fai_file = args.fai
species_file = args.species_file
comp_species = args.species
window_size = args.w
min_den = args.m

if min_den < 2:
    raise ValueError("-m must be at least 2")

# Get SLL and SP accessions
# Store in a list so their indices can be used in results matrices
sll = []
sp = []
slc = []
sg = []
sc = []
with open(species_file, "r") as f:
    for line in f:
        species, acc = line.rstrip().split("\t")
        if species == "SLL":
            sll.append(acc)
        elif species == "SP":
            sp.append(acc)
        elif species == "SLC":
            slc.append(acc)
        elif species == "GAL":
            sg.append(acc)
        elif species == "CHE":
            sc.append(acc)
        else:
            raise ValueError("Illegal species in species file. Can only be SLL, SLC, SP, GAL, or CHE")

# Associate the list of samples with a species
comp_species_dict = {"SLC": slc, "SP": sp, "GAL": sg, "CHE": sc}

# Get the chromosome sizes
chr_lens = dict()
with open(fai_file, "r") as f:
    for line in f:
        header, length, x, y, z = line.rstrip().split("\t")
        chr_lens[header] = int(length)

# Iterate through the VCF file and get the distances
# Assumes only one chromosome at a time
acc_vec = None
n_windows = chr_lens[chr_used] // window_size
distances = np.zeros((len(sll), n_windows))
comp_max_accs = np.zeros((len(sll), n_windows), dtype=np.int32)
current_window = 0
supp_matrix = []

with open(vcf_file, "r") as f:
    for line in f:
        line = line.rstrip()
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                acc_vec = [i.replace(".ont.s", "") for i in line.split("\t")[9:]]
                assert set(acc_vec) == set(sll + sp + slc + sg + sc)
        else:
            fields = line.split("\t")
            tags = fields[7].split(";")
            start = int(fields[1])
            widx = start // window_size

            # Get the support vector
            supp_vec = None
            for j in tags:
                if j.startswith("SUPP_VEC="):
                    supp_vec = [int(i) for i in j[9:]]
            if supp_vec is None:
                raise ValueError("Missing 'SUPP_VEC' field")

            # Build the support matrix for this window or start a new matrix for a new window
            if widx == current_window:
                supp_matrix.append(supp_vec)
            else:
                # We have moved on to the next window. Save the distances for the finished window
                sm = np.asarray(supp_matrix)

                # For now, I will calculate the distances in a 'for' loop. Perhaps vectorize in the future
                # Iterate over the SLLs
                t_distances = [[] for i in range(len(sll))]
                for i in range(len(sll)):
                    this_acc = sll[i]
                    supp_idx = acc_vec.index(this_acc)
                    this_vec = sm[:, supp_idx]

                    # Iterate over the comp species:
                    for comp_acc in comp_species_dict[comp_species]:
                        comp_supp_idx = acc_vec.index(comp_acc)
                        this_comp_vec = sm[:, comp_supp_idx]
                        if np.count_nonzero(this_comp_vec) >= min_den and np.count_nonzero(this_vec) >= min_den:
                            # Get the Jaccard distance
                            num = np.count_nonzero(np.logical_and(this_vec, this_comp_vec))  # Intersection
                            den = np.count_nonzero(np.logical_or(this_vec, this_comp_vec))  # Union
                            t_distances[i].append(num / den)
                        else:
                            t_distances[i].append(-1)

                # Find which comp sample gave the max
                t_distances_argmax = np.asarray([np.argmax(i) for i in t_distances])
                comp_max_accs[:, current_window] = t_distances_argmax
                # Get the max % shared SVs between a given SLL and each comp species sample
                t_distances = np.asarray([np.max(i) for i in t_distances])

                distances[:, current_window] = t_distances

                # Now that we have calculated the distances for the finished window, start the next one
                if widx == n_windows:
                    break
                current_window = widx
                supp_matrix = []
                supp_matrix.append(supp_vec)

# Write out the comparison species that gave the max
with open("comp_matrix." + comp_species + ".txt", "w") as f:
    f.write("Sample\t" + "\t".join( [str(i*window_size) for i in range(n_windows)]) + "\n")
    for i in range(len(sll)):
        f.write(sll[i] + "\t" + "\t".join( [comp_species_dict[comp_species][j] for j in list(comp_max_accs[i, :])] ) + "\n")

# Print the output matrix
print("Sample\t" + "\t".join( [str(i*window_size) for i in range(n_windows)] ))
for i in range(len(sll)):
    print(sll[i] + "\t" + "\t".join( [str(j) for j in list(distances[i, :])] ).replace("-1.0", "NA"))