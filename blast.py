from Bio.Align import substitution_matrices
from Bio.Data.IUPACData import protein_letters
from Bio import SeqIO
import itertools
import aho_corasick


def blosum_62():
    return substitution_matrices.load("BLOSUM62")


# ================== 1st step: Generate Neighbourhood ===================== #
def generate_neighbourhood(
    all_kmers: list, possible_alphabets: list, neighbourhood_threshold_T
):
    """generate neighbourhood"""
    all_neighbourhood = []
    for kmer in all_kmers:
        for alphabet in possible_alphabets:
            current_score = 0
            for i in range(kmer_size):
                current_score += blosum_62()[kmer[i]][alphabet[i]]

            if current_score >= neighbourhood_threshold_T:
                all_neighbourhood.append(alphabet)

    return all_neighbourhood


def generate_kmers(query, k):
    """generate all possible k-mers from the query (of proteins)"""
    all_kmers = []
    for i in range(len(query) - k + 1):
        all_kmers.append(query[i : i + k])

    return all_kmers


def generate_possible_alphabets(alphabet, k):
    """generate all possible alphabets of size k (of proteins)"""
    possible_alphabets = []
    for x in itertools.product(alphabet, repeat=k):
        possible_alphabets.append("".join(x))

    return possible_alphabets


# ========================================================================== #


# ================== 2nd step: Search for seeds ============================ #
def search_potential_seeds(fasta_file, all_neighbourhood):
    all_potential_seeds = {}
    total_count = 0
    aho_trie = aho_corasick.Trie(all_neighbourhood)

    for record in SeqIO.parse(fasta_file, "fasta"):
        dictionary = {}

        # print(f"Protein ID: {record.id}")
        # print(f"Protein Sequence: {record.seq}\n")
        aho_corasick.aho_corasick(record.seq, aho_trie, dictionary)

        # inserting to 'hitPos' dict
        if dictionary:  # dont store if seed is not found in the sequence
            all_potential_seeds[record.id] = dictionary

        for item in dictionary:
            total_count += len(dictionary[item])

    return all_potential_seeds, total_count


# ========================================================================== #


def initialize_output_dictionary(sequence_list_patterns):
    """
    Init output directory with empty list
    Purpose: store 'count' and 'location' of the patterns
    """
    dictionary = {}
    for pattern in sequence_list_patterns:
        if pattern not in dictionary:
            dictionary[pattern] = []
    return dictionary


if __name__ == "__main__":
    # temp inputs, default
    file_name = "proteins.fasta"
    query = "AVEKQLAEP"
    kmer_size = 3
    neighbourhood_threshold_T = 14
    extension_threshold_eT = 15
    hsp_threshold = ""

    possible_alphabets = generate_possible_alphabets(protein_letters, kmer_size)

    # all_kmers = generate_kmers(query, kmer_size)

    # all_neighbourhood = generate_neighbourhood(
    #     all_kmers, possible_alphabets, neighbourhood_threshold_T
    # )

    all_neighbourhood = ["VEK", "EKQ", "KQL", "AEP"]

    print(f"all_neighbourhood: {all_neighbourhood}")

    all_potential_seeds, total_count = search_potential_seeds(
        file_name, all_neighbourhood
    )

    # print(f"all_potential_seeds: {all_potential_seeds}")

    print(f"total_count potential seeds: {total_count}")
