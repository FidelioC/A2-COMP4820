from Bio.Align import substitution_matrices
from Bio.Data.IUPACData import protein_letters
import itertools


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

if __name__ == "__main__":
    # temp inputs, default
    file_name = ""
    query = "AVEKQLAEP"
    kmer_size = 3
    neighbourhood_threshold_T = 14
    extension_threshold_eT = 15
    hsp_threshold = ""

    possible_alphabets = generate_possible_alphabets(protein_letters, kmer_size)

    all_kmers = generate_kmers(query, kmer_size)
    print(
        generate_neighbourhood(all_kmers, possible_alphabets, neighbourhood_threshold_T)
    )
