from Bio.Align import substitution_matrices
from Bio.Data.IUPACData import protein_letters
from Bio import SeqIO
import itertools
import aho_corasick


# ================== 1st step: Generate Neighbourhood ===================== #
def generate_neighbourhood(
    all_kmers: list,
    possible_alphabets: list,
    neighbourhood_threshold_T,
    kmer_size,
    blosum,
):
    """generate neighbourhood"""
    all_neighbourhood = []
    for kmer in all_kmers:
        for alphabet in possible_alphabets:
            current_score = 0
            for i in range(kmer_size):
                current_score += blosum[kmer[i]][alphabet[i]]

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
def search_potential_seeds(fasta_file, aho_trie):
    """search all potential seeds using aho corasick algorithm"""
    all_potential_seeds = {}
    total_count = 0

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


# ================== 3rd step: Extending Seeds ============================= #
def extend_seeds(
    extension_threshold_eT,
    query_sequence,
    protein_sequence,
    query_start_index,
    db_start_index,
    blosum,
):
    """Extend seeds:
    - Extend left or right for each seed, calculate the total score
    - Extend until you either reach the end of the query,
    - or the whole score of the aligned sequence and query drops below the maximum score by more than the threshold eT
    """
    (
        right_query_sequence,
        right_record_sequence,
        score_right,
        query_index_right,
        db_index_right,
    ) = extend_right(
        query_sequence,
        protein_sequence,
        query_start_index,
        db_start_index,
        extension_threshold_eT,
        blosum,
    )

    (
        left_query_sequence,
        left_record_sequence,
        score_left,
        query_index_left,
        db_index_left,
    ) = extend_left(
        query_sequence,
        protein_sequence,
        query_start_index,
        db_start_index,
        extension_threshold_eT,
        blosum,
    )

    record_sequence = left_record_sequence + right_record_sequence
    query_sequence = left_query_sequence + right_query_sequence
    total_score = score_right + score_left
    return (
        record_sequence,
        query_sequence,
        total_score,
        (query_index_left, query_index_right),
        (db_index_left, db_index_right),
    )


def extend_right(
    query_sequence,
    protein_sequence,
    query_start_index,
    db_start_index,
    extension_threshold_eT,
    blosum,
):
    """extend right"""
    blosum_result = 0
    curr_query_sequence = ""
    curr_record_sequence = ""
    max = blosum_result
    prev_max = 0

    query_index = query_start_index
    db_index = db_start_index

    while (prev_max - extension_threshold_eT <= max) and (
        query_index < len(query_sequence) and db_index < len(protein_sequence)
    ):
        curr_query_sequence += query_sequence[query_index]
        curr_record_sequence += protein_sequence[db_index]
        # print(
        #     f"curr_sequence {curr_sequence}, query_index: {query_index}, db_index: {db_index}, query_char: {query_sequence[query_index]}, protein_char: {protein_sequence[db_index]}"
        # )

        blosum_result = blosum[query_sequence[query_index]][protein_sequence[db_index]]
        # print(f"max: {max}")
        # print(f"blosum_result: {blosum_result}")

        prev_max = max
        max += blosum_result

        query_index += 1
        db_index += 1

    return curr_query_sequence, curr_record_sequence, max, query_index - 1, db_index - 1


def extend_left(
    query_sequence,
    protein_sequence,
    query_start_index,
    db_start_index,
    extension_threshold_eT,
    blosum,
):
    """extend_left"""
    blosum_result = 0
    curr_query_sequence = ""
    curr_record_sequence = ""
    max = blosum_result
    prev_max = 0

    query_index = query_start_index - 1
    db_index = db_start_index - 1

    while (prev_max - extension_threshold_eT <= max) and (
        query_index >= 0 and db_index >= 0
    ):
        curr_query_sequence = query_sequence[query_index] + curr_query_sequence
        curr_record_sequence = protein_sequence[db_index] + curr_record_sequence
        # print(
        #     f"curr_sequence {curr_sequence}, query_index: {query_index}, db_index: {db_index}, query_char: {query_sequence[query_index]}, protein_char: {protein_sequence[db_index]}"
        # )

        blosum_result = blosum[query_sequence[query_index]][protein_sequence[db_index]]
        # print(f"max: {max}")
        # print(f"blosum_result: {blosum_result}")

        prev_max = max
        max += blosum_result

        query_index -= 1
        db_index -= 1

    return curr_query_sequence, curr_record_sequence, max, query_index + 1, db_index + 1


# ========================================================================== #


# ================== 4th step: Calculate HSPs ============================== #
def calculate_HSPs(
    fasta_file,
    all_potential_seeds,
    aho_trie,
    query,
    extension_threshold_eT,
    hsp_threshold,
    blosum,
):
    query_dictionary = {}
    query_dictionary = aho_corasick.aho_corasick(query, aho_trie, query_dictionary)
    print(query_dictionary)

    # iterate over every potential seeds to find HSP
    for record in SeqIO.parse(fasta_file, "fasta"):
        record_id = record.id
        if record_id in all_potential_seeds:
            seeds_dict = all_potential_seeds[record_id]
            # print(record.id)
            for seed, db_indexes in seeds_dict.items():
                record_seq = record.seq
                # print(seed)
                for db_index in db_indexes:
                    # print(db_index)
                    (
                        record_sequence,
                        query_sequence,
                        total_score,
                        query_indices,
                        db_indicies,
                    ) = extend_seeds(
                        extension_threshold_eT,
                        query,
                        record_seq,
                        query_dictionary[seed][0],
                        db_index,
                        blosum,
                    )

                    if total_score >= hsp_threshold:
                        print(record.description)
                        print(
                            record_sequence,
                            query_sequence,
                            total_score,
                            query_indices,
                            db_indicies,
                        )
                    total_score = 0


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


def main():
    # temp inputs, default
    file_name = "proteins.fasta"
    query = "AVEKQLAEP"
    kmer_size = 3
    neighbourhood_threshold_T = 14
    extension_threshold_eT = 5
    hsp_threshold = 36

    blosum = substitution_matrices.load("BLOSUM62")

    possible_alphabets = generate_possible_alphabets(protein_letters, kmer_size)
    # print(f"possible_alphabets: {possible_alphabets}")

    all_kmers = generate_kmers(query, kmer_size)
    # print(f"all_kmers: {all_kmers}")

    all_neighbourhood = generate_neighbourhood(
        all_kmers, possible_alphabets, neighbourhood_threshold_T, kmer_size, blosum
    )

    # all_neighbourhood = ["VEK", "EKQ", "KQL", "AEP"]

    print(f"all_neighbourhood: {all_neighbourhood}")

    aho_trie = aho_corasick.Trie(all_neighbourhood)

    all_potential_seeds, total_count = search_potential_seeds(file_name, aho_trie)

    # print(f"all_potential_seeds: {all_potential_seeds}")

    print(f"total_count potential seeds: {total_count}")

    # parameters for extend seeds
    # protein_sequence = "MQYFLCLADEKNVTRAARRLNIVQPALSMQIAKLEVELGQRLFDRSVQGMTLTSAGEALVRLTAPIVRDAEYARQEMAQIGGRISGRVAVGLITSVAQSTMASSSATVARRYPEIILSACEGYTETLVDWVNSGQLDFALINVPRRRTPLAAHHIMDEEMVFACRKDGPIRPAAKLRFDHIANFDLVLPSKRHGLRLILDEHAAALGIDLRPRLELDTLPALCDVIATTDFATVLPTIALRQSLASGTTRAHRFDAQRIVRSIAWVHHPRRAVSVAAKAVLDVISHDLAQAAAVAKQLAEPGSGGAASSSRKQRRKTGKTIS"
    # query_start_index = 3
    # db_start_index = 295

    # print(
    #     extend_seeds(
    #         extension_threshold_eT,
    #         query,
    #         protein_sequence,
    #         query_start_index,
    #         db_start_index,
    #     )
    # )

    calculate_HSPs(
        file_name,
        all_potential_seeds,
        aho_trie,
        query,
        extension_threshold_eT,
        hsp_threshold,
        blosum,
    )


if __name__ == "__main__":
    main()
