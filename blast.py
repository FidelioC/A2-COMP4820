from Bio.Align import substitution_matrices
from Bio.Data.IUPACData import protein_letters
from Bio import SeqIO
import itertools
import aho_corasick
import click


# ================== 1st step: Generate Neighbourhood ===================== #
def generate_neighbourhood(
    all_kmers: list,
    possible_alphabets: list,
    neighbourhood_threshold_T,
    kmer_size,
    blosum,
):
    """generate neighbourhood"""
    all_neighbourhood = {}
    query_index = 0
    for kmer in all_kmers:
        # print(f"kmer: {kmer}")
        for alphabet in possible_alphabets:

            current_score = 0
            for i in range(kmer_size):
                current_score += blosum[kmer[i]][alphabet[i]]

            if current_score >= neighbourhood_threshold_T:
                # print(
                #     f"kmer: {kmer} alphabets: {alphabet}, current_score: {current_score}"
                # )
                if query_index not in all_neighbourhood:
                    all_neighbourhood[query_index] = [alphabet]
                else:
                    all_neighbourhood[query_index].append(alphabet)

        query_index += 1

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
        [query_index_left, query_index_right],
        [db_index_left, db_index_right],
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
    all_neighbourhood,
    query,
    extension_threshold_eT,
    hsp_threshold,
    blosum,
):
    output_dict = {}
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
                        db_indices,
                    ) = extend_seeds(
                        extension_threshold_eT,
                        query,
                        record_seq,
                        find_key_dict_by_value(all_neighbourhood, seed),
                        db_index,
                        blosum,
                    )

                    if total_score >= hsp_threshold:
                        record_desc = record.description
                        output_result = {
                            "Score": total_score,
                            "Record sequence": record_sequence + str(db_indices),
                            "Query sequence": query_sequence + str(query_indices),
                        }
                        if record_desc not in output_dict:
                            output_dict[record_desc] = output_result
                    total_score = 0

    return output_dict


def find_key_dict_by_value(dictionary, target_val):
    for key, value_list in dictionary.items():
        if target_val in value_list:
            return key
    return None


# ========================================================================== #
def print_output(all_kmers, all_neighbourhood_list, output_dict, total_count):
    print("Output\n======")
    print("\nNeighbourhood\n-------------")
    print(f"k-mers: {all_kmers}")
    print(f"neighbourhood: {all_neighbourhood_list}")
    print("\nHSPs\n----")
    index = 1
    for key, value in output_dict.items():
        print(f"{index}. {key}.")
        for item_key, item_value in value.items():
            print(f"   * {item_key}: {item_value}")
        index += 1
    print("\nSummary\n-------")
    print(f"* {total_count} seeds were considered exceeding T")
    print(f"* {len(output_dict)} hits were found >= S")


@click.command()
@click.option("--database", required=True)
@click.option("--query", required=True)
@click.option("--kmer-size", default=3)
@click.option("--neighbourhood-threshold", default=14)
@click.option("--extension-threshold", default=5)
@click.option("--hsp-threshold", default=36)
def commands_processing(
    database,
    query,
    kmer_size,
    neighbourhood_threshold,
    extension_threshold,
    hsp_threshold,
):
    main(
        database,
        query,
        kmer_size,
        neighbourhood_threshold,
        extension_threshold,
        hsp_threshold,
    )


def main(
    file_name,
    query,
    kmer_size,
    neighbourhood_threshold_T,
    extension_threshold_eT,
    hsp_threshold,
):
    # init
    blosum = substitution_matrices.load("BLOSUM62")

    # step 1: generate alphabets, kmers, neighbourhood
    possible_alphabets = generate_possible_alphabets(protein_letters, kmer_size)
    all_kmers = generate_kmers(query, kmer_size)
    all_neighbourhood = generate_neighbourhood(
        all_kmers, possible_alphabets, neighbourhood_threshold_T, kmer_size, blosum
    )

    # Combine all values into a single list to create trie
    all_neighbourhood_list = sum(all_neighbourhood.values(), [])
    aho_trie = aho_corasick.Trie(all_neighbourhood_list)

    # step 2: search all potential seeds
    all_potential_seeds, total_count = search_potential_seeds(file_name, aho_trie)

    # step 3: extend seeds and calculate HSP
    output_dict = calculate_HSPs(
        file_name,
        all_potential_seeds,
        all_neighbourhood,
        query,
        extension_threshold_eT,
        hsp_threshold,
        blosum,
    )

    print_output(all_kmers, all_neighbourhood_list, output_dict, total_count)


if __name__ == "__main__":
    commands_processing()
