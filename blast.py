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
def extend_seeds(fasta_file, all_potential_seeds, extension_threshold_eT, aho_trie):
    """Extend seeds:
    - Extend left or right for each seed, calculate the total score
    - Extend until you either reach the end of the query,
    - or the whole score of the aligned sequence and query drops below the maximum score by more than the threshold eT
    """
    query_dictionary = {}
    query_dictionary = aho_corasick.aho_corasick(query, aho_trie, query_dictionary)

    protein_sequence = "MKQLLVRSGKVFLQDVPAPVAGPKNVLVRVQRSCVSVGTEMAGVKMSGLPLYRRALKQPHHVKRVLQLMRDQGVARVYKQVKGKLDAGLPTGYSAAGTVIAVGSEVDGIAVGDRVACAGAGVANHAEVIDVPVNLCVPVPQQVSFDAAATVTLGAIAMQGVRRAQPTLGETVVVIGLGILGQITAQLLTANGCRVIGTDVDNKRIATALENGLDHGINPNDGNLVDSIIKLTDGFGADVAIITAASASSDILAQAFQSCRKKARVVIVGDVGLNMARSDIYTKELDVLISCSYGPGRYDPVYEDEGGDYPLAYVRWTENRNMGEYLRLLAAGRVRLDNMLQEPYPIDRAEEAYGRLAGEGEKPLLVLLQYPHREEAVRSVLQIAPPKPVDGRIKIAVVGAGSFAQGMHLPNLKRLGDKFDLRSVVSRTGLSARTAAERFGFSTASTDFQAVLDDPQVDLVLIATRHDLHAEMTLAALKADKHVFVEKPTSMTEEGLDAIEAFYADNPNGPLLMTGFNRRFAPAVTAAREAIKGRLSPMIVNYRMNAGYIPSDHWVHGPHGGGRNIGEACHIYDLFNALTGSQPVEVQARSIVPASGHWRRDDNFVATVRYADGSLCTLTYTSLGSKEFPKERFDIFVDGKVLVLDDYKRLEVTGAKGGWKGLTIEKGQLEELVALAEAFKPGGEWPISLADQLSATRVSFAVEKQLAE"
    right_sequence, score_right, query_index_right, db_index_right = extend_right(
        query, protein_sequence, 1, 701, extension_threshold_eT
    )
    left_sequence, score_left, query_index_left, db_index_left = extend_left(
        query, protein_sequence, 1, 701, extension_threshold_eT
    )

    return (
        left_sequence + right_sequence,
        score_right + score_left,
        (query_index_left, query_index_right),
        (db_index_left, db_index_right),
    )


def extend_right(
    query_sequence,
    protein_sequence,
    query_start_index,
    db_start_index,
    extension_threshold_eT,
):
    """extend right"""
    blosum_result = 0
    curr_sequence = ""
    max = blosum_result
    prev_max = 0

    query_index = query_start_index
    db_index = db_start_index

    while (prev_max - extension_threshold_eT <= max) and (
        query_index < len(query) and db_index < len(protein_sequence)
    ):
        curr_sequence += query_sequence[query_index]
        # print(
        #     f"curr_sequence {curr_sequence}, query_index: {query_index}, db_index: {db_index}, query_char: {query_sequence[query_index]}, protein_char: {protein_sequence[db_index]}"
        # )

        blosum_result = blosum_62()[query_sequence[query_index]][
            protein_sequence[db_index]
        ]
        # print(f"max: {max}")
        # print(f"blosum_result: {blosum_result}")

        prev_max = max
        max += blosum_result

        query_index += 1
        db_index += 1

    return curr_sequence, max, query_index - 1, db_index - 1


def extend_left(
    query_sequence,
    protein_sequence,
    query_start_index,
    db_start_index,
    extension_threshold_eT,
):
    """extend_left"""
    blosum_result = 0
    curr_sequence = ""
    max = blosum_result
    prev_max = 0

    query_index = query_start_index - 1
    db_index = db_start_index - 1

    while (prev_max - extension_threshold_eT <= max) and (
        query_index >= 0 and db_index >= 0
    ):
        curr_sequence += query_sequence[query_index]
        # print(
        #     f"curr_sequence {curr_sequence}, query_index: {query_index}, db_index: {db_index}, query_char: {query_sequence[query_index]}, protein_char: {protein_sequence[db_index]}"
        # )

        blosum_result = blosum_62()[query_sequence[query_index]][
            protein_sequence[db_index]
        ]
        # print(f"max: {max}")
        # print(f"blosum_result: {blosum_result}")

        prev_max = max
        max += blosum_result

        query_index -= 1
        db_index -= 1

    return curr_sequence, max, query_index + 1, db_index + 1


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
    extension_threshold_eT = 5
    hsp_threshold = ""

    # possible_alphabets = generate_possible_alphabets(protein_letters, kmer_size)

    # all_kmers = generate_kmers(query, kmer_size)

    # all_neighbourhood = generate_neighbourhood(
    #     all_kmers, possible_alphabets, neighbourhood_threshold_T
    # )

    all_neighbourhood = ["VEK", "EKQ", "KQL", "AEP"]

    # print(f"all_neighbourhood: {all_neighbourhood}")

    aho_trie = aho_corasick.Trie(all_neighbourhood)

    # all_potential_seeds, total_count = search_potential_seeds(
    #     file_name, aho_trie
    # )

    # print(f"all_potential_seeds: {all_potential_seeds}")

    # print(f"total_count potential seeds: {total_count}")
    query_dictionary = {}
    query_dictionary = aho_corasick.aho_corasick(query, aho_trie, query_dictionary)
    print(query_dictionary)

    print(extend_seeds(file_name, "", extension_threshold_eT, aho_trie))
