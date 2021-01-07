import sys
import Bio.SeqIO
from tqdm import tqdm


def get_data(seq, k):
    seq_dict = {}
    l = len(str(seq[0].seq))
    #every sequence, count str with length k in them. Made a dict
    for j in tqdm(range(len(seq))): 
        i = 0
        seq_ = str(seq[j].seq)
        while i <= l-k :
            chunk = seq_[i:i+k]
            if chunk in seq_dict.keys() :
                seq_dict[chunk] += 1
            else: 
                seq_dict[chunk] = 1
            i += 1
    return seq_dict

def find_missing_keys_in_dict(from_dict, to_dict, file_name_1, file_name_2):
    print("Looking at sub-reads from {} file that miss in the {} file".format(file_name_1, file_name_2))
    
    chosen_chunks = {}

    for chunk, count in tqdm(from_dict.items()):
        if chunk not in to_dict and count > 5:
            chosen_chunks[chunk] = count
    
    return chosen_chunks

def combine_found_snps(arr):
    to = {}
    from_ = {}
    
    for el1 in arr:
        for el2 in arr:
            if el1 != el2 and el1[1:] == el2[:-1]:
                assert(el1 not in from_.keys())
                from_[el1] = el2
                to[el2] = el1
    
    res = []
    for el in set(from_.keys()).difference(set(to.keys())):
        cur = el
        part = el
        while(cur in from_.keys()):
            cur = from_[cur]
            part += cur[-1]
        res.append(part)
    combine(set(chosen_chunks_1.keys()))
    return res


def list_of_snps(argv1, argv2):
    dict_argv1 = get_data(argv1, 30)
    dict_argv2 = get_data(argv2, 30)

    chunks_missing_in_second_dict = find_missing_keys_in_dict(dict_argv1, dict_argv2, "first", "second")
    chunks_missing_in_first_dict = find_missing_keys_in_dict(dict_argv2, dict_argv1, "second", "first")

    found_snps = combine_found_snps(set(chunks_missing_in_second_dict.keys()))
    original_without_snps = combine_found_snps(set(chunks_missing_in_first_dict.keys()))

    print("The found snps are:")
    print(found_snps)
    print("The orignial is:")
    print(original_without_snps)


if __name__ == "__main__":

    try:
        if (sys.argv[1] == sys.argv[2]):
            print("You have given the same file twice")
            sys.exit(2)
        
        print("Reading the first file:")
        first_file = tqdm(Bio.SeqIO.parse(sys.argv[1],"fasta"))
        first_file = list(first_file)
        print("Reading the second file:")
        second_file = tqdm(Bio.SeqIO.parse(sys.argv[2],"fasta"))
        second_file = list(second_file)

        list_of_snps(first_file, second_file)

    except IndexError:
        print("Two files are required for execution")
        sys.exit(1)
