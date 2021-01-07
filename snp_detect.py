import sys
import Bio.SeqIO
from tqdm import tqdm


def get_data(fasta_reads, k):
    fasta_reads_dict = {}
    l = len(str(fasta_reads[0].seq))
    #every sequence, count str with length k in them. Made a dict
    for j in tqdm(range(len(fasta_reads))): 
        i = 0
        seq_ = str(fasta_reads[j].seq)
        while i <= l-k :
            chunk = seq_[i:i+k]
            if chunk in fasta_reads_dict.keys() :
                fasta_reads_dict[chunk] += 1
            else: 
                fasta_reads_dict[chunk] = 1
            i += 1
    return fasta_reads_dict

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
    return res

def levenshtein_distance(a, b):
    d = np.zeros((len(a), len(b)))

    for i in range(len(a)):
        d[i, 0] = i

    for j in range(len(b)):
        d[0, j] = j

    for j in range(1, len(b)):
        for i in range(1, len(a)):
            if a[i] == b[j]:
                cost = 0
            else:
                cost = 1

            d[i, j] = min(d[i - 1, j] + 1, 
                          d[i, j - 1] + 1, 
                          d[i - 1, j - 1] + cost)

    return int(d[-1][-1])

def list_of_snps(argv1, argv2):
    dict_argv1 = get_data(argv1, 30)
    dict_argv2 = get_data(argv2, 30)

    chunks_missing_in_second_dict = find_missing_keys_in_dict(dict_argv1, dict_argv2, "first", "second")
    chunks_missing_in_first_dict = find_missing_keys_in_dict(dict_argv2, dict_argv1, "second", "first")

    found_snps = combine_found_snps(set(chunks_missing_in_second_dict.keys()))
    original_without_snps = combine_found_snps(set(chunks_missing_in_first_dict.keys()))

    for original in original_without_snps:
	    for found in found_snps:
	        if levenshtein_distance(original, found) < 10:
	            print("The found snps are: ", end='')
	            print(found)
	            print("The orignial is:    ", end='')
	            print(original)


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
