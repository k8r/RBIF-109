# Given forward and reverse primers and a DNA sequence, return the predicted
# PCR product: the subsequence from the forward primer to the reverse primer
# (inclusive), accounting for the reverse complement of the reverse primer.
def pcr_product(forward_primer, reverse_primer, sequence):
  # find the index of the forward primer (handle not finding it)
  # make the reverse complement of the reverse primer
  # find the index of the reverse complement (handle not finding it)
  # return the string between the two locations plus the reverse complement at the end
  return sequence

# Given a list of sequences, return a tuple containing the maximum count
# and the list of k-mers that occur with that count.
def most_frequent_kmer(sequences, k):
  dict_count = {}
  for seq_index in range(len(sequences)):
    curr_seq = sequences[seq_index]
    if k > len(curr_seq):
        print(f"k is bigger than the length of the sequence: {curr_seq}. This sequence won't be considered.")
        continue
    
    for i in range((len(curr_seq) - k) + 1):
        kmer = curr_seq[i:i+k]
        if kmer not in dict_count:
          dict_count[kmer] = 1
        else:
          dict_count[kmer] += 1

  first_max_kmer = max(dict_count.items(), key=lambda x: x[1])
  all_max_kmers = [(k) for k, v in dict_count.items() if v == first_max_kmer[1]] # check for any ties
  return (first_max_kmer[1], all_max_kmers)

# Run all functions
print("\n------------------------2) IN SILICO PCR-------------------------")
pcr_product("TCCGA", "GGATC", "AAAATCCGAAATCGGTACCGCTAGGGGTCCAAAAAAAA")

print("\n---------3) MOST FREQUENT K-MERS IN A LIST OF SEQUENCES----------")
sequences = ["ATTTTTTT", "ATGT", "GCTA", "CCAAAAAA", "AAA", "AAA", "TTT", "CA"]
kmer_freq = most_frequent_kmer(sequences, k = 3)
formatted = ", ".join(f"{k}" for k in kmer_freq[1]) + f" ({kmer_freq[0]} occurences)"
print(f"The most frequent kmer(s) are: {formatted}")
print("-----------------------------------------------------------------\n")