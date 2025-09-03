# Given a list of sequences, return the kmer that appears most frequently in all
def most_frequent_kmer(sequence, k):
  length=len(sequence)
  if k>length:
    return "k is bigger than the length of the sequence"
  else:
    if k==length:
      return sequence
    else:
      # The dictionary to store the count of kmers
      dic_count={}
      for i in range((length-k)+1):
        kmer=sequence[i:i+k]
        if kmer not in dic_count:
          dic_count[kmer] =1
        else:
          dic_count[kmer] += 1
  value_key_list=[(value,key) for key, value in dic_count.items()]
  #return max(value_key_list)[1] # Returns only the most frequently identified k-mer
  return value_key_list
KmerFreq = most_frequent_kmer("ATGCTGCCGTAATGCCGATCAACGTCGGACTATGC", k = 3)
print(KmerFreq)