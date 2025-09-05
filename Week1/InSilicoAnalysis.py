import random
import argparse

CODON_TABLE = {
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y"
}

# Given two sequences of the same size, return the edit distance
def edit_distance(seq1, seq2):
   # Check that the two sequences are the same length
  if len(seq1) != len(seq2):
    raise ValueError("Sequences must be the same length.")
  
  distance = 0
  for i in range(len(seq1)):
    if seq1[i] != seq2[i]:
        distance += 1
  return distance

# Given forward and reverse primers and a DNA sequence, return the predicted
# PCR product: the subsequence from the forward primer to the reverse primer
# (inclusive), accounting for the reverse complement of the reverse primer.
# Assumptions: 
# - The forward primer's index <  the reverse primer's index.
# - If a primer sequence occurs more than once in the DNA, only the first occurrence is considered.
def pcr_product(forward_primer, reverse_primer, sequence):
  forward_primer_pos = sequence.find(forward_primer)
  if forward_primer_pos == -1: raise ValueError(f"Unable to locate the forward primer ({forward_primer}) in the given sequence ({sequence})")

  reverse_complement_of_reverse_primer = reverse_complement(reverse_primer)

  reverse_primer_pos = sequence.find(reverse_complement_of_reverse_primer)
  if reverse_primer_pos == -1: raise ValueError(f"Unable to locate the reverse primer ({reverse_primer}) in the given sequence ({sequence})")

  return sequence[forward_primer_pos:reverse_primer_pos+len(reverse_primer):1]

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

# Translate the given DNA coding sequence directly into its amino acid sequence 
# (skipping transcription since we're working from the coding strand)
# Assumptions: 
# - No introns need to be removed before translation
# - The start codon is the standard ATG (DNA)
# - The stop codons are the standard TAA, TAG, TGA (DNA)
# - Only one reading frame is considered, defined by the first start codon
def dna_to_protein(seq):
  start_codon = "ATG"
  start_index = seq.find(start_codon)
  if start_index == -1: raise ValueError("Missing the start codon")

  stop_codons = ["TAA", "TAG", "TGA"]
  seq_after_start = seq[start_index::1]
  codons = [seq_after_start[i:i+3:1] for i in range(0, len(seq_after_start), 3)]

  stop_codon_locations = [codons.index(stop_codon) for stop_codon in stop_codons if stop_codon in codons]
  if len(stop_codon_locations) == 0: raise ValueError("Missing a stop codon")
  first_stop_codon = min(stop_codon_locations)

  amino_acids = [CODON_TABLE[codons[i]] for i in range(first_stop_codon)]
  return "".join(amino_acids)

def reverse_complement(sequence):
  # Make a table to change the nucleotide characters to their complements, and then reverse the resulting str
  complement_table = str.maketrans("ATCGatcg", "TAGCtagc")
  return sequence.translate(complement_table)[::-1]

def random_sequence(length):
  return ''.join(random.choice("ATGC") for i in range(length))

# Runs all the analysis functions with randomly generated sequences
def run_all():
  seq1 = random_sequence(20)
  seq2 = random_sequence(20)
  print("\n-------1) EDIT DISTANCE------------------------------------------")
  try:
    print(f"The Edit Distance of {seq1} and {seq2} is {edit_distance(seq1, seq2)}")
  except Exception as e:
    print("Error:", e)
  print("-----------------------------------------------------------------\n")

  print("\n-------2) IN SILICO PCR------------------------------------------")
  sequence = random_sequence(40)
  forward_primer = sequence[random.randint(0, 10):random.randint(12, 17):1]
  reverse_primer = reverse_complement(sequence[random.randint(22, 27):random.randint(35, 39):1])
  try:
    product = pcr_product(forward_primer, reverse_primer, sequence)
    print(f"The PCR product of {sequence} using primers, {forward_primer} and {reverse_primer} is {product}")
  except Exception as e:
    print("Error:", e)
  print("-----------------------------------------------------------------\n")

  print("\n-------3) MOST FREQUENT K-MERS IN A LIST OF SEQUENCES------------")
  sequences = [random_sequence(5), random_sequence(5), random_sequence(5)]
  kmer_freq = most_frequent_kmer(sequences, k = random.randint(2, 4))
  formatted_sequences = ", ".join(sequences)
  formatted_kmers = ", ".join(f"{k}" for k in kmer_freq[1]) + f" ({kmer_freq[0]} occurences)"
  print(f"In the sequences {formatted_sequences}, the most frequent kmer(s) are: {formatted_kmers}")
  print("-----------------------------------------------------------------\n")

  print("\n-------4) TRANSLATION--------------------------------------------")
  for i in range(10):
    curr_seq = ''.join(random.choice("ATGC") for i in range(200))
    try:
      curr_protein = dna_to_protein(curr_seq)
      print(f"{curr_seq} codes for {curr_protein}")
    except Exception as e:
      print(f"Unable to translate {curr_seq}: {e}")
  print("-----------------------------------------------------------------\n")
   

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="In silico analysis tools")
    subparsers = parser.add_subparsers(dest="command", required=False)

    parser_edit = subparsers.add_parser("edit_distance", help="Calculate edit distance between two sequences")
    parser_edit.add_argument("sequence1", help="First DNA sequence")
    parser_edit.add_argument("sequence2", help="Second DNA sequence")

    parser_pcr = subparsers.add_parser("pcr", help="Return the predicted PCR product")
    parser_pcr.add_argument("forward", help="The forward primer")
    parser_pcr.add_argument("reverse", help="The reverse primer")
    parser_pcr.add_argument("sequence", help="The template DNA")

    parser_kmer = subparsers.add_parser("most_frequent_kmer", help="Find the most frequent k-mer in a list of sequences")
    parser_kmer.add_argument("sequences", nargs="+", help="List of DNA sequences (space-separated)")
    parser_kmer.add_argument("k", type=int, help="Length of the k-mer (integer)")

    parser_translator = subparsers.add_parser("translate", help="Translate the given DNA to a protein")
    parser_translator.add_argument("sequence", help="DNA sequence (coding strand)")


    args = parser.parse_args()

    if args.command == None:
       run_all()
    elif args.command == "edit_distance":
      print(edit_distance(args.sequence1, args.sequence2))
    elif args.command == "pcr":
      print(pcr_product(args.forward, args.reverse, args.sequence))
    elif args.command == "most_frequent_kmer":
      print(most_frequent_kmer(args.sequences, args.k))
    elif args.command == "translate":
      print(dna_to_protein(args.sequence))
