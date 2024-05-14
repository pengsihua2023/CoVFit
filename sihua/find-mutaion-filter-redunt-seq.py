from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

def find_variants(alignment):
    reference = alignment[0]
    variants = defaultdict(list)

    for record in alignment[1:]:
        for i in range(len(reference.seq)):
            if reference.seq[i] != record.seq[i]:
                variants[record.id].append((i+1, reference.seq[i], record.seq[i]))
    
    return variants

def remove_redundant_sequences(alignment):
    unique_sequences = {}
    for record in alignment:
        seq_str = str(record.seq)
        if seq_str not in unique_sequences:
            unique_sequences[seq_str] = record
    return list(unique_sequences.values())

def write_variants_to_file(variants, filename):
    with open(filename, "w") as f:
        for seq_id, changes in variants.items():
            f.write(f">{seq_id}\n")
            for change in changes:
                f.write(f"Position {change[0]}: {change[1]} -> {change[2]}\n")
            f.write("\n")

def main():
    alignment = AlignIO.read("USA_0401_0514_aligned.aln", "clustal")

    # Find variants
    variants = find_variants(alignment)

    # Write variants to file
    write_variants_to_file(variants, "USA_0401-0514_variants.txt")

    # Remove redundant sequences
    unique_sequences = remove_redundant_sequences(alignment)

    # Write unique sequences to a new file
    with open("unique_sequences_USA_04-1_0514.fasta", "w") as output_handle:
        SeqIO.write(unique_sequences, output_handle, "fasta")

if __name__ == "__main__":
    main()

