## 去除冗余序列

```
from Bio import SeqIO

# 输入文件名
input_file = "2020-1-3.fasta"
output_file = "2020-1-3_unique.fasta"

# 使用集合来存储唯一序列
unique_sequences = set()

# 读取FASTA文件
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        seq_str = str(record.seq)
        if seq_str not in unique_sequences:
            unique_sequences.add(seq_str)
            # 将唯一序列写入新的FASTA文件
            SeqIO.write(record, outfile, "fasta")

print(f"冗余序列已删除，唯一序列已保存到 {output_file}")

```

## 安装多序列比对软件
```
conda install -c bioconda clustalo  
```
检查安装是否成功
```
clustalo -h
```
## 多序列比对
```
clustalo -i sequences.fasta -o aligned_sequences.aln --outfmt=clu
```

## 解析比对结果
为了生成包含变异的序列文件并删除冗余的序列。  
识别变异位点。  
创建包含变异的序列文件。  
删除相同冗余的序列。  

```
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


```
