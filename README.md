# retrieve-seq

Retrieve the sequence of every ORF from a genome or retrieve sequences relative to the first, last or both nucleotides of all ORFs (e.g. from -100 to +50 relative to the start codon of every ORF).

## Usage

**Input:** FASTA file and .txt file
**Output:** multiFASTA file

**Input .txt file format example**

```
ORF-1 13092 13769
ORF-2 13738 12953
ORF-3 54562 54678
...
```

Run the script:

```
python3 retrieve-seq.py file_name.FASTA file_name.txt [-start/-stop/-startstop rp_from rp_to]
[optional commands]
```

If you don't use the optional commands, the script will retrieve the sequence of every ORF.
Using the optional commands you can retrieve a sequence relative to that nt (e.g. -start -100 +150 will retrieve the sequence from -100 to 150 relative to the first nt of every ORF)
