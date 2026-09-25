# Singleton host-filtering correction

## Issue

The original `remove_host_sing` rule filtered single-end reads with:

```bash
samtools view -b -f 12 -F 256
```

FLAG 12 combines:

- `4`: read unmapped
- `8`: mate unmapped

The singleton reads are passed to Bowtie2 with `-U`, so they are unpaired and do not have a mate. Requiring FLAG 8 therefore removed all singleton reads after host filtering.

The rule was corrected to:

```bash
samtools view -b -f 4 -F 256
```

which retains single-end reads that do not map to the Bos taurus reference.

## Verification

The previously generated cleaned singleton files were checked with:

```bash
column -t -s $'\t' tables/read_stats/seqkit_cleaned_reads_sing.tsv | head
```

Example result:

```text
file                                              num_seqs  sum_len
snakestream/reads_clean/3550_sing_clean.fastq.gz  0         0
snakestream/reads_clean/3551_sing_clean.fastq.gz  0         0
snakestream/reads_clean/3552_sing_clean.fastq.gz  0         0
```

The trimmed singleton files produced by BBDuk were not empty. They were checked with:

```bash
apptainer exec \
  --bind /scratch:/scratch \
  /archive/extra/daniel.remondini2/bioinf_db/containers/seqkit/2.13.0/seqkit_2.13.0.sif \
  seqkit stats -T \
  /scratch/applicata/alessandro.fuschi2/milk-sewage-wgs/data/reads_trim/*_sing.fastq.gz
```

Representative results:

```text
sample   num_seqs
3550     897
3551     1087
3552     1209
3558     4914
4425     4179
4426     4273
```

Across the 66 rows visible in the terminal output, there were 155,112 singleton reads in total, with an average of about 2,350 reads per sample and a maximum of 4,914 reads.

## Consequence for existing assemblies

MEGAHIT was configured to use singleton reads through `-r`, but the cleaned singleton files were empty. Existing assemblies were therefore generated from the correctly host-filtered paired-end reads without the singleton contribution.

Because the singleton fraction was very small compared with the paired-end read sets (hundreds of MB to >1 GB per mate file), the existing assemblies were retained rather than recomputed solely for this issue.

Future runs use the corrected `-f 4` filter.
