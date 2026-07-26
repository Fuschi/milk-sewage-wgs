# Pipeline workspace

Questa cartella contiene tutti i file pesanti e rigenerabili della pipeline ed è
ignorata da Git.

I FASTQ grezzi devono essere copiati in:

```text
pipeline/workspace/reads_raw/{sample_id}_R1.fastq.gz
pipeline/workspace/reads_raw/{sample_id}_R2.fastq.gz
```

Durante l'esecuzione Snakemake creerà qui anche:

```text
reads_trim/
reads_clean/
host/
qc/
logs/
.snakemake/
```

Le sole tabelle leggere destinate a Git vengono scritte in `data/raw/`.
