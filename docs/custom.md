### Preprocess GTF file

```bash
less GCF_009035845.1_ASM903584v1_genomic.sorted.gtf | grep -v "#" | awk -F"\t" '{if($3 == "gene"){ start = $4; end = $5; split($9, a, "; "); for (i in a){ split(a[i], b, " "); if (b[1] == "gene_id"){ split(b[2], c, "\""); gene = c[2]; } } print $1"\t"$4"\t"$5"\t"gene } }' > genes.bed
```

### Test

```bash
docker run -itv /data/GL/nextflow/namxle/samples:/workspace/samples namxle/ubuntu-tools:22.04-1.9 bash

seqtk seq -1 samples/SRR2098616.fastq.gz | bgzip -@ 10  > SRR2098616_R1.fastq.gz
seqtk seq -2 samples/SRR2098616.fastq.gz | bgzip -@ 10 > SRR2098616_R2.fastq.gz
```
