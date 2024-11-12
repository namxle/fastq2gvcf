## fastq2gvcf

### Run

```bash
rm -rf /data/GL/nextflow/namxle/results && \
nextflow run main.nf \
-profile docker \
--fastq /data/GL/nextflow/namxle/samples/SRR2098616.fastq.gz \
--sample_id SRR2098616 \
--sample_name SRR2098616 \
--outdir /data/GL/nextflow/namxle/results \
--docker_registry namxle
```

### Preprocess GTF file

```bash
less GCF_009035845.1_ASM903584v1_genomic.sorted.gtf | grep -v "#" | awk -F"\t" '{if($3 == "gene"){ start = $4; end = $5; split($9, a, "; "); for (i in a){ split(a[i], b, " "); if (b[1] == "gene_id"){ split(b[2], c, "\""); gene = c[2]; } } print $1"\t"$4"\t"$5"\t"gene } }' > genes.bed
```

### Docker

```bash
docker build -t namxle/gatk:4.6.1.0 -f docker/gatk.docker .

docker build -t namxle/ubuntu-tools:22.04-1.9 -f docker/ubuntu-tools.docker .
```

### Test

```bash
docker run -itv /data/GL/nextflow/namxle/samples:/workspace/samples namxle/ubuntu-tools:22.04-1.9 bash

seqtk seq -1 samples/SRR2098616.fastq.gz | bgzip -@ 10  > SRR2098616_R1.fastq.gz
seqtk seq -2 samples/SRR2098616.fastq.gz | bgzip -@ 10 > SRR2098616_R2.fastq.gz
```