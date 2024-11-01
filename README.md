## fastq2gvcf

### Run

```bash
nextflow run main.nf \
-profile docker \
--fastq /data/GL/nextflow/namxle/samples/SRR2098616.fastq.gz \
--fasta /data/GL/nextflow/namxle/samples/SRR2098616.fasta.gz \
--sample_id SRR2098616 \
--sample_name SRR2098616 \
--outdir /data/GL/nextflow/namxle/results \
--docker_registry namxle
```

### Docker

```bash
docker build -t namxle/gatk:4.6.1.0 -f docker/gatk.docker .

docker build -t namxle/ubuntu-samtools:22.04-1.9 -f docker/ubuntu-samtools.docker .
```