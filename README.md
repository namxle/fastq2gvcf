## fastq2gvcf

### Run

```bash
nextflow run main.nf \
-profile docker \
--fastq /data/GL/nextflow/namxle/samples/SRR2098616.fastq.gz \
--sample_id SRR2098616 \
--sample_name SRR2098616 \
--outdir /data/GL/nextflow/namxle/results \
--docker_registry namxle
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