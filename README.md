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