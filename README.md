## fastq2gvcf

### Prequisite

1. Install [nextflow](https://www.nextflow.io/docs/latest/install.html).
2. Install [docker](https://docs.docker.com/engine/install).
3. Build docker images:

```bash
docker build -t namxle/gatk:4.6.1.0 -f docker/gatk.docker .

docker build -t namxle/ubuntu-tools:22.04-1.9 -f docker/ubuntu-tools.docker .
```

4. Download `SRR2098616.fastq.gz` to `examples` folder. From [here](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2098616&display=download).

### Run

Run example below:

```bash
rm -rf results && \
nextflow run main.nf \
-profile test \
--fastq examples/SRR2098616.fastq.gz \
--sample_id SRR2098616 \
--sample_name SRR2098616 \
--outdir results \
--docker_registry namxle
```
