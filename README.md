# Blast example

A basic Blast pipeline using Nextflow, written with [static types](https://docs.seqera.io/nextflow/static-typing).

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A526.10-brightgreen.svg)](http://nextflow.io)

## Get started

Install Nextflow

```bash
curl https://get.nextflow.io | bash
```

Run the pipeline:

```bash
nextflow run nextflow-io/blast-example -profile docker
```

The matching sequences are printed to the console and saved to `results/result.txt`.

## Parameters

| Parameter     | Default                  | Description                                  |
|---------------|--------------------------|----------------------------------------------|
| `--query`     | `data/sample.fa`         | Query sequences (FASTA)                      |
| `--db`        | `blast-db/pdb/tiny`      | BLAST database (directory and database name) |
| `--out`       | `result.txt`             | Name of the output file                      |
| `--chunkSize` | `100`                    | Number of sequences per query chunk          |

## Dependencies

* Nextflow 26.10 or later
* Java 17 or later
* Docker
