#!/usr/bin/env nextflow

/*
 *
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */


nextflow.enable.types = true

/*
 * Pipeline parameters. Each parameter can be overridden on the command line (e.g. `--query query.fa`)
 */
params {
    // Query sequences (FASTA)
    query: Path = "${projectDir}/data/sample.fa"

    // BLAST database (directory and database name)
    db: String = "${projectDir}/blast-db/pdb/tiny"

    // Name of the output file with the matching sequences
    out: String = "result.txt"

    // Number of sequences per query chunk
    chunkSize: Integer = 100
}


// Entry workflow
workflow {
    main:
    db_dir = file(params.db).parent
    db_name = file(params.db).name

    // Split the query file into chunks of `chunkSize` sequences
    ch_fasta = channel.of(params.query)
        .flatMap { fasta -> fasta.splitFasta(by: params.chunkSize, file: true) }

    // Run BLAST on each chunk and emit the top hits
    ch_hits = blast(ch_fasta, db_dir, db_name)

    // Extract the sequences of the top hits from the database
    ch_sequences = extract(ch_hits, db_dir, db_name)

    // Merge the matching sequences into a single file and print it
    val_result = merge_sequences(ch_sequences.collect(), params.out)
    val_result.view { result -> "matching sequences:\n${result.text}" }

    publish:
    result = val_result
}


// Pipeline outputs (published to the `results` directory by default)
output {
    result {
    }
}


// Search each query chunk against the database and keep the IDs of the top 10 hits
process blast {
    input:
    query: Path
    db_dir: Path
    db_name: String

    output:
    file('top_hits')

    script:
    """
    blastp -db ${db_dir}/${db_name} -query ${query} -outfmt 6 > blast_result
    cat blast_result | head -n 10 | cut -f 2 > top_hits
    """
}


// Fetch the sequences of the top hits from the database
process extract {
    input:
    top_hits: Path
    db_dir: Path
    db_name: String

    output:
    file('sequences')

    script:
    """
    blastdbcmd -db ${db_dir}/${db_name} -entry_batch ${top_hits} | head -n 10 > sequences
    """
}


// Concatenate the matching sequences into a single file
process merge_sequences {
    input:
    sequences: Bag<Path>
    name: String

    stage:
    stageAs sequences, 'sequences_*'

    output:
    file(name)

    script:
    """
    cat sequences_* > ${name}
    """
}
