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


/*
 * Defines the pipeline parameters (with a default value for each one)
 * Each parameter can be specified on the command line (e.g. `--query query.fa`)
 */
params.query = "$projectDir/data/sample.fa"
params.db = "$projectDir/blast-db/pdb/tiny"
params.out = "result.txt"
params.chunkSize = 100


workflow {
    db_dir = file(params.db).parent
    db_name = file(params.db).name

    /*
     * Create a channel emitting the given query fasta file(s).
     * Split the file into chunks containing as many sequences as defined by the parameter 'chunkSize'.
     * Finally, assign the resulting channel to the variable 'ch_fasta'
     */
    ch_fasta = channel.fromPath(params.query)
        .splitFasta(by: params.chunkSize, file: true)

    /*
     * Execute a BLAST job for each chunk emitted by the 'ch_fasta' channel
     * and emit the resulting BLAST matches.
     */
    ch_hits = blast(ch_fasta, db_dir, db_name)

    /*
     * Each time a file emitted by the 'blast' process, an extract job is executed,
     * producing a file containing the matching sequences.
     */
    ch_sequences = extract(ch_hits, db_dir, db_name)

    /*
     * Collect all the sequences files into a single file
     * and print the resulting file contents when complete.
     */
    ch_sequences
        .collectFile(name: params.out)
        .view { file -> "matching sequences:\n ${file.text}" }
}


process blast {
    input:
    path 'query.fa'
    path db_dir
    val db_name

    output:
    path 'top_hits'

    script:
    """
    blastp -db $db_dir/$db_name -query query.fa -outfmt 6 > blast_result
    cat blast_result | head -n 10 | cut -f 2 > top_hits
    """
}


process extract {
    input:
    path 'top_hits'
    path db_dir
    val db_name

    output:
    path 'sequences'

    script:
    """
    blastdbcmd -db $db_dir/$db_name -entry_batch top_hits | head -n 10 > sequences
    """
}
