/*
 * Alignment with Cellranger
 */

include {CELLRANGER_COUNT} from "../../modules/nf-core/modules/cellranger/count/main.nf"
include {CELLRANGER_AGGREGATE_REPORTS} from "../../modules/local/cellranger_aggregate_reports.nf"

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_PIPELINE {
    take:
        cellranger_reference
        ch_fastq

    main:
        ch_versions = Channel.empty()

        // Obtain read counts
        CELLRANGER_COUNT (
            ch_fastq.map{ meta, reads -> [meta + ["gem": meta.id, "samples": [meta.id]], reads]},
            cellranger_reference
        )
        ch_versions = ch_versions.mix(CELLRANGER_COUNT.out.versions)

        CELLRANGER_AGGREGATE_REPORTS (
            CELLRANGER_COUNT.out.metrics_summary.collect(),
            CELLRANGER_COUNT.out.meta.collect()
        )

    emit:
        ch_versions
        cellranger_out  = CELLRANGER_COUNT.out.outs
        aggregated_report = CELLRANGER_AGGREGATE_REPORTS.out.csv
}
