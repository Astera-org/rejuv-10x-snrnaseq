process CELLRANGER_AGGREGATE_REPORTS {
    input:
    path(metric_reports), stageAs: 'metrics_summary??.csv'
    val(meta)

    output:
    path "aggregated_metrics.csv"       , emit: csv

    script:
    def sample_names = meta.collect{el -> el.gem}.join(",")

    """
    sed "1q;d" ${metric_reports[0]} | paste -d"," <(echo Sample) - > aggregated_metrics.csv &&\
    cat $metric_reports | sed -n 'n;p' | paste -d"," <(echo ${sample_names} | sed 's/,/\\n/g') - >> aggregated_metrics.csv
    """
}
