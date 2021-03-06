SELECT
--gene_expression_analysis_results.*,
gene_expression_analysis_results.container,
gene_expression_analysis.analysis_accession,
feature_id,
gene_symbol,
adj_p_val,
ave_expr,
log_fc,
p_value,
statistic,
gene_expression_analysis.arm_name AS cohort,
gene_expression_analysis.coefficient
FROM
gene_expression_analysis,
gene_expression_analysis_results
WHERE
gene_expression_analysis_results.container IN (SELECT Container FROM study.participant)
AND gene_expression_analysis.arm_name IN (SELECT name FROM study.cohort_membership)
AND gene_expression_analysis_results.analysis_accession = gene_expression_analysis.analysis_accession
AND gene_expression_analysis_results.container = gene_expression_analysis.container
