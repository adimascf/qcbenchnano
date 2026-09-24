rule plot_main_figure_assembly:
	input:
		quast_csv=TABLES / "assess/assembly/metrics/combo_quast_compiled_metrics.csv",
		missed_csv=TABLES / "assess/assembly/metrics/combo_assembly_missed_contigs.csv",
		contam_csv=TABLES / "assess/assembly/metrics/combo_contaminant_summary_count.csv"
	log:
		LOGS / "assess/assembly/plot_main_figure_assembly_{model}.log"
	resources:
		mem="16GiB",
		runtime="20m"
	conda:
		ENVS / "generate_figure_python.yaml"
	output:
		figure=FIGURES / "assess/assembly/metrics/combo_main_figure_assembly_{model}.png"
	script:
		"../scripts/plot_main_figure_assembly.py"


rule plot_main_figure_call:
	input:
		csv=TABLES / "assess/call/metrics/combo_variant_summary.csv"
	log:
		LOGS / "assess/call/plot_main_figure_call_{model}.log"
	resources:
		mem="16GiB",
		runtime="20m"
	conda:
		ENVS / "generate_figure_python.yaml"
	output:
		figure=FIGURES / "assess/call/metrics/combo_main_figure_call_{model}.png"
	script:
		"../scripts/plot_main_figure_call.py"


rule plot_supp_figure_call:
	input:
		csv=TABLES / "assess/call/metrics/combo_variant_summary.csv"
	log:
		LOGS / "assess/call/plot_supp_figure_call_{model}.log"
	resources:
		mem="16GiB",
		runtime="20m"
	conda:
		ENVS / "generate_figure_python.yaml"
	output:
		figure=FIGURES / "assess/call/metrics/combo_supp_figure_call_{model}.png"
	script:
		"../scripts/plot_supp_figure_call.py"


rule compile_assembly_metrics_master:
	input:
		quast=rules.compile_quast_metrics.output.csv,
		contaminants=rules.aggregate_assembly_contam.output.summary,
		missed_contigs=rules.plot_missed_contig.output.table
	log:
		LOGS / "assess/assembly/compile_assembly_metrics_master.log"
	resources:
		mem="8GiB",
		runtime="10m"
	conda:
		ENVS / "generate_figure_python.yaml"
	output:
		master=TABLES / "assess/assembly/metrics/combo_assembly_metrics_master.csv"
	script:
		"../scripts/generate_master_csv.py"


rule plot_assembly_aggregate_score:
	input:
		master_csv=rules.compile_assembly_metrics_master.output.master,
		survey_csv=DATA / "microbial-qc-survey.csv"
	log:
		LOGS / "assess/assembly/plot_assembly_aggregate_score.log"
	resources:
		mem="16GiB",
		runtime="20m"
	conda:
		ENVS / "generate_figure_python.yaml"
	output:
		fig_perdepth_sup=FIGURES / "assess/assembly/metrics/combo_assembly_aggregate_score_perdepth_sup.png",
		fig_perdepth_hac=FIGURES / "assess/assembly/metrics/combo_assembly_aggregate_score_perdepth_hac.png",
		scores_csv_perdepth=TABLES / "assess/assembly/metrics/combo_assembly_survey_scores_perdepth.csv",
		summary_csv_perdepth=TABLES / "assess/assembly/metrics/combo_assembly_survey_ranking_summary_perdepth.csv"
	script:
		"../scripts/plot_assembly_aggregate_score.py"


rule summarise_assembly_metrics:
	input:
		master_csv=rules.compile_assembly_metrics_master.output.master
	log:
		LOGS / "assess/assembly/summarise_assembly_metrics_{model}.log"
	resources:
		mem="8GiB",
		runtime="10m"
	conda:
		ENVS / "generate_figure_python.yaml"
	output:
		csv=TABLES / "assess/assembly/metrics/combo_assembly_summary_metrics_{model}.csv"
	script:
		"../scripts/summarise_assembly_metrics.py"


