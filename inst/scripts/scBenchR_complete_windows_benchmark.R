###############################################################
# scBenchR COMPLETE WINDOWS MULTISESSION BENCHMARK PIPELINE
# Author: Bacem Saada
#
# End-to-end Seurat v5 benchmark for SMALL, LARGE and MERGED
# objects. The script performs repeated FindAllMarkers runs,
# statistical summaries, reproducibility checks, source-data
# export, and publication-ready figure generation.
###############################################################

###############################
# 0. USER CONFIGURATION
###############################
config <- list(
  project_dir = "D:/scRNAseq",
  output_dir = "D:/scRNAseq/revision_benchmark/scBenchR_complete_run",
  objects = c(
    SMALL  = "D:/scRNAseq/output/seu_small_step3.rds",
    LARGE  = "D:/scRNAseq/output/seu_large_step3.rds",
    MERGED = paste0("D:/scRNAseq/output/integrated/",
                    "seu_merged_integrated_step3_RNA_joined.rds")
  ),
  workers = c(1, 2, 4, 6, 8),
  replicates = 5,
  assay = "RNA",
  test_use = "wilcox",
  only_positive = TRUE,
  seed = 123,
  bootstrap_repetitions = 2000,
  confidence_level = 0.95,
  dpi = 600
)

###############################
# 1. PACKAGES
###############################
required_packages <- c(
  "Seurat", "future", "parallelly", "peakRAM", "ggplot2",
  "dplyr", "tidyr", "readr", "purrr", "patchwork",
  "scales", "viridis", "rstatix", "lobstr", "sessioninfo"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}
invisible(lapply(required_packages, install_if_missing))

suppressPackageStartupMessages({
  library(Seurat); library(future); library(parallelly)
  library(peakRAM); library(ggplot2); library(dplyr)
  library(tidyr); library(readr); library(purrr)
  library(patchwork); library(scales); library(viridis)
  library(rstatix); library(lobstr); library(sessioninfo)
})

###############################
# 2. DIRECTORIES AND LOGGING
###############################
dirs <- list(
  root = config$output_dir,
  environment = file.path(config$output_dir, "00_environment"),
  characteristics = file.path(config$output_dir, "01_dataset_characteristics"),
  raw = file.path(config$output_dir, "02_raw_results"),
  statistics = file.path(config$output_dir, "03_statistics"),
  figures = file.path(config$output_dir, "04_main_figures"),
  supplementary = file.path(config$output_dir, "05_supplementary_figures"),
  tables = file.path(config$output_dir, "06_tables"),
  logs = file.path(config$output_dir, "07_logs"),
  manuscript = file.path(config$output_dir, "08_manuscript_ready"),
  source_data = file.path(config$output_dir, "09_figure_source_data")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
log_file <- file.path(dirs$logs, "scBenchR_complete_pipeline.log")

write_log <- function(message) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

safe_write_csv <- function(x, path) readr::write_csv(x, path, na = "")
set.seed(config$seed)

###############################
# 3. ENVIRONMENT CAPTURE
###############################
environment_summary <- tibble(
  Parameter = c("Operating system", "R version", "Seurat version",
                "SeuratObject version", "future version",
                "Available logical cores", "Workers", "Replicates"),
  Value = c(paste(Sys.info()[["sysname"]], Sys.info()[["release"]]),
            R.version.string, as.character(packageVersion("Seurat")),
            as.character(packageVersion("SeuratObject")),
            as.character(packageVersion("future")),
            parallelly::availableCores(),
            paste(config$workers, collapse = ", "), config$replicates)
)
safe_write_csv(environment_summary,
               file.path(dirs$environment, "system_information.csv"))
capture.output(sessioninfo::session_info(),
               file = file.path(dirs$environment, "session_info.txt"))
saveRDS(config, file.path(dirs$environment, "analysis_config.rds"))

###############################
# 4. LOAD AND VALIDATE OBJECTS
###############################
load_object <- function(name, path) {
  if (!file.exists(path)) stop("Missing object: ", path)
  obj <- readRDS(path)
  if (!inherits(obj, "Seurat")) stop(name, " is not a Seurat object")
  if (!(config$assay %in% Assays(obj))) stop(name, " lacks RNA assay")
  if (name == "MERGED") {
    DefaultAssay(obj) <- config$assay
    layers <- Layers(obj[[config$assay]])
    if (!all(c("counts", "data") %in% layers)) {
      stop("MERGED must contain joined RNA layers: counts and data")
    }
  }
  if (length(levels(Idents(obj))) < 2) stop(name, " has fewer than two clusters")
  write_log(paste(name, ncol(obj), "cells", length(levels(Idents(obj))), "clusters"))
  obj
}
datasets <- imap(config$objects, load_object)

###############################
# 5. DATASET CHARACTERISTICS
###############################
dataset_characteristics <- imap_dfr(datasets, function(obj, name) {
  cs <- as.numeric(table(Idents(obj)))
  tibble(
    Dataset = name,
    Cells = ncol(obj),
    Features = nrow(obj[[config$assay]]),
    Clusters = length(levels(Idents(obj))),
    Smallest_cluster = min(cs),
    Median_cluster_size = median(cs),
    Largest_cluster = max(cs),
    Default_assay = DefaultAssay(obj),
    Available_assays = paste(Assays(obj), collapse = ", "),
    Object_size_MB = round(as.numeric(lobstr::obj_size(obj)) / 1024^2, 2)
  )
})
safe_write_csv(dataset_characteristics,
               file.path(dirs$tables, "Table_S1_dataset_characteristics.csv"))

###############################
# 6. BENCHMARK FUNCTIONS
###############################
run_markers <- function(obj) {
  FindAllMarkers(
    object = obj,
    assay = config$assay,
    only.pos = config$only_positive,
    test.use = config$test_use,
    verbose = FALSE
  )
}

run_one <- function(obj, dataset, workers, replicate) {
  write_log(paste("START", dataset, "workers", workers, "replicate", replicate))
  tryCatch({
    if (workers == 1) plan(sequential) else plan(multisession, workers = workers)
    gc()
    start <- Sys.time()
    mem <- peakRAM({ markers <- run_markers(obj) })
    end <- Sys.time()
    out <- tibble(
      Dataset = dataset,
      Workers = workers,
      Replicate = replicate,
      Start_time = as.character(start),
      End_time = as.character(end),
      Runtime_seconds = as.numeric(difftime(end, start, units = "secs")),
      Peak_RAM_MB = mem$Peak_RAM_Used_MiB,
      Markers_detected = nrow(markers),
      Unique_marker_genes = n_distinct(markers$gene),
      Status = "SUCCESS",
      Error = ""
    )
    rm(markers); gc(); out
  }, error = function(e) {
    tibble(Dataset = dataset, Workers = workers, Replicate = replicate,
           Start_time = NA_character_, End_time = as.character(Sys.time()),
           Runtime_seconds = NA_real_, Peak_RAM_MB = NA_real_,
           Markers_detected = NA_real_, Unique_marker_genes = NA_real_,
           Status = "FAILED", Error = conditionMessage(e))
  })
}

###############################
# 7. CHECKPOINTED BENCHMARK GRID
###############################
checkpoint <- file.path(dirs$raw, "marker_runtime_checkpoint.csv")
results <- if (file.exists(checkpoint)) {
  read_csv(checkpoint, show_col_types = FALSE)
} else tibble()

for (dataset in names(datasets)) {
  for (workers in config$workers) {
    for (replicate in seq_len(config$replicates)) {
      done <- nrow(results) > 0 && any(
        results$Dataset == dataset & results$Workers == workers &
          results$Replicate == replicate
      )
      if (done) next
      results <- bind_rows(results,
                           run_one(datasets[[dataset]], dataset, workers, replicate)) %>%
        arrange(factor(Dataset, levels = names(datasets)), Workers, Replicate)
      safe_write_csv(results, checkpoint)
      safe_write_csv(filter(results, Status == "FAILED"),
                     file.path(dirs$raw, "errors_checkpoint.csv"))
    }
  }
}
plan(sequential)
safe_write_csv(results, file.path(dirs$raw, "marker_runtime_final.csv"))
safe_write_csv(filter(results, Status == "FAILED"),
               file.path(dirs$raw, "errors_final.csv"))

successful <- filter(results, Status == "SUCCESS")
expected_runs <- length(datasets) * length(config$workers) * config$replicates
if (nrow(successful) != expected_runs) {
  warning("Expected ", expected_runs, " successful runs; found ", nrow(successful))
}

###############################
# 8. BOOTSTRAP CI AND RUNTIME SUMMARY
###############################
bootstrap_median_ci <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(c(NA_real_, NA_real_))
  medians <- replicate(config$bootstrap_repetitions,
                       median(sample(x, length(x), replace = TRUE)))
  alpha <- 1 - config$confidence_level
  quantile(medians, c(alpha / 2, 1 - alpha / 2), names = FALSE)
}

runtime_basic <- successful %>%
  group_by(Dataset, Workers) %>%
  summarise(
    N = n(),
    Mean_runtime_seconds = mean(Runtime_seconds),
    Median_runtime_seconds = median(Runtime_seconds),
    SD_runtime_seconds = sd(Runtime_seconds),
    IQR_runtime_seconds = IQR(Runtime_seconds),
    Minimum_runtime_seconds = min(Runtime_seconds),
    Maximum_runtime_seconds = max(Runtime_seconds),
    Range_runtime_seconds = max(Runtime_seconds) - min(Runtime_seconds),
    CV_percent = 100 * sd(Runtime_seconds) / mean(Runtime_seconds),
    Mean_runtime_minutes = mean(Runtime_seconds) / 60,
    Median_runtime_minutes = median(Runtime_seconds) / 60,
    .groups = "drop"
  )

runtime_ci <- successful %>%
  group_by(Dataset, Workers) %>%
  summarise(CI = list(bootstrap_median_ci(Runtime_seconds)), .groups = "drop") %>%
  mutate(Median_CI_lower_seconds = map_dbl(CI, 1),
         Median_CI_upper_seconds = map_dbl(CI, 2)) %>%
  select(-CI)

baseline <- runtime_basic %>%
  filter(Workers == 1) %>%
  select(Dataset, Baseline_median_seconds = Median_runtime_seconds)

runtime_summary <- runtime_basic %>%
  left_join(runtime_ci, by = c("Dataset", "Workers")) %>%
  left_join(baseline, by = "Dataset") %>%
  mutate(
    Speedup = Baseline_median_seconds / Median_runtime_seconds,
    Parallel_efficiency_percent = 100 * Speedup / Workers,
    Runtime_reduction_percent = 100 *
      (Baseline_median_seconds - Median_runtime_seconds) / Baseline_median_seconds
  ) %>%
  group_by(Dataset) %>%
  arrange(Workers, .by_group = TRUE) %>%
  mutate(
    Previous_median_seconds = lag(Median_runtime_seconds),
    Incremental_runtime_reduction_percent = 100 *
      (Previous_median_seconds - Median_runtime_seconds) / Previous_median_seconds,
    Incremental_speedup_gain = Speedup - lag(Speedup),
    Additional_workers = Workers - lag(Workers),
    Speedup_gain_per_added_worker = Incremental_speedup_gain / Additional_workers
  ) %>% ungroup()

safe_write_csv(runtime_summary,
               file.path(dirs$statistics, "runtime_summary_complete.csv"))
safe_write_csv(runtime_summary,
               file.path(dirs$tables, "Table_S3_runtime_summary.csv"))

###############################
# 9. STATISTICAL TESTS
###############################
stats_input <- successful %>%
  mutate(Worker_factor = factor(Workers, levels = config$workers))
kruskal <- stats_input %>% group_by(Dataset) %>%
  kruskal_test(Runtime_seconds ~ Worker_factor) %>% ungroup()
pairwise <- stats_input %>% group_by(Dataset) %>%
  pairwise_wilcox_test(Runtime_seconds ~ Worker_factor,
                       p.adjust.method = "BH", exact = FALSE) %>% ungroup()
safe_write_csv(kruskal,
               file.path(dirs$statistics, "kruskal_wallis_worker_effect.csv"))
safe_write_csv(pairwise,
               file.path(dirs$tables, "Table_S5_pairwise_worker_tests.csv"))

###############################
# 10. REPRODUCIBILITY, MEMORY, FASTEST CONFIGURATION
###############################
marker_summary <- successful %>%
  group_by(Dataset, Workers) %>%
  summarise(
    Runs = n(),
    Unique_marker_counts = n_distinct(Markers_detected),
    Unique_gene_counts = n_distinct(Unique_marker_genes),
    Marker_count = first(Markers_detected),
    Unique_marker_gene_count = first(Unique_marker_genes),
    Median_runtime_seconds = median(Runtime_seconds),
    Markers_per_second = Marker_count / Median_runtime_seconds,
    Markers_per_minute = 60 * Marker_count / Median_runtime_seconds,
    .groups = "drop"
  )

memory_summary <- successful %>%
  group_by(Dataset, Workers) %>%
  summarise(
    N = n(),
    Median_parent_peak_RAM_MB = median(Peak_RAM_MB, na.rm = TRUE),
    Mean_parent_peak_RAM_MB = mean(Peak_RAM_MB, na.rm = TRUE),
    SD_parent_peak_RAM_MB = sd(Peak_RAM_MB, na.rm = TRUE),
    IQR_parent_peak_RAM_MB = IQR(Peak_RAM_MB, na.rm = TRUE),
    .groups = "drop"
  )

fastest <- runtime_summary %>% group_by(Dataset) %>%
  slice_min(Median_runtime_seconds, n = 1, with_ties = FALSE) %>% ungroup()

safe_write_csv(successful, file.path(dirs$tables, "Table_S2_raw_runtime_benchmark.csv"))
safe_write_csv(runtime_summary %>% select(Dataset, Workers, Median_runtime_seconds,
  Speedup, Parallel_efficiency_percent, Runtime_reduction_percent,
  Incremental_runtime_reduction_percent, Speedup_gain_per_added_worker),
  file.path(dirs$tables, "Table_S4_parallel_performance.csv"))
safe_write_csv(fastest, file.path(dirs$tables, "Table_S6_fastest_tested_configuration.csv"))
safe_write_csv(marker_summary,
               file.path(dirs$tables, "Table_S7_marker_reproducibility_throughput.csv"))
safe_write_csv(memory_summary,
               file.path(dirs$tables, "Table_S8_parent_process_memory.csv"))

###############################
# 11. FIGURE HELPERS
###############################
publication_theme <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(plot.title = element_text(face = "bold", size = base_size + 2),
          plot.subtitle = element_text(size = base_size, margin = margin(b = 8)),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(colour = "black"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom",
          legend.title = element_text(face = "bold"))
}

save_figure <- function(plot, stem, width, height, directory = dirs$figures) {
  ggsave(file.path(directory, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = config$dpi, bg = "white")
  ggsave(file.path(directory, paste0(stem, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(file.path(directory, paste0(stem, ".tiff")), plot,
         width = width, height = height, dpi = config$dpi,
         compression = "lzw", bg = "white")
}

###############################
# 12. MAIN FIGURES 1-5
###############################
# Figure 1
nodes <- tibble(
  Node = c("Public PDAC\nscRNA-seq datasets", "Quality control\nand preprocessing",
           "RPCA integration\nand clustering", "Joined RNA-assay\nmarker detection",
           "Worker grid\n1, 2, 4, 6, 8", "Five repeated\nexecutions",
           "Runtime, speedup,\nefficiency, stability"),
  X = c(1, 2.75, 4.5, 6.25, 8, 9.75, 11.5),
  W = c(1.48, 1.48, 1.48, 1.58, 1.35, 1.40, 1.58), H = 0.42
) %>% mutate(xmin = X-W/2, xmax = X+W/2, ymin = -H/2, ymax = H/2)
arrows <- tibble(x = head(nodes$xmax, -1)+0.07,
                 xend = tail(nodes$xmin, -1)-0.07, y = 0, yend = 0)
f1a <- ggplot() + geom_segment(data=arrows, aes(x,xend=xend,y=y,yend=yend),
  arrow=grid::arrow(length=grid::unit(.105,"inches"), type="closed")) +
  geom_rect(data=nodes,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="white") +
  geom_text(data=nodes,aes(X,0,label=Node),size=3.3) + theme_void()
f1b <- ggplot(dataset_characteristics,aes(Dataset,Cells,fill=Dataset))+
  geom_col(show.legend=FALSE)+geom_text(aes(label=comma(Cells)),vjust=-.4)+
  scale_fill_viridis_d()+scale_y_continuous(labels=comma,expand=expansion(mult=c(0,.15)))+
  labs(title="B  Cell counts",x=NULL,y="Cells")+publication_theme()
f1c <- ggplot(dataset_characteristics,aes(Dataset,Clusters,fill=Dataset))+
  geom_col(show.legend=FALSE)+geom_text(aes(label=Clusters),vjust=-.4)+
  scale_fill_viridis_d()+labs(title="C  Cluster counts",x=NULL,y="Clusters")+publication_theme()
fig1 <- f1a/(f1b|f1c)+plot_layout(heights=c(.48,1.52))+
  plot_annotation(title="Figure 1. Benchmark design and dataset overview")
save_figure(fig1,"Fig1_benchmark_design_dataset_overview",14,8)

# Figure 2
f2a <- successful %>% mutate(Workers=factor(Workers,levels=config$workers)) %>%
  ggplot(aes(Workers,Runtime_seconds/60,fill=Workers))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=.1)+facet_wrap(~Dataset,scales="free_y",nrow=1)+
  scale_fill_viridis_d()+labs(title="A  Replicate-level runtime distributions",
  x="Workers",y="Runtime (minutes)")+publication_theme()+theme(legend.position="none")
f2b <- ggplot(runtime_summary,aes(Workers,Median_runtime_seconds/60,colour=Dataset,group=Dataset))+
  geom_ribbon(aes(ymin=Median_CI_lower_seconds/60,ymax=Median_CI_upper_seconds/60,fill=Dataset),alpha=.14,colour=NA)+
  geom_line()+geom_point()+scale_x_continuous(breaks=config$workers)+
  scale_colour_viridis_d()+scale_fill_viridis_d()+
  labs(title="B  Median runtime with bootstrap 95% confidence intervals",
       x="Workers",y="Median runtime (minutes)")+publication_theme()
fig2 <- f2a/f2b+plot_annotation(title="Figure 2. Repeated runtime benchmarking")
save_figure(fig2,"Fig2_repeated_runtime_benchmarking",13,10)

# Figure 3
f3a <- ggplot(runtime_summary,aes(Workers,Speedup,colour=Dataset,group=Dataset))+
  geom_abline(slope=1,intercept=0,linetype="dashed")+geom_line()+geom_point()+
  scale_x_continuous(breaks=config$workers)+scale_colour_viridis_d()+
  labs(title="A  Observed speedup",x="Workers",y="Speedup relative to one worker")+publication_theme()
f3b <- ggplot(runtime_summary,aes(Workers,Parallel_efficiency_percent,colour=Dataset,group=Dataset))+
  geom_hline(yintercept=100,linetype="dashed")+geom_line()+geom_point()+
  scale_x_continuous(breaks=config$workers)+scale_y_continuous(labels=label_percent(scale=1))+
  scale_colour_viridis_d()+labs(title="B  Parallel efficiency",x="Workers",y="Parallel efficiency")+publication_theme()
f3c <- runtime_summary %>% filter(Workers>1) %>%
  ggplot(aes(factor(Workers),Runtime_reduction_percent,fill=Dataset))+
  geom_col(position="dodge")+scale_fill_viridis_d()+scale_y_continuous(labels=label_percent(scale=1))+
  labs(title="C  Runtime reduction relative to one worker",x="Workers",y="Runtime reduction")+publication_theme()
fig3 <- (f3a|f3b)/f3c+plot_annotation(title="Figure 3. Speedup, efficiency, and runtime reduction")
save_figure(fig3,"Fig3_speedup_efficiency_runtime_reduction",13,10)

# Figure 4
f4a <- runtime_summary %>% filter(Workers>1) %>%
  ggplot(aes(Workers,Incremental_runtime_reduction_percent,colour=Dataset,group=Dataset))+
  geom_line()+geom_point()+scale_colour_viridis_d()+
  labs(title="A  Incremental runtime reduction",x="Workers",y="Additional runtime reduction")+publication_theme()
f4b <- runtime_summary %>% filter(Workers>1) %>%
  ggplot(aes(Workers,Speedup_gain_per_added_worker,colour=Dataset,group=Dataset))+
  geom_line()+geom_point()+scale_colour_viridis_d()+
  labs(title="B  Speedup gained per additional worker",x="Workers",y="Incremental speedup")+publication_theme()
f4c <- fastest %>% mutate(label=paste0(Workers," workers\n",round(Median_runtime_minutes,1)," min")) %>%
  ggplot(aes(Dataset,Median_runtime_minutes,fill=Dataset))+geom_col(show.legend=FALSE)+
  geom_text(aes(label=label),vjust=-.3)+scale_fill_viridis_d()+
  scale_y_continuous(expand=expansion(mult=c(0,.3)))+
  labs(title="C  Fastest tested configuration",x=NULL,y="Median runtime (minutes)")+publication_theme()
fig4 <- (f4a|f4b)/f4c+plot_annotation(title="Figure 4. Diminishing returns and worker selection")
save_figure(fig4,"Fig4_diminishing_returns_fastest_tested",13,10)

# Figure 5
marker_dataset <- marker_summary %>% group_by(Dataset) %>%
  summarise(Marker_count=first(Marker_count),.groups="drop")
f5a <- ggplot(marker_dataset,aes(Dataset,Marker_count,fill=Dataset))+geom_col(show.legend=FALSE)+
  geom_text(aes(label=comma(Marker_count)),vjust=-.4)+scale_fill_viridis_d()+
  scale_y_continuous(labels=comma,expand=expansion(mult=c(0,.15)))+
  labs(title="A  Positive cluster-marker associations",x=NULL,y="Associations")+publication_theme()
f5b <- ggplot(marker_summary,aes(Workers,Markers_per_second,colour=Dataset,group=Dataset))+
  geom_line()+geom_point()+scale_colour_viridis_d()+
  labs(title="B  Computational throughput",x="Workers",y="Associations per second")+publication_theme()
f5c <- successful %>% mutate(Workers=factor(Workers,levels=config$workers)) %>%
  ggplot(aes(Workers,Markers_detected,colour=Dataset,group=interaction(Dataset,Replicate)))+
  geom_line(alpha=.45)+geom_point()+facet_wrap(~Dataset,scales="free_y",nrow=1)+
  scale_colour_viridis_d()+scale_y_continuous(labels=comma)+
  labs(title="C  Output reproducibility",x="Workers",y="Associations")+
  publication_theme()+theme(legend.position="none")
fig5 <- (f5a|f5b)/f5c+plot_annotation(title="Figure 5. Marker-output reproducibility and throughput")
save_figure(fig5,"Fig5_marker_reproducibility_throughput",13,10)

###############################
# 13. SUPPLEMENTARY FIGURES S1-S8
###############################
# S1 replicate trajectories
s1 <- ggplot(successful,aes(Replicate,Runtime_seconds/60,colour=factor(Workers),group=Workers))+
  geom_line()+geom_point()+facet_wrap(~Dataset,scales="free_y",nrow=1)+
  scale_colour_viridis_d()+labs(title="Supplementary Figure S1. Runtime across replicates",
  y="Runtime (minutes)",colour="Workers")+publication_theme()
save_figure(s1,"FigS1_replicate_runtime_trajectories",13,5,dirs$supplementary)
# S2 CV
s2 <- ggplot(runtime_summary,aes(Workers,CV_percent,colour=Dataset,group=Dataset))+
  geom_line()+geom_point()+scale_colour_viridis_d()+labs(title="Supplementary Figure S2. Runtime variability",y="CV (%)")+publication_theme()
save_figure(s2,"FigS2_runtime_coefficient_variation",8,6,dirs$supplementary)
# S3 runtime heatmap
s3 <- ggplot(runtime_summary,aes(factor(Workers),Dataset,fill=Median_runtime_minutes))+
  geom_tile(colour="white")+geom_text(aes(label=round(Median_runtime_minutes,1)))+
  scale_fill_viridis_c()+labs(title="Supplementary Figure S3. Median runtime heatmap",x="Workers",fill="Minutes")+publication_theme()
save_figure(s3,"FigS3_median_runtime_heatmap",8,5,dirs$supplementary)
# S4 speedup heatmap
s4 <- ggplot(runtime_summary,aes(factor(Workers),Dataset,fill=Speedup))+
  geom_tile(colour="white")+geom_text(aes(label=paste0(round(Speedup,2),"×")))+
  scale_fill_viridis_c()+labs(title="Supplementary Figure S4. Speedup heatmap",x="Workers")+publication_theme()
save_figure(s4,"FigS4_speedup_heatmap",8,5,dirs$supplementary)
# S5 parent-process memory
s5 <- ggplot(memory_summary,aes(Workers,Median_parent_peak_RAM_MB,colour=Dataset,group=Dataset))+
  geom_line()+geom_point()+scale_colour_viridis_d()+
  labs(title="Supplementary Figure S5. Parent-process peak RAM",
       subtitle="Exploratory only; does not sum child-process memory",y="Peak RAM (MB)")+publication_theme()
save_figure(s5,"FigS5_parent_process_peak_RAM",8,6,dirs$supplementary)
# S6 observed range
s6 <- ggplot(runtime_summary,aes(Workers,Mean_runtime_seconds/60,colour=Dataset,group=Dataset))+
  geom_errorbar(aes(ymin=Minimum_runtime_seconds/60,ymax=Maximum_runtime_seconds/60),width=.15)+
  geom_line()+geom_point()+scale_colour_viridis_d()+
  labs(title="Supplementary Figure S6. Mean runtime and observed range",y="Runtime (minutes)")+publication_theme()
save_figure(s6,"FigS6_runtime_observed_range",8,6,dirs$supplementary)
# S7 cluster sizes
cluster_long <- dataset_characteristics %>%
  select(Dataset,Smallest_cluster,Median_cluster_size,Largest_cluster) %>%
  pivot_longer(-Dataset,names_to="Metric",values_to="Cells")
s7 <- ggplot(cluster_long,aes(Dataset,Cells,fill=Metric))+geom_col(position="dodge")+
  scale_fill_viridis_d()+labs(title="Supplementary Figure S7. Cluster-size characteristics")+publication_theme()
save_figure(s7,"FigS7_cluster_size_characteristics",8,6,dirs$supplementary)
# S8 descriptive burden
serial <- runtime_summary %>% filter(Workers==1) %>% select(Dataset,Serial=Median_runtime_minutes)
desc <- dataset_characteristics %>% left_join(serial,by="Dataset") %>%
  select(Dataset,Cells,Features,Clusters,Serial) %>% pivot_longer(-Dataset) %>%
  group_by(name) %>% mutate(Relative=value/value[Dataset=="SMALL"]) %>% ungroup()
s8 <- ggplot(desc,aes(Dataset,Relative,fill=Dataset))+geom_col(show.legend=FALSE)+
  facet_wrap(~name,scales="free_y",nrow=1)+geom_hline(yintercept=1,linetype="dashed")+
  scale_fill_viridis_d()+labs(title="Supplementary Figure S8. Descriptive characteristics and serial burden")+publication_theme()
save_figure(s8,"FigS8_descriptive_dataset_characteristics_serial_burden",13,5.5,dirs$supplementary)

###############################
# 14. SOURCE DATA AND MANUSCRIPT TEXT
###############################
safe_write_csv(successful,file.path(dirs$source_data,"SourceData_Fig2_runtime_distributions.csv"))
safe_write_csv(runtime_summary,file.path(dirs$source_data,"SourceData_Fig2_Fig3_Fig4_runtime_summary.csv"))
safe_write_csv(dataset_characteristics,file.path(dirs$source_data,"SourceData_Fig1_FigS7_FigS8_dataset_characteristics.csv"))
safe_write_csv(marker_summary,file.path(dirs$source_data,"SourceData_Fig5_marker_reproducibility_throughput.csv"))
safe_write_csv(memory_summary,file.path(dirs$source_data,"SourceData_FigS5_parent_process_memory.csv"))

writeLines(paste(
  "Marker-detection performance was evaluated using Seurat's FindAllMarkers",
  "function with the Wilcoxon rank-sum test and positive markers only.",
  "Each dataset-worker configuration was repeated five times. Speedup was",
  "calculated relative to the median one-worker runtime, and parallel",
  "efficiency was defined as speedup divided by worker count."
), file.path(dirs$manuscript,"Methods_insert_complete_benchmark.txt"))

writeLines(paste(
  "Peak memory values collected using peakRAM primarily reflected the parent",
  "R process and did not reliably aggregate memory across multisession child",
  "processes. They are presented only as exploratory diagnostics."
), file.path(dirs$manuscript,"Memory_measurement_limitation.txt"))

###############################
# 15. FINAL VALIDATION AND INVENTORY
###############################
final_report <- tibble(
  Item = c("Successful benchmark runs","Datasets","Workers tested",
           "Replicates per configuration","Main figures",
           "Supplementary figures","Supplementary tables"),
  Result = c(nrow(successful),paste(names(datasets),collapse=", "),
             paste(config$workers,collapse=", "),config$replicates,5,8,8)
)
safe_write_csv(final_report,file.path(dirs$manuscript,"Final_analysis_validation.csv"))
files <- list.files(config$output_dir,recursive=TRUE,full.names=TRUE)
safe_write_csv(tibble(File=basename(files),Full_path=files,
  Extension=tools::file_ext(files),Size_KB=round(file.info(files)$size/1024,2)),
  file.path(dirs$manuscript,"Generated_output_inventory.csv"))
capture.output(sessionInfo(),file=file.path(dirs$manuscript,"Session_information_complete_pipeline.txt"))
write_log("scBenchR complete pipeline finished successfully")
cat("\nscBenchR COMPLETE PIPELINE FINISHED\nOutput:",config$output_dir,"\n")
