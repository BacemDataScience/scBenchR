
<!-- README.md is generated from README.Rmd. Do not edit README.md directly. -->

# scBenchR

scBenchR is an R package for system-aware benchmarking of Seurat v5
single-cell RNA-seq workflows on Windows systems. It provides
reproducible utilities to quantify runtime performance, repeated runtime performance, parallel scaling, computational efficiency, and marker-output reproducibility under Windows multisession execution.

The package is intended for method development, performance evaluation,
and publication-ready reporting of scRNA-seq analyses.

# Scope

- Reproducible benchmarking with deterministic defaults
- Explicit handling of Windows multisession parallelism
- Separation of serial and parallel workflow components
- CSV and figure-based outputs suitable for manuscripts
- Repeated benchmarking across configurable worker counts
- Runtime variability and statistical performance assessment
- Marker-output reproducibility across repeated executions

Linux and HPC benchmarking are intentionally out of scope.

# Workflow overview

``` mermaid
flowchart TD
  A[Seurat v5 object] --> B[Quality control]
  B --> C[Step 3: Serial pipeline]
  C --> D[Optional RPCA integration]
  C --> E[Step 4: Marker benchmarking]
  E --> F[Worker scaling analysis]
  F --> G[Tables and figures]
```

# Installation

``` r
# install.packages("pak")
pak::pak("BacemDataScience/scBenchR")
```

# Example execution

``` r
library(scBenchR)

timing_example <- data.frame(
  step        = c("serial", "parallel"),
  elapsed_sec = c(120, 55)
)

timing_example
#>       step elapsed_sec
#> 1   serial         120
#> 2 parallel          55
```

# Typical scBenchR workflow

``` r
library(Seurat)
library(scBenchR)

# Pre-existing Seurat v5 object
# seu <- ...

bench_res <- scbenchr_run_step4(
  seu      = seu,
  assay   = "RNA",
  layer   = "data",
  workers = c(1, 2, 4,6 , 8),
  tag     = "example_dataset"
)

bench_res$timing
bench_res$markers
```

# Outputs

| Output type     | Description                                  |
|-----------------|----------------------------------------------|
| Timing tables   | Per-step and per-worker runtime summaries    |
| Marker tables   | Cluster-wise differential expression results |
| Scaling metrics | Runtime, speedup, parallel efficiency, runtime reduction, and variability summaries         |
| Figures         | Publication-quality plots (PDF/PNG)          |

# Windows-specific considerations

- Core pipelines are executed serially for stability
- Parallelization is evaluated using future::multisession
- Scaling analysis quantifies speedup, parallel efficiency, runtime reduction, and diminishing returns

# Citation

If this package is used in academic work, please cite:

Saada, B. scBenchR: System-aware benchmarking of Seurat v5 scRNA-seq
workflows on Windows systems. Manuscript in preparation.

# License

MIT License © Bacem Saada
