library(scBenchR)
library(Seurat)

base_dir <- "D:/scRNAseq"
out_dir  <- file.path(base_dir, "output")

scbenchr_init_outdirs(out_dir)
scbenchr_write_session_info(out_dir)

# 1) Load
dat <- scbenchr_load_pdac(base_dir = base_dir)
seu_small <- dat$seu_small
seu_large <- dat$seu_large

saveRDS(seu_small, file.path(out_dir, "seu_small_raw.rds"))
saveRDS(seu_large, file.path(out_dir, "seu_large_raw.rds"))

# 2) QC
seu_small <- scbenchr_add_percent_mt(seu_small)
seu_large <- scbenchr_add_percent_mt(seu_large)

scbenchr_qc_summary(seu_small, "SMALL (before)")
scbenchr_qc_summary(seu_large, "LARGE (before)")

seu_small_f <- scbenchr_qc_filter(seu_small, min_features=200, max_features=6000, min_counts=500, max_mt=20)
seu_large_f <- scbenchr_qc_filter(seu_large, min_features=300, max_features=5000, min_counts=800, max_mt=15)

scbenchr_qc_summary(seu_small_f, "SMALL (after)")
scbenchr_qc_summary(seu_large_f, "LARGE (after)")

saveRDS(seu_small_f, file.path(out_dir, "seu_small_filtered.rds"))
saveRDS(seu_large_f, file.path(out_dir, "seu_large_filtered.rds"))

# 3) Step3 serial
res_small <- scbenchr_run_serial_step3(seu_small_f, "SMALL", npcs=30, resolution=0.6)
res_large <- scbenchr_run_serial_step3(seu_large_f, "LARGE", npcs=30, resolution=0.6)

saveRDS(res_small$seu, file.path(out_dir, "seu_small_step3.rds"))
saveRDS(res_large$seu, file.path(out_dir, "seu_large_step3.rds"))

write.csv(res_small$timing, file.path(out_dir, "timing_step3_small_serial.csv"), row.names=FALSE)
write.csv(res_large$timing, file.path(out_dir, "timing_step3_large_serial.csv"), row.names=FALSE)

tim3 <- rbind(
  transform(res_small$timing, dataset="SMALL"),
  transform(res_large$timing, dataset="LARGE")
)
write.csv(tim3, file.path(out_dir, "timing_step3_all.csv"), row.names=FALSE)

# 3b) Integration
seu_merged_integrated <- scbenchr_integrate_rpca(seu_small_f, seu_large_f, npcs=30, resolution=0.6)
saveRDS(seu_merged_integrated, file.path(out_dir, "integrated", "seu_merged_integrated_step3.rds"))
scbenchr_export_merged_metadata(seu_merged_integrated, out_dir)

# 4) Step4 markers + worker grid
scbenchr_run_step4(res_small$seu, res_large$seu, seu_merged_integrated, out_dir)
scbenchr_worker_grid_benchmark(res_small$seu, res_large$seu, out_dir)

message("DONE: outputs written to ", out_dir)
