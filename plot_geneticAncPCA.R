# Required packages
library(data.table)
library(ggplot2)

# Load data --------------------------------------------------------------------
pca_dir <- "/lustre07/scratch/chanalex/CARTaGENE_HGDP-1KG/PCA_projection/12_reference_PCA"
meta_dir <- "/lustre06/project/6061810/shared/HGDP_1KG"

ref <- fread(file.path(pca_dir, "HGDP_1KG.reference_self_projected.sscore"))

# Clean column names in case PLINK writes #IID
names(ref) <- sub("^#", "", names(ref))

# Rename projected score columns if needed
if ("PC1_AVG" %in% names(ref)) setnames(ref, "PC1_AVG", "PC1")
if ("PC2_AVG" %in% names(ref)) setnames(ref, "PC2_AVG", "PC2")

meta <- fread(file.path(meta_dir, "hgdp_tgp_meta.tsv"))

# Keep and rename useful metadata columns
meta_plot <- meta[, .(
  IID = s,
  project = hgdp_tgp_meta.Project,
  study_region = hgdp_tgp_meta.Study.region,
  population = hgdp_tgp_meta.Population,
  genetic_region = hgdp_tgp_meta.Genetic.region
)]

# -----------------------
# Merge metadata onto PCA scores
ref_meta <- merge(ref, meta_plot, by = "IID", all.x = TRUE)

# Check that merge worked (both tables should show all FALSE)
table(is.na(ref_meta$genetic_region))
table(is.na(ref_meta$population))
table(ref_meta$project, useNA = "ifany")

# Create project labels for figure (in column "project_shape") for plotting
ref_meta[, project_shape := fifelse(
  grepl("HGDP", project, ignore.case = TRUE),
  "HGDP",
  fifelse(
    grepl("TGP|1000", project, ignore.case = TRUE),
    "1000G",
    project
  )
)]

# Calculate variance explained by PCs ------------------------------------------
# Read variance explained from PLINK eigenvalues
eigenval <- scan(file.path(pca_dir, "HGDP_1KG.reference_PCA.eigenval"))

eigenvec <- fread(file.path(pca_dir, "HGDP_1KG.reference_PCA.eigenvec"))
names(eigenvec) <- sub("^#", "", names(eigenvec))

n_train <- nrow(eigenvec)

pc_var <- eigenval / n_train * 100

pc1_lab <- paste0("PC1 (", round(pc_var[1], 2), "%)")
pc2_lab <- paste0("PC2 (", round(pc_var[2], 2), "%)")

# Genetic ancestry PCA coloured by broad genetic region ------------------------
ggplot(ref_meta, aes(
  x = PC1,
  y = PC2,
  fill  = genetic_region,
  color = genetic_region,
  shape = project_shape
)) +
  geom_point(
    size = 1.8,
    alpha = 0.55,      # lower = more transparent fill
    stroke = 0.25      # border thickness
  ) +
  scale_shape_manual(
    name = "Project",
    values = c(
      "1000G" = 16,  # circle
      "HGDP"  = 17   # triangle
    ),
    labels = c(
      "1000G" = "1000 Genomes Project",
      "HGDP"  = "Human Genome Diversity Project"
    )
  ) +
  # Modify legends
  guides(
    shape = guide_legend(
      order = 1,
      override.aes = list(size = 3.5, alpha = 1)
    ),
    color = guide_legend(
      order = 2,
      override.aes = list(size = 3.5, alpha = 1)
    ),
    fill = "none"
  ) +
  theme_bw(base_size = 12) +
  coord_equal() +
  labs(
    x = pc1_lab,
    y = pc2_lab,
    shape = "Projects",
    color = "Genetic regions",
    title = "HGDP-1000G reference PCA",
    subtitle = "Reference samples coloured by genetic region"
  )

# ggplot(ref_meta, aes(x = PC1, y = PC2, color = genetic_region)) +
#   geom_point(size = 0.8, alpha = 0.75) +
#   theme_bw(base_size = 12) +
#   coord_equal() +
#   labs(
#     x = "PC1",
#     y = "PC2",
#     color = "Genetic region",
#     title = "HGDP-1000G reference PCA",
#     subtitle = "Reference samples coloured by genetic region"
#   )

# Genetic ancestry PCA coloured by population label
ggplot(ref_meta, aes(x = PC1, y = PC2, color = population)) +
  geom_point(size = 0.7, alpha = 0.75) +
  theme_bw(base_size = 12) +
  coord_equal() +
  labs(
    x = "PC1",
    y = "PC2",
    color = "Population",
    title = "HGDP-1000G reference PCA",
    subtitle = "Reference samples coloured by population"
  )

ggsave(
  filename = file.path(pca_dir, "HGDP_1KG.reference_self_projected_PC1_PC2.png"),
  width = 7,
  height = 5,
  dpi = 300
)