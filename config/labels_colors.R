celltype_order <- c(
  "germ",
  "degenerated",
  "granulosa",
  "stroma",
  "endothelial",
  "mural",
  "immune",
  "erythroid"
)

celltype_colors <- c(
  germ = "#1B9E77",
  degenerated = "#2EC4B6",
  granulosa = "#D95F02",
  stroma = "#7570B3",
  endothelial = "#E7298A",
  mural = "#66A61E",
  immune = "#E6AB02",
  erythroid = "#A6761D"
)

subcluster_label_map <- c(
  germ_1 = "Meiotic-entry germ cells",
  germ_2 = "Stalled meiotic germ cells",
  germ_3 = "Mitotic oogonia",
  germ_4 = "Pachytene/diplotene germ cells",
  germ_5 = "Leptotene/zygotene germ cells",
  germ_7 = "Primordial follicle oocytes",
  germ_8 = "Degenerating germ cells",
  germ_0 = "Atresia-stressed degenerating follicle cells",
  germ_6 = "Clearance-associated degenerating follicle cells",

  granulosa_0 = "Supportive pre-granulosa",
  granulosa_1 = "Signaling granulosa",
  granulosa_2 = "Primordial follicle granulosa",
  granulosa_3 = "Morphogenetic granulosa RELN+",
  granulosa_4 = "Matrix-remodeling granulosa",
  granulosa_5 = "Epithelial-like granulosa",
  granulosa_6 = "Proliferative granulosa progenitors",
  granulosa_7 = "Morphogenetic granulosa SOX5+",
  granulosa_8 = "Stress-activated granulosa",

  stroma_0 = "Cortical stroma",
  stroma_1 = "Medullary stroma",
  stroma_2 = "Signaling stroma",
  stroma_3 = "Proliferative stromal progenitors",
  stroma_4 = "Perineural stroma",

  endothelial_0 = "Angiogenic endothelial",

  mural_0 = "Pericytes",
  mural_1 = "Contractile VSMC",

  erythroid_0 = "Late erythroid",
  erythroid_1 = "Early erythroid",

  immune_0 = "Tissue macrophages",
  immune_1 = "NK T cells"
)
