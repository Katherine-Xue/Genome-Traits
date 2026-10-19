### ratio ----
library(tidyr)
library(dplyr)
library(ggpubr)


# =========================
# 1️⃣ 数据准备
# =========================
df1 <- df0[,c(1:12,46:63,73:87)]

df1$trophic <- df1$Biosynthesis_growth / df1$Cell_maintenance
df1$trade <- df1$Biosynthesis_growth / df1$Stress_tolerance
df1$usage <- df1$Resource_acquisition / df1$Biosynthesis_growth
df1$maintenance <- df1$Resource_acquisition / df1$Cell_maintenance

traits <- c("GC", "geneDensity", "Genome_Size")
functions <- c("trophic","trade","usage","maintenance")

# =========================
# 2️⃣ long 数据
# =========================
df_long <- df1 %>%
  pivot_longer(
    cols = all_of(traits),
    names_to = "Traits",
    values_to = "Values"
  ) %>%
  pivot_longer(
    cols = all_of(functions),
    names_to = "Functions",
    values_to = "Abundance"
  )

# =========================
# 3️⃣ 统计
# =========================
df_stat <- df_long %>%
  group_by(Functions, Traits) %>%
  summarise(
    r = cor(Values, Abundance, use = "complete.obs"),
    r2 = r^2,
    p.value = cor.test(Values, Abundance)$p.value,
    .groups = "drop"
  )

# =========================
# 4️⃣ 合并 + 控制变量
# =========================
df_plot <- df_long %>%
  left_join(df_stat, by = c("Functions", "Traits")) %>%
  mutate(
    sig_p = ifelse(p.value < 0.05, "sig", "ns"),
    sig_r2 = ifelse(r2 >= 0.1, "high", "low"),
    
    color_plot = ifelse(sig_p == "sig" & sig_r2 == "high",
                        as.character(Traits), "ns"),
    
    linetype_plot = ifelse(sig_p == "sig", "dashed", "solid")
  )

# =========================
# ❗关键修复：最后一步统一 factor
# =========================
df_plot <- df_plot %>%
  mutate(
    Traits = factor(Traits,
                    levels = c("GC","Genome_Size","geneDensity")),
    
    Functions = factor(Functions,
                       levels = c("trophic","trade","usage","maintenance"))
  )

# =========================
# 5️⃣ label 数据（必须在 df_plot 后）
# =========================
label_df <- df_plot %>%
  distinct(Functions, Traits, r, r2, p.value) %>%
  mutate(
    p_label = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.01  ~ "p < 0.01",
      p.value < 0.05  ~ "p < 0.05",
      TRUE ~ paste0("p = ", signif(p.value, 2))   # ✅ 非显著显示原值
    ),
    
    label = paste0(
      "r = ", round(r, 2),
      " ; R² = ", round(r2, 2),
      " ; ", p_label
    )
  )

# =========================
# 6️⃣ 颜色
# =========================
trait_levels <- levels(df_plot$Traits)
default_cols <- scales::hue_pal()(length(trait_levels))
color_map <- setNames(default_cols, trait_levels)
color_map <- c(color_map, ns = "grey40")

# =========================
# 7️⃣ 作图
# =========================
ggplot(df_plot, aes(x = Values, y = Abundance)) +
  
  geom_point(aes(color = color_plot), alpha = 0.6, size = 1) +
  
  geom_smooth(
    aes(color = color_plot, linetype = linetype_plot),
    method = "lm",
    se = FALSE,
    linewidth = 1
  ) +
  
  # ✅ label 修复（关键：带 facet 变量）
  geom_text(
    data = label_df,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.1,
    vjust = 1.3,
    size = 3,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(values = color_map) +
  
  facet_grid(
    Functions ~ Traits,
    scales = "free",
    labeller = labeller(
      Traits = c(
        "Genome_Size" = "Genome Size",
        "geneDensity" = "Gene Density",
        "GC" = "GC"
      ),
      Functions = c(
        "trophic" = "Growth vs Maintenance",
        "trade" = "Growth vs Stress",
        "usage" = "Acquisition vs Growth",
        "maintenance" = "Acquisition vs Maintenance"
      )
    )
  ) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0, "lines"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )


###
colnames(df_plot)
df8 <- unique(df_plot[,c(42,44,46:50)])
df9 <- unique(df_plot[,c(43,45,47:49)])

### PLS-PM ----
library(plspm)

colnames(df1)

env_vars <- c("TK","TP","pH","EC","CEC",
              "clay","sand","TN","TC","CN",
              "SOC","MAP","TRA","MAT","PVA","ADM")

trait_vars <- c("GC","geneDensity","Genome_Size")

func_vars <- c("Biosynthesis_growth",
               "Cell_maintenance",
               "Energy_metabolism",
               "Resource_acquisition",
               "Stress_tolerance")

pls_data <- df1[, c(env_vars, trait_vars, func_vars)]

path_matrix <- rbind(
  Environment = c(0, 0, 0),
  Traits      = c(1, 0, 0),
  Functions   = c(1, 1, 0)
)

colnames(path_matrix) <- rownames(path_matrix)

blocks <- list(
  env_vars,
  trait_vars,
  func_vars
)

modes <- c("A","A","A")

pls_model <- plspm(
  pls_data,
  path_matrix,
  blocks,
  modes = modes,
  scaled = TRUE
)

summary(pls_model)

pls_model$path_coefs
pls_model$inner_model
plot(pls_model)



#### SEM ----
setwd("D:/OneDrive - The University of Sydney (Staff)/0_Projects/0_MetaDiverse/SEM")

library(bnlearn)
library(Rgraphviz)
library(dplyr)

env_vars <- c(
  "TK","TP","TN","pH","EC","CEC",
  "clay","SOC",
  "MAP","TRA","MAT","PVA","ADM")

trait_vars <- c(
  "GC",
  "geneDensity",
  "Genome_Size")

func_vars <- c(
  "Biosynthesis_growth",
  "Cell_maintenance",
  "Energy_metabolism",
  "Resource_acquisition",
  "Stress_tolerance")

vars <- c(
  env_vars,
  trait_vars,
  func_vars)

df_bn <- df1[, vars]

df_bn <- na.omit(df_bn)

df_bn <- as.data.frame(scale(df_bn))

# -----------------------------------------
# traits/functions CANNOT affect environment
# -----------------------------------------

blacklist1 <- expand.grid(
  from = c(trait_vars, func_vars),
  to = env_vars
)

# -----------------------------------------
# functions CANNOT affect traits
# -----------------------------------------

blacklist2 <- expand.grid(
  from = func_vars,
  to = trait_vars
)

# -----------------------------------------
# functions CANNOT affect functions
# -----------------------------------------

blacklist3 <- expand.grid(
  from = func_vars,
  to = func_vars
)

blacklist3 <- blacklist3[
  blacklist3$from != blacklist3$to,]

# -----------------------------------------
# combine
# -----------------------------------------

blacklist <- rbind(
  blacklist1,
  blacklist2,
  blacklist3)

bn <- hc(
  df_bn,
  blacklist = blacklist,
  score = "bic-g"
)

boot <- boot.strength(
  df_bn,
  R = 1000,
  algorithm = "hc",
  algorithm.args = list(
    blacklist = blacklist
  )
)

avg_bn <- averaged.network(
  boot,
  threshold = 0.8
)

strength <- arc.strength(
  avg_bn,
  data = df_bn,
  criterion = "bic-g"
)

avg_bn
strength

#plot(avg_bn)

edges <- strength %>%
  
  filter(
    strength < -5
  )  #%>%
  # 
  # filter(
  #   
  #   # Environment → Traits
  #   (
  #     from %in% env_vars &
  #       to %in% trait_vars
  #   ) |
  #     
  #     # Traits → Functions
  #     (
  #       from %in% trait_vars &
  #         to %in% func_vars
  #     ) |
  #     
  #     # Environment → Functions
  #     (
  #       from %in% env_vars &
  #         to %in% func_vars
  #     )
  # )


edges


### plot

library(tidygraph)
library(ggraph)
library(igraph)
library(dplyr)
library(ggplot2)
library(grid)

# =========================================================
# 1️⃣ edge table
# =========================================================

edges <- data.frame(
  
  from = c(
    "TK","TK","TP","TN","TN","TN","TN",
    "pH",
    "CEC","CEC","CEC","CEC","CEC","CEC",
    "clay","clay","clay","clay",
    "SOC","SOC",
    "MAP",
    "TRA","TRA","TRA","TRA","TRA","TRA","TRA","TRA",
    "MAT","MAT","MAT","MAT",
    "PVA","PVA","PVA","PVA","PVA",
    "GC","GC","GC","GC",
    "Genome_Size","Genome_Size","Genome_Size","Genome_Size","Genome_Size"
  ),
  
  to = c(
    "EC",
    "Cell_maintenance",
    "Genome_Size",
    "TP",
    "Cell_maintenance",
    "Energy_metabolism",
    "Stress_tolerance",
    "Genome_Size",
    "TK",
    "pH",
    "PVA",
    "Genome_Size",
    "Energy_metabolism",
    "Stress_tolerance",
    "TK",
    "CEC",
    "TRA",
    "PVA",
    "TN",
    "Cell_maintenance",
    "ADM",
    "PVA",
    "ADM",
    "geneDensity",
    "Biosynthesis_growth",
    "Cell_maintenance",
    "Energy_metabolism",
    "Resource_acquisition",
    "Stress_tolerance",
    "pH",
    "PVA",
    "ADM",
    "GC",
    "ADM",
    "Biosynthesis_growth",
    "Cell_maintenance",
    "Energy_metabolism",
    "Stress_tolerance",
    "geneDensity",
    "Biosynthesis_growth",
    "Cell_maintenance",
    "Energy_metabolism",
    "GC",
    "Biosynthesis_growth",
    "Cell_maintenance",
    "Energy_metabolism",
    "Resource_acquisition"
  ),
  
  strength = c(
    -12.153597,
    -12.335005,
    -15.701174,
    -101.026698,
    -12.189352,
    -17.544787,
    -8.674870,
    -10.893624,
    -15.076110,
    -38.141609,
    -17.448520,
    -17.955637,
    -27.240819,
    -54.483284,
    -22.132166,
    -122.208580,
    -41.542909,
    -30.727658,
    -230.442581,
    -6.571958,
    -104.782187,
    -31.621975,
    -26.668452,
    -5.739660,
    -31.824560,
    -12.575716,
    -31.301993,
    -54.604213,
    -14.493633,
    -14.198671,
    -195.093344,
    -140.227052,
    -84.716629,
    -9.514626,
    -11.155521,
    -25.862740,
    -10.424573,
    -54.192774,
    -44.998110,
    -27.110432,
    -44.426572,
    -29.029527,
    -19.424335,
    -55.438861,
    -62.820450,
    -45.443544,
    -33.535861
  )
)



# =========================================================
# 2️⃣ groups
# =========================================================

env_vars <- c(
  "TP","TK","TN","EC","SOC","CEC",
  "pH","clay","MAT","TRA","MAP",
  "PVA","ADM"
)

trait_vars <- c(
  "GC",
  "Genome_Size",
  "geneDensity"
)

func_vars <- c(
  "Biosynthesis_growth",
  "Cell_maintenance",
  "Energy_metabolism",
  "Resource_acquisition",
  "Stress_tolerance"
)

# =========================================================
# 3️⃣ edge metadata
# =========================================================

edges <- edges %>%
  
  mutate(
    
    abs_strength = abs(strength),
    
    # =====================================================
    # edge color groups
    # =====================================================
    
    edge_group = case_when(
      
      # ---------------------------------------------------
      # Environment ↔ Environment
      # ---------------------------------------------------
      
      from %in% env_vars &
        to %in% env_vars
      ~ "Environment",
      
      # ---------------------------------------------------
      # Trait ↔ Trait
      # ---------------------------------------------------
      
      from %in% trait_vars &
        to %in% trait_vars
      ~ "Traits",
      
      # ---------------------------------------------------
      # Environment → Trait
      # ---------------------------------------------------
      
      from %in% env_vars &
        to %in% trait_vars
      ~ "Env_to_Trait",
      
      # ---------------------------------------------------
      # Environment → Function
      # ---------------------------------------------------
      
      from %in% env_vars &
        to %in% func_vars
      ~ "Env_to_Function",
      
      # ---------------------------------------------------
      # Trait → Function
      # ---------------------------------------------------
      
      from %in% trait_vars &
        to %in% func_vars
      ~ "Trait_to_Function",
      
      TRUE
      ~ "Other"
    ),
    
    # =====================================================
    # line type
    # =====================================================
    
    edge_type = case_when(
      
      # only within-group are dotted
      
      from %in% env_vars &
        to %in% env_vars
      ~ "within",
      
      from %in% trait_vars &
        to %in% trait_vars
      ~ "within",
      
      # Environment → Trait should be SOLID
      
      TRUE
      ~ "between"
    )
  )

# =========================================================
# 5️⃣ fixed node positions
# =========================================================

nodes <- nodes %>%
  
  mutate(
    
    x = case_when(
      
      group == "Environment" ~ 0,
      
      group == "Traits" ~ 6,
      
      TRUE ~ 12
    )
  )

nodes$y <- NA

# ---------------------------------------------------------
# Environment order
# ---------------------------------------------------------

env_order <- c(
  "TP","TK","TN","EC","SOC","CEC",
  "pH","clay","MAT","TRA","MAP",
  "PVA","ADM"
)

env_y <- seq(
  12,
  -12,
  length.out = length(env_order)
)

for(i in seq_along(env_order)){
  
  nodes$y[
    nodes$name == env_order[i]
  ] <- env_y[i]
}

# ---------------------------------------------------------
# Trait order
# ---------------------------------------------------------

trait_order <- c(
  "GC",
  "Genome_Size",
  "geneDensity"
)

trait_y <- c(4, 0, -4)

for(i in seq_along(trait_order)){
  
  nodes$y[
    nodes$name == trait_order[i]
  ] <- trait_y[i]
}

# ---------------------------------------------------------
# Function order
# ---------------------------------------------------------

func_order <- c(
  "Biosynthesis_growth",
  "Cell_maintenance",
  "Energy_metabolism",
  "Resource_acquisition",
  "Stress_tolerance"
)

func_y <- seq(
  8,
  -8,
  length.out = length(func_order)
)

for(i in seq_along(func_order)){
  
  nodes$y[
    nodes$name == func_order[i]
  ] <- func_y[i]
}

# =========================================================
# 6️⃣ graph object
# =========================================================

graph <- tbl_graph(
  nodes = nodes,
  edges = edges,
  directed = TRUE
)

# =========================================================
# 7️⃣ colors
# =========================================================

node_cols <- c(
  "Environment" = "#00BFC4",
  "Traits" = "#C77CFF",
  "Functions" = "#7CAE00"
)

edge_cols <- c(
  
  # ------------------------------------------------------
  # within Environment
  # ------------------------------------------------------
  
  "Environment" = "#007439",
  
  # ------------------------------------------------------
  # within Traits
  # ------------------------------------------------------
  
  "Traits" = "#A65B00",
  
  # ------------------------------------------------------
  # Environment → Trait
  # ------------------------------------------------------
  
  "Env_to_Trait" = "#00BFC4",
  
  # ------------------------------------------------------
  # Environment → Function
  # ------------------------------------------------------
  
  "Env_to_Function" = "gray70",
  
  # ------------------------------------------------------
  # Trait → Function
  # ------------------------------------------------------
  
  "Trait_to_Function" = "#C77CFF",
  
  "Other" = "gray70"
)

# =========================================================
# 8️⃣ plot
# =========================================================

ggraph(
  graph,
  layout = "manual",
  x = nodes$x,
  y = nodes$y
) +
  
  # ---------------------------------------------------------
# edges
# ---------------------------------------------------------

geom_edge_arc(
  
  aes(
    width = abs_strength,
    color = edge_group,
    linetype = edge_type
  ),
  
  alpha = 0.9,
  
  strength = 0.14,
  
  arrow = arrow(
    length = unit(3.5, "mm"),
    type = "closed"
  ),
  
  end_cap = circle(3.5, "mm")
) +
  
  # ---------------------------------------------------------
# nodes
# ---------------------------------------------------------

geom_node_point(
  
  aes(fill = group),
  
  shape = 21,
  
  size = 8,
  
  color = "gray50",
  
  stroke = 0
) +
  
  # ---------------------------------------------------------
# labels
# ---------------------------------------------------------

geom_node_text(
  
  aes(label = name),
  
  family = "sans",
  
  size = 3
) +
  
  # ---------------------------------------------------------
# scales
# ---------------------------------------------------------

scale_fill_manual(
  values = node_cols
) +
  
  scale_edge_color_manual(
    values = edge_cols
  ) +
  
  scale_edge_width(
    range = c(0.6, 2)
  ) +
  
  # dotted internal edges
  
  scale_edge_linetype_manual(
    
    values = c(
      
      "within" = "solid",
      
      "between" = "solid"
    )
  ) +
  
  # ---------------------------------------------------------
# coordinate
# ---------------------------------------------------------

coord_cartesian(
  xlim = c(-2, 15),
  ylim = c(-14, 14),
  clip = "off"
) +
  
  # ---------------------------------------------------------
# theme
# ---------------------------------------------------------

theme_void(
  base_family = "sans"
) +
  
  theme(
    
    legend.position = "right",
    
    plot.margin = margin(
      20, 80, 20, 40
    ),
    
    plot.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    
    panel.background = element_rect(
      fill = "transparent",
      colour = NA
    )
  )


