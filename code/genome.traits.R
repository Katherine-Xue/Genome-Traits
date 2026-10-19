library(readxl)

df0 <- readRDS("sam.v4.4.RDS")
df0$Transect <- factor(df0$Transect, levels = c("NS", "EW","West"))

file <- "all.sample.kingdom.GC.xls"

df.GC <- read.delim(
  file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(df.GC)[-1] <- gsub("^d__", "GC_", colnames(df.GC)[-1])
colnames(df.GC)

df0 <- merge(df0,df.GC,by.x="SequenceID",by.y="sample", all=TRUE)

#write.csv(df0, file = "sam.KingdomeGC.csv", row.names = FALSE)

colnames(df0)
df1 <- df0[,c(1:12,46:63,73:87)]

### box plot genome traits ----
colnames(df1)
vars_to_plot <- colnames(df1)[c(31,35,34)]
vars_to_plot
# 转换为长格式
df_long <- df1 %>%
  pivot_longer(
    cols = all_of(vars_to_plot),
    names_to = "Variable",
    values_to = "Value"
    )

df_long$Variable <- factor(df_long$Variable, levels = c("GC","Genome_Size","geneDensity"))
df_long$LandUse <- factor(df_long$Landuse,levels = c("Na","Ag"))

###plot
library(ggpubr)
library(rstatix)
library(multcompView)
library(ggplot2)
### sig
letters_df <- df_long %>%
  group_by(Variable) %>%
  dunn_test(Value ~ cluster, p.adjust.method = "BH") %>%
  mutate(comparison = paste(group1, group2, sep = "-")) %>%
  group_by(Variable) %>%
  summarise(
    Letters = list(multcompLetters(setNames(p.adj, comparison))$Letters),
    .groups = "drop"
  ) %>%
  unnest_longer(Letters) %>%
  mutate(cluster = names(Letters),
         letter = Letters) %>%
  select(Variable, cluster, letter)

### position
label_df <- df_long %>%
  group_by(Variable, cluster) %>%
  summarise(y = max(Value, na.rm = TRUE), .groups = "drop") %>%
  left_join(letters_df, by = c("Variable","cluster")) %>%
  mutate(y = y * 1.08)

### cluster ----
pd <- position_dodge(width = 0.6)

ggplot(df_long, #%>% filter(!(Variable == "SOC" & Value > 10)), 
       aes(x = cluster, y = Value, colour = cluster, fill = cluster,group = cluster)) +
  geom_boxplot(position = pd, width = 0.6, outlier.shape = NA, alpha = 0.3) +
  geom_point(aes(shape = Landuse),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
             size = 2, alpha = 0.8) +
  facet_wrap(~Variable, scales = "free", ncol = 5) +
  geom_text(data = label_df,
            aes(x = cluster, y = y, label = letter),
            colour = "black",
            size = 4)+
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 11),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    y = "Value",
    colour = "Cluster"
  )

#### smooth
colnames(df1)

 ggplot(data = df1 , aes(x = MAT, y=geneDensity, colour = `Landuse`)) + #%>% filter(geneDensity<1700)
  geom_point( position = "dodge",alpha=0.6)+
   geom_smooth(method = "loess", se = FALSE,size=0.5, linetype = "dashed") + #aes(linetype = Landuse),
  theme_bw()+
  scale_y_continuous()+ #labels = scientific
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  labs(y = "Gene Density (genes/Mbp)",x="MAT (0C)")

### latitude ----
 library(ggpmisc)
 colnames(df0)
 
 ggplot(df0, aes(x =Y_Latitude,  y=`GC_Bacteria`,group = Landuse)) +
  # facet_grid(~Transect,scales="free")+
   geom_point(aes(colour = MAT,shape=Landuse),  position = "dodge",alpha=0.8,size=1.8)+
     geom_smooth(aes(linetype = Landuse),method = "lm", se = FALSE, size = 0.5) +
   stat_poly_eq(
     aes(label = paste(..eq.label.., "*\",\"*",
                       ..rr.label.., "*\",\"*",
                       ..p.value.label..)),
     formula = y ~ x,
     parse = TRUE,lineheight = 2,
     size = 3
   ) +
   scale_colour_viridis_c(
     option = "plasma",      # or "viridis", "magma", "inferno", "cividis"
     name = "MAT (°C)"
   ) +
   theme_bw(base_size = 15)+
   theme(axis.text.x = element_text(angle = 0, hjust = 1))+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.line = element_line(colour = "black"))+
   labs(y = "GC (%)",x="Latitude")

 #### soil type ----
 colnames(df1)
 
 vars <- c("Genome_Size","GC","geneDensity")
 
 df_summary <- df1 %>%
   group_by(SoilType) %>%
   summarise(
     mean_val = mean(geneDensity, na.rm = TRUE),
     median_val = median(geneDensity, na.rm = TRUE),
     se_val = sd(geneDensity, na.rm = TRUE) / sqrt(n()),
     .groups = "drop"
   ) %>%
   arrange(mean_val)
 
 df2 <- df1 %>%
   mutate(SoilType = factor(SoilType, levels = df_summary$SoilType))

 # plot boxplot
  ggplot(df2, aes(x = SoilType, y = geneDensity, colour = SoilType)) +
   geom_jitter(width = 0.2, size = 0.3, alpha = 0.6) +   # small points
   #  geom_boxplot(alpha = 0)+
   stat_summary(fun = mean, geom = "point", size = 2) +  # mean points
   stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, size = 0.8) + # error bars
   coord_flip() +
   #scale_y_continuous(labels = scales::label_scientific())+
   labs(x = "Soil Type", y = "Gene Density (genes/Mb)") +
   theme_bw(base_size = 15) +
   theme(legend.position = "none")
 
  ###
 colnames(df0)
 
 df_summary <- df0 %>%
   group_by(SoilType) %>%
   summarise(
     mean_val = mean(GC, na.rm = TRUE),
     median_val = median(GC, na.rm = TRUE),
     se_val = sd(GC, na.rm = TRUE) / sqrt(n()),
     .groups = "drop"
   ) %>%
   arrange(mean_val)
 
 df2 <- df1 %>%
   mutate(SoilType = factor(SoilType, levels = df_summary$SoilType))
 
 ggplot(df2, aes(x = SoilType, y = GC, colour = SoilType)) +
   geom_jitter(width = 0.2, size = 0.3, alpha = 0.6) +   # small points
   #  geom_boxplot(alpha = 0)+
   stat_summary(fun = mean, geom = "point", size = 2) +  # mean points
   stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, size = 0.8) + # error bars
   coord_flip() +
   #scale_y_continuous(labels = scales::label_scientific())+
   labs(x = "Soil Type", y = "GC content (%)") +
   theme_bw(base_size = 15) +
   theme(legend.position = "none")
 #### function percentage ----
 
 library(dplyr)
 library(tidyr)
 library(ggplot2)
 
 df1_pct <- df1 %>%
   mutate(
     total_func = Biosynthesis_growth + Cell_maintenance + Energy_metabolism +
       Resource_acquisition + Stress_tolerance,
     
     Biosynthesis_growth_pct  = Biosynthesis_growth  / total_func * 100,
     Cell_maintenance_pct     = Cell_maintenance     / total_func * 100,
     Energy_metabolism_pct    = Energy_metabolism    / total_func * 100,
     Resource_acquisition_pct = Resource_acquisition / total_func * 100,
     Stress_tolerance_pct     = Stress_tolerance     / total_func * 100
   )
 
 df_long <- df1_pct %>%
   select(Transect, Buffer,
          Biosynthesis_growth_pct,
          Cell_maintenance_pct,
          Energy_metabolism_pct,
          Resource_acquisition_pct,
          Stress_tolerance_pct) %>%
   pivot_longer(
     cols = ends_with("_pct"),
     names_to = "Function",
     values_to = "Percentage"
   )
 
 df_long$Transect <- factor(df_long$Transect, levels = c("NS", "EW", "West"))
 
 ggplot(df_long, aes(x = Buffer, y = Percentage + 0.01, color = Function)) +
     geom_smooth(method = "loess", se = TRUE, linewidth = 1, alpha = 0) +
     geom_point(alpha = 1, size = 1) +
    facet_wrap(~Transect, scales = "free_x") +
    scale_color_brewer(palette = "Set2") +
   scale_y_log10(breaks = c(0,1,5,6,10,13,20,40, 50))+
      theme_bw(base_size = 13) +
       theme(
     strip.background = element_rect(fill = "grey90", color = NA),
     strip.text = element_text(face = "bold"),
     legend.position = "right",
     panel.grid.minor = element_blank()
   ) +
      labs(
     x = "Buffer",
     y = "Functional composition (%)",
     color = "Functions"
   )
 
 ####SHAP ----
library(SHAPforxgboost)
library(xgboost)
library(ggplot2)
library(corrplot)
library(caret)
 
 colnames(df1)
 df1 <- df0[,c(1:12,46:63,73:87)]
 
 df6 <- df1[,c(31,35,34,13:30,33)]
 colnames(df6)
 #df6[df6 == 0] <- NA
 colnames(df6)[1:3] <- c("GC","GenomeSize","GeneDensity")
 colnames(df6) <- make.names(colnames(df6))
 
 var <-  df6[,c(1:3,22)]
 env <- df6[,c(5:13,15,16,18:21)]
 
 ### filter var
 colnames(env)
 cor_mat <- cor(env, method = "pearson", use = "pairwise.complete.obs")
 
 corrplot(
   cor_mat,
   method = "color",
   type = "upper",
   tl.cex = 0.8,
   col = colorRampPalette(c("blue", "white", "red"))(200)
 )
 
 cor_cut <- 0.85 
 
 remove_idx <- findCorrelation(
   cor_mat,
   cutoff = cor_cut,
   exact = TRUE)
 
 colnames(env)[remove_idx]
 
 env <- env[, -remove_idx]
 colnames(env)
 
 ###
 y_var <-  colnames(var)
 y_var
 
 ##### Shap normalize ----
 library(scales)
 env_scaled <- as.data.frame(lapply(env, function(x) rescale(x, to = c(0, 1))) )
 #var <- as.data.frame(lapply(var, function(x) rescale(x, to = c(0, 1))) )
 
 ##### Shap plot ----
 # 删除目标变量列并转换为矩阵
 dataX <- as.matrix(env_scaled)
 
 # 2. 参数
 params <- list(
   objective = "reg:squarederror",
   eta = 0.01,
   max_depth = 10,
   gamma = 0.01,
   subsample = 0.98,
   colsample_bytree = 0.86
 )
 
 shap_long.all <- data.frame()
 shap_long <- list()
 p <- list()
 g <- list()
 
 for (i in 1:ncol(var)) {
   
   y_var <- colnames(var)[i]
   
   dtrain <- xgb.DMatrix(
     data = dataX,
     label = var[[i]] )
   
   mod <- xgb.train(
     params = params,
     data = dtrain,
     nrounds = 300,
     verbose = 0
   )
   
   shap_values <- shap.values(xgb_model = mod, X_train = dataX)
   
   shap_long1 <- shap.prep(xgb_model = mod, X_train = dataX)
   shap_long1$group <- y_var
   
   shap_long[[y_var]] <- shap_long1
   shap_long.all <- rbind(shap_long.all, shap_long1)
   
   x_range <- max(abs(shap_long1$value)) * 1.1
   
   p[[y_var]] <- shap.plot.summary(
     shap_long1,
     x_bound = x_range,
     dilute = 5
   ) +
     ggtitle(y_var) +
     theme(plot.title = element_text(hjust = 0.5))
   
   plot_data <- shap.prep.stack.data(
     shap_contrib = shap_values$shap_score,
     top_n = 5,
     n_groups = 1
   )
   
   g[[y_var]] <- shap.plot.force_plot(
     plot_data,
     zoom_in = FALSE
   ) +
     ggtitle(y_var) +
     theme(plot.title = element_text(hjust = 0.5))
 }
 
 
 ### plot
 library(gridExtra)
 
  grid.arrange(grobs = c(p[[4]],g[[4]]))

###
 dfs1 <- unique(shap_long.all[,c(2,6,7)])
 dfs1$Var <- "Soil"
dfs1$Var[dfs1$variable %in% c("ADM","MAT","PVA","TRA")] <- "Climate"
 
#### plot bar and donut ----
library(dplyr)
library(ggplot2)
library(ggrepel)
library(grid)

# ---------------------------
# 1. 数据处理
# ---------------------------
df2 <- dfs1 %>%
  group_by(group) %>%
  arrange(desc(mean_value), .by_group = TRUE) %>%
  mutate(
    variable = factor(variable, levels = unique(variable)),
    perc = mean_value / sum(mean_value) * 100
  ) %>%
  ungroup()

# ---------------------------
# 3. 画图函数
# ---------------------------
plot_fun <- function(d){
  
  d <- d %>%
    arrange(mean_value) %>%
    mutate(variable = factor(variable, levels = variable))
  ### colour
  vars <- d$variable
  pal <- colorRampPalette(rev(c("#08306B", "#4292C6", "#FEE0D2", "#99000D")))(length(vars))
  names(pal) <- vars
  
  # -------- BAR --------
  p_bar <- ggplot(d, aes(x = variable, y = mean_value, fill = variable)) +
    geom_col(width = 0.8) +
    coord_flip(clip = "off") +   # ⭐关键：不裁切
    scale_fill_manual(values = pal) +
    labs(
      title = unique(d$group),
      x = NULL,
      y = "Mean |SHAP| value"
    ) +
    theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
      panel.grid = element_blank(),
      axis.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      plot.margin = ggplot2::margin(10, 80, 10, 10)   # ⭐给右侧留空间
    )
  
  # -------- DONUT 数据 --------
  d_donut <- d %>%
    arrange(desc(perc)) %>%
    mutate(
      ymax = cumsum(perc),
      ymin = lag(ymax, default = 0),
      label_pos = (ymax + ymin) / 2
    )
  
  # -------- DONUT --------
  donut <- ggplot(d_donut, aes(ymax = ymax, ymin = ymin, xmax = 4.5, xmin = 3.2, fill = variable)) +
    geom_rect(color = "white") +
    coord_polar(theta = "y", start = -pi/2) +
    xlim(c(2, 6)) +   # ⭐扩展范围防止裁切
    scale_fill_manual(values = pal) +
    theme_void()+theme(legend.position = "none")
  
  # -------- 百分比标签（外侧 + 引线）--------
  donut <- donut +
    geom_segment(
      data = d_donut,
      aes(
        x = 4.5,
        xend = 5,
        y = label_pos,
        yend = label_pos
      ),
      inherit.aes = FALSE
    ) +
    geom_text(
      data = d_donut,
      aes(
        x = 5.5,
        y = label_pos,
        label = paste0(round(perc, 1), "%")
      ),
      size = 3,
      inherit.aes = FALSE
    )
  
  # -------- 合并（嵌入右下角）--------
  p_bar +
    annotation_custom(
      grob = ggplotGrob(donut),
      xmin = 0,
      xmax = length(d$variable) * 0.4,
      ymin = max(d$mean_value) * 0.3,
      ymax = max(d$mean_value) * 1.1
    )
}

# ---------------------------
# 4. 按 group 出图
# ---------------------------
plot_list <- df2 %>%
  group_split(group) %>%
  lapply(plot_fun)

plot_list <- lapply(plot_list, function(p){
  p + theme(
    plot.margin = ggplot2::margin(2, 20, 2, 20)  # 默认一般是 5.5，调小就更紧凑
  )
})

# 查看某一张
library(gridExtra)

grid.arrange(grobs = c(p[4],plot_list[c(1)]), ncol = 2)

#### response ----

colnames(df0)
df1 <- df0[,c(1:12,46:63,73:87)]
colnames(df1)

df6 <- df1[,c(1:12,31,35,34,28,30,20,27)]
colnames(df6)
###
df_long.1 <- df6 %>%
  pivot_longer(
    cols = all_of(c(13,14,15)),
    names_to = "Traits",
    values_to = "Value"
  )
###
colnames(df_long.1)

df_long.2 <- df_long.1 %>%
  pivot_longer(
    cols = all_of(c(13,14,15,16)),
    names_to = "Variables",
    values_to = "Var_response"
  )
###
colnames(df_long.2)

ggplot(df_long.2, aes(x = Var_response, y = Value, colour = Traits,shape=Landuse)) +
   geom_point(position = "dodge", alpha = 1,size=1) +
   geom_smooth(aes(linetype=Landuse),method = "loess", se = FALSE, size = 0.5) +
   facet_grid(Variables~Traits, scales = "free") +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
   labs(x = "Environmental Variable", y = "Genome Trait Value")

colnames(df6)

ggplot(df6, aes(x = ADM, y = GC, colour = Landuse,shape=Landuse)) +
  geom_point(position = "dodge", alpha = 1,size=1) +
  geom_smooth(aes(linetype=Landuse),method = "loess", se = FALSE, size = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


### partial dependence plot ----
library(dplyr)
library(randomForest)
library(pdp)
library(ggplot2)
library(tidyr)

colnames(df6)

var <-  df6[,c(1:3)]
y_vars <- colnames(var)
y_vars
env_vars <- colnames(env)
env_vars

###### calculate pdp ----
compute_pdp_all <- function(df, y_vars, env_vars) {
  
  df <- df %>% drop_na()
  
  pd_list <- list()
  
  for (y in y_vars) {
    
    df_model <- df %>% select(all_of(c(y, env_vars)))
    
    # RF model
    mod <- randomForest(
      as.formula(paste(y, "~ .")),
      data = df_model,
      ntree = 500
    )
    
    for (x in env_vars) {
      
      pd <- pdp::partial(
        object = mod,
        pred.var = x,
        train = df_model,
        grid.resolution = 100
      )
      
      colnames(pd) <- c("x", "yhat")
      
      pd$env_var <- x
      pd$y_var <- y
      
      pd_list[[paste(y, x, sep = "_")]] <- pd
    }
  }
  
  bind_rows(pd_list)
}

pdp_df <- compute_pdp_all(df6, y_vars, env_vars)

###### pdp plot ---
y_vars

target_list <- list(
  GC = c("ADM", "MAT"),
  Genome.Size = c("clay", "MAT"),
  Gene.Density = c("ADM", "TRA")
)

col_map <- c(
  "GC" = "#0080FF",
  "Genome.Size" = "#FF9700",
  "Gene.Density" = "#00A087"
)

p_list <- list()

for (y in names(target_list)) {
  
  for (x in target_list[[y]]) {
    
    df_plot <- pdp_df %>%
      filter(y_var == !!y, env_var == !!x)
    
    p <- ggplot(df_plot, aes(x = x, y = yhat)) +
      geom_line(colour = "grey60", alpha = 0.5) +
      geom_smooth(method = "loess", span = 0.4,
        se = TRUE,alpha = 0.5,level = 0.99, 
        colour = col_map[[y]],
        size = 1
      ) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      labs(
        title = paste(y),
        x = x,
        y = "yhat"
      )
    
    p_list[[paste(y, x, sep = "_")]] <- p
  }
}

###
grid.arrange(grobs = p_list[c(1,3,5,2,4,6)], ncol = 3)

### radar ----
colnames(df1)

library(fmsb)
library(dplyr)
library(scales)

traits <- c("Biosynthesis_growth", "Cell_maintenance", 
            "Energy_metabolism", "Resource_acquisition", 
            "Stress_tolerance")

# 1️⃣ 按 cluster 求平均
df_trait <- df1 %>%
  group_by(cluster) %>%
  summarise(across(all_of(traits), ~ round(mean(.x, na.rm = TRUE), 4))) %>%
  as.data.frame()

# 设置行名
rownames(df_trait) <- paste("cluster", df_trait$cluster, sep = "_")
df_trait$cluster <- NULL

# 2️⃣ 标准化到 0–1（推荐）
#df_scaled <- as.data.frame(apply(df_trait, 2, rescale))
#df_scaled <- as.data.frame(rescale(as.matrix(df_trait)))

df_scaled <- as.data.frame(as.matrix(df_trait))

max(df_scaled, na.rm = TRUE)
min(df_scaled, na.rm = TRUE)

# 3️⃣ 添加 max / min 行（fmsb 必须）
df_plot <- rbind(rep(0.0945, length(traits)),
                 rep(0, length(traits)),
                 df_scaled)

colnames(df_plot) <- traits

# 4️⃣ 颜色（自动生成，适配 cluster 数量）
n <- nrow(df_trait)
colors_border <- scales::hue_pal()(n)
colors_in <- scales::alpha(colors_border, 0)

df_plot <- df_plot[,c(1,2,5,4,3)]

# 5️⃣ 画图（完全按你给的模板风格）
radarchart(
  df_plot,
  axistype = 1,
  pcol = colors_border,
  pfcol = colors_in,
  plwd = 2,
  plty = 1,
  cglcol = "grey",
  cglty = 1,
  axislabcol = "grey",
  caxislabels = c("0", "", "", "", "0.0945"),
  seg = 4,
  cglwd = 0.8,
  vlcex = 0.9
)

# 6️⃣ 图例（去掉前两行 max/min）
legend(
  x = 0.8, y = 0.2,
  legend = rownames(df_plot[-c(1,2),]),
  bty = "n",
  pch = 20,
  col = colors_border,
  text.col = "grey30",
  cex = 1,
  pt.cex = 2.5,
  y.intersp = 0.4,   # 👈 行间距（默认 ~1）
  x.intersp = 0.5    # 👈 点和文字距离
)

### plot_all

par(mfrow = c(1,4))  # 4个 cluster

for(i in 1:nrow(df_trait)) {
  
  df_single <- rbind(
    rep(0.1, length(traits)),   # max
    rep(0, length(traits)),     # min
    df_trait[i, ]
  )
  
  radarchart(
    df_single,
    axistype = 1,
    pcol = colors_border[i],
    cglcol = "grey",
    axislabcol = "grey",
    pfcol = scales::alpha(colors_border[i], 0.2),
    plwd = 3,
    cglty = 1,
   # title = rownames(df_trait)[i],
    caxislabels = seq(0,0.1, 0.02)
  )
  title(main = rownames(df_trait)[i], line = 1.2)
}

### function ratio
traits_no_self <- c("Biosynthesis_growth",
                    "Cell_maintenance",
                    "Energy_metabolism",
                    "Resource_acquisition",
                    "Stress_tolerance")
df2 <- df1 %>%
  mutate(across(all_of(traits_no_self), 
                ~ round(.x / Stress_tolerance, 2)))

### correlation ----
colnames(df1)

traits <- c("GC", "geneDensity", "Genome_Size")

functions <- c("Biosynthesis_growth",
               "Cell_maintenance",
               "Energy_metabolism",
               "Resource_acquisition",
               "Stress_tolerance")

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

###plot
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(broom)

# =========================
# 1️⃣ 计算每个 facet 的 p 和 R²
# =========================
df_stat <- df_long %>%
  group_by(Functions, Traits) %>%
  do({
    m <- lm(Abundance ~ Values, data = .)
    tidy_m <- tidy(m)
    glance_m <- glance(m)
    
    data.frame(
      p.value = tidy_m$p.value[tidy_m$term == "Values"],
      r2 = glance_m$r.squared
    )
  }) %>%
  ungroup()

# =========================
# 2️⃣ 合并 + 生成控制变量
# =========================
df_plot <- df_long %>%
  left_join(df_stat, by = c("Functions", "Traits")) %>%
  mutate(
    # p 控制显著性
    sig_p = ifelse(p.value < 0.05, "sig", "ns"),
    
    # R² 控制强度
    sig_r2 = ifelse(r2 >= 0.1, "high", "low"),
    
    # 🎨 颜色逻辑
    color_plot = ifelse(sig_p == "sig" & sig_r2 == "high",
                        as.character(Functions), "ns"),
    
    # 📏 线型逻辑
    linetype_plot = ifelse(sig_p == "sig", "dashed", "solid")
  )

# =========================
# 3️⃣ 颜色设置（默认 + 深灰）
# =========================
trait_levels <- sort(unique(df_plot$Functions))
default_cols <- scales::hue_pal()(length(trait_levels))
color_map <- setNames(default_cols, trait_levels)

# ns 用深灰
color_map <- c(color_map, ns = "grey40")

# =========================
# 4️⃣ 作图
# =========================
unique(df_plot$Traits)
df_plot$Traits <- factor(df_plot$Traits, levels = c("GC","Genome_Size","geneDensity"))

ggplot(df_plot, aes(x = Values, y = Abundance)) +
  
 # geom_point(aes(color = color_plot), alpha = 0.6,size=0.1) +
  
  # ① 显著 + 高 R² → 线性
  geom_smooth(
    data = subset(df_plot, sig_p == "sig" & sig_r2 == "high"),
    aes(color = color_plot),
    method = "lm",
    se = FALSE,
    linewidth = 1
  ) +
  
  # ② 不显著 or 低 R² → loess
  geom_smooth(
    data = subset(df_plot, !(sig_p == "sig" & sig_r2 == "high")),
    aes(color = color_plot),
    method = "gam",
    se = FALSE,
    linewidth = 0.6,
    linetype = "longdash"   # 可选：统一实线
  )+
    scale_color_manual(values = color_map) +
  #scale_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed")) +
  
  facet_grid(Functions ~ Traits, scales = "free") +
  
  theme_bw() +
  theme( panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 1)
  )+theme(panel.spacing = unit(0, "lines"))+
  theme(
    axis.title = element_blank(),   # 去标题
    axis.text = element_blank(),    # 去刻度文字
    axis.ticks = element_blank()    # 去刻度线
  )
   #theme(legend.position = "none")

