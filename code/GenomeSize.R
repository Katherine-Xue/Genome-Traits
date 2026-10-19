
library(ggplot2)
library(scales)
library(ggpmisc)
library(dplyr)
library(tidyr)
library(corrplot)


setwd("0_MetaDiverse/MicroCensus")

df0 <- readRDS("sam.v4.4.RDS")
df0$Transect <- factor(df0$Transect, levels = c("NS", "EW","West"))

colnames(df0)
df1 <- df0[,c(1:12,46:63,73:87)]

# df0 <- readRDS("D:/OneDrive - The University of Sydney (Staff)/0_Projects/0_MetaDiverse/sam.v4.RDS")
# colnames(df0)
# 
# ### update names
# #colnames(df0)[c(2,17)] <- c("Genome Size (Mbp)", "GC content (%)")
# 
# colnames(df0)[c(35,38:40,62:69)] <- c(
#   "pH","clay","sand","silt","Iron","TS","TK","TP","TN","TC","CN","SOC")
# 
# colnames(df0)[c(74:78)] <- c("MAP","TRA","MAT","ARV","ADM")
# 
# 
# ### subset ----
# 
# df1 <- df0[,c(1:10,17,18:28,35:39,62:79)]
# colnames(df1)
# 
# # df1$average_genome_size <- df1$average_genome_size/1e6 # convert to Mbp
# # colnames(df1)[colnames(df1) == "average_genome_size"] <- "GenomeSize"
# df1$`Cells(log10)` <- log10(df1$cells)
# 
# cols <- c("Iron","TK", "TP", "TS")
# 
# df1[,cols] <- as.data.frame(lapply(df1[,cols], function(x) {
#   x <- gsub("^<.*", "0", x)   # replace strings starting with "<"
#   as.numeric(x)               # convert to numeric
# }))


### models ----
library(piecewiseSEM)
library(lme4)

colnames(df1)

sem_mod <- psem(
  lmer(SOC ~ MAT + PTA+MAP+ARV +pH+ (1 | LandUse), data = df1),
lmer(GenomeSize ~ SOC +`PH.H2O.` + (1 | LandUse), data = df1),
  lmer(`GC_mean` ~ SOC + `PH.H2O.`+Clim_MAT + (1 | LandUse), data = df1),
  lmer(Cell.log ~ GenomeSize + `GC_mean` + (1 | LandUse), data = df1)
)

summary(sem_mod)


### RDA ----
library(vegan)
library(compositions)
colnames(df1)

Y <- df1[, c("DNA_Repair", "Metabolism", "Stress_Response", "Transporters")]
X <- df1[, c(2,11,23:27,28,29:35,40:45)]
colnames(X)

Y_clr <- clr(Y + 1e-6)

rda_res <- rda(Y_clr ~ ., data = X)
anova(rda_res, by = "term")

anova_term <- anova(rda_res, by = "term")

sig_vars <- intersect(
  rownames(anova_term)[
    !is.na(anova_term$`Pr(>F)`) &
      anova_term$`Pr(>F)` < 0.05
  ],
  colnames(X))

envfit_res <- envfit(rda_res, X[, sig_vars], permutations = 999)

###
library(tidyverse)

site_scores <- scores(rda_res, display = "sites", scaling = 2) %>%
  as.data.frame()
site_scores <- cbind(df1, site_scores)

env_scores <- scores(envfit_res, display = "vectors") %>%
  as.data.frame() %>%
  rownames_to_column("Variable")

ggplot() +
  geom_point(data = site_scores %>% filter(RDA2>-1),
             aes(x=RDA1, y=RDA2,colour = MAT, shape = Landuse),
             size = 2, alpha = 0.7) +
  geom_segment(data = env_scores,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "firebrick") +
  geom_text(data = env_scores,
            aes(RDA1, RDA2, label = Variable),
            colour = "firebrick",
            hjust = 1.1, vjust = 1.1) +
  scale_colour_viridis_c(option = "D")+
  coord_fixed(ratio = 1)+
  theme_bw() +theme(aspect.ratio = 1)+
  labs(x = "RDA1", y = "RDA2")

#### remove multicollinearity ----
library(corrplot)
library(caret)

colnames(df1)
Y <- df1[, c("GC", "Genome_Size", "geneDensity")]
X <- df1[, c(13:22,24:30)]
colnames(X)
# 环境变量（只放连续型变量）
env <- X   # 你的 X 已经是环境矩阵
colnames(env)
# 响应变量（RDA 用）

cor_mat <- cor(env, method = "pearson", use = "pairwise.complete.obs")

corrplot(
  cor_mat,
  method = "color",
  type = "upper",
  tl.cex = 0.8,
  col = colorRampPalette(c("blue", "white", "red"))(200)
)

cor_cut <- 0.9 

remove_idx <- findCorrelation(
  cor_mat,
  cutoff = cor_cut,
  exact = TRUE
)

env_nocor <- env[, -remove_idx]

colnames(env_nocor)

###
anova_df <- as.data.frame(anova(rda_res, by = "term")) %>%
  rownames_to_column("Variable") %>%
  filter(`Pr(>F)` < 0.05)

ggplot(anova_df,
       aes(x = reorder(Variable, F),
           y = F)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_bw() +
  labs(x = "Environmental factors",
       y = "F value")

#### genome by soiltype ----

df2 <- df1 %>%
  mutate(SoilType = reorder(SoilType, average_genome_size, FUN = mean))

# plot boxplot
ggplot(df2, aes(x = SoilType, 
                y = average_genome_size, 
                fill = SoilType)) +
  geom_boxplot(alpha = 0.7) +
  labs(x = "Soil Type", y = "Average bacterial/archaeal genome size") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

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


###
df2 <-  df0 %>%
  pivot_longer(
    cols = c(4,14),
    names_to = "Methods",
    values_to = "Size"
  ) %>%
  relocate(Methods, Size, .after = SequenceID)

ggplot(df2, aes(x = SoilType, y = Size/1000000, fill = Methods)) +
  geom_boxplot(alpha = 0.7) 

#### pairwise correlation ----
library(corrplot)
library(Hmisc)
library(reshape2)

colnames(df1)

df6 <- df1[,c(31,35,34,32,33,13:30)]
colnames(df6)
df6[df6 == 0] <- NA
colnames(df6)[1:3] <- c("GC","Genome Size","Gene Density")


#cor_mat <- cor(df6, use = "pairwise.complete.obs", method = "pearson")
colnames(df6)


#### person correlation plot ----
# 计算相关性
library(Hmisc)
library(reshape2)

cor_res <- rcorr(as.matrix(df6), type = "pearson")
cor_mat <- cor_res$r
p_mat   <- cor_res$P

# 转为长格式
cor_df <- melt(cor_mat, varnames = c("Var1", "Var2"), value.name = "cor")
p_df   <- melt(p_mat, varnames = c("Var1", "Var2"), value.name = "p")

plot_df <- left_join(cor_df, p_df, by = c("Var1","Var2"))

# 添加显著性星号
plot_df$signif <- cut(plot_df$p,
                      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                      labels = c("***", "**", "*", ""))

# X/Y 变量筛选
colnames(df6)

y_vars <-  colnames(df6)[c(1:3)]
x_vars <- colnames(df6)[c(5:13,15,16,18:21)]

plot_df2 <- plot_df %>%
  filter(Var1 %in% y_vars & Var2 %in% x_vars)

plot_df2$Var1 <- factor(plot_df2$Var1, levels = y_vars)
plot_df2$Var2 <- factor(plot_df2$Var2, levels = rev(x_vars))


# 计算点大小（绝对相关性）
plot_df2$size_val <- abs(plot_df2$cor)

# 绘图
ggplot(plot_df2, aes(x = Var1, y = Var2)) +
  geom_point(aes(color = cor, size = size_val)) +
  geom_text(aes(label = signif),
            color = "black",
            size = 4,
            vjust = 0.5,   # 正中
            hjust = 0.5) + # 正中
  scale_color_gradient2(low = "blue", high = "red", mid = "#FFFFF0",
                        midpoint = 0, limit = c(-1, 1)) +
  scale_size(range = c(1, 10)) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 0, hjust =0.8),
        panel.grid = element_blank(),
        legend.title = element_text(size = 10),
        legend.text  = element_text(size = 10)) +
  labs(x = "",
       y = "",
       color = "Correlation",
       size = "|Correlation|")

### PCA analysis ---
library(factoextra)
colnames(df6)

df7 <- df6[,c(-4,-14,-17)]
colnames(df7)
#colnames(df7)[1:3] <- c("Genome Size","GC","Cells")

#df7[is.na(df7)] <- 0

pca_res <- prcomp(df7, scale. = TRUE) 

#plot(pca_res, type = "l")
var_groups <- c(
  rep("Biology", 3),
  rep("Soil", 11),
  rep("Climate", 4))


fviz_pca_var(
  pca_res,
 col.var = var_groups,     # 按分组上色
  legend.title = "Property", # 图例标题
  repel = TRUE)

df8 <- cbind(df1, pca_res$x)
colnames(df8)

library(ggplot2)

ggplot(df8, aes(x = PC1, y = PC2,
                #size = average_genome_size/1000000,
                color = MAT,
                shape = SoilType)) +
  geom_point(alpha = 0.8) +
  # scale_size_continuous(name = "Genome size (Mb)",
  #                       range = c(2, 7)) +
 # scale_color_continuous(name = "Soil Type") +
  scale_color_viridis_c(
    option = "D",
    name = "MAT (°C)"
  ) +
  #scale_shape_discrete(name = "SoilType") +
  scale_shape_manual(
    name = "Soil type",
    values = c(20, 18, 17, 8, 7, 4, 25, 0)
  ) +
  coord_fixed(ratio = 1)+
  theme_bw(base_size = 13) +
  labs(x = "PC1 (46.3%)", y = "PC2 (42%)"
     #  title = "PCA Scatter Plot: PC1 vs PC2"
     ) +
  theme(aspect.ratio = 1,
   # plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )


#### Random Forest---- 

library(randomForest)
library(missForest)

colnames(df6)

y_imp <-df6[,3]
                
x_imp <-df6[, c(3:21)]

#x_imp<- missForest(x_imp)$ximp

x_imp[is.na(x_imp)] <- 0


#set.seed(123)
rf_model <- randomForest(
  x = x_imp,     # predictors
  y = y_imp,        # response (column 1)
  importance = TRUE,
  mtry=5,
  nodesize=5,
  ntree = 500)

rf_model

varImpPlot(rf_model)

# Variable importance
vip <- data.frame(importance(rf_model))

# Visualize importance
varImpPlot(rf_model)

###
row.names(vip)
vip$properties <- "soil"
vip$properties[c(1:2)] <- 'biology'
vip$properties[c(14:17)] <- 'climate'

vip$IncNodePurity <- round(vip$IncNodePurity,1)
vip$X.IncMSE <- round(vip$X.IncMSE, 1)

vip$variable <- rownames(vip)
vip$variable <- factor(vip$variable, 
                       levels = vip$variable[order(vip$X.IncMSE)])

vip$properties <- factor(vip$properties, 
                         levels = c("climate","soil","biology")) #"biology"

ggplot(vip, aes(x = variable, y = X.IncMSE, fill = properties)) +
  geom_col(width = 0.6) +
  coord_flip() +
  labs(x = "Variable", y = "Importance (%IncMSE)", title = "GC content",fill = "") +
  #theme_bw(base_size = 14)+
  theme_minimal(base_size = 16) +
  #scale_y_continuous(limits = c(5, 21), oob = scales::oob_squish) +
  scale_fill_brewer(palette = "Set2")  # 你也可以换成 "Dark2" 等

### 折线图
#p2 <- list()
nn <- colnames(df6)[8]

p2[[4]] <- ggplot(vip, aes(x = variable, y = X.IncMSE)) +
  geom_line(group = 1, linewidth = 1, color = "grey50") +
  geom_point(aes(fill = properties),
             shape = 21, size = 3, color = "black") +
  coord_flip() +
  labs(
    x = "Variable",
    y = "Importance (%IncMSE)",
    title = nn,
    fill = ""
  ) +
  theme_minimal(base_size = 16) +
  scale_fill_brewer(palette = "Set2")

library(cowplot)

plot_grid(plotlist = p2[c(2, 1, 4, 3)], ncol = 4)


###
vip$variable <- factor(vip$variable, 
                       levels = vip$variable[order(vip$IncNodePurity)])

ggplot(vip, aes(x = variable, y = IncNodePurity, fill = properties)) +
  geom_col(width = 0.6) +
  coord_flip() +
  scale_y_continuous(
    breaks = pretty_breaks(n = 3)
  )+
  labs(x = "Variable", y = "Importance (IncNodePurity)", fill = "Properties") +
  #theme_bw(base_size = 14)+
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set2")  # 你也可以换成 "Dark2" 等


#### SHAP----
library(SHAPforxgboost)
library(xgboost)
library(ggplot2)

colnames(df7)
env <- df7[,c(5:18)]
var <- df7[,1:3]

y_var <-  colnames(var)
y_var
# 删除目标变量列并转换为矩阵
dataX <- as.matrix(env)
colnames(dataX) <- colnames(env)

# hyperparameter tuning results
library(xgboost)

# 1. 构建 DMatrix（关键）
dtrain <- xgb.DMatrix(
  data = dataX,
  label = var[[y_var]]
)

# 2. 参数
params <- list(
  objective = "reg:squarederror",
  eta = 0.01,
  max_depth = 10,
  gamma = 0.01,
  subsample = 0.98,
  colsample_bytree = 0.86
)

# 3. 训练模型（稳定）
mod <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 300,
  verbose = 0
)

# To store SHAP values and plots
shap_long.all <- data.frame()
shap_long <- list()
p <- list()
g <- list()

for (i in 1:3) {
  y_var <- colnames(var)[i]
  y_var
  
  dtrain <- xgb.DMatrix(
    data = dataX,
    label = var[[y_var]]
  )
  
  mod <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 300,
    verbose = 0
  )
  
  shap_values <- shap.values(xgb_model = mod, X_train = dataX)
  shap_values$mean_shap_score
  
    # To prepare the long-format data:
  shap_long1 <- shap.prep(xgb_model = mod, X_train = dataX)
  shap_long1$group <-y_var 
  shap_long[[y_var]] <- shap_long1
  
  p[[y_var]] <- shap.plot.summary(shap_long1,x_bound=0.5,dilute = 5)+
    ggtitle(y_var)+
    theme(plot.title = element_text(hjust = 0.5))
  
  shap_long.all <- rbind(shap_long.all,shap_long1)
  
  plot_data <- shap.prep.stack.data(shap_contrib = shap_values$shap_score, top_n = 5, n_groups = 1)
  
  g[[y_var]] <- shap.plot.force_plot(plot_data,zoom_in = FALSE)+
    ggtitle(y_var)+
    theme(plot.title = element_text(hjust = 0.5))
  }

####
library(gridExtra)
grid.arrange(grobs = p, ncol = 3)

#g <- g[c(1,4,2,5,3,6)]
names(g)
grid.arrange(grobs = g, ncol = 1)


#### linear regression plot ----
library(mgcv)
library(purrr)

colnames(df1)
df6 <- df1[,c(1,13:22,2,11,44,3,7:12,23:43)]
colnames(df6)
df1$Transect <- factor(df1$Transect, levels = c("NS", "West","EW"))

ggplot(data = df1 %>% filter(cells<1.5e10), aes(y = `CN`, #%>% filter(EC < 400)
                             x = `Buffer`, 
                             colour = Landuse)) +
  facet_wrap(~Transect, scales = "free") +
  geom_point(position = "dodge", alpha = 0.6) +
  #geom_smooth(method = "loess", se = FALSE, size = 0.5,linetype = "dashed") +
 geom_smooth(method = "gam",formula = y ~ s(x, k = 4),se = FALSE, size = 0.5) +
 # geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "dashed") +
 # stat_poly_eq(
 #   aes(label = paste(..eq.label.., "*\",\"*",
 #                     ..rr.label.., "*\",\"*",
 #                     ..p.value.label..)),
 #   formula = y ~ x,
 #   parse = TRUE,
 #   size = 3
 # ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  labs(y = "C:N",x="Buffer") #+scale_y_continuous(labels = scientific)

### multiple regression all----
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpmisc)

vars <- c("ADM","PVA","TRA","CEC","TS","TK","TP","TN","pH","SOC","CN","clay")

df_long <- df6 %>%
 # filter(cells < 1.5e10) %>%
  pivot_longer(cols = all_of(vars),
               names_to = "Variable",
               values_to = "Value")
colnames(df_long)

ggplot(df_long, aes(x = Value,
                    y = geneDensity,
                    colour = Landuse)) +
  geom_point(alpha = 0.6) +
  
  # GAM（可选）
  geom_smooth(method = "gam",
              formula = y ~ s(x, k = 4),
              se = FALSE,
              size = 0.5) +
  
  # 线性拟合
  geom_smooth(method = "lm",
              se = FALSE,
              size = 0.5,
              linetype = "dashed") +
  
  # 方程 + R² + p值
  stat_poly_eq(
    aes(label = paste(..eq.label.., "*\",\"*",
                      ..rr.label.., "*\",\"*",
                      ..p.value.label..)),
    formula = y ~ x,
    parse = TRUE,
    size = 3
  ) +
  
  facet_wrap(~Variable, ncol = 3, scales = "free_x") +
  
  theme_bw() +
  theme(
   # strip.background = element_blank(),   # 去掉facet标题背景
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  
  labs(y = "Gene density (genes/Mb)", x = NULL)

###GC-size ----
library(ggpmisc)
library(ggnewscale)

colnames(df0)

ggplot(data = df0%>% filter(`GC` < 70), 
       aes(y = `geneDensity`, x=`GC`,shape = `Landuse`)) +
  geom_point(aes(colour=Y_Latitude),alpha=0.6)+
  scale_color_gradient(low = "#0080FF", high = "#FF9700")+
  #facet_grid(~Transect,scales="free")+
  # scale_fill_manual(values = c("#0077b6","#dda15e"))+
  #geom_text(data=df3,aes(x = SoilType,y=Freq,label = Freq), vjust = -0.5, color = "black")+
  #geom_smooth(method = "loess", se = FALSE,size=0.5) +
  new_scale_color() + 
  geom_smooth(aes(linetype = Landuse,colour=Landuse),method = "lm", se = FALSE, size = 0.5) + 
  #scale_colour_manual(values = c("#0080FF", "#FF9700")) +
  scale_linetype_manual(values = c("dashed", "twodash")) +
  stat_poly_eq(
    aes(colour=Landuse,label = paste(..eq.label.., "*\",\"*", 
                      ..rr.label.., "*\",\"*", 
                      ..p.value.label..)),
    formula = y ~ x,
    parse = TRUE,lineheight = 2,
    size = 4
  ) +
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  #scale_y_continuous(labels = scientific)+
  labs(x = "GC (%)",y="Gene Density (genes/Mb)")
  #scale_color_gradient(low = "#0080FF", high = "#FF9700")

colnames(df1)
ggplot(data = df1%>% filter(CN<60), 
       aes(y = `cells`, x=`CN`,shape = `Landuse`,colour=Landuse)) +
  geom_point(position = "dodge",alpha=0.6)+
  geom_smooth(method = "loess", se = FALSE,size=0.5,span = 0.8) +
  geom_smooth(aes(colour=Landuse),method = "lm", se = FALSE, size = 0.5,linetype = "dashed") +
  stat_poly_eq(
    aes(colour=Landuse,label = paste(..eq.label.., "*\",\"*",
                                     ..rr.label.., "*\",\"*",
                                     ..p.value.label..)),
    formula = y ~ x,
    parse = TRUE,lineheight = 2,
    size = 4
  ) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  #scale_y_continuous(labels = scientific)+
  labs(y = "Gene Density (genes/Mb)") #,x="pH"
  # stat_ellipse(
  #   aes(group = Landuse, colour = Landuse),
  #   type = "norm",     # 正态置信椭圆
  #   level = 0.95,      # 95% 置信区间
  #   linewidth = 0.3
  # )


#### GAM（Generalized Additive Model）----
library(mgcv)
library(dplyr)
library(ggplot2)
library(purrr)

# ----------------------------
# 1. 数据筛选
# ----------------------------
df_gam <- df6 %>%
  filter(cells < 1.5e10, pH > 4.5)

# ----------------------------
# 2. 按 Landuse 拟合 GAM 并提取统计量
# ----------------------------
gam_stats <- df_gam %>%
  group_by(Landuse) %>%
  nest() %>%
  mutate(
    gam_mod = map(data, ~ gam(cells ~ s(pH, k = 4), data = .x)),
    edf = map_dbl(gam_mod, ~ summary(.x)$s.table[1, "edf"]),
    p_value = map_dbl(gam_mod, ~ summary(.x)$s.table[1, "p-value"])
  ) %>%
  ungroup()

# ----------------------------
# 3. 生成注释文本和放置位置
# ----------------------------
label_df <- df_gam %>%
  group_by(Landuse) %>%
  summarise(
    x = min(pH) + 0.05 * diff(range(pH)),
    y = max(cells) * 0.95,
    .groups = "drop"
  ) %>%
  left_join(
    gam_stats %>% dplyr::select(Landuse, edf, p_value),
    by = "Landuse"
  ) %>%
  mutate(
    label = paste0(
      "edf = ", round(edf, 2),
      "\n", "p = ", format.pval(p_value, digits = 2)
    )
  )

# ----------------------------
# 4. 作图：散点 + GAM 曲线 + 统计注释
# ----------------------------
ggplot(
  df_gam,
  aes(x = pH, y = cells, colour = Landuse)
) +
  geom_point(alpha = 0.6) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 4),
    se = FALSE,
    size = 0.6,
    linetype = "dashed"
  ) +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label, colour = Landuse),
    hjust = 0,
    vjust = 1,
    show.legend = FALSE,
    size = 3.5
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "Soil pH",
    y = "Microbial cells (g⁻¹ soil)",
    colour = "Land use"
  )
