
library(ggplot2)
library(ggpmisc)
library(dplyr)
library(scales)

setwd("0_Projects/0_MetaDiverse/fastQ")

df1 <- read.table("GC_content_results_20250905_1255.txt",
                  header = TRUE,    # TRUE if first row is column names
                  sep = "\t",       # tab separator
                  stringsAsFactors = FALSE)
colnames(df1) <- c("SequenceID", "GC_R1", "GC_R2", "GC")

###
### read sam file
sam <- readRDS("0_MetaDiverse/sam.v2.rds")
colnames(sam)

# Merge data frames
df3 <- merge(df1, sam, by.x = "SequenceID", by.y = "SequenceID", all.x = TRUE)
#colnames(df3)[1] <- "SequenceID"

###

df3 <- readRDS("0_MetaDiverse/sam_v4.RDS")
df3$Transect <- factor(df3$Transect, levels = c("NS", "EW","West"))
colnames(df3)

p1 <- ggplot(data = df3, aes(x = Buffer, y=average_bacterial_archaeal_genome_size, colour = `LandUse`)) +
  geom_point( position = "dodge",alpha=0.6)+
  facet_grid(~Transect,scales="free")+
    geom_smooth(method = "loess", se = FALSE,size=0.5) +
  theme_bw()+
  scale_y_continuous(labels = scientific)+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  labs(y = "Genome size",x="Buffer")
p1

colnames(df3)
p2 <- ggplot(df3 %>% filter(C.N < 60), aes(x = C.N, y=GC, colour = `LandUse`)) +
  geom_point( position = "dodge",alpha=0.6)+
    geom_smooth(method = "lm", se = FALSE, size = 0.5, linetype = "dashed") +
    stat_poly_eq(
    aes(label = paste(..eq.label.., "*\",\"*", 
                      ..rr.label.., "*\",\"*", 
                      ..p.value.label..)),
    formula = y ~ x,
    parse = TRUE,lineheight = 2,
    size = 3
  ) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  labs(y = "GC ratio",x="C:N ratio")
p2

###
colnames(df3)
#df3[,43] <- as.numeric(gsub("^<.*", "0", df3[,43]))
colnames(df1)

p2 <- ggplot(data = df1%>% filter(GC < 70), aes(y = Genome_Size, x=GC,shape = `Landuse`)) +
  geom_point( data = df1,aes(colour=ADM),position = "dodge",alpha=0.6)+
    geom_smooth(aes(linetype = Landuse),method = "lm", se = FALSE, size = 0.5, colour="grey36") +
  scale_linetype_manual(values = c( "dashed", "twodash"))+
  stat_poly_eq(
    aes(label = paste(..eq.label.., "*\",\"*", 
                      ..rr.label.., "*\",\"*", 
                      ..p.value.label..)),
    formula = y ~ x,
    parse = TRUE,lineheight = 2,
    size = 3
  ) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_y_continuous(labels = scientific)+
  labs(x = "GC ratio",y="Genome size")+
  scale_color_gradient(low = "#0080FF", high = "#FF9700")

  #scale_color_manual(values = c("#FF9700","#0080FF"))
p2

###



# Install and load the package
library(randomForest)
library(data.table)

df3 <- readRDS("sam.v4.4.RDS")
df3$Transect <- factor(df3$Transect, levels = c("NS", "West","EW"))

colnames(df3)
colnames(df3)[c(5,8,23)] <- c("Genome_Size","GC_ratio","pH")

### density
###
library(scales)

peak_df <- df3 %>%
  group_by(Landuse) %>%
  summarise(
    density_data = list(density(Genome_Size)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    peak_x = round(density_data$x[which.max(density_data$y)], 2),
    peak_y = max(density_data$y)
  ) %>%
  ungroup()

ggplot(df3, aes(x = Genome_Size, fill = Landuse, colour = Landuse)) +
  geom_density(alpha = 0.4, linewidth = 1)+
#---- 添加虚线 ----
geom_vline(
  data = peak_df,
  aes(xintercept = peak_x, colour = Landuse),
  linetype = "dashed",
  linewidth = 0.8
) +
  
  #---- 添加峰值标签 ----
geom_text(
  data = peak_df,
  aes(x = peak_x, y = peak_y, label = paste0("", round(peak_x, 2),"Mbp")),
  hjust = 0,
  vjust = 0,
  size = 4
) +
  
  #---- X轴科学计数法、减少标签 ----
scale_x_continuous(
    breaks = pretty_breaks(n = 3)
) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1)
  ) +
  labs(x = "Gene Density (genes/Mb)",
       y = "Density",
       fill = "Land Use",
       colour = "Land Use")




#########
colnames(df3)
p2 <- ggplot(df3, aes(y = metagenome_size, x=Rain_Var, colour = `LandUse`)) +
  geom_point( position = "dodge",alpha=0.6)+
    geom_smooth(method = "loess", se = FALSE,size=0.5) +
    theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  labs(y = "Genome size") #,x="C:N ratio"
p2



#####
df4 <- df3[, c(1,5,8,23,25,26,28,35,40,42,49,53,54,55:58,65:69)]

df4[,8:10] <- as.data.frame(lapply(df4[,8:10], function(x) {
  x <- gsub("^<.*", "0", x)   # replace strings starting with "<"
  as.numeric(x)               # convert to numeric
}))

colSums(is.na(df4))
colnames(df4)


df5 <- df4 [,-1]
colnames(df5)

df5$PDD <- df5$Clim_PTA*df5$Rain_Var

df7<- df5[, -c(2,6,9,11,12,14,17,18,20,22)]
colnames(df7)
rf_model <- randomForest(`Genome_Size` ~ ., data = df7, importance = TRUE, ntree = 300)
varImpPlot(rf_model)
rf_model
# Variable importance
vip <- data.frame(importance(rf_model))

# Visualize importance
varImpPlot(rf_model)

###
row.names(vip)
vip$properties <- "soil"
#vip$properties[c(1)] <- 'biology'
vip$properties[c(10:11)] <- 'climate'

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
  labs(x = "Variable", y = "Importance (%IncMSE)", fill = "Properties") +
  #theme_bw(base_size = 14)+
  theme_minimal(base_size = 14) +
  scale_y_continuous(limits = c(5, 21), oob = scales::oob_squish) +
  scale_fill_brewer(palette = "Set2")  # 你也可以换成 "Dark2" 等

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

### PCA
# 假设 df 是你的数据框，且都是数值变量
pca_res <- prcomp(df7, scale. = TRUE) 
colnames(df7)
plot(pca_res, type = "l")
var_groups <- c("Biology","Biology",
  "Soil", "Soil", "Soil", "Soil", "Soil", "Soil","Soil", "Soil", "Soil",
  "Climate", "Climate")

fviz_pca_var(
  pca_res,
  col.var = var_groups,     # 按分组上色
  legend.title = "Property" # 图例标题
)

df8 <- cbind(df3[,c(1,8,9:17,59,62)], pca_res$x)
colnames(df8)

library(ggplot2)

ggplot(df8, aes(x = PC1, y = PC2,
                size = GC_ratio,
                color = SoilType,
                shape = LandUse)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(name = "GC ratio") +
  scale_color_discrete(name = "Soil Type") +
  scale_shape_discrete(name = "Transect") +
  theme_bw(base_size = 14) +
  labs(x = "PC1", y = "PC2",
       title = "PCA Scatter Plot: PC1 vs PC2") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

#### sedv

library(dplyr)
library(ggplot2)

df3 <- readRDS("sam_v4.RDS")
df3$Transect <- factor(df3$Transect, levels = c("NS", "West","EW"))

colnames(df3)[c(5,8,23)] <- c("Genome_Size","GC_ratio","pH")

df4 <- df3[,c(5,8,62)]
#df4$Genome_Size <- df4$Genome_Size/1000000
colnames(df4)

# 计算每个 SoilType 的均值和标准差
df_summary <- df4 %>%
  group_by(SoilType) %>%
  summarise(
    mean_GS = mean(Genome_Size, na.rm = TRUE),
    sd_GS = sd(Genome_Size, na.rm = TRUE),
    mean_GC = mean(GC_ratio, na.rm = TRUE),
    sd_GC = sd(GC_ratio, na.rm = TRUE)
  )

# 作图
ggplot() +
  # 原始点
  geom_point(data = df4,
             aes(x = Genome_Size, y = GC_ratio,
                 color = SoilType, shape = SoilType),
             alpha = 0.6, size = 0.5) +
  
  # 每组均值点
  geom_point(data = df_summary,
             aes(x = mean_GS, y = mean_GC,
                 color = SoilType, shape = SoilType),
             size = 1, stroke = 1) +
  
  # 水平误差线（Genome_Size 的 SD）
  geom_errorbarh(data = df_summary,
                 aes(y = mean_GC,
                     xmin = mean_GS - sd_GS,
                     xmax = mean_GS + sd_GS,
                     color = SoilType),
                 height = 0.2, size = 0.3) +
  
  # 垂直误差线（GC_ratio 的 SD）
  geom_errorbar(data = df_summary,
                aes(x = mean_GS,
                    ymin = mean_GC - sd_GC,
                    ymax = mean_GC + sd_GC,
                    color = SoilType),
                width = 0.01, size = 0.3) +
  scale_shape_manual(values = c(3,8,7,10,16,17,18,15)) +
  labs(x = "Genome Size (MB)",
       y = "GC Ratio",
       color = "Soil Type",
       shape = "Soil Type") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black")
  )

