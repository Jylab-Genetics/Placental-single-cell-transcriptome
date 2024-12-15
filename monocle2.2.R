
########################*************************************************************************
########################*基因可视化************************************************************
library(dplyr)
library(monocle)
library(ggplot2)
library(Seurat)
library(reshape2)

colnames(pData(dog))#cds是已经做好拟时的monocle2对象
#将需要展示的基因表达量添加到cds， 这里我使用的是log2处理，也可以使用log10进行标准化
genes <- c('CEP152',
           'CEPT1',
           'GPATCH8',
           'GTF2E2',
           'IFT74',
           'PPP3CB',
           'TMEM218',
           'ZNF182','CDH2','CLDN3','EPAS1')

genes_exp <- list()

for(i in 1:length(genes)){
  A <- log2(exprs(human)[genes[i],]+1)
  A <- as.data.frame(A)
  genes_exp[[i]] <- A
}

gene_exp <- do.call(cbind, genes_exp)
colnames(gene_exp) <- genes
#将上述几个基因的拟时表达添加到monocle
pData(dog) = cbind(pData(dog), gene_exp)


#提取作图数据，只需要游基因表达和拟时即可
data <- pData(dog)

#选择需要的列即可，我这里的origin.ident就是分组
data$species<-'dog'
head(data)

data<-data %>% select("species","subclass2","Pseudotime",
                      'CEP152',
                      'CEPT1',
                      'GPATCH8',
                      'GTF2E2',
                      'IFT74',
                      'PPP3CB',
                      'TMEM218',
                      'ZNF182','CDH2','CLDN3','EPAS1')
#ggplot作图

# data<-data %>% select("species","celltypes_Abbreviation1","Pseudotime",'SEC23A',
#                       'GRN',
#                       'BNC1',
#                       'SLC43A3',
#                       'FUT11',
#                       'LTBP1',
#                       'SLC6A2',
#                       'SEC61A1',
#                       'SLC16A3',
#                       'PDGFRB',
#                       'TOR1B',
#                       'TM9SF2',
#                       'CPEB4')
#使用分屏，应该就是文献中的办法
#首先将data宽数据转化为长数据
head(data)
data_long_m<-melt(data, id.vars = c("species","subclass2","Pseudotime"), #需保留的不参与聚合的变量列名
                  measure.vars = 4:14,#选择需要转化的列
                  variable.name = c('gene'),#聚合变量的新列名
                  value.name = 'value')#聚合值的新列名
colnames(data_long_m)
head(data_long_m)
#除了dog
#其他物种在这里要提取相应的亚群
dog_long_m <- subset(data_long_m, subclass3 %in% c('CTB','STB'))

# 生成图形并将背景设为透明
ggplot(dog_long_m, aes(x = Pseudotime, y = value, color = gene)) +
  geom_smooth(aes(fill = gene)) + # 平滑的填充
  xlab('pseudotime') + 
  ylab('Relative expression') +
  theme(
    axis.text = element_text(color = 'black', size = 12),
    axis.title = element_text(color = 'black', size = 14),
    strip.text = element_text(color = 'black', size = 14), # 分面标题
    panel.background = element_rect(fill = "transparent", color = NA), # 设置面板背景透明
    plot.background = element_rect(fill = "transparent", color = NA), # 设置绘图区域背景透明
    panel.grid.major = element_blank(), # 移除主要网格线
    panel.grid.minor = element_blank(), # 移除次要网格线
    panel.border = element_blank(), # 移除面板边框
    axis.line = element_line(color = 'black') # 保留 x 和 y 轴
  ) +
  scale_color_manual(name = NULL, values = c("#089E86","#3D5387","#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00")) + # 修改颜色
  scale_fill_manual(name = NULL, values = c("#089E86","#3D5387","#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00")) # 修改颜色

ggsave("EPAS1_exp_dog_stb.pdf")

#绘制平滑曲线并关闭阴影
ggplot(dog_long_m, aes(x = Pseudotime, y = value, color = gene)) +
  geom_smooth(aes(fill = gene), method = "loess", se = FALSE) + # 使用geom_smooth并关闭阴影
  xlab('pseudotime') + 
  ylab('Relative expression') +
  theme(
    axis.text = element_text(color = 'black', size = 12),
    axis.title = element_text(color = 'black', size = 14),
    strip.text = element_text(color = 'black', size = 14), # 分面标题
    panel.background = element_rect(fill = "transparent", color = NA), # 设置面板背景透明
    plot.background = element_rect(fill = "transparent", color = NA), # 设置绘图区域背景透明
    panel.grid.major = element_blank(), # 移除主要网格线
    panel.grid.minor = element_blank(), # 移除次要网格线
    panel.border = element_blank(), # 移除面板边框
    axis.line = element_line(color = 'black') # 保留 x 和 y 轴
  ) +
  scale_color_manual(name = NULL, values = c("#089E86","#3D5387","#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00")) + # 修改颜色
  scale_fill_manual(name = NULL, values = c("#089E86","#3D5387","#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00")) # 修改颜色


#**********************************************************************************
#************获取相应通路****************************************************
##########
#代谢通路随拟时序改变
setwd("H:\\result5\\代谢和KEGG通路改变\\STB代谢")
rm(list = ls())
library(KEGGREST)
library(stringr)
##提取KEGG所有的通路
hsa_path <- keggLink("pathway","hsa")
#KEGG共收集8443个相关通路的基因
length(unique(names(hsa_path)))
#涉及通路353个
length(unique(hsa_path))
#代谢相关的通路在KEGG中以“00……”开头，提取”00……“开头通路共84个
pathway= unique(hsa_path)[grepl('hsa00',unique(hsa_path))];length(pathway)
#提取这84个通路的基因
hsa_info <- lapply(pathway, keggGet)


#提取通路和基因的信息
meta= unique(hsa_path)[grepl('hsa00',unique(hsa_path))]
hsa_info <- lapply(meta, keggGet)  

nm=unlist(lapply( hsa_info , function(x) x[[1]]$NAME))
library(stringr)
genes = unlist(lapply( hsa_info , function(x) {
  g = x[[1]]$GENE
  paste(str_split(g[seq(2,length(g),by=2)],';',simplify = T)[,1],collapse =';')
}))
df =data.frame(
  hsa= meta,
  nm=nm,
  genes =genes
)

head(df)

# #分别把基因名字和基因的ID提取出来，合成一个数据框
# gene_symbol = unlist(lapply( hsa_info , function(x) {
#   g = x[[1]]$GENE
#   str_split(g[seq(2,length(g),by=2)],';',simplify = T)[,1]
# }))
# 
# gene_ID = unlist(lapply( hsa_info , function(x) {
#   g = x[[1]]$GENE
#   paste0("hsa:",g[seq(1,length(g),by=2)])
# }))
# 
# genelist <- data.frame(gene_ID,gene_symbol)
# 
# #去重复,提取84个代谢通路中1728个基因
# genelist <- genelist[!duplicated(genelist$gene_symbol),]
# length(genelist$gene_symbol)
# 
# head(genelist)
# write.csv(genelist,"genelist.csv")
#保存为Rdata，随时可以加载使用
save(df,file = "metabolism.Rdata")

write.csv(df,"metabolism_dog.csv")


########***********************************************
#######         通路随拟时表达变化
#Rscript pathway_diff_trajectory.R



#*********************************************************

#########**************************************************************
#########**************************************************************
##########结果可视化
library(dplyr)
library(monocle)
library(ggplot2)
library(Seurat)
library(reshape2)

pathways_pseu<-read.csv("dog_combined_pathways_data.csv") #全部同源后KEGG通路的数据
metabolism<-read.csv("metabolism_dog.csv")
metabolismID<-metabolism$ID
# 提取 dog 数据框中 pathway 列值存在于 metabolismID 的行
# filtered_pp <- pathways_pseu %>%
#   filter(pathway %in% metabolismID)
# head(filtered_pp)

 filtered_pp <- pathways_pseu %>%
   filter(pathway %in% c('hsa00380','hsa00600','hsa00630','hsa00982'))
 head(filtered_pp)

 #,'hsa00380','hsa00600','hsa00630','hsa00982'

#根据分化线路纤薄类型在筛选
#CTB-STB
filtered_celltype <- filtered_pp[filtered_pp$celltypes_merge %in% c('dog_CTB','dog_STB'), ]

table(filtered_celltype$pathway)

#作图



ggplot(filtered_celltype, aes(x = pseudotime, y = value, color = pathway)) +
  geom_smooth(aes(fill = pathway),show.legend = FALSE) + # 平滑的填充
  xlab('pseudotime') + 
  ylab('Relative expression') +
  facet_wrap(~ pathway, scales = "free_y") + # 使用 pathway 进行分面
  theme(axis.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 12),
        strip.text = element_text(color = 'black', size = 14)) 

# +scale_color_manual(name=NULL, values = c("#089E86","#3D5387"))+#修改颜色
#scale_fill_manual(name=NULL, values = c("#089E86","#3D5387"))#修改颜色

ggsave('dog_metabolism_pseu.pdf')

#绘制在一张图上
# 生成77种不同的颜色
library(RColorBrewer)
color_palette <- brewer.pal(8, "Dark2")
extended_colors <- colorRampPalette(color_palette)(77)
filtered_pp$pathway<-as.factor(filtered_pp$pathway)
filtered_pp$orig.ident<-as.factor(filtered_pp$orig.ident)

ggplot(filtered_celltype, aes(x = pseudotime, y = value, color = pathway)) +
  geom_smooth(aes(fill = pathway)) + # 平滑的填充
  xlab('pseudotime') + 
  ylab('Relative expression') +
  theme(
    axis.text = element_text(color = 'black', size = 12),
    axis.title = element_text(color = 'black', size = 14),
    strip.text = element_text(color = 'black', size = 14), # 分面标题
    panel.background = element_rect(fill = "transparent", color = NA), # 设置面板背景透明
    plot.background = element_rect(fill = "transparent", color = NA), # 设置绘图区域背景透明
    panel.grid.major = element_blank(), # 移除主要网格线
    panel.grid.minor = element_blank(), # 移除次要网格线
    panel.border = element_blank(), # 移除面板边框
    axis.line = element_line(color = 'black') # 保留 x 和 y 轴
  ) +
  scale_color_manual(name = NULL, values = c("#089E86","#3D5387","#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00")) + # 修改颜色
  scale_fill_manual(name = NULL, values = c("#089E86","#3D5387","#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00")) # 修改颜色


# Note: 如果你有更多或更少的pathway，可以相应调整生成颜色的数量。

