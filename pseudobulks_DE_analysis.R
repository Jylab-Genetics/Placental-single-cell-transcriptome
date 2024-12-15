#******************************************************
##################单细胞差异基因
#******************************************************
library(ggrepel)
library(Seurat)
#**************************
#1、FindMarkers/FindAllmarkers 分析
#*****************************
#概念：
#● FindMarkers/FindAllmarkers是直接在单细胞RNA测序数据的细胞级别进行差异表达分析。
##它使用统计检验方法（如Wilcoxon秩和检验、MAST或t检验）来比较不同细胞群体之间的基因表达差异。
#优点：
#● 保留了单细胞数据的分辨率，能够捕捉到细胞群体内的异质性。
#● 可以分析细胞亚群之间的差异，适合于细胞类型复杂的研究。
#缺点：
#● 由于单细胞数据的稀疏性和高噪声，差异表达分析的结果可能不如pseudobulk分析稳定。
#● 分析结果依赖于选择的细胞群体，因此需要谨慎进行群体定义和数据预处理。

#****************************************
#* 单个样本
#****************************************
markers <- FindMarkers(scRNA_tumor,ident.1 = 'acitve', only.pos = FALSE, 
                       min.pct = 0.25, logfc.threshold = 0.25)


#****************************************
#循环两两组比较
#****************************************
# 获取所有物种对的唯一组合
species_combinations <- combn(unique(placenta$species), 2)
# 循环遍历每一对物种
for (i in 1:ncol(species_combinations)) {
  species1 <- species_combinations[1, i]
  species2 <- species_combinations[2, i]
  # 找到两个物种之间的差异基因
  DEGs <- FindMarkers(placenta,
                      ident.1 = species1,
                      ident.2 = species2,
                      group.by = "species",
                      logfc.threshold = 0,
                      min.pct = 0)
  # 生成文件名
  filename <- paste("Comparison_", species1, "_vs_", species2, ".csv", sep="")
  # 将结果写入 CSV 文件
  write.csv(DEGs, file=filename, row.names=TRUE)
}



#****************************************
#*作图
#****************************************
#常规火山图的代码：
#增加一列显著基因，待会标名字
markers$sign <- ifelse(markers$p_val_adj < 0.005 & abs(markers$avg_log2FC) > 2,
                       rownames(markers),NA)
#当然也可以自定义的，随机的
k <- c("TP53","CD34","CD68")
markers <- markers %>%
  mutate(Sign = ifelse(rownames(markers) %in% k, rownames(markers), NA))

#自定义阈值
log2FC = 0.585
padj = 0.05

colnames(markers)
library(ggplot2)
library(ggrepel) # 这个R包用于添加基因名称
ggplot(data = markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha=0.4, size=3.5, aes(color=change)) +
  #scale_color_manual(values=c("blue", "grey","red"))+
  scale_color_manual('change',labels=c(paste0("down(",table(markers$change)[[1]],')'),
                                       'ns',
                                       paste0("up(",table(markers$change)[[3]],')' )),
                     values=c("blue", "grey","red" ))+
  geom_vline(xintercept=c(-log2FC,log2FC),lty=4,col="black",linewidth=0.8) +
  geom_hline(yintercept = -log10(padj),lty=4,col="black",linewidth=0.8) +
  geom_label_repel(aes(label = sign ), fontface="bold", 
                   color="black", box.padding=unit(0.35, "lines"),
                   point.padding=unit(0.5, "lines"), segment.colour = "grey50") +
  theme_bw()


#新的火山图的代码：
# 单细胞火山图展示方式-new
# Findmarkers参数设置为0！
markers <- FindMarkers(scRNA_tumor,ident.1 = 'acitve', 
                       only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
# 

markers <- markers %>%
  mutate(Difference = pct.1 - pct.2) %>% 
  rownames_to_column("gene")

#画图
ggplot(markers, aes(x=Difference, y=avg_log2FC, color = change)) + 
  geom_point(size=1.2) + 
  scale_color_manual('change',labels=c(paste0("down(",table(markers$change)[[1]],')'),
                                       'ns',
                                       paste0("up(",table(markers$change)[[3]],')' )),
                     values=c("#0d1b46", "grey","tomato" ))+
  geom_label_repel(data=subset(markers, 
                               avg_log2FC >= 1 & Difference >= 0.2 & p_val_adj <= 0.05), 
                   aes(label=gene),  #添加label
                   color="black", #设置label中标签的颜色
                   segment.colour = "black",#设置label框的颜色
                   label.padding = 0.2, 
                   segment.size = 0.3,  #框的大小
                   size=4,
                   max.overlaps = 20) +  # 增加 max.overlaps 参数
  geom_label_repel(data=subset(markers, 
                               avg_log2FC <= -1 & Difference <= -0.1 & p_val_adj <= 0.05), 
                   aes(label=gene), 
                   label.padding = 0.2, 
                   color="black",
                   segment.colour = "black",
                   segment.size = 0.3, size=4,
                   max.overlaps = 20) + 
  geom_vline(xintercept = 0,linetype = 2) +
  geom_hline(yintercept = 0,linetype = 2) +
  labs(x="△ Percentage difference",y="Log-Fold change") + 
  theme_bw()
ggsave("volcanoplot_new.pdf",width = 9,height = 7)




#**************************************************************************
#####################pseudobulks单细胞差异基因分析########################## 
##***************************************************
##*Pseudobulk 分析
#除了FindMarkers/FindAllmarkers这种方法以外，pseudobulks是另一种单细胞差异基因分析的方法
#Pseudobulk分析将单细胞RNA测序数据中的细胞按特定的条件（如样本、群体、时间点等）聚合为“伪散装”样本，然后对这些聚合样本进行差异表达分析


#异同点总结
#相同点：
#● 两者都用于识别在不同条件或群体之间存在差异表达的基因。
#● 都需要预处理和标准化单细胞RNA测序数据。
#不同点：
#● 数据处理方式：pseudobulk分析将细胞聚合为群体进行分析，而FindMarkers分析直接在单细胞层面上进行。
#● 解析水平：pseudobulk分析关注的是群体间的差异，FindMarkers分析则关注细胞群体内及群体间的异质性。
#● 工具和方法：pseudobulk分析可以使用传统的bulk RNA-seq分析工具，而FindMarkers/FindAllmarkers通常依赖于专门为单细胞数据设计的统计检验方法。


#适用场景
#● Pseudobulk分析：适用于样本数较多且希望降低单细胞数据噪声的研究，或希望利用传统bulk RNA-seq分析工具进行下游分析的场景。
#● FindMarkers/FindAllmarkers分析：适用于细胞异质性较高的研究，或希望深入探索特定细胞亚群差异的场景。


#分析流程
#1、加载R包及读取数据
rm(list=ls())
library(qs)
library(ggplot2)
library(DESeq2)
library(Seurat)
library(tidyverse)
library(dplyr)

scRNA = readRDS('./Discoid_placenta_Invasive_tro_down1800_inter.rds') 
dim(scRNA)
#[1]  2000 10800
DefaultAssay(scRNA)<-"RNA"
dim(scRNA)
#[1] 11068 10800


#2、数据处理
table(scRNA$species)
# gunpig  human macaca  mouse rabbit    rat 
# 1800   1800   1800   1800   1800   1800 
table(scRNA$sample)
# 
# E140-M        E140-S            GP           GP1    mouse_12.5 
# 905           895           512          1288           702 
# mouse_14.5           PN1           PN2        rabbit Rat_19_5_1_PD 
# 1098           139          1661          1800           459 
# Rat_19_5_2_PD Rat_19_5_3_PD 
# 907           434 

scRNA$sample.id<-rownames(scRNA@meta.data)

#split(colnames(scRNA), scRNA$sample.id)：将所有细胞的列名按 sample.id 进行分组。
#split 函数返回一个列表，其中每个元素包含属于同一 sample.id 的细胞列名。
#bs 是一个列表，列表的每个元素代表一个物种，并包含所有属于该物种的细胞名。
bs = split(colnames(scRNA),scRNA$sample)

ct = do.call(
  cbind,lapply(names(bs), function(x){ 
    # x=names(bs)[[1]]
    kp =colnames(scRNA) %in% bs[[x]]
    rowSums(as.matrix(scRNA@assays$RNA@counts[, kp]  )) #SeuratV5 rowSums(as.matrix(scRNA@assays$RNA@layers$counts[, kp]  ))
  })
)
colnames(ct) <- names(bs)
head(ct)
rownames(ct) <- rownames(scRNA)
#lapply(names(bs), function(x){...}): 对于每一个ID（即 names(bs) 的每一个元素），执行函数体内部的操作。
#kp =colnames(scRNA) %in% bs[[x]]: 这一步通过布尔索引筛选出当前样本ID对应的所有细胞列。
#rowSums(as.matrix(scRNA@assayscounts[, kp])): 对选定的细胞列(不同组)中的基因表达矩阵进行行求和，得到每个基因在该样本中的总表达量。
#这里需要思考一下，我们使用的kp，这里的kp其实代表的是bs中的ID，所以按照这个数据而言，分别是对不同物种组的数据的基因表达矩阵进行行求和。
#最终通过 cbind 函数将所有样本的基因表达总和结果列绑定（即按列组合），生成矩阵 ct，其中每一列对应一个样本，每一行对应一个基因。

#不过此时需要注意的是，Seurat V5  时候ct表格中没有行名，也就是没有基因名，因此我们需要把scRNA的行名加上去。也可以在差异分析前加上去
#
head(ct)
#          gunpig human    macaca mouse rabbit rat
# A4GALT       0   382 2005.0333   411     45 859
# AADACL3      0   301    0.0000     0      4   0
# AADAT        5    10  317.8071     0    106   0
# AAGAB      178   145  257.9245   120   1048  10
# AAK1       415    95  452.6811  1097    840  45
# AAMDC      564  2030    0.0000   166    419 107


#3、获取分组信息
# 获取分组信息
phe = scRNA@meta.data[,c('sample','species')]
phe = unique(scRNA@meta.data[,c('sample','species')]) 
phe[1:5,1:2]
phe
#                                     sample species sample_adjusted
# AAATCCCTGTAAAG-1_1                    PN1   human             PN1
# AAACATTGGCTGTA-1_2                    PN2   human             PN2
# E140-M_AAACCCACAAGATCCT-29         E140-M  macaca          macaca
# E140-S_AAACGCTAGAGGATCC-30         E140-S  macaca          macaca
# Rat_E191_AAACGCTTCGGCTTGG-1 Rat_19_5_1_PD     rat   Rat_19_5_1_PD
# Rat_E192_AAACGAATCCGTATAG-1 Rat_19_5_2_PD     rat   Rat_19_5_2_PD
# Rat_E193_AAACCCAGTGCCTACG-1 Rat_19_5_3_PD     rat   Rat_19_5_3_PD
# AAACCCATCACCGCTT_3             mouse_12.5   mouse      mouse_12.5
# AAACCCAGTCATGCAT_4             mouse_14.5   mouse      mouse_14.5
# CELL57_N3_1                        rabbit  rabbit          rabbit
# CELL20_N1_1_1                          GP  gunpig              GP
# CELL1_N1_1_2                          GP1  gunpig             GP1

group_list <- phe[match(names(bs), phe$sample), 'species']
group_list

#[1] "gunpig" "human"  "macaca" "mouse"  "rabbit" "rat"
#第一行代码从 scRNA 对象的 meta.data 中提取两列数据：sample.id（样本ID）和 tissue.type（组织类型）。
#meta.data 是存储每个细胞对应的元数据信息的表格。提取后的结果 phe 是一个数据框，其中包含每个细胞的样本ID和对应的组织类型。
#第二行代码使用 unique 函数对刚才提取的数据进行去重操作。
#unique 函数会移除数据框中重复的行，因此生成的 phe 数据框会包含每个样本ID唯一对应的一行记录，即每个样本ID对应的组织类型。
#这样处理后，phe 数据框的每一行代表一个样本，而不是一个细胞。

# 接下来的group_list代码是匹配样本ID并提取对应的组织类型：
# names(bs): 这个部分提取的是之前创建的列表 bs 中的样本ID（样本的列名）。
# match (names(bs), phe $ sample.id) : match函数用于在 phe$sample.id 中查找与 names(bs) 相匹配的样本ID，返回匹配到的位置索引。
# 简单来说，它会告诉你每个 bs 列表中的样本ID在 phe 数据框中的位置。
# phe[...]: 这里使用这些位置索引来从 phe 数据框中提取相应行的 tissue.type 列，
# 最终得到的 group_list 是一个向量，包含了 bs 中样本ID对应的物种。


#4、过滤数据
# 赋值并对每一行的
exprSet = ct
dim(exprSet)
exprSet=exprSet[apply(exprSet,1, function(x) sum(x>1) > 1),]
dim(exprSet)  
table(group_list)

#建议group_list设置成因子模式
group_list <- factor(group_list,levels = c("human","macaca","rat","mouse","gunpig","rabbit"))


# apply(exprSet, 1, function(x) sum(x > 1) > 1) 这段代码对 exprSet 矩阵的每一行（即每个基因）进行操作：
# apply(exprSet, 1, function(x) ...)：apply 函数在矩阵的每一行（1 表示行操作）上应用给定的函数。sum(x > 1) > 1：对于每个基因（每行），
# 计算在多少个样本（列）中该基因的表达量大于1，如果该数量大于1（即至少在两个样本中有表达量大于1），则保留该基因。

#5、接下来就是差异分析和绘图，除了DESeq2分析，差异分析方法有很多(deseq2，edgeR，limma_voom)，选择自己适合的，可以参考
#https://mp.weixin.qq.com/s/dCgLjzOGKeSpZlpFMzIR9g
# DESeq2分析
colData <- data.frame(row.names=colnames(exprSet),group_list=group_list)
exprSet_int <- apply(exprSet, 2, as.integer) #这一步不需要哦，但如果报错矩阵不是整数  则需要加上
#rownames(exprSet) <- rownames(scRNA) #如果这步上面没有做，可以再加上去  Seurat V5 运行这个 加上基因行名
rownames(exprSet_int) <- rownames(exprSet)
#
#

#
dds <- DESeqDataSetFromMatrix(countData = exprSet_int,
                              colData = colData,
                              design = ~ group_list)
dds <- DESeq(dds) 
table(group_list)
# group_list
# human macaca    rat  mouse gunpig rabbit 
# 2      2      3      2      2      1

#****************************************
#############选择自己想分析的样本单个分析
#*****************************************
#注意这点，谁在前面相当于谁VS谁，比如group_list)[1]为human，group_list)[2]为macaca, 则得到的差异基因为正，表示human上调，其他与相同
# res <- results(dds, 
#                contrast=c("group_list",
#                           levels(group_list)[1],        
#                           levels(group_list)[2]))
# resOrdered <- res[order(res$padj),]
# head(resOrdered)
# DEG =as.data.frame(resOrdered)
# DEG_deseq2 = na.omit(DEG)
# write.csv(DEG_deseq2,paste0(levels(group_list)[1],"_VS_",levels(group_list)[2],"_DESeq2.csv"))


#****************************************
#######写个循环依次获取 获取human与其他物种的比较
#****************************************
# 注意：这个代码假设你的 group_list 是因子类型，并且包含 "human" 和其他物种名称。
# 获取human与其他物种的比较
for (species in levels(group_list)[-1]) {  # 排除human本身
  # 计算差异分析
  res <- results(dds, 
                 contrast=c("group_list",
                            "human",        
                            species))
  # 排序结果
  resOrdered <- res[order(res$padj),]
  DEG <- as.data.frame(resOrdered)
  DEG_deseq2 <- na.omit(DEG)
  # 保存结果到CSV文件
  output_filename <- paste0("DESeq2_results/human_VS_", species, "_DESeq2.csv")
  write.csv(DEG_deseq2, output_filename)
  print(paste("Results saved for human vs", species))
}


#****************************************
##             作图
#****************************************
#添加上下调信息
DEG_deseq2 <- DEG_deseq2 %>%
  mutate(Type = if_else(padj > 0.05, "stable",
                        if_else(abs(log2FoldChange) < 1, "stable",
                                if_else(log2FoldChange >= 1, "up", "down")))) %>%
  arrange(desc(abs(log2FoldChange))) %>% rownames_to_column("Symbol")
# ggplot绘图
ggplot(DEG_deseq2, aes(log2FoldChange,-log10(padj))) +
  geom_point(size = 3.5, alpha = 0.8,
             aes(color = Type),show.legend = T)  +
  scale_color_manual(values = c("#00468B", "gray", "#E64B35")) +
  ylim(0, 15) +
  xlim(-10, 10) +
  labs(x = "Log2(fold change)", y = "-log10(padj)") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = 'black',lwd=0.8) + 
  geom_vline(xintercept = c(-1, 1), linetype = 2, color = 'black',lwd=0.8)+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave("volcano.pdf",width = 9,height = 7)





