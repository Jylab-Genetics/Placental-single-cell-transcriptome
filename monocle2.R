library(dplyr)
library(Seurat)
library(dplyr)
library(Matrix)
library(patchwork)
library(monocle)
library(tidyverse)
library(ggplot2)
library(viridis)
#library(RColor)
library(reshape2)
##确定物种
species<-"cow"
#
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

#cow
#BNC1,BNC2, BNC3,  UNC1, UNC2, UNC3,UNC4

#c("#7B68EE","#7B68EE","#20B2AA", "#FFA500","#FFA500","#87CEEB","#FFA500")

#导入cds对象
cds<-readRDS("./dog_monocleFindAllMarkers_cds.RDS")
colnames(pData(cds))
plot_cell_trajectory(cds)+ scale_color_manual(values = colour)
plot_cell_trajectory(cds,color_by = "subclass3",show_branch_points = T)+ scale_color_manual(values = colour)
plot_cell_trajectory(cds)+ scale_color_manual(values = colour)
plot_cell_trajectory(cds, color_by = "Pseudotime") 
p1=plot_cell_trajectory(cds,color_by = "subclass3",show_branch_points = T)+ scale_color_manual(values = colour)
p4=plot_cell_trajectory(cds, color_by = "Pseudotime") +scale_color_viridis(10,option = "H")
#cow
#p1=plot_cell_trajectory(cds,color_by = "subclass3",show_branch_points = T)+ scale_color_manual(values = c("#5F9EA0","#5F9EA0","#DC143C", "#4682B4","#4682B4","#87CEEB","#4682B4"))


p2=plot_cell_trajectory(cds, color_by = "State") 

data_df<-t(reducedDimS(cds)) %>% 
  as.data.frame() %>%
  select_('Component 1' = 1, 'Component 2' = 2) %>%
  rownames_to_column("Cells") %>%mutate(pData(cds)$State,
                                        pData(cds)$Pseudotime,
                                        pData(cds)$sample,
                                        pData(cds)$subclass3)

colnames(data_df) <- c("cells","Component_1","Component_2","State",
                       "Pseudotime","orig.ident","subclass3")
reduced_dim_coords <- reducedDimK(cds)
ca_space_df <- Matrix::t(reduced_dim_coords) %>% 
  as.data.frame() %>% 
  select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>% 
  mutate(sample_name = rownames(.), sample_state = rownames(.)) 

dp_mst <- minSpanningTree(cds)
edge_df <- dp_mst %>% 
  igraph::as_data_frame() %>%
  select_(source = "from", target = "to") %>%
  left_join(ca_space_df %>% select_(source="sample_name",source_prin_graph_dim_1="prin_graph_dim_1",source_prin_graph_dim_2="prin_graph_dim_2"), 
            by = "source") %>%
  left_join(ca_space_df %>% select_(target="sample_name",target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"),
            by = "target")
p3=ggplot() +geom_point(data = data_df, aes(x = Component_1, y = Component_2,color =Pseudotime)) +scale_color_viridis()+theme_bw()+
  theme(panel.background = element_blank(),panel.border = element_blank(),panel.grid = element_blank(),axis.ticks.length = unit(0.8, "lines"),
        axis.ticks = element_blank(),axis.line = element_blank(),axis.title = element_text(size=15),)+
  
  
  geom_segment(aes_string(x="source_prin_graph_dim_1",y="source_prin_graph_dim_2",
                          xend="target_prin_graph_dim_1",yend="target_prin_graph_dim_2"),
               size=0.75, linetype="solid", na.rm=TRUE, data=edge_df)

pdf(paste0(species,'_plot_cell_trajectory.pdf'),width=20,height=5)
p1|p2|p3|p4
dev.off()

#===============================================================================
#                           基因随拟时表达变化
#===============================================================================
colnames(pData(cds))#cds是已经做好拟时的monocle2对象
#将需要展示的基因表达量添加到cds， 这里我使用的是log2处理，也可以使用log10进行标准化
genes <- rownames(cds)
genes_exp <- list()

for(i in 1:length(genes)){
  A <- log2(exprs(cds)[genes[i],]+1)
  A <- as.data.frame(A)
  genes_exp[[i]] <- A
}

gene_exp <- do.call(cbind, genes_exp)
colnames(gene_exp) <- genes
#将上述几个基因的拟时表达添加到monocle
pData(cds) = cbind(pData(cds), gene_exp)
#提取作图数据，只需要游基因表达和拟时即可
data <- pData(cds)
data[1:4,1:24]
# 提取指定列
selected_columns <- data[, c("subclass3", "Pseudotime", "State")]
# 提取第 24 列及其后的所有列  提取所有基因表达列

#***************************************注意此处根据自己数据结构修改24:ncol(data)******************************************************************
additional_columns <- data[, 24:ncol(data)] #**********************************  

# 合并提取的列，得到基因随时序的表达量
result <- cbind(selected_columns, additional_columns)
# 查看结果
result[1:4,1:10]  

write.csv(result,paste0(species, '_Pseudotime_exp.csv'))

#应该先根据state选择属于分化线路的细胞，然后根据细胞类型去除因数据造成的少量偏差细胞
##提取CTB-EVT分化线路细胞和基因,根据时序图确定State
#******************************************************************************
#*#****************************************************************************
CTB2EVT_state <- result[result$State %in% c("2", "1"), ]
#筛选细胞
CTB2EVT <- CTB2EVT_state[CTB2EVT_state$subclass3 %in% c("cow_CTB", "cow_EVT"), ]
CTB2EVT[1:4,1:10]
write.csv(CTB2EVT,paste0(species,'CTB2EVTPseudotime_exp.csv'))
#提取CTB-STB分化线路细胞和基因
CTB2STB_state <- result[result$State %in% c("2", "3"), ]
CTB2STB <- CTB2STB_state[CTB2STB_state$subclass3 %in% c("Dog_CTB", "Dog_STB"), ]
CTB2STB[1:4,1:10]
write.csv(CTB2STB,paste0(species,"_CTB2STBPseudotime_exp.csv"))

###########然后每个阶段分化的细胞山脊图
library(monocle)
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)

#test=readRDS("test_monocle.rds")
#plotdf=pData(test)

ggplot(CTB2STB, aes(x=Pseudotime,y=subclass3,fill=subclass3))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(10,24),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )
ggsave(paste0(species,"_CTB2STB_density.pdf"))

ggplot(CTB2STB, aes(x=Pseudotime,y=subclass3,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(10,24),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )
ggsave(paste0(species,"_CTB2STB_PseudotimeDensity.pdf"))

ggplot(pdata, aes(x = Pseudotime, y = subclass2, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 1) +
  geom_vline(xintercept = c(10, 24), linetype = 2) +
  scale_fill_gradientn(name = "Pseudotime", colors = viridis::viridis(10, option = "H")) +
  scale_y_discrete("") +
  theme_minimal() +
  theme(panel.grid = element_blank())

#*******************************************************
##############提取伪细胞坐标，自己作图
#*******************************************************

##########然后画每个分支的分化线路图
plotdf2=as.data.frame(t(cds@reducedDimS))
colnames(plotdf2)=c("component1","component2")
plotdf2$Pseudotime=cds$Pseudotime
head(plotdf2)
#             component1 component2 Pseudotime
# CELL18_N1_1   2.997840   4.634530   10.40013
# CELL20_N2_1  -6.099351  -2.041395   13.04920
# CELL27_N2_1  -7.902289  -1.659195   14.86363
# CELL29_N1_1  -7.942900  -1.620100   14.90600
# CELL30_N3_1   2.208798   7.529382   13.35936
# CELL46_N1_1   2.856301   7.171360   12.98451
write.csv(plotdf2,paste0(species,"_Pseudotime_coordinate.csv"))
max_pseudotime <- max(plotdf2$Pseudotime)
max_pseudotime
phase1 = 5
phase2 = 10
phase3 = 15.37778
plotdf2$phase = "phase2"
plotdf2[plotdf2$Pseudotime<=5,"phase"] = "phase1"
plotdf2[plotdf2$Pseudotime>5 & plotdf2$Pseudotime<10,"phase"] = "phase2"
plotdf2[plotdf2$Pseudotime>=10,"phase"] = "phase3"
#######
plotdf2%>%ggplot(aes(component1,component2,color=phase))+
  geom_point()+
  theme_minimal() +scale_color_manual(values = colour)+  
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
ggsave(paste0(species,"_phase.pdf"))


#取与上面提取单个分支的相同细胞，如CTB_TO_EVT
selected_plotdf <- plotdf2[rownames(CTB2EVT), ]
selected_plotdf%>%ggplot(aes(component1,component2,color=Pseudotime))+
  geom_point()+
  theme_minimal()
ggsave(paste0(species,"__CTB2STB_Pseudotime_branch.pdf"))

