#Marker gene heatmap(Fig1E)
#all marker(Contains top genes)
library(ggplot2)
library("viridis")

markers <- c(
  # 第一部分
  "EGFR", "ITGB6", "PARP1", "KRT18", "XAGE3", "SMAGP", "ADIRF", "MIR205HG", "VGLL1", "CDO1", 
  "FXYD3", "ISYNA1", "LDHB", "SMS", "EFEMP1", "HLA-G", "PAPPA", "TAC3", "AOC1", "PRG2", 
  "EBI3", "QSOX1", "HPGD", "CSH1", "PLAC8", "HTRA1", "ISM2", "CGA", "CYP19A1", "INSL4", 
  "ERVW-1", "ERVV-1", "HOPX", "HSPB8", "PRKCZ", "FAR2", "CSF3R", "SLC13A4", "ADHFE1", "TACC2", 
  "LGALS13", "GCM1", "CLDN3", "CKB", "C1QBP", "TGM2", "UCHL1", "TYMS", "STMN1", "H2AFX", 
  "CKS1B", "KIAA0101", "TBL1X", "RANBP1", "BMP7", "DNMT1", "LAIR2", "ITGA5", "TNFSF10", 
  "CX3CL1", "HS3ST1", "DEFB1", "CD59", "SRI", "MFAP2", "ISM2", "GJA1", "SERPINB6", "MFAP5", 
  # 第二部分
  "PSG7", "CSHL1", "PSG3", "GDF15", "PSG4", "CSH2", "LGALS14", "CRYAB", "HPGD", "UBE2A", 
  "PLAC1", "GH1", "PLAC2", "CDH1", "MKI67", "KRT8","LGALS2", "LIFR", "HSPA8", "TAF1D", "TRRAP", "STK26", "ADGRG6", "TFPI", "MMP11",  
  "KRT7", "PRLRP2", "LEP", "NOTUM", "CTNNB1", "LOC106027897", "FST", "LOC100724235", "LOC100725464", 
  "GPM6A", "UTRN", "DDX3X", "TBX3", "LAMA1", "TEAD3", "GATA3", "LAMA5", "SAT1", "PLCB4", 
  "KCNMA1", "TSC22D1", "FAP", "NEB", "DST", "ALDH1A2", "SORBS2", "PRL7B1", "PRL5A1", "PRL4A1", 
  "PRL2A1", "PRL5A2", "PRL2C1",  "MMP15", "TPBPA", "GALNT6", "PRL7B1", "PRL2C1", 
  "LOC684107", "HSD17B2", "C1QTNF6", "NCAM1", "PLAC8", "TIAM2", "RORA", "NUAK1", "GNA14", 
  "ALDH1A3", "TIAM1", "TCF4", "ALDH1A3", "ANO6", "PHACTR1", "OLR1", "RGS12", "LEPR", "NOS1AP", 
  # 第三部分
  "CTSQ", "PODXL", "APBA1", "PRKCE", "ST3GAL4", "CTDSPL", "CDH23", "PDE10A", "PRKCH", "GLIS1", 
  "STRA6", "ROR2", "LGR5", "DAB1", "CASKIN1", "PLEKHA6", "KANK1", "TACC2", "IQGAP2", "SNAP91", 
  "PRL8A9", "MITF", "SLCO2A1", "1700025G04RIK", "ST6GAL1", "CEACAM14", "ARHGAP39", "MGAT4A", 
  "CEACAM5", "SRGAP1", "GM42418", "MICAL3", "TNFRSF9", "FLT1", "GATA3", "LGALS3", "GATA2", 
  "LOC100358868", "SRD5A2", "RALGPS2", "MAP9", "JMJD1C", "PTPRK", "PLAU", "LOC103348889", 
  "MMP2", "MMP14", "DIO2", "SERPINE1", "HSPG2", "RAC1", "NGFR", "FBLN2", "GJD4", "IL21R", 
  "SEMA7A", "NTRK1", "CRISPLD2", "LRP8", "TFRC", "DUSP9", "TP63", "ANGPT4", "TFEC", "SAT1", 
  # 第四部分
  "GRB10", "LOC108178888", "FOXRED1", "LOC100346323", "SNAP91", "STON2", "PHF8", "FAM13A", 
  "PLET1", "PPAG3", "PAG6", "ERBB3", "PLBD1", "LOC110260671", "LOC110259250", "LOC100524786", 
  "LOC100524130", "SOX5", "LRP2", "SMPDL3A", "DAB2", "LGMN", "EPB41L3", "DEFB1", "LOC100511891", 
  "CTSD", "PAPSS2", "CHST8", "SLC45A3", "CHST3", "ST3GAL5", "PECR", "ST3GAL3", "FGL2", "SLC15A1", 
  "PAG10", "PAG3", "PRP3", "PAG5", "PAG4", "PAG16", "PRP8", "PAG20", "PAG15", "PAG18", "PRP6", 
  "PRP14", "PRP1", "PRP4", "PRP3", "PAG2", "PAG8", "TEAD1", "LOC101907852", "LOC112448071", 
  "B4GALNT3", "LOC508098", "PLEKHD1", "CTNND2", "HSD3B1", "CYP11A1", "PHLDA2", "RAI14", "PAG6", 
  "CITED1", "PAG-6", "TGFBI", "CPXM1", "FIBIN", "ANGPTL2", "PPAG3", "GATM", "LGALS16", "LGALS15", 
  "PAGE4", "PAG-8", "ERV3-1", "PEG3", "LOC102171766", "LOC102177258", "GULP1", "CCL17", "SRSF7", 
  "FMR1", "TFAP2A", "PEG10", "CCNE1", "ITGA6", "KRT19", "LOC477365", "LOC102153162", "LOC611209", 
  "PAPPA2", "ADM", "SERPINE2", "MMP1", "MMP12", "ANGPT2", "C35H6orf62", "ESM1", "VASN", "OLR1", 
  "RAB27B", "LOC482182", "RRAS2", "MFSD2A", "PTHLH", "CAR1", "SMPD1", "CEMIP", "JAG2", "LDLR", 
  "GPT2", "PGF", "EGFR", "INHBA", "NPPC", "PECAM1", "VWF", "EGFL7", "KDR", "CD93", "FCGR2B", 
  "RAMP2", "PLVAP", "LDB2", "COL3A1", "COL1A1", "LUM", "COL1A2", "DCN", "COL15A1", "CLU", "POSTN", 
  "PINLYP", "GUCA2B", "OXT", "MUC1", "EPCAM", "MUC4", "CDH16", "PDPN", "CTSL", "C1QC", "LYVE1", 
  "CTSC", "CD68", "CXCL14", "CD52", "CD83", "CD86", "GZMA", "CD69", "CD52", "CD3G", "PTPRC", 
  "CD3E", "TBXAS1", "ITGA2B", "CD226", "PRKCQ", "CD163", "CD14", "CD209", "C1QA", "C1QB", "AIF1", 
  "MRC1", "NKG7", "CD3D", "GZMA", "GZMB", "GZMK", "TMEM163", "RUNX3", "DOCK2", "DOCK10", "DSCAM", 
  "FAM178B", "FYB", "S100A9", "S100A8", "IFIT1BL", "G0S2", "IL1B", "LRG1", "LYZ2", "CCL6", 
  "CENPF", "TOP2A", "SMC4", "HMGB2", "CTSC", "B2M", "SPP1"
)

markers1 <- as.data.frame(markers)
markerdata <- ScaleData(tro, features = as.character(unique(markers1$markers)), assay = "RNA")

markerdata$celltype_merge <- factor(x=markerdata$celltype_merge,
                                    levels = c("Cytotrophoblast(human)","Extravillous_trophoblast(human)","Syncytiotrophoblast(human)", "CTB(macaca)", "EVT(macaca)", 'STB(macaca)',
                                               "CTB(gp)", "EVT(gp)",'STB(gp)',"Invasive(rat)","GC(mouse)", "S-TGC(mouse)","SynT(mouse)",'SpT(mouse)',
                                               "CTB(rabbit)", "EVT(rabbit)", "STB(rabbit)","pUNC1", "pUNC2","pUNC3", "pUNC4", "pUNC5",
                                               "BNC1(cow)", "BNC2(cow)", "BNC3(cow)","UNC1(cow)", "UNC2(cow)", "UNC3(cow)", "UNC4(cow)",
                                               "BNC1(goat)", "BNC2(goat)", "BNC3(goat)", "BNC4(goat)", "BNC5(goat)", "BNC6(goat)","UNC(goat)",
                                               "CTB(dog)", "EVT(dog)", "STB(dog)", 
                                               "Endothelial_cell(human)", "Endo(macaca)","Endo(gp)",  "Endo(rat)",'Endo(mouse)',"Endo(rabbit)","Endo(pig)",  "Endo(cow)","Endo(goat)","Endo(dog)",
                                               "Stromal_cell(human)",'Stro(macaca)', 'Stro(gp)',"Stro(rat)","Stro(mouse)","Mes(mouse)", "Stro(rabbit)","Stro(pig)", "Stro(cow)", "Stro(goat)","Stro(dog)",
                                               'Epi(rat)','Epi(rabbit)','Epi(pig)','Epi(cow)','Epi(goat)','Epi(dog)',
                                               "Dendritic_cell(human)","T_cell(human)","Tcells(macaca)","Tcells(rat)",'Tcells(mouse)',
                                               'Macrophage(Hofbauer cell)(human)',"Mac(macaca)", "Mac(gp)", "Mac(rat)",'Mac(mouse)',"Mac(rabbit)",'Mac(goat)','Mac(cow)',"Mac(dog)",
                                               "NKcells(macaca)", "NKcells(rat)", "NKcells(mouse)",
                                               "Bcells(rat)", 'Neu(goat)'))


cols <- viridis(100)[c(20, 50, 100)]
pdf(file = "heatmap_output_FINAL.pdf", width = 10, height =25)  # Adjust width and height as needed

# Generate the heatmap and save it to the PDF
DoHeatmap(markerdata,
          features = as.character(unique(markers1$markers)),
          group.by = "celltype_merge",
          assay = 'RNA',draw.lines = FALSE,
          group.colors = c("#a6cee3","#94c2dc","#82b7d6","#257cb2",'#398aac','#4c97a7',
                           '#a5d981','#95d074','#84c868','#42a737','#eebc6a','#fdbb68','#fdb259','#fdaa4b',
                           '#fe8524','#fd8838','#fd8c4c','#f16464','#ed5354','#ea4344','#e73233','#e42123',
                           '#dd3841','#da4c59','#d76072','#d4748a','#d088a3','#cd9cbb','#cab0d3',
                           '#8c66af','#7f57a7','#73489f','#6f4399','#825d99','#967699','#a99099',
                           '#eedb80','#e4c571','#d9af63',
                           '#5fa0ca','#5fa5a1','#74c05c','#35a02c', '#fe992e','#fc9060','#f78685','#bea4ce','#bda999','#cf9a54', #Endo
                           '#4d95c4','#72b29c','#63b84f','#4fa435','#fda23c','#fb9374','#f47575', '#b294c7','#f7f599','#bb6e36', #Stro
                           '#6aa83e','#fb9788','#e12429','#a585bf','#b15928','#d6bcc0','#b15928',#Epi
                           '#70acd0','#3b8abe','#99cd91','#d4b861','#fe8002',#T
                           '#2a7fb7','#acdb8b','#52af43','#9fb04f','#fe8110','#fa9696','#9876b7','#bb7784','#f8f18f','#c58445',#,mac
                           '#86c096','#b9b458','#b9b458','#84ac47','#d0c299' )) +
  scale_fill_gradientn(colors = c("#355E8DFF", "#218F8DFF", "#FDE725FF"))+theme(
    legend.position = "none")
#


# Close the PDF device
dev.off()


##########heatmap或者ComplexHeatmap做图。
data<-GetAssayData(markerdata,slot="scale.data")
#整理亚群注释信息
celltype_info<-sort(markerdata$celltype_merge)
table(celltype_info)

celltype_info<-factor(celltype_info,levels=c("Cytotrophoblast(human)","Extravillous_trophoblast(human)","Syncytiotrophoblast(human)", "CTB(macaca)", "EVT(macaca)", 'STB(macaca)',
                                             "CTB(gp)", "EVT(gp)",'STB(gp)',"Invasive(rat)","GC(mouse)", "S-TGC(mouse)","SynT(mouse)",'SpT(mouse)',
                                             "CTB(rabbit)", "EVT(rabbit)", "STB(rabbit)","pUNC1", "pUNC2","pUNC3", "pUNC4", "pUNC5",
                                             "BNC1(cow)", "BNC2(cow)", "BNC3(cow)","UNC1(cow)", "UNC2(cow)", "UNC3(cow)", "UNC4(cow)",
                                             "BNC1(goat)", "BNC2(goat)", "BNC3(goat)", "BNC4(goat)", "BNC5(goat)", "BNC6(goat)","UNC(goat)",
                                             "CTB(dog)", "EVT(dog)", "STB(dog)", 
                                             "Endothelial_cell(human)", "Endo(macaca)","Endo(gp)",  "Endo(rat)",'Endo(mouse)',"Endo(rabbit)","Endo(pig)",  "Endo(cow)","Endo(goat)","Endo(dog)",
                                             "Stromal_cell(human)",'Stro(macaca)', 'Stro(gp)',"Stro(rat)","Stro(mouse)","Mes(mouse)", "Stro(rabbit)","Stro(pig)", "Stro(cow)", "Stro(goat)","Stro(dog)",
                                             'Epi(rat)','Epi(rabbit)','Epi(pig)','Epi(cow)','Epi(goat)','Epi(dog)',
                                             "Dendritic_cell(human)","T_cell(human)","Tcells(macaca)","Tcells(rat)",'Tcells(mouse)',
                                             'Macrophage(Hofbauer cell)(human)',"Mac(macaca)", "Mac(gp)", "Mac(rat)",'Mac(mouse)',"Mac(rabbit)",'Mac(goat)','Mac(cow)',"Mac(dog)",
                                             "NKcells(macaca)", "NKcells(rat)", "NKcells(mouse)",
                                             "Bcells(rat)", 'Neu(goat)'))

#将data中由变量g指定和行和celltype_info中的细胞类型对应的列转换为一个矩阵，并存储在mat变量中
mat<-as.matrix(data[markers,names(celltype_info)])


#设置列标签的颜色
#将列标签的col设置为和umap一样的配色，并且与亚群的数量统一起来
#使用HeatmapAnnotation整理需要添加的列标签信息
#在Heatmap中使用top_annotation添加上美化后的列标签信息
#设置列标签颜色，并整理top_anno标签注释信息
col<-c("#a6cee3","#94c2dc","#82b7d6","#257cb2",'#398aac','#4c97a7',
       '#a5d981','#95d074','#84c868','#42a737','#eebc6a','#fdbb68','#fdb259','#fdaa4b',
       '#fe8524','#fd8838','#fd8c4c','#f16464','#ed5354','#ea4344','#e73233','#e42123',
       '#dd3841','#da4c59','#d76072','#d4748a','#d088a3','#cd9cbb','#cab0d3',
       '#8c66af','#7f57a7','#73489f','#6f4399','#825d99','#967699','#a99099',
       '#eedb80','#e4c571','#d9af63',
       '#5fa0ca','#5fa5a1','#74c05c','#35a02c', '#fe992e','#fc9060','#f78685','#bea4ce','#bda999','#cf9a54', #Endo
       '#4d95c4','#72b29c','#63b84f','#4fa435','#fda23c','#fb9374','#f47575', '#b294c7','#f7f599','#bb6e36', #Stro
       '#6aa83e','#fb9788','#e12429','#a585bf','#b15928','#d6bcc0','#b15928',#Epi
       '#70acd0','#3b8abe','#99cd91','#d4b861','#fe8002',#T
       '#2a7fb7','#acdb8b','#52af43','#9fb04f','#fe8110','#fa9696','#9876b7','#bb7784','#f8f18f','#c58445',#,mac
       '#86c096','#b9b458','#b9b458','#84ac47','#d0c299' )

names(col)<-levels(celltype_info)


top_anno<-HeatmapAnnotation(
  cluster=anno_block(gp=gpar(fill=col),#设置填充色
                       labels=levels(celltype_info),
                       labels_gp=gpar(cex=0.5,col="white")))#设置字体

#只展示特定基因

gene <- c(
  'EGFR', 'PARP1', 'HLA-G', 'PAPPA', 'CGA', 'INSL4', 'ERVV-1', 'GCM1',
  'C1QBP', 'LAIR2', 'ITGA5', 'PSG7', 'PSG3', 'CSH2', 'CDH1', 'MKI67',
  'KRT8', 'MMP11', 'LGALS2', 'PRLRP2', 'NOTUM', 'TBX3', 'TEAD3', 'GATA3',
  'PRL7B1', 'PRL2A1', 'NCAM1', 'PLAC8', 'NOS1AP', 'PODXL', 'STRA6', 'ROR2',
  'PRL8A9', 'MITF', 'FLT1', 'LGALS3', 'GATA2', 'DIO2', 'RAC1', 'TFRC', 'TP63',
  'PLET1', 'PAG6', 'LRP2', 'PAG3', 'PAG4', 'PAG2', 'PAG8', 'PAG6', 'CITED1',
  'GATM', 'LGALS16', 'LGALS15', 'PAGE4', 'PAG-8', 'ERV3-1', 'TFAP2A', 'PEG10',
  'PAPPA2', 'SERPINE2', 'MFSD2A', 'PTHLH', 'CAR1', 'PECAM1', 'VWF', 'COL1A1',
  'LUM', 'DCN', 'MUC1', 'EPCAM', 'CD52', 'CD83', 'GZMA', 'PTPRC', 'CD163',
  'MRC1', 'NKG7', 'GZMB', 'S100A9', 'S100A8', 'LRG1', 'CENPF', 'SPP1'
)

# 找到 gene 列表中不在 mat_scaled 行名中的基因
missing_genes <- gene[!gene %in% rownames(mat)]

gene_pos<-which(rownames(mat)%in%gene)

row_anno<-rowAnnotation(gene=anno_mark(at=gene_pos,labels=gene))




pdf(file = "heatmap_output_fixed.pdf", width = 15, height = 25)

# 绘制热图并保存
draw(
  Heatmap(
    mat,
    name = "Expression",
    row_order = 1:nrow(mat),
    col = colorRamp2(c(-1, 0, 2), c("#355E8DFF", "#218F8DFF", "#FDE725FF")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    border = FALSE,  # 去除边框，避免间隔线
    row_names_gp = gpar(fontface = 'italic', fontsize = 4),
   column_split = celltype_info,
   top_annotation = top_anno,  # 热图上边增加注释
    right_annotation = row_anno,
    column_gap = unit(0, "mm"),  # 移除列之间的间隔
  )
)

# 关闭 PDF 设备
dev.off()

# 提示完成保存
cat("Heatmap has been saved to 'heatmap_output_fixed.pdf'\n")
