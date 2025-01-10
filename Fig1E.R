merge$celltype_merge <- factor(x = merge$celltype_merge,
                               levels = c("CTB(human)",
                                          "EVT(human)",
                                          "STB(human)",
                                          "Stro(human)",
                                          "Endo(human)",
                                          "Mac(human)",
                                          "DC(human)",
                                          "Tcell(human)",
                                          "CTB(macaca)",
                                          "EVT(macaca)",
                                          "STB(macaca)",
                                          "Stro(macaca)",
                                          "Endo(macaca)",
                                          "Mac(macaca)",
                                          "NKcells(macaca)",
                                          "Tcells(macaca)",
                                          "CTB(guinea pig)",
                                          "EVT(guinea pig)",
                                          "STB(guinea pig)",
                                          "Stro(guinea pig)",
                                          "Endo(guinea pig)",
                                          "Epi(guinea pig)",
                                          "Mono(guinea pig)",
                                          "Mac(guinea pig)",
                                          
                                          "Invasive(rat)",
                                          "Stro(rat)",
                                          "Endo(rat)",
                                          "Epi(rat)",
                                          "Mac(rat)",
                                          "Tcells(rat)",
                                          "NKcells(rat)",
                                          "Bcells(rat)",
                                          
                                          "SpT(mouse)",
                                          "GC(mouse)",
                                          "S-TGC(mouse)",
                                          "SynTI(mouse)",
                                          "SynTII(mouse)",
                                          "Stro(mouse)",
                                          "Endo(mouse)",
                                          "Mac(mouse)",
                                          "NKcells(mouse)",
                                          "Tcells(mouse)",
                                          
                                          "CTB(rabbit)",
                                          "EVT(rabbit)",
                                          "STB(rabbit)",
                                          "Stro(rabbit)",
                                          "Endo(rabbit)",
                                          "Epi(rabbit)",
                                          "Mac(rabbit)",
                                          
                                          "UNC1(pig)",
                                          "UNC2(pig)",
                                          "UNC3(pig)",
                                          "UNC4(pig)",
                                          "UNC5(pig)",
                                          "Stro(pig)",
                                          "Endo(pig)",
                                          "Epi(pig)",
                                          "Mac(pig)",
                                         
                                          "UNC1(cow)",
                                          "UNC2(cow)",
                                          "UNC3(cow)",
                                          "UNC4(cow)",
                                          "BNC1(cow)",
                                          "BNC2(cow)",
                                          "BNC3(cow)",
                                          "Stro(cow)",
                                          "Endo(cow)",
                                          'Epi(cow)',
                                          "Mac(cow)",
                                          "Leu(cow)",
                                          "Neu(cow)",
                                         
                                          "UNC1(goat)",
                                          "UNC2(goat)",
                                          "BNC1(goat)",
                                          "BNC2(goat)",
                                          "BNC3(goat)",
                                          "Stro(goat)",
                                          "Endo(goat)",
                                          "Epi(goat)",
                                          "Mac(goat)",
                                          "Neu(goat)",
                                          
                                          
                                          "CTB(dog)",
                                          "EVT(dog)",
                                          "STB(dog)",
                                          "Stro(dog)",
                                          "Endo(dog)",
                                          "Epi(dog)",
                                          "Mac(dog)",
                                          "Tcells(dog)"))
DimPlot(merge,group.by="celltype_merge",reduction="tsne",theme_use = "theme_blank",pt.size=0.02,pt.alpha=0.5,raster=FALSE)

markers <- c(
  "PEG10", "PRDX5", "XAGE3", "SMAGP", "ADIRF", "MIR205HG", "VGLL1", "CDO1", "FXYD3", "ISYNA1",
  "LDHB", "SMS", "EFEMP1", "HLA-G", "PAPPA", "TAC3", "AOC1", "PRG2", "EBI3", "QSOX1",
  "HPGD", "CSH1", "PLAC8", "HTRA1", "ISM2", "CGA", "CYP19A1", "INSL4", "ERVW-1", "ERVV-1",
  "HOPX", "HSPB8", "PRKCZ", "FAR2", "CSF3R", "SLC13A4", "ADHFE1", "TACC2", "LGALS13", "GCM1",
  "CLDN3", "CKB", "C1QBP", "TGM2", "UCHL1", "TYMS", "STMN1", "H2AFX", "CKS1B", "KIAA0101",
  "TBL1X", "RANBP1", "BMP7", "DNMT1", "LAIR2", "ITGA5", "TNFSF10", "CX3CL1", "HS3ST1", "DEFB1",
  "CD59", "SRI", "MFAP2", "ISM2", "GJA1", "SERPINB6", "MFAP5", "PSG7", "CSHL1", "PSG3",
  "GDF15", "PSG4", "CSH2", "LGALS14", "CRYAB", "HPGD", "UBE2A", "PLAC1", "GH1", "PLAC2",
  "MKI67", "KRT8", "CDH1","TAF1D", "TRRAP", "LGALS2", "PHLDA2", "YBX1", "MIF", "TPX2", "PARP1",
  "LEP", "RPS19", "RPS3A", "RPS7", "MMP11", "PHLDB2", "TSPAN5", "GPM6A", "FST", "NOTUM",
  "VEGFA", "PRLRP2", "LIFR", "TPK1", "EPB41L1", "VEGFA", "NXPE2", "TBX3", "TFAP2C", "LAMA1",
  "TEAD3", "GATA3", "TSC22D1", "FAP", "NEB", "KCNMA1", "ATL1", "ALDH1A2", "PLCB4", "ERGIC1",
  "TFPI", "LAMA5", "PRL7B1", "PRL5A1", "PRL4A1", "PRL2A1", "PRL5A2", "PRL2C1","MMP15",
  "TPBPA", "GALNT6", "PRL7B1", "PRL2C1", "LOC684107", "HSD17B2", "C1QTNF6", "NCAM1", "PLAC8",
  "TIAM2", "RORA", "NUAK1", "GNA14", "ALDH1A3", "TIAM1", "TCF4", "ALDH1A3", "ANO6", "PHACTR1",
  "OLR1", "RGS12", "LEPR", "NOS1AP", "CTSQ", "PODXL", "APBA1", "PRKCE", "ST3GAL4", "CTDSPL",
  "CDH23", "PDE10A", "PRKCH", "GLIS1", "DAB1", "CASKIN1", "PLEKHA6", "STRA6", "KANK1", "TACC2",
  "IQGAP2", "SNAP91", "RBM20", "RALGPS1", "MPP7", "ARL15", "PSMD9", "TGFA", "ROR2",  
  "ZNRF3", "CMSS1", "TACO1", "GM47283", "PRL8A9", "MITF", "SLCO2A1",
  "1700025G04RIK", "ST6GAL1", "CEACAM14", "ARHGAP39", "MGAT4A", "CEACAM5", "SRGAP1", "GM42418",
  "MICAL3", "TNFRSF9", "GATA2", "CD9", "GRHL1", "PTHLH", "ELF3", "TFAP2A", "TGFB1", "PLAU",
  "WSB1", "RALGPS2", "ANXA1", "MMP2", "MMP14", "DIO2", "EFNB1", "SERPINE1", "HSPG2", "RAC1",
  "VCAN", "TWSG1", "AQP3", "ADA", "FBLN2", "ALDH1A3", "GJA1", "LRP8", "TFRC", "DUSP9", "TP63",
  "SRSF7", "SAT1", "RLN1", "CLK1", "RESF1", "SRSF2", "RSRP1", "PJA1", "ZRANB2", "AKR1B10", "LGALS1",
  "FLT1", "ANO4", "TCF7L2", "SOX5", "GRN", "PLET1",
  "PLBD1", "CITED1", "C1QTNF4",  "FOLR1",
  "CTBS", "PPAG3", "FOLR1", "PAG6", "CITED1",
  "DAB2", "LRP2", "KRT7",
  "CHST8", "PAPSS2", "CHST3", "ST3GAL3", "SLC45A3", "ST3GAL5", "SLC15A1",
  "PAG10", "PAG3",
  "PRP3", "PAG5", "PAG4", "PAG16", "PRP8", "PAG20", "PAG15", "PAG18", "PRP14", "PRP1", "PRP4", "PRP3",
  "PAG2", "PAG8", "TEAD1", "LOC101907852", "LOC112448071", "B4GALNT3", "LOC508098", "PLEKHD1",
  "CTNND2", "HSD3B1", "CYP11A1", "RAI14", "PAG-3", "PAG-11", "PRP6", "PRP1", "PAG-6", "ADA", "ANGPT1",
  "WIPF1", "GULP1", "PEG3", "CYP19A1", "PAG-8", "ERV3-1", "PAGE4", "CCL17", 
  "SOLD1", "FGD4", "EPS8", "ARL15", "ECI2", "ITGB6", "ITGA6", "KRT19", "CCNE1", "RBFOX1", "GPM6A",
  "FMR1", "AGMAT", "WSB1", "LAMA1", "PAPPA2", "MMP12", "CITED2", "ADM", "ARAP2", "EXT1", "ST14",
  "NPY", "NPC1", "FHDC1", "RRAS2", "LUZP1", "PPP2R2B", "MFSD2A", "PTHLH", "CAR1", "EGFR", "JAG2",
  "ADAM18", "PGF", "SFRP2", "NPPC", "LRP2", "TFPI2", "SCD", "CEMIP", "PECAM1", "VWF", "EGFL7", "KDR",
  "FCGR2B", "RAMP2", "PLVAP", "LDB2", "COL3A1", "COL1A1", "LUM", "COL1A2", "DCN", "COL15A1",
  "CLU", "POSTN", "COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A3", 'MMRN2','WDR72',"GUCA2B",  "MUC1",
  "EPCAM",'WFDC2', "MUC4", "CDH16", "PDPN", "CTSL", "C1QC", "LYVE1", "CTSC", "CD68", "CXCL14", "CD52", "CD83",
  "CD86", "GZMA", "CD69", "CD52", "CD3G", "CXCL8", "CD4", "PTPRC", "CD3E", "TBXAS1", "ITGA2B", "CD226",
  "PRKCQ",'CXCL8','CD4', "CD163", 'CD74','RBPJ',"CD14", "CD209", "C1QA", "C1QB", "AIF1", "MRC1", "NKG7", "CD3D", "GZMA", "GZMB",
  "GZMK", "TMEM163", "RUNX3", "DOCK2", "DOCK10", "DSCAM", "FAM178B", "FYB",  "IFIT1BL",
  "G0S2", "IL1B", "LRG1", "LYZ2", "CCL6", "CENPF", "TOP2A", "SMC4", "HMGB2", "CTSC", "B2M", "SPP1",
  "S100A8",  "S100A9", "S100A12","CXCL8", "CXCL2")





markers <- as.data.frame(markers)
markerdata <- ScaleData(merge, features = as.character(unique(markers$markers)), assay = "RNA")




markerdata$celltype_merge <- factor(x=markerdata$celltype_merge,
                                    levels = c("CTB(human)","EVT(human)","STB(human)", "CTB(macaca)", "EVT(macaca)", 'STB(macaca)',
                                               "CTB(guinea pig)", "EVT(guinea pig)",'STB(guinea pig)',"Invasive(rat)","GC(mouse)", "S-TGC(mouse)","SynTI(mouse)","SynTII(mouse)",'SpT(mouse)',
                                               "CTB(rabbit)", "EVT(rabbit)", "STB(rabbit)", "UNC1(pig)", "UNC2(pig)","UNC3(pig)","UNC4(pig)", "UNC5(pig)","BNC1(cow)", "BNC2(cow)", "BNC3(cow)","UNC1(cow)", "UNC2(cow)", "UNC3(cow)", "UNC4(cow)",
                                               "BNC1(goat)","BNC2(goat)", "BNC3(goat)","UNC1(goat)", "UNC2(goat)",
                                               "CTB(dog)", "EVT(dog)", "STB(dog)", 
                                               "Endo(human)", "Endo(macaca)","Endo(guinea pig)","Endo(rat)",'Endo(mouse)',"Endo(rabbit)","Endo(pig)",  "Endo(cow)","Endo(goat)","Endo(dog)",
                                               "Stro(human)",'Stro(macaca)', 'Stro(guinea pig)',"Stro(rat)","Stro(mouse)", "Stro(rabbit)","Stro(pig)", "Stro(cow)", "Stro(goat)","Stro(dog)",
                                               "Epi(guinea pig)",'Epi(rat)','Epi(rabbit)','Epi(pig)','Epi(cow)','Epi(goat)','Epi(dog)',
                                               "DC(human)","Tcell(human)","Tcells(macaca)","Tcells(rat)",'Tcells(mouse)',"Tcells(dog)",
                                               'Mac(human)',"Mac(macaca)", "Mono(guinea pig)", "Mac(guinea pig)", "Mac(rat)",'Mac(mouse)',"Mac(rabbit)","Mac(pig)",'Mac(cow)','Mac(goat)',"Mac(dog)",
                                               "NKcells(macaca)", "NKcells(rat)", "NKcells(mouse)", "Bcells(rat)", "Neu(cow)",'Neu(goat)',"Leu(cow)"))





library(ggplot2)




library("viridis")
cols <- viridis(100)[c(30, 50, 100)]
cols
#[1] "#355E8DFF" "#218F8DFF" "#FDE725FF"
# Open a PDF device
pdf(file = "heatmap_output_FINAL3.pdf", width = 15, height =20)  # Adjust width and height as needed

# Generate the heatmap and save it to the PDF
DoHeatmap(markerdata,
          features = as.character(unique(markers$markers)),
          group.by = "celltype_merge",
          assay = 'RNA',draw.lines = FALSE,
          group.colors = c('#a4c7db','#93bad4','#85afce',
                           '#2f6fa7','#3a7ca4','#49889f',
                           '#a9d486','#99cd7b','#8ac670',
                           '#39993a','#f0b362','#efaa56','#f0a34b','#ef9b40','#f1bb6f',
                           '#ec8434','#ec8842','#ed8a51',
                           '#e88684','#e57675','#e16967','#dd5b5a','#da4e4d',
                           '#c787a1','#c69ab7','#c5add0','#cc464c','#cb5460','#c96475','#c7758a',
                           '#785291','#8b6992','#9c8295','#704796','#673a8f',
                           '#eae488','#e0d07b','#d7bc6d',
                           '#6599c0','#68a499','#6cb65a','#5fa03e','#ed8d2f','#ed8f71','#d43433','#ac8ebe','#c1b396','#c49454',
                           '#75a3c6','#57969c','#7bbe65','#4a9b37','#ee9535','#ed8d60','#d74141','#b89dc6','#af9a96','#cea961',
                           '#5eae4f','#78a547','#ee9381','#d12d29','#a080b6','#d2ca96','#ba8047',
                           '#4983b2','#3979ac','#96ca8b','#a8ad55','#eb8229','#a65a2f',
                           '#558eb9','#76b193','#4ca644','#409f3a','#90a94e','#ec8726','#ee9794','#ce3739','#9371ae','#e3e095','#af6c3b',
                           '#88be8f','#c0b25e','#eb7e21','#d8b665','#7c569e','#f2f499','#8863a5')) +
  scale_fill_gradientn(colors = c("#355E8DFF", "#218F8DFF", "#FDE725FF"))+theme(
    legend.position = "none")
#

# Close the PDF device
dev.off()





#Heatmap
gene <- c(
  "CDH1","TP63", "FLT1", "ANO4", "TCF7L2", "SOX5","FOLR1", "GRN", "PLET1", "PLBD1", "CITED1",
  "FOLR1", "CTBS",  "PPAG3",  "PAG6", "CITED1",
  "DAB2", "LRP2", "KRT7","ITGB6", "PPARG", "STRA6", "CHST8", "SGK1", "TFAP2C", "GATA2","GATA3","YBX3","CITED2",
  "CHST3", "ST3GAL3", "SLC45A3",
  "ST3GAL5", "SLC15A1", "SLC15A1", "SLC16A10", "SLC1A3", "SLC20A1", "SLC22A3",
  "SLC22A5", "SLC23A2", "SLC25A21", "SLC25A3", "SLC25A3", "SLC25A6", "SLC25A6",
  "SLC2A2", "SLC2A2", "SLC39A11", "SLC40A1", "SLC45A3", "SLC4A7", "SLC7A4",
  "SLC7A8", "SLC7A8", "SLC9A3", "SOX5", "CTNNA1",  "DAB2"
  
)

#计算平均表达量
gene_cell_exp <- AverageExpression(asthma,
                                   features = unique(gene),
                                   group.by = 'celltype_rep',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)

#complexheatmap作图
library(ComplexHeatmap)
#顶部细胞类型注释
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'class'



top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('UNC5'="#9ECABE",
                                                  'UNC4'="#F6F5B4",
                                                  'UNC3'="#2F528F",
                                                  "UNC2"="#E3AD68",
                                                  "UNC1"="#ACD45E")))#颜色设置
#数据标准化缩放一下
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
# 打开 PDF 设备并设置文件名和尺寸
pdf(file = "/home/hkli/database/TGH_analyse/Cow_dev_placenta/PIG_REP2/PIG_TR_monocle/monocle2/pig_gene_Heatmap.pdf", width = 3.43, height = 7)  # 根据需要调整宽度和高度

# 绘制热图
Heatmap(marker_exp,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        column_title = NULL,
        heatmap_legend_param = list(title = ""),
        col = colorRampPalette(c("#355E8DFF", "#218F8DFF", "#FDE725FF"))(100),
        border = 'black',
        rect_gp = gpar(col = "black", lwd = 0.5),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        column_order = c("UNC1", "UNC2", "UNC3", "UNC4", "UNC5"), # 设置列顺序
        top_annotation = top_anno)
dev.off()






#DotPlot
genes <- c(
  "MITF", "TEAD1", "GATA3", "TP63", "PEG10", "TCF7L2", "GATA2", "TFAP2A",
  "PAG-8", "ERV3-1", "PHLDA2", "PEG3", "PAGE4", "CDH1", "PAG-2", "KRT7", "CTNNB1",
  "PLET1", "PLBD1", "PAG6", "PPAG3",
  "TFPI", "PAG-3", "PAG-11", "PRP6", "PRP1", "PAG-6",
  "PAG-3", "PAG-6", "PRP6",
  "CDH18", "GPC6", "CDH12", "FSTL5",
  "KDR", "PECAM1", "VWF",
  "GATM", "SERPINE2", "LGALS15",
  "EPCAM", "KRT18", "WFDC2",
  "TGFBI", "CPXM1", "ANGPTL2",
  "WFDC1", "COL3A1", "ACTA2",
  "CENPF", "CKAP2", "TOP2A"
)

goat$celltype_harmony <- factor(x = goat$celltype_harmony,
                               levels = c(
                                 "UNC1", "UNC2", "BNC1", "BNC2", "BNC3", 
                                 "Endo1", "Endo2", "GECs", "LECs", "Stro", 
                                 "SMC", "Neu","Mac"
                               ))


DotPlot(goat,features = unique(genes),group.by="celltype_harmony")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  theme_bw()+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))






































