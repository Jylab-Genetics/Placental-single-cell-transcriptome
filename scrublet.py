#####################################
##双细胞预测-scrublet （10x）
#scrublet批量双细胞预测代码
############################################
import sys
print(sys.path)
%matplotlib inline
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import gzip
import io
#修改路径
print(os.getcwd())#显示当前路径
input_dir ='H:/cellranger_V6_matrix_bam/TP_C4_2/filtered_feature_bc_matrix/'
#导入表达矩阵和基因列表
#读入10X的scRNA-seq矩阵，读入raw counts矩阵为scipy sparse矩阵，cells作为行，genes作为列：
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
#读入基因
filename = input_dir + '/features.tsv.gz'
with gzip.open(filename, 'rb') as f:
    with io.TextIOWrapper(f, encoding='utf-8') as file:
        gene_names = [line.strip()for line in file]
        genes = np.array(gene_names)
#
print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))
#####
#expected_doublet_rate值看下面链接选择0.05-0.1
#
#Multiplet Rate (%)  # of Cells Loaded # of Cells Recovered
#~0.4%                ~870                 ~500
#~0.8%               ~1700                ~1000
#~1.6%                ~3500                ~2000
#~2.3%                 ~5300                 ~3000
#~3.1%                 ~7000                 ~4000
#~3.9%                 ~8700                 ~5000
#~4.6%                  ~10500                 ~6000
#~5.4%                   ~12200                ~7000
#~6.1%                   ~14000                 ~8000
#~6.9%                 ~15700                 ~9000
#~7.6%                      ~17400            ~10000

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.05)#看上面链接                  
#计算doublet score
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
#################
#scrub.call_doublets(threshold=0.25)
scrub.plot_histogram();

print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
print('Done.')

#UMAP可视化
scrub.plot_embedding('UMAP', order_points=True)
# doublets占比
print (scrub.detected_doublet_rate_)
# 0.0015981735159817352

#输出结果
barcode =np.array(pd.read_csv(input_dir + '/barcodes.tsv.gz',header=None,index_col=None))
a=np.array([barcode[:,0],doublet_scores,predicted_doublets])
data = pd. DataFrame({'barcode': a[0, :], 'score': a[1, :], 'prediction': a[2, :]}) 
data.to_csv('./TP_C4_2_doublet_prediction.csv', index=False,header=True)
#
data.head()




