import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def exp_pca():
    # 加载exp数据，输出维度
    exp_data = np.loadtxt('PRAD/data/RNA_exp.csv', dtype=str, delimiter=',')
    print("exp数据数据维度： " + str(exp_data.shape))

    # 加载基因名，输出维度
    genes_id = exp_data[:, 0]
    genes_id = genes_id.reshape(-1, 1)
    print("基因个数 ： " + str(genes_id.shape))

    # 加载数据部分，输出维度
    exp_data = np.delete(exp_data, 0, axis=1)
    exp_data = exp_data.astype(float)  # 把数据部分转换数据类型
    print("exp数据部分维度“ " + str(exp_data.shape))
    # print(exp_data.dtype)


    # 建模
    pca = PCA(n_components = 10)

    exp_pca = pca.fit_transform(exp_data)
    print(exp_pca.shape)
    print(pca.explained_variance_ratio_)
    print("方差贡献率： " + str(sum(pca.explained_variance_ratio_)))

    print(exp_pca.shape)
    exp_pca_data = np.concatenate((genes_id, exp_pca), axis=1)  # 将降维后的数据部分与基因名合并
    print(exp_pca_data.shape)
    np.savetxt('output/exp_pca.csv', exp_pca_data, delimiter=',', fmt='%s')


exp_pca()
