import numpy as np
import math


def hit_rate2():
    sort_result = []
    with open("BRCA/contrast/p3_4_InterOmics.txt") as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            sort_result.append(line)
    print("list长度：" + str(len(sort_result)))

    # cancer_id = np.loadtxt("data/breast_cancer_genes.txt", delimiter='\t', dtype=str)
    cancer_id = np.loadtxt("data/breast_cancer_genes.txt", delimiter='\t', dtype=str)

    group = list(range(10, 110, 10))
    print(group)
    print(len(group))

    rate = []
    for i in range(50):
        current_modul = sort_result[i]
        print(current_modul)
        count = 0
        for item in current_modul:
            if item in cancer_id:
                count += 1
        current_rate = count / len(current_modul)
        rate.append(current_rate)

    print(rate)
    # np.savetxt("GBM/contrast/p7_4_hit_rate.csv", rate, delimiter=',', fmt='%s')

    av_rate = []
    for i in group:
        # print(i)
        # print(rate[0:i])
        av_rate.append(sum(rate[0:i])/i)

    print(av_rate)


# hit_rate2()

def save_list(filename, file_list):
    # output = open('E:/6.1/data/ppi_sp.xls','w',encoding='gbk')
    output = open(filename, 'w', encoding='gbk')
    for i in range(len(file_list)):
        for j in range(len(file_list[i])):
            output.write(str(file_list[i][j]))  # write函数不能写int类型的参数，所以使用str()转化
            output.write('\t')  # 相当于Tab一下，换一个单元格
        output.write('\n')  # 写完一行立马换行
    output.close()


# cancer_list = []
# filename = "E:/5.15/cancer_label.txt"
# with open(filename) as f:
#     lines = f.readlines()
#     for line in lines:
#         line = line.split()
#         cancer_list.append(line[0])

# modul_list = []
# data = np.loadtxt("E:/5.15/out/ModulOmics_ILP2.txt",dtype=str,skiprows=1,delimiter="\t")
# for i in range(0,len(data)):
#     for j in range (0,2):
#         count = 0
#         if data[i][j] in cancer_list:
#             count = count + 1
#     if count != 0:
#         datai=np.insert(data[i],0,i)
#         datai=np.append(datai,count/2)
#         modul_list.append(datai.tolist())
#
# np.savetxt('E:/5.15/result/2.csv',np.array(modul_list),delimiter=',',fmt='%s')

# cancer_list = []
# filename = "data/canonical_drivers_Breast cancer.txt"
# with open(filename) as f:
#     lines = f.readlines()
#     for line in lines:
#         line = line.split()
#         cancer_list.append(line[0])



# filepath = 'data/out/ModulOmics_ILP'
# savepath = 'result/modul_out_ILP_'

# filepath = 'simIN_out/modul/ModulOmics'
# savepath = 'result/modul_out_'

# analyse_modul(filepath, savepath)

def hit_1():
    sort_result = np.loadtxt('BRCA/out/p9/p9_sort_2_name.txt', delimiter='\t', dtype=str)
    cancer_id = np.loadtxt('data/breast_cancer_genes.txt', delimiter='\t', dtype=str)
    # cancer_id = np.loadtxt('data/cancer_list.txt', delimiter='\t', dtype=str)
    print("数据维度：" + str(sort_result.shape))

    count_list = []
    count = 0

    for i in range(200):
        # if sort_result[i][0] in cancer_id or sort_result[i][1] in cancer_id or sort_result[i][2] in cancer_id:
        if sort_result[i][0] in cancer_id or sort_result[i][1] in cancer_id:
            count = count + 1
        count_list.append(count)
    print(count)
    print(count_list)
    print(len(count_list))
    # np.savetxt("result/mReg_K2_hit1.csv", count_list, delimiter=',', fmt='%s')
    # np.savetxt("result/p12_K2_hit1.csv", count_list, delimiter=',', fmt='%s')

# hit_1()


def hit_2():
    sort_result = np.loadtxt('LUNG/out/p3/p3_sort_3_name.txt', delimiter='\t', dtype=str)  # 排序好的文件
    cancer_id = np.loadtxt('data/GBM_cancer_genes.txt', delimiter='\t', dtype=str)  # 已知列表
    print(sort_result.shape)

    count_list = []
    count = 0

    for i in range(200):
        hit_count = 0
        for j in range(0, 3):
            # for j in range(0,3):
            if sort_result[i][j] in cancer_id:
                hit_count = hit_count + 1
        if hit_count >=  2:
            # if hit_count >= 3:
            # if hit_count >= 2:
            count = count + 1
        count_list.append(count)

    print(count)
    print(count_list)
    print(len(count_list))
    # np.savetxt("result/p12_K5_hit4.csv",count_list,delimiter=',',fmt='%s')


# hit_2()



def cancer_hit_count():
    sort_result = np.loadtxt('data/out/mutExpPpi_sort5_name.txt', delimiter='\t', dtype=str)
    cancer_id = np.loadtxt('data/canonical_drivers_Breast cancer.txt', delimiter='\t', dtype=str)

    cancer_set = set()
    sort_set = set()
    for i in range(0, 200):
        for j in range(0, 5):
            sort_set.add(sort_result[i][j])
            if sort_result[i][j] in cancer_id:
                cancer_set.add(sort_result[i][j])
    print(str(len(cancer_set)))
    print(str(len(sort_set)))
    print(cancer_set)


# cancer_hit_count()




def load_modul_list():
    filename = 'PRAD/out/p18/com_p18.txt'
    modul_list = []
    score = np.loadtxt('PRAD/out/p18/score_p18.txt', delimiter='\t', dtype=str)
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            modul_list.append(line)
    print(len(modul_list))
    print(modul_list)

    score_1 = []  # 直接均值
    score_2 = []  # 考虑size
    for i in range(len(modul_list)):
        list_size = len(modul_list[i])
        # print(list_size)
        list_edges = len(modul_list[i]) * (len(modul_list[i]) - 1) / 2
        list_score = float(score[i])

        current_score1 = list_score / list_edges
        x = math.log(1 + 1 / (1 - current_score1 ** 2))
        s = (list_size / 2) ** (1 / 2)
        z = 4 * math.exp(s) / (1 + 4 * math.exp(s))
        # current_score2 = (list_size/2)**2*0.5*x
        current_score2 = z * x

        score_1.append(current_score1)
        score_2.append(current_score2)
        # score_2.append((list_num ** (1 / 4)) * list_score / (2 * list_edges))
        modul_list[i].append(str(current_score1))

        # modul_list[i].append(str(x))
        #
        # modul_list[i].append(str(z))

        modul_list[i].append(str(current_score2))
        # modul_list[i].append((list_num ** (1 /4)) * list_score / (2 * list_edges))
        # print(list_num)

    score1 = np.array(score_1)
    s1 = (-score1).argsort()
    score2 = np.array(score_2)
    s2 = (-score2).argsort()
    modul_list_score1 = []
    modul_list_score2 = []
    for i in s1:
        for j in range(len(modul_list)):
            if j == i:
                modul_list_num = modul_list[j]
                modul_list_score1.append(modul_list_num)
    for i in s2:
        for j in range(len(modul_list)):
            if j == i:
                modul_list_num = modul_list[j]
                modul_list_score2.append(modul_list_num)

    print(score_1)
    print(score_2)
    print(len(score_1))
    print(len(score_2))

    save_list('PRAD/out/p18/p18_modul_4_sort1.txt', modul_list_score1)
    save_list('PRAD/out/p18/p18_modul_4_sort2.txt', modul_list_score2)


load_modul_list()


def hit_rate():
    sort_result = []
    with open("PRAD/contrast/p18_4_InterOmics.txt") as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            sort_result.append(line)
    print("list长度：" + str(len(sort_result)))

    cancer_id = np.loadtxt("data/PRAD_cancer_genes.txt", delimiter='\t', dtype=str)

    group = list(range(5, 55, 5))
    print(len(group))

    rate = []
    for i in group:
        current_set = set()
        current_modul = sort_result[:i]
        print(current_modul)
        for m in current_modul:
            for n in m:
                current_set.add(n)
        print(current_set)
        count = 0
        for item in current_set:
            if item in cancer_id:
                count += 1
        current_rate = count / len(current_set)
        rate.append(current_rate)

    print(rate)
    np.savetxt("PRAD/contrast/p18_4_hit_rate.csv", rate, delimiter=',', fmt='%s')


hit_rate()



def gene_AUCn():
    sort_result = []
    with open("PRAD/contrast/hotnet.txt") as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            sort_result.append(line)
    print("list长度：" + str(len(sort_result)))

    # print(sort_result)
    # sort_result = sort_result[0:100]

    print(len(sort_result))

    sort_modul = []

    nodeRank = np.loadtxt("PRAD/data/coReg/nodePageRank.csv", delimiter=',', dtype = str)
    cancer_id = np.loadtxt('data/PRAD_cancer_genes.txt', delimiter='\t', dtype=str)

    for i in range(len(sort_result)):
        res = []
        for item in nodeRank[:,0]:
            if item in sort_result[i]:
                res.append(item)
        sort_modul.append(res)

    print(sort_modul)

    count_set = dict()
    cancer_set = dict()
    for i in range(len(sort_modul)):
        curScore = 1 / (i + 1) ** (1 / 2)
        for j in sort_modul[i]:
            if j in cancer_set:
                count_set[j] = count_set.get(j) + 1
                # cancer_set[j] = cancer_set.get(j) + curScore
            else:
                count_set[j] = 1
                cancer_set[j] = curScore

    print(count_set)
    print(len(cancer_set))

    print(cancer_set)
    print(len(cancer_set))

    # for i in cancer_set:
    #     cancer_set[i] = cancer_set.get(i)/count_set.get(i)

    gene_order = sorted(cancer_set.items(), key=lambda x: x[1], reverse=True)
    print(gene_order)
    T = 100
    print("T的长度:",T)
    print(gene_order[0][0])

    gene_order = gene_order[0:100]


    count_list = []
    auc_list = []
    aim = 0
    # count = 0
    for i in range(len(gene_order)):
        if gene_order[i][0] in cancer_id:
            aim = aim + 1
        else:
            count = aim
            count_list.append(count)

    print(count_list)
    print(len(count_list))
    count_list = np.array(count_list)
    for n in range(len(count_list)):
        auc = 1 / ((n + 1) * T) * np.sum(count_list[:n + 1])
        # auc = 1 / (200 * 98) * np.sum(count_list[:n + 1])
        auc_list.append(auc)
    print(auc_list)
    print(len(auc_list))
    # np.savetxt("data/contrast/p16_InteOmics_aucn_4.csv", auc_list, delimiter=',', fmt='%s')
    np.savetxt("PRAD/contrast/hotnet2_aucn_T.csv", auc_list, delimiter=',', fmt='%s')


# gene_AUCn()
