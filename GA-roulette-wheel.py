# -*- coding: utf-8 -*-
"""
Created on Sun Apr  8 01:52:05 2018
@author: lovedoglion
Genetic Algorithm
roulette wheel selection
"""

import random
import math
import matplotlib.pyplot as plt
#畫圖的圖表

string_length = 8  # 基因長度
pop_gene_num = 6  # 母群體數量
pool = []  # 交配池
best_gene = []  # 最佳基因
itera = 50  # 迭代次數
mutation_rate = 0.01  # 突變率
mutation_num=0
#母群數量
gene_1 = []
gene_2 = []
gene = []

adept = []  # 適應值
Random = 0  # 輪盤法亂數

#基因範圍判斷
gene_1_U = 255
gene_1_L = 0
gene_2_U = 149
gene_2_L = 120

#繪圖參數
pltX = []  # x軸為迭代次數
pltY = []  # y軸為平均適應值
"""
轉換
"""
#實際對應數值換算成十進位
def floTurnTen(x):
    xx = (x-(-3.0))*255/15.1
    return round(xx)  # 無條件進位
#十進位轉換成字串基因碼
def turnStrGene(x):
    # 轉換二進位後取序列2之後的數值，在前填滿0達到指定字串長度
    return bin(x)[2:].zfill(string_length) 
#二進位的字串轉回十進位
def bin_Int(x):
    return int(x, 2) # 二進制的數值x用int轉回十進位
#十進位換算成真實浮點數
def tenTurnflo(x):
    xx = 15.1*x/255-3.0
    return round(float(xx), 1) # 取到小數點後1位


"""
基因的操作
"""
#指定二進位位置突變  要重寫
def _invert_at(s, index):
    if s[index]=='0':
        return s[:index]+'1'+s[index+1:]
    else:
        return s[:index]+'0'+s[index+1:]

#大到小排序，會變更原順序
def order(x):
    temp = 0
    for i in range(6):
        for j in range(6):
            if x[j] < x[i]:
                temp = x[j]
                x[j] = x[i]
                x[i] = temp
    return x

#解壓縮
def zipReturn(pool):
    x = []
    adept[:], x[:] = zip(*pool) # 解壓縮出適應值與基因陣列
    print('adept,x', adept, x)
    print('x', x)
    gene_1[:], gene_2[:] = zip(*x) #  得出基因1與基因2兩條陣列
    print('gene_1,gene_2', gene_1, gene_2)
    return gene_1, gene_2

#判斷大小，輸入int
def range_gene_1(x):
    if x >= gene_1_L and x <= gene_1_U:
        return 1
    else:
        return 0
def range_gene_2(x):
    if x >= gene_2_L and x <= gene_2_U:
        return 1
    else:
        return 0

#刪掉前6個基因
def delList(x):
    if itera != 50:
       del x[0:6]

"""
# 適應函數Adaptation function
"""
def Adaptation(x1, x2):
    f = 21.5+x1*math.sin(4*math.pi*x1)+x2*math.sin(20*math.pi*x2)
    return f

"""
產生初代基因
"""
x = ''
for i in range(pop_gene_num):
    s1 = random.randint(gene_1_L, gene_1_U) # 在限制內隨機產生基因母體
    s2 = random.randint(gene_2_U, gene_2_U)
    gene_1.append(turnStrGene(s1)) # 產生初代基因陣列
    gene_2.append(turnStrGene(s2))
print('初代基因組1', gene_1)
print('初代基因組2', gene_2)
#初代基因組
x = list(zip(gene_1, gene_2)) # 產生初代基因組
print('初代基因組', x)
print('初代基因組排序', order(x))

"""
繪圖函數
"""
def plotData(plt, data):
    x = [p[0] for p in data]
    y = [p[1] for p in data]
    plt.plot(x, y)

"""
迭代開始
"""
while(itera > 0):
    #適應值總和
    adeptSUM = 0
    #計算實際適應值float
    for i in range(pop_gene_num):
        x1 = tenTurnflo(bin_Int(gene_1[i]))
        x2 = tenTurnflo(bin_Int(gene_2[i]))
        adept.append(Adaptation(x1, x2)) # 算出每組適應值
        #print('adept['+str(i)+']'+str(adept[i]))
        adeptSUM += adept[i] # 計算輪盤法的分母
    #print(SUM)
    delList(adept)  # 刪掉前6個基因
    pltY.append(adeptSUM/pop_gene_num)
    #初代基因跟適應值綁定
    gene = list(zip(adept, x))
    print('初代基因跟適應值', gene)
    print('排序適應值', order(adept))
    print('排序適應值+基因', order(gene)) # 適應值被排序由大到小 要存成新的
    #留適應值最大的基因組
    if best_gene < list(gene[0]):
        best_gene = list(gene[0])
    #輪盤法 隨機式
    #計算適應值機率，以及製造輪盤
    adeptR = []  # 適應值機率
    plusAdeptR = []  # 輪盤
    R = 0 # 該適應值所在的機率區間
    for i in range(pop_gene_num):
        adeptR.append(adept[i]/adeptSUM) # 按照已排序的基因適應值順序計算
        R += adeptR[i]
        plusAdeptR.append(R) # 對應基因順序製造機率輪盤
        #print('adeptR['+str(i)+']'+str(adeptR[i]))
    print('基因出現機率順序', adeptR)
    print('輪盤', plusAdeptR)
    #取6組亂數加判斷到交配池
    for i in range(pop_gene_num):
        Random = random.random()
        print('隨機亂數', Random) 
        if Random <= plusAdeptR[i]: # 判斷隨機亂數在哪個區間
            pool.append(gene[0])
        elif Random > plusAdeptR[0] and Random <= plusAdeptR[1]:
            pool.append(gene[1])
        elif Random > plusAdeptR[1] and Random <= plusAdeptR[2]:
            pool.append(gene[2])
        elif Random > plusAdeptR[2] and Random <= plusAdeptR[3]:
            pool.append(gene[3])
        elif Random > plusAdeptR[3] and Random <= plusAdeptR[4]:
            pool.append(gene[4])
        else:
            pool.append(gene[5])
    #刪掉前6個
    delList(pool)
    print('交配池', pool)

    #交配 偶數情況
    i = 0
    cut = 0
    x = []
    gene1 = []
    gene2 = []
    while(cut < pop_gene_num):
        adept[:], x[:] = zip(*pool) # 分出適應值與基因組2個陣列
        #print('adept,x',adept,x)
        #print('x',x)
        gene_1[:], gene_2[:] = zip(*x) # 分出兩個基因組序列
        #print('gene_1,gene_2',gene_1,gene_2)
        #隨機點交配
        dot = random.randint(0, 7)
        gene1.append(gene_1[cut][:dot]+gene_1[cut+1][dot:])
        gene1.append(gene_1[cut+1][:dot]+gene_1[cut][dot:])

        gene2.append(gene_2[cut][:dot]+gene_2[cut+1][dot:])
        gene2.append(gene_2[cut+1][:dot]+gene_2[cut][dot:])
        #兩條複製完畢
        cut += 2

    #交配完看是否在範圍裡，有才放回並替代掉母體
    for j in range(6):
        if range_gene_1(bin_Int(gene1[j])) == 1:
            gene_1[j] = gene1[j]  # 將子代輸回母體
        if range_gene_2(bin_Int(gene2[j])) == 1:
            gene_2[j] = gene2[j]

    #突變，pool突變部分基因
    for i in range(pop_gene_num):
        rand = random.random()
        if rand <= mutation_rate:
            mutation_num+=1
            s1 = ''
            s2 = ''            
            #突變，輸入為2位元，輸出字串
            h = random.randint(0, 7)  # 隨機位置突變
            s1=_invert_at(gene_1[i], h) # 將基因碼突變
            s2=_invert_at(gene_2[i], h)
            #突變完看是否在範圍裡        
            if range_gene_1(bin_Int(s1)) == 1:
                gene_1[i] = s1
            if range_gene_2(bin_Int(s2)) == 1:
                gene_2[i] = s2

    print('gene_1', gene_1)
    print('gene_2', gene_2)
    pltX.append(50-itera)  # 圖的X軸為次數
    itera -= 1

print('best_gene', best_gene)
print('mutation_num=', mutation_num)

"""
繪圖
"""
#將X軸
geneChart = list(zip(pltX, pltY))
plotData(plt, geneChart)  # X軸
plt.show()