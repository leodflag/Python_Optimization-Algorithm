# Optimization-Algorithm
This is a project that implement optimization algorithms. For example : gene algorithm.

## 目標
    Max f (x1, x2 ) = 21.5 + x1 sin(4πx1) + x2 sin(20πx2 ) 
## 條件限制
    −3.0 ≤ x1 ≤ 12.1,  4.1 ≤ x2 ≤ 5.8 
## 虛擬碼
    a. 產生x1, x2的各6個十進位母體基因
    b. 將十進位母體基因轉換成二進位
    c. 計算適應值
    d. 綁定母體基因、以及用母體基因算出的適應值
    e. 用輪盤法挑選基因
    f. 奇數基因與偶數基因交配產生兩條子代基因
    g. 判斷子代基因範圍，範圍內的才接收，範圍外的仍用母體基因
    h. 判斷6個子代基因是否需突變
    i. 突變後仍在範圍內才將子代基因加到母體並取代母體
    j. 重複c-i ，迭代50次
## 設定
    交配機率設1，突變機率設0.1 
## 函式
    str[start:end] # 從start索引開始取，在end前結束，字串可相+
    list(zip(A,B)) # 將2個陣列按序列結合成1個zip，再用list轉成陣列
    A[:],B[:]=zip(*x) # 將x分開成2個陣列
    bin(x)[2:].zfill(num) # 用bin轉成二位元，去除前兩個後，用zfill在數字前補0
