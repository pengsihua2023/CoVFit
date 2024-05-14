import pandas as pd

# 文件名列表
filenames = ['Book1.txt', 'Book2.txt', 'Book3.txt']

# 读取第一个文件初始化mean_dataframe
mean_dataframe = pd.read_csv(filenames[0], sep='\t', header=None)

# 遍历剩余的文件名列表，读取每个文件到DataFrame，并与mean_dataframe相加
for filename in filenames[1:]:
    df = pd.read_csv(filename, sep='\t', header=None)
    mean_dataframe += df

# 计算平均值
mean_dataframe /= len(filenames)

# 将计算得到的平均值DataFrame保存为TSV文件
mean_dataframe.to_csv('cell_mean_value-sihua2.txt', sep='\t', header=None, index=False)
