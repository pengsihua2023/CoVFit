## PyMOL比较三维结构
```
# 加载结构文件
load output/wildtype.pdb, wildtype
load output/mutant.pdb, mutant

# 对齐结构
align mutant, wildtype

# 可视化差异
show cartoon, wildtype
show cartoon, mutant
color blue, wildtype
color red, mutant


```
