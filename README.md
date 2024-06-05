# Discrete-Gaussian-Distribution
离散高斯分布采样器

**使用方法**  
`make`  
`./example`  
将会生成一个output.txt文件存储一秒内所有的采样结果

**内存检测**
安装：  
`sudo apt update`
`sudo apt install valgrind`  
内存占用检测：  
`valgrind --tool=massif ./example`  
内存使用情况检测：  
`valgrind --leak-check=full --show-leak-kinds=all --undef-value-errors=no ./example`


**包含的算法：**  

1. CDT反演法 $c=0, sigma=0.75$  
2. Knuth_Yao法 $c = 0, sigma = 4$
3. 查表法 $c = 0, sigma = 0.75$  
4. 拒绝采样法 $c = 0, sigma = 0.75$  
5. 高斯卷积法 $c = 0, sigma = 1024$  
6. 伯努利拒绝采样法 $c \in [0,1), sigma = 1.5$
7. 伯努利拒绝采样法 $c \in [0,1), sigma \in (0.8,1.6)$
8. Karney算法
9. 实时计算的拒绝采样法 $c \in [0,1)$
10. AVX2并行优化的高斯卷积法 $c = 0, sigma = 1024$  
