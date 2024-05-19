# 卡方正态性检验
# 分析工具集合

import numpy as np
from math import ceil, sqrt, exp
from scipy.stats import chisquare, moment, skew, kurtosis
# 数据分析用的
from matplotlib import pyplot as plt
import pandas as pd

# 读取数据，至少要有五个数据，赛题检验样本数为十万以上
with open(r'output.txt', 'r') as file:
    data_array = [int(i) for i in file.read().split()]

# 正态分布参数
mu = 0
sigma = 4

# 尾截率
tau = 14
# Minimal size of a bucket for the chi-squared test (must be >= 5)
# 如果某一类别样本数小于这个值，会导致卡方检验结果不准确
# 将样本数小于这个值的类别合并到相邻的类别中，直到所有类别样本数都大于这个值
# 合并为新类别，新类别的样本数 = 合并前类别的样本数之和
chi2_bucket = 10
# Minimal p-value
pmin = 0.001


#——————————————————————————————————————————————#
# 计算离散高斯分布的概率分布，输出一个字典，键值对 {x:p(x)} ，x属于[-zmax, zmax)
def gaussian(x, mu, sigma):
    """
    Gaussian function of center mu and "standard deviation" sigma.
    """
    return exp(- ((x - mu) ** 2) / (2 * (sigma ** 2)))


def make_gaussian_pdt(mu, sigma):
    """
    Make the probability distribution table (PDT) of a discrete Gaussian.
    The output is a dictionary.
    """
    # The distribution is restricted to [-zmax, zmax).
    zmax = ceil(tau * sigma)
    pdt = dict()
    for z in range(int(np.floor(mu)) - zmax, int(ceil(mu)) + zmax):
        pdt[z] = gaussian(z, mu, sigma)
    gauss_sum = sum(pdt.values())
    for z in pdt:
        pdt[z] /= gauss_sum
    return pdt

exp_histogram = make_gaussian_pdt(mu, sigma)
#——————————————————————————————————————————————#


zmax = ceil(tau * sigma)
nsamples = len(data_array)
histogram = dict() # 把data_array的数据放入一个字典中,键值对为{x: 样本中x的数量}
exp_histogram1 = dict() # 创建一个新字典，防止在原数据上修改，键值对为 {x: x的样本数的期望值}
outlier = 0
outs = np.array([])
Values = np.array(range(int(np.floor(mu)) - zmax, int(ceil(mu)) + zmax))
# Initialize histogram
for z in Values:
    histogram[z] = 0
    exp_histogram1[z] = 0
for z in data_array:
    # Detect and count outliers (samples not in [-zmax, zmax))
    if z not in histogram:
        outlier += 1
        outs = [outs, z] # 把出界值加入outs，方便之后输出
    # Fill histogram according to the samples
    else:
        histogram[z] += 1

observed_nums = np.array(list(histogram.values()))

for z in Values:
    exp_histogram1[z] = exp_histogram[z] * nsamples # 期望样本数 = 期望概率 * 样本数

# 将字典的值转化为列表
expected_freq_raw = list(exp_histogram1.values())
observed_freq_raw = list(histogram.values())

# 将样本数小于chi2_bucket的类别合并到相邻的类别中，直到所有类别样本数都大于chi2_bucket
# 初始化
expected_freq = expected_freq_raw
observed_freq = observed_freq_raw

z = 0
while z < len(expected_freq) - 1:
    # 如果当前类别的观察频数小于阈值
    if observed_freq[z] < chi2_bucket:
        # 将当前类别和下一个类别合并
        expected_freq[z] += expected_freq[z + 1]
        observed_freq[z] += observed_freq[z + 1]
        # 删除下一个类别
        expected_freq.pop(z + 1)
        observed_freq.pop(z + 1)
    else:
        # 如果当前类别的观察频数不小于阈值，继续检查下一个类别
        z += 1

# 如果最后一个类别的观察频数小于阈值，将其与前一个类别合并
if observed_freq[-1] < chi2_bucket / nsamples:
    expected_freq[-2] += expected_freq[-1]
    observed_freq[-2] += observed_freq[-1]
    expected_freq.pop(-1)
    observed_freq.pop(-1)

# ValueError: For each axis slice, the sum of the observed frequencies must agree with the sum of the expected frequencies to a relative tolerance of 1e-08
# 实在不行,再加上下面一句吧
# expected_freq = expected_freq * (observed_freq.sum() / expected_freq.sum())

# 卡方检验
chi2_stat, p_value = chisquare(observed_freq, f_exp=expected_freq)
# 其他统计量计算
mean = sum(data_array) / len(data_array)
stdev = sqrt(moment(data_array, 2))
skewness = skew(data_array)
kurtosisV = kurtosis(data_array)

# 输出结果
print(f"样本数：{nsamples}")
print(f"无效样本数(必须为0)：{outlier}")
if outlier != 0:
    print("无效样本值：",outs)
print(f"卡方统计量(越小越接近): {chi2_stat}")
print(f"P值: {p_value}")
if p_value > pmin and outlier == 0:
    print("样本有效")
else:
    print("样本无效")
print(f"均值mu(中心): {mu} || {mean}   ")
print(f"标准差sigma: {sigma} || {stdev}")
print(f"偏度(期望为0): {skewness}")
print(f"峰度(期望为0): {kurtosisV}")

exp_ = np.array(list(exp_histogram1.values()))
obs_ = np.array(list(histogram.values()))
freq_diff = np.array(expected_freq) - np.array(observed_freq) # 直接相减是正确的,因为一一对应。通过样本数放大
freq_diff_raw = exp_ - obs_

# 使用csv文件存储分析数据
# 创建一个DataFrame
df1 = pd.DataFrame({
    'Values': Values,
    'obs_nums': observed_nums,
    'exp_freq_row': exp_,
    'obs_freq_row': obs_,
    'Freq_Diff_Raw': freq_diff_raw
})

df2 = pd.DataFrame({
    'exp_freq': expected_freq,
    'obs_freq': observed_freq,
    'Freq_Diff': freq_diff
})

# 导出为CSV文件
df1.to_csv('data1_raw.csv', index=False)
df2.to_csv('data2_compute.csv', index=False)

# 绘制图形
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.bar(Values, obs_, label='Observed', alpha=0.6)
plt.plot(Values, exp_, label='Expected', color='red')
plt.title('Frequency Distribution')
plt.legend()

plt.subplot(1, 2, 2)
plt.bar(Values, freq_diff_raw, label='Difference', alpha=0.6)
plt.title('Frequency Difference')
plt.legend()

plt.tight_layout()
plt.show()
