from saga import UnivariateSamples

#初始化变量
mu =  0
sigma = 1024

# 打开txt文件进行读取
with open(r'output.txt', 'r') as file:
    # 将文件的全部内容读取为一个字符串
    data = file.read()
# 将字符串分割为字符串列表，以空白字符分隔
data_list = data.split()
# 将字符串列表转换为整数列表
data_array = [int(i) for i in data_list]

sub_array_step = data_array[::10] #截取一部分

a = []

us1 = UnivariateSamples(mu, sigma, data_array)
print(us1)

