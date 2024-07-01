# Discrete-Gaussian-Distribution
离散高斯分布采样器

**本项目为2024年密码学挑战赛赛题二参赛代码，现公开，不会再更新。**

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


使用pip安装ELFFile库：  
`pip3 install pyelftools`  
然后运行testdata.py文件

*说明*  
在一个可执行文件中，.text、.data、.bss 和 .rodata 段代表了程序的不同部分：  
.text段：包含程序的机器代码，也就是编译后的二进制指令。这是程序的实际“代码”部分。  
.data段：包含已初始化的全局变量和静态变量。这些变量在程序开始执行前已经被赋予了初始值。  
.bss段：包含未初始化的全局变量和静态变量。在程序开始执行时，这些变量会被自动初始化为零。  
.rodata段：包含只读数据，如字符串常量和其他常量数据。这部分数据在程序运行时不会被修改。  


**使用 avstack.pl 工具**  
avstack.pl 是一个Perl脚本，用于分析GCC编译器生成的 .su 文件，这些文件包含了函数的堆栈使用信息。  
脚本原网站：https://dlbeer.co.nz/oss/avstack.html  
要使用这个工具，你需要：  
1. 在编译时加入`-fstack-usage`选项  
一旦你有了 .su 文件，你可以运行 avstack.pl 脚本，传递所有 .o 文件作为参数。
.su 文件假定位于与其对应的 .o 文件相同的目录中。avstack.pl 脚本将读取所有 .su 文件，并解析 .o 文件来构建调用图，
然后计算每个函数的最大堆栈使用量。
2. 下载文件后使用 chmod 命令给文件添加执行权限： 
`chmod +x avstack.pl`  
3. avr-objdump 是AVR工具链的一部分，通常用于处理AVR微控制器的对象文件。
如果你在处理非AVR架构的代码，你可能需要使用适合你架构的 objdump 工具。
如果你的代码是为x86或x86_64架构编译的，
你应该使用 objdump 而不是 avr-objdump。你可以尝试以下步骤来解决这个问题：  
3.1. 安装binutils：确保你的系统上安装了 binutils，因为 objdump 包含在这个包中。你可以使用以下命令安装：  
  `sudo apt-get install binutils`  
3.2. 修改avstack.pl：  
   编辑 avstack.pl 脚本，将所有的 avr-objdump 调用替换为 objdump。这可以通过简单的文本替换来完成。
5. 在命令行中使用 avstack.pl 的示例：  
`./avstack.pl *.o`


**包含的算法：**  

1. CDT反演法 $c=0, sigma=0.75$  
2. Knuth_Yao法 $c = 0, sigma = 4$
3. 查表法 $c = 0, sigma = 0.75$  
4. 拒绝采样法 $c = 0, sigma = 0.75$  
5. 高斯卷积法 $c = 0, sigma = 1024$  
6. 伯努利拒绝采样法 $c \in [0,1), sigma = 1.5$
7. 伯努利拒绝采样法 $c \in [0,1), sigma \in (0.8,1.6)$
8. Karney算法（未完成）
9. 实时计算的拒绝采样法 $c \in [0,1)$
10. AVX2并行优化的高斯卷积法 $c = 0, sigma = 1024$  


**@感谢**
好友龙神的代码测试  
好友鹏的论文写作
