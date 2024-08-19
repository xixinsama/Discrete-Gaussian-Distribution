// 生成双精度浮点数的离散高斯分布函数
#ifndef GAUSS_H
#define GAUSS_H

// 定义链表节点结构体
typedef struct Node {
	double data;
	struct Node* next;
} Dnode;

typedef struct
{
	int xmin;
	int xmax;
	Dnode* result;
}data;

// 计算整数高斯函数exp(-x^2/2*sigma^2)
double gaussexp(int x, double sigma);

// 生成双精度浮点数的离散高斯分布函数
// 返回值为概率总和和概率值，输入中心（均值）center，标准差sigma
double gaussianDistribution(data* s, double center, double sigma);

double HalfGaussianDistribution(data* x, double sigma);

// 计算指数函数exp(-k^2/2)
double expk(int k);

// 生成指数函数的概率分布
double expkDistribution(data* s, int theta);

// 将结果累加CDF函数，将累加结果放入新链表中
Dnode* accumulateCDF(data* s);
Dnode* accumulateCDF_n(data* s);

void print_list(Dnode* head);

void free_list(Dnode* head);

#endif // GAUSS_H
