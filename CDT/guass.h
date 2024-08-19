// ����˫���ȸ���������ɢ��˹�ֲ�����
#ifndef GAUSS_H
#define GAUSS_H

// ��������ڵ�ṹ��
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

// ����������˹����exp(-x^2/2*sigma^2)
double gaussexp(int x, double sigma);

// ����˫���ȸ���������ɢ��˹�ֲ�����
// ����ֵΪ�����ܺͺ͸���ֵ���������ģ���ֵ��center����׼��sigma
double gaussianDistribution(data* s, double center, double sigma);

double HalfGaussianDistribution(data* x, double sigma);

// ����ָ������exp(-k^2/2)
double expk(int k);

// ����ָ�������ĸ��ʷֲ�
double expkDistribution(data* s, int theta);

// ������ۼ�CDF���������ۼӽ��������������
Dnode* accumulateCDF(data* s);
Dnode* accumulateCDF_n(data* s);

void print_list(Dnode* head);

void free_list(Dnode* head);

#endif // GAUSS_H
