#include "guass.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// 创建新节点
static Dnode* create_node(double data) {
	Dnode* new_node = (Dnode*)malloc(sizeof(Dnode));
	if (new_node == NULL) {
		printf("内存分配失败\n");
		return NULL;
	}
	new_node->data = data;
	new_node->next = NULL;
	return new_node;
}

// 添加节点到链表末尾
static void append_node(Dnode** head, double data) {
	Dnode* new_node = create_node(data);
	if (*head == NULL) {
		*head = new_node;
		return;
	}
	Dnode* temp = *head;
	// 找到链表末尾
	while (temp->next != NULL) {
		temp = temp->next;
	}
	temp->next = new_node;
}

// 打印链表
void print_list(Dnode* head) {
	Dnode* temp = head;
	while (temp != NULL) {
		printf("%.30lf \n", temp->data);
		temp = temp->next;
	}
	printf("NULL\n");
}

// 释放链表内存
void free_list(Dnode* head) {
	Dnode* temp;
	while (head != NULL) {
		temp = head;
		head = head->next;
		free(temp);
	}
}

double gaussexp(int x, double sigma) {
	return exp(-1 * x * x / (2 * sigma * sigma));
}

double gaussianDistribution(data* x, double center, double sigma) {
	double result = 0.0;
	int tau = 14; //尾部截断参数
	int zmax = (int)ceil(tau * sigma);
	int xmax = (int)ceil(center) + zmax;
	int xmin = (int)floor(center) - zmax;
	int xrange = xmax - xmin;
	double* temp = (double*)malloc(xrange * sizeof(double));

	for (int i = 0; i < xrange; i++)
	{
		temp[i] = gaussexp(xmin + i, sigma);
		result += temp[i];
	}

	x->xmin = xmin;
	x->xmax = xmax;
	x->result = create_node(temp[0] / result); //初始化第一个节点为具体数值
	for (int i = 1; i < xrange; i++)
	{
		append_node(&x->result, temp[i] / result);
	}

	free(temp);

	return result;
}

// 半高斯分布Z+
double HalfGaussianDistribution(data* x, double sigma) {
	double result = 0.0;
	int tau = 20; //尾部截断参数
	int zmax = (int)ceil(tau * sigma);
	double* temp = (double*)malloc(zmax * sizeof(double));

	for (int i = 0; i < zmax; i++)
	{
		temp[i] = gaussexp(i, sigma);
		result += temp[i];
	}

	x->xmin = 0;
	x->xmax = zmax;
	x->result = create_node(temp[0] / result); //初始化第一个节点为具体数值
	for (int i = 1; i < zmax; i++)
	{
		append_node(&x->result, temp[i] / result);
	}

	free(temp);

	return result;
}

double expk(int k) {
	return exp(-1.0 * k * k / 2);
}

double expkDistribution(data* s, int theta) {
	double result = 0.0;
	int xmax = 1;
	double eps = pow(10, -theta);
	while (expk(xmax) > eps)
	{
		xmax++;
	}
	double* temp = (double*)malloc(xmax * sizeof(double));

	for (int i = 0; i < xmax; i++)
	{
		temp[i] = expk(i + 1);
		result += temp[i];
	}

	s->xmin = 1;
	s->xmax = xmax;
	s->result = NULL;
	for (int i = 0; i < xmax; i++)
	{
		append_node(&s->result, temp[i] / result);
	}

	free(temp);

	return result;
}

// 对称形式下的累加概率密度函数
Dnode* accumulateCDF(data* s) {
	Dnode* temp = s->result;
	Dnode* result = NULL;
	// 先找到链表中值最大的节点
	while (temp != NULL && temp->next != NULL) {
		if (temp->data > temp->next->data) {
			result = create_node(temp->data);
			temp = temp->next;
			break;
		}
		else {
			temp = temp->next;
		}
	}
	Dnode* temp1 = result;
	// 从最大值开始向后累加
	while (temp != NULL && temp1 != NULL) {
		append_node(&result, temp1->data + 2 * temp->data);
		temp1 = temp1->next;
		temp = temp->next;
	}

	return result;
}

// 累加概率密度函数
Dnode* accumulateCDF_n(data* s) {
	Dnode* temp = s->result;
	Dnode* result = create_node(temp->data);
	Dnode* head = result;
	while (result != NULL && temp->next != NULL) {
		append_node(&result, result->data + temp->next->data);
		temp = temp->next;
		result = result->next;
	}
	return head;
}