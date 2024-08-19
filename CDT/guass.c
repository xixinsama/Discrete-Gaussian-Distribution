#include "guass.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// �����½ڵ�
static Dnode* create_node(double data) {
	Dnode* new_node = (Dnode*)malloc(sizeof(Dnode));
	if (new_node == NULL) {
		printf("�ڴ����ʧ��\n");
		return NULL;
	}
	new_node->data = data;
	new_node->next = NULL;
	return new_node;
}

// ��ӽڵ㵽����ĩβ
static void append_node(Dnode** head, double data) {
	Dnode* new_node = create_node(data);
	if (*head == NULL) {
		*head = new_node;
		return;
	}
	Dnode* temp = *head;
	// �ҵ�����ĩβ
	while (temp->next != NULL) {
		temp = temp->next;
	}
	temp->next = new_node;
}

// ��ӡ����
void print_list(Dnode* head) {
	Dnode* temp = head;
	while (temp != NULL) {
		printf("%.30lf \n", temp->data);
		temp = temp->next;
	}
	printf("NULL\n");
}

// �ͷ������ڴ�
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
	int tau = 14; //β���ضϲ���
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
	x->result = create_node(temp[0] / result); //��ʼ����һ���ڵ�Ϊ������ֵ
	for (int i = 1; i < xrange; i++)
	{
		append_node(&x->result, temp[i] / result);
	}

	free(temp);

	return result;
}

// ���˹�ֲ�Z+
double HalfGaussianDistribution(data* x, double sigma) {
	double result = 0.0;
	int tau = 20; //β���ضϲ���
	int zmax = (int)ceil(tau * sigma);
	double* temp = (double*)malloc(zmax * sizeof(double));

	for (int i = 0; i < zmax; i++)
	{
		temp[i] = gaussexp(i, sigma);
		result += temp[i];
	}

	x->xmin = 0;
	x->xmax = zmax;
	x->result = create_node(temp[0] / result); //��ʼ����һ���ڵ�Ϊ������ֵ
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

// �Գ���ʽ�µ��ۼӸ����ܶȺ���
Dnode* accumulateCDF(data* s) {
	Dnode* temp = s->result;
	Dnode* result = NULL;
	// ���ҵ�������ֵ���Ľڵ�
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
	// �����ֵ��ʼ����ۼ�
	while (temp != NULL && temp1 != NULL) {
		append_node(&result, temp1->data + 2 * temp->data);
		temp1 = temp1->next;
		temp = temp->next;
	}

	return result;
}

// �ۼӸ����ܶȺ���
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