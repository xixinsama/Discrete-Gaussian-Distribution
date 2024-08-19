#include <stdio.h>
#include "n2b.h"
#include "guass.h"

int main() {
	// 计算整数高斯分布
	data s;
	double sum = gaussianDistribution(&s, 0.0, 1.5);
	print_list(s.result); //NULL是链表结束标志

	int theta = 15;
	Dnode* point = s.result; //链表头指针
	for (int i = 0; i < s.xmax - s.xmin; i++) {
		uint192_t* result = double2int(point->data, theta);
		if (result != NULL) {
			printf("%llu\n", result->bits[0]);
			//printf("bits 1: %llu\n", result->bits[1]);
			//printf("bits 2: %llu\n", result->bits[2]);
		}
		else {
			printf("Invalid theta value\n");
		}
		point = point->next;
	}

	// 计算累加CDF
	Dnode* cdf = accumulateCDF(s.result);
	print_list(cdf);

	point = cdf;
	// uint32_t half0, half1;
	while (point != NULL) {
		uint192_t* resultCDF = double2int(point->data, theta);
		if (resultCDF != NULL) {
			printf("%lluu,\n", resultCDF->bits[0]);
			//printf("bits 1: %llu\n", result->bits[1]);
			//printf("bits 2: %llu\n", result->bits[2]);
			// 分割为两个32位数据，需要根据theta值来判断
			//half0 = (uint32_t)(resultCDF->bits[0] << 34 >> 34);
			//half1 = (uint32_t)((resultCDF->bits[0] >> 30));
			//printf("half0: %lu half1: %lu\n", half0, half1);
		}
		else {
			printf("Invalid theta value\n");
		}
		point = point->next;
	}

	// 计算累加CDF
	Dnode* cdt = accumulateCDF_n(s.result);
	print_list(cdt);

	point = cdt;
	// uint32_t half0, half1;
	while (point != NULL) {
		uint192_t* resultCDF = double2int(point->data, theta);
		if (resultCDF != NULL) {
			printf("%llu\n", resultCDF->bits[0]);
			//printf("bits 1: %llu\n", result->bits[1]);
			//printf("bits 2: %llu\n", result->bits[2]);
			// 分割为两个32位数据，需要根据theta值来判断
			//half0 = (uint32_t)(resultCDF->bits[0] << 34 >> 34);
			//half1 = (uint32_t)((resultCDF->bits[0] >> 30));
			//printf("half0: %lu half1: %lu\n", half0, half1);
		}
		else {
			printf("Invalid theta value\n");
		}
		point = point->next;
	}

	// 计算半高斯分布
	data s1;
	double sum1 = HalfGaussianDistribution(&s1, 1.6);
	print_list(s1.result);

	// 计算累加CDF
	Dnode* cdf1 = accumulateCDF_n(&s1);
	print_list(cdf1);

	point = cdf1;
	while (point != NULL) {
		uint192_t* resultCDF = double2int(point->data, theta);
		if (resultCDF != NULL) {
			printf("%lluu,\n", resultCDF->bits[0]);
			//printf("bits 1: %llu\n", result->bits[1]);
			//printf("bits 2: %llu\n", result->bits[2]);
			// 分割为两个32位数据，需要根据theta值来判断
			//half0 = (uint32_t)(resultCDF->bits[0] << 34 >> 34);
			//half1 = (uint32_t)((resultCDF->bits[0] >> 30));
			//printf("half0: %lu half1: %lu\n", half0, half1);
		}
		else {
			printf("Invalid theta value\n");
		}
		point = point->next;
	}

	// 计算指数分布
	//data s1;
	//double sum1 = expkDistribution(&s1, 30);
	//print_list(s1.result);

	//// 计算累加CDF
	//Dnode* expcdf = accumulateCDF_n(&s1);
	//print_list(expcdf);

	//Dnode* pt = expcdf;
	//// uint32_t half0, half1;
	//while (pt != NULL) {
	//	uint192_t* resultCDF = double2int(pt->data, 30);
	//	if (resultCDF != NULL) {
	//		printf("%llu\n", resultCDF->bits[0]);
	//	}
	//	else {
	//		printf("Invalid theta value\n");
	//	}
	//	pt = pt->next;
	//}

	//释放链表内存
	free_list(s.result);
	free_list(cdf);
	free_list(cdt);

	free_list(s1.result);
	free_list(cdf1);


	//free_list(s1.result);
	//free_list(expcdf);

	return 0;
}
