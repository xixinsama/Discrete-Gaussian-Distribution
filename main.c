#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include "sampler.h"
#include <time.h>
#include <stdlib.h>


int main() {
	sampler_shake256_context rng;
	sampler_context sc;

	sampler_shake256_init(&rng);
	sampler_shake256_inject(&rng, (const void*)"test sampler", 12);
	sampler_shake256_flip(&rng);
	Zf(prng_init)(&sc.p, &rng);
	
	sc.center = 0;
	sc.sigma = 1;
	// 定义函数指针
	int (*func_ptr)();

	int choice;
	printf("1. CDT反演法c=0, sigma=0.75 \n2. Knuth_Yao法c = 0, sigma = 4 \n3. 查表法c = 0, sigma = 0.75 \n4. 拒绝采样法c = 0, sigma = 0.75 \n5. 高斯卷积法c = 0, sigma = 1024 \n6. 伯努利拒绝采样法c = [0,1), sigma = 1.5 \n7. 伯努利拒绝采样法c = [0,1), sigma = (0.8,1.6) \n8. Karney算法 \n9. 实时计算的拒绝采样法 \n10. 并行优化的高斯卷积法 \n输入序号：\n  ");
	scanf("%d", &choice);

	// 根据用户的选择，设置函数指针
	if (choice == 1) {
		func_ptr = sampler_1_CDT;
	}
	else if (choice == 2) {
		func_ptr = sampler_1_KY;
	}
	else if (choice == 3) {
		func_ptr = sampler_1_LUT;
	}
	else if (choice == 4) {
		func_ptr = sampler_1_Reject;
	}
	else if (choice == 5) {
		func_ptr = sampler_2;
	}
	else if (choice == 6) {
		sc.sigma = 1.5; // 伯努利拒绝采样法
		printf("输入标准差[0.8, 1.6]：");
		scanf("%lf", &sc.sigma);
		func_ptr = sampler_3;
		printf("center = %f\n", sc.center);
		printf("sigma = %f\n", sc.sigma);
	}
	else if (choice == 7) {
		printf("输入中心[0, 1)：");
		scanf("%lf", &sc.center);
		printf("输入标准差[0.8, 1.6]：");
		scanf("%lf", &sc.sigma);
		func_ptr = sampler_4;
		printf("center = %f\n", sc.center);
		printf("sigma = %f\n", sc.sigma);
	}
	else if (choice == 8) {
		printf("输入中心[0, 1)：");
		scanf("%lf", &sc.center);
		printf("输入标准差[0.8, 1.6]：");
		scanf("%lf", &sc.sigma);
		func_ptr = sampler_karney; // Karney算法
		printf("center = %f\n", sc.center);
		printf("sigma = %f\n", sc.sigma);
	}
	else if (choice == 9) {
		printf("输入中心[0, 1)：");
		scanf("%lf", &sc.center);
		printf("输入标准差[0.8, 1.6]：");
		scanf("%lf", &sc.sigma);
		func_ptr = sampler_5; // Karney算法
		printf("center = %f\n", sc.center);
		printf("sigma = %f\n", sc.sigma);
	}
	else if (choice == 10) {
		func_ptr = sampler_2_Vector;
	}
	else {
		printf("无效的选择\n");
		return 0;
	}

	//要去测试一秒内采样结果的准确性，请把结果放入output.txt文件中
	clock_t start_t, end_t;

	// 初始化采样计数器
	int sample_count = 0;
	int max_samples = 10000000; // 假设最多采样10000000次
	int* sample_results = (int*)malloc(max_samples * sizeof(int)); // 用于存储采样结果

	// 检查内存分配是否成功
	if (sample_results == NULL) {
		printf("Memory allocation failed\n");
		return 1;
	}

	// 记录开始时间
	start_t = clock();
	// 采样循环
	while (1) {
		// 调用采样器函数
		int sample_result = func_ptr(&sc);

		// 将采样结果添加到列表中
		if (sample_count < max_samples) {
			*(sample_results + sample_count) = sample_result; //sample_results[sample_count] = sample_result;
			sample_count++;
		}

		// 检查是否已经过了1秒或达到最大采样次数
		end_t = clock();
		if (end_t - start_t > CLOCKS_PER_SEC || sample_count >= max_samples) {
			break;
		}
	}

	// 计算程序执行时间
	double cpu_time_used = ((double)(end_t - start_t)) / CLOCKS_PER_SEC;
	printf("程序执行耗费的CPU时间：%f秒\n", cpu_time_used);
	// 输出采样数量和结果列表
	printf("采样数量：%d\n", sample_count);

	// 打开文件
	FILE* file = fopen("output.txt", "w");
	if (file == NULL) {
		printf("无法打开文件\n");
		return 1;
	}

	for (int i = 0; i < sample_count; i++) {
		fprintf(file, "%d ", sample_results[i]);
	}
	// 关闭文件
	fclose(file);

	// 释放动态分配的内存
	free(sample_results);

	return 0;
}