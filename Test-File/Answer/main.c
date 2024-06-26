#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include "sampler.h"

int main() {
    sampler_shake256_context rng;

    printf("Test sampler: \n");
    fflush(stdout);

    sampler_shake256_init(&rng);
    sampler_shake256_inject(&rng, (const void *)"test sampler", 12);
    sampler_shake256_flip(&rng);

    int z;
    double mu = 0.99625;
    double sigma = 1.53984;
    
    z = sampler_4(&rng, sigma, mu); //改变为直接输入随机数指针

    printf("z: %d\n", z);
    
    clock_t start_t, end_t;
	// 初始化采样计数器
	int sample_count = 0;
	int max_samples = 100000;
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
		int sample_result = sampler_4(&rng, sigma, mu);

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
    printf("预计一秒内采样数：%d\n", (int)(sample_count/cpu_time_used));

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
