#include "sampler.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main() {
    sampler_shake256_context shake_rng;
    prng rng;

    // 初始化 SHAKE256 上下文
    i_shake256_init(&shake_rng);
    i_shake256_inject(&shake_rng, (const void*)"test sampler", 12); // 注入种子，不同环境下生成的随机数相同
    i_shake256_flip(&shake_rng);

    uint8_t random_bytes[4];
    // 基于 SHAKE256 算法生成随机数
    i_shake256_extract(&shake_rng, random_bytes, 4);
    uint64_t random_number = 0;
    for (int i = 0; i < 4; i++) {
        random_number |= ((uint64_t)random_bytes[i] << (8 * i));
    }
    printf("Random number: %llu\n", random_number);

    // 初始化 PRNG 上下文
    // 防止侧信道攻击，使用 SHAKE256 上下文作为种子
    prng_init(&rng, &shake_rng); //以 SHAKE256 上下文为种子初始化 PRNG
    prng_refill(&rng);

    uint8_t random_bytes2[4];
    // 从 PRNG 中获取随机数, 按字节提取，理论上更快
    prng_get_bytes(&rng, random_bytes2, 4);
    uint64_t random_number2 = 0;
    for (int i = 0; i < 4; i++) {
        random_number2 |= ((uint64_t)random_bytes2[i] << (8 * i));
    }
    printf("Random number2: %llu\n", random_number2);

    /*——————————————————————————————————————————————————————————*/
    
    sampler_context sc;
    sc.p = rng;
    sc.center = 0;
    sc.sigma = 1;
    // 定义函数指针
    int (*func_ptr)();

    int choice;
    printf(" 1. CDT查表法c=0, sigma=0.75 \n \
2. Knuth_Yao法c = 0, sigma = 4 \n \
3. 查表法c = 0, sigma = 0.75 \n \
4. 拒绝采样法c = 0, sigma = 0.75 \n \
5. 并行优化的CDT查表法（没有） \n \
以上为情形1(6000k) \n \
6. AVX2整数向量并行优化的高斯卷积法(1300k) \n \
7. 高斯卷积法c = 0, sigma = 1024 \n \
8. 改进Karney算法，对称中心，固定标准差（没有）(2600k) \n \
以上为情形2 \n \
9. 伯努利拒绝采样法c = [0,1), sigma = 1.5 \n \
以上为情形3(2900k) \n \
10. 伯努利拒绝采样法c = [0,1), sigma = (0.8,1.6)(2600k) \n \
11. Karney算法，任意中心，任意标准差>1（没有） \n \
12. 实时计算的拒绝采样法 \n \
以上为情形4 \n \
输入序号：\n  ");
    scanf_s("%d", &choice);
    // 11和12不需要预计算表，占用空间=0
    // 编译器优化，在vs直接测，第二次在虚拟机上测
    // 根据用户的选择，设置函数指针
    if (choice == 1) {
        func_ptr = sampler_1;
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
        func_ptr = sampler_1_vector;
    }
    else if (choice == 6) {
        func_ptr = sampler_2;
    }
    else if (choice == 7) {
        func_ptr = sampler_2_ori;
    }
    else if (choice == 8) {
        func_ptr = sampler_2_karney;
    }
    else if (choice == 9) {
        func_ptr = sampler_3;
    }
    else if (choice == 10) {
        func_ptr = sampler_4;
    }
	else if (choice == 11) {
		func_ptr = sampler_4_karney;
	}
	else if (choice == 12) {
		func_ptr = sampler_5;
	}
    else {
        printf("无效的选择\n");
        return 0;
    }

    // 获取输入参数
    if (choice == 9) {
        printf_s("输入中心值：");
        scanf_s("%lf", &sc.center);
    }
	else if (choice == 10 || choice == 11 || choice == 12) {
	    printf_s("输入中心值：");
		scanf_s("%lf", &sc.center);
		printf_s("输入标准差：");
		scanf_s("%lf", &sc.sigma);
	}
    else if (choice == 1 || choice == 3 || choice == 4 || choice == 5) {
        sc.center = 0;
		sc.sigma = 0.75;
	}
	else if (choice == 2) {
		sc.center = 0;
		sc.sigma = 4;
	}
	else if (choice == 6 || choice == 7 || choice == 8) {
		sc.center = 0;
		sc.sigma = 1024;
	}
	/*——————————————————————————————————————————————————————————*/

    clock_t start_t, end_t;
    int sample_count = 0; 	// 初始化采样计数器
    int max_samples = 100000; // 最大采样次数
    int *sample_results = (int *)malloc(max_samples * sizeof(int)); // 用于存储采样结果

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
        sample_results[sample_count] = sample_result;
        sample_count++;

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
    printf("预计一秒内采样数：%d\n", (int)(sample_count / cpu_time_used));

    // 打开文件
    FILE* file = fopen("output.txt", "w");
    if (file == NULL) {
        printf("无法打开文件\n");
        return 1;
    }

    // 将center, sigma写入文件
    fprintf(file, "%lf\n", sc.center);
    fprintf(file, " %lf\n", sc.sigma);
    if (sample_results != NULL && sample_count > 0) {
        for (int i = 0; i < sample_count; i++) {
            fprintf(file, "%d ", sample_results[i]);
        }
    } else {
        printf("Invalid sample results or count\n");
    }

    // 关闭文件
    fclose(file);

    // 释放动态分配的内存
    free(sample_results);
    

    return 0;
}
