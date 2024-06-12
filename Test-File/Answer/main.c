#include <stdio.h>
#include <malloc.h>
#include "sampler.h"

int main() {
    sampler_shake256_context rng;

    printf("Test sampler: \n");
    fflush(stdout);

    sampler_shake256_init(&rng);
    sampler_shake256_inject(&rng, (const void *)"test sampler", 12);
    sampler_shake256_flip(&rng);

    int z;
    double mu = 0.5;
    double sigma = 1.5;
    
    struct mallinfo before = mallinfo();
    z = sampler_3(&rng, mu); //改变为直接输入随机数指针
    struct mallinfo after = mallinfo();

    // 计算差异
    int allocated_diff = after.uordblks - before.uordblks;
    printf("该函数的堆内存分配: %d bytes\n", allocated_diff);

    printf("z: %d\n", z);
    
    return 0;
}
