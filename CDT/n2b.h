// 各种数据类型的转换
#ifndef N2B
#define N2B

#include <stdint.h>

typedef union {
    double d;
    uint64_t bits;
} DoubleBits;

typedef struct {
	uint64_t bits[3];
} uint192_t;

// 将十进制四精度浮点数转化为二进制字符串
void ldouble2bin(long double f, char* binstr);

// 将十进制双精度浮点数转化为二进制字符串
void double2bin(double d, char* binstr);

// 将二进制字符串转换为十进制整数
uint64_t bin2int(char* binstr, int strlens);

// 将0-1之间的浮点数转换为整数
// 如果浮点数没有那么高精度，则低位全为1，出现无符号最大值
// 如果大于等于1，则返回值不符合转换规则，有结果但无用
uint192_t* double2int(double num, int theta);


#endif N2B