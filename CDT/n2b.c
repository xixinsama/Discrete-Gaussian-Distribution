#include "n2b.h"
#include <math.h>
#include <stdlib.h>

// 将十进制四精度浮点数转化为二进制字符串
void ldouble2bin(long double f, char* binstr) {
    // 将四精度浮点数转换为无符号长整型
    unsigned long long* ptr = (unsigned long long*) &f;
    unsigned long long bits = *ptr;

    // 将二进制位转换为字符串
    for (int i = 0; i < 128; i++) {
        binstr[127 - i] = (bits & 1) ? '1' : '0';
        bits >>= 1;
    }
    binstr[128] = '\0';
}

// 将十进制双精度浮点数转化为二进制字符串
void double2bin(double d, char* binstr) {
    uint64_t* p = (uint64_t*) &d;
    uint64_t value = *p;
    for (int i = 0; i < 64; i++) {
        binstr[63 - i] = (value & 1) ? '1' : '0';
        value >>= 1;
    }
    binstr[64] = '\0';
}

// 将二进制字符串转换为十进制整数
uint64_t bin2int(char* binstr, int strlens) {
    uint64_t result = 0;
    for (int i = 0; i < strlens - 1; i++) {
        result <<= 1;
        result += binstr[i] - '0';
    }
    return result;
}

// 将0-1之间的浮点数转换为整数
// theta为精度
// 2^theta * result
// 将如果超出64位，则将超出的部分放在bits[1]中，如果超出128位，则将超出的部分放在bits[2]中
uint192_t* double2int(double num, int theta) {
    uint192_t* result = (uint192_t*)malloc(sizeof(uint192_t));
    result->bits[0] = 0;
    result->bits[1] = 0;
    result->bits[2] = 0;
    if (theta > 192 || theta < 0) {
        return NULL;
    }
    else if(theta > 128) {
		for (int i = theta - 1; i >= 0; i--) {
			num *= 2;
			if (num > 1) {
				num -= 1;
				if (i >= 128) {
					result->bits[2] |= 1 << (i - 128);
				}
				else if (i >= 64) {
					result->bits[1] |= 1 << (i - 64);
				}
				else {
					result->bits[0] |= 1 << i;
				}
			}
		}
    }
    else if(theta > 64) {
        for (int i = theta - 1; i >= 0; i--) {
            num *= 2;
            if (num > 1) {
                num -= 1;
                if (i >= 64) {
                    result->bits[1] |= 1 << (i - 64);
                }
                else {
                    result->bits[0] |= 1 << i;
                }
            }
        }
	}
	else {
       double constValue = pow(2, theta);
       result->bits[0] = (uint64_t)(num * constValue);
	}
    return result;
}
