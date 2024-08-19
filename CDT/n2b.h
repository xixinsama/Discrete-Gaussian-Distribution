// �����������͵�ת��
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

// ��ʮ�����ľ��ȸ�����ת��Ϊ�������ַ���
void ldouble2bin(long double f, char* binstr);

// ��ʮ����˫���ȸ�����ת��Ϊ�������ַ���
void double2bin(double d, char* binstr);

// ���������ַ���ת��Ϊʮ��������
uint64_t bin2int(char* binstr, int strlens);

// ��0-1֮��ĸ�����ת��Ϊ����
// ���������û����ô�߾��ȣ����λȫΪ1�������޷������ֵ
// ������ڵ���1���򷵻�ֵ������ת�������н��������
uint192_t* double2int(double num, int theta);


#endif N2B