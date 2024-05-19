#include "sampler.h"

/*
 *
 * ********************************Attention!!!************************************************
 * ctx : Pointer to the sampler_context structure, including PRNG or other.
 * ********************************************************************************************
 * Please read the following precautions carefully!!!!
 * 1. In supplied folder "example/sampler", it offers a falcon's sampler as an instance for giving the codes about pointer "ctx".
 * 2. In folder "example/sampler", it includes the implementation of SHAKE256 in shake.c and PRNG in rng.c, please just use it.
 * 3. In the head file "example/sampler/sampler.h", it gives the definition of "ctx" in the struct "sampler_context", please just use it and don't change it.
 * 4. In the main file "example/sampler/main.c", it provides the procedure of how to initialize the pointer ctx, you should learn it and use it.
 * 5. The API of sampler function in the "example/sampler" is similar to ours. Pointer ctx is the same and fpr also denotes the type of "double", but we use "double" here.
 *
 *
 *In conclusion, you should use the method of defining and initializing for pointer ctx in given "example/sampler",
 *please don't change it and just use it. You should focus on the four samplers themselves. Pointer ctx should be the same for everyone.
 *
 *
 * */

 /*******************Pre-Computation Table********************/
  //c = 0, sigma = 0.75, 28bit
static const uint32_t DistForSampler_1_CDT[] = {
142782702u,
260182150u,
268339469u,
268435265u,
268435455u
};

// 32bit sigma = 0.75
static const uint32_t DistForSampler_1_KY[] = {
21012607352u,
1878391165u,
130517100u,
1532744u,
3042u,
1u
};

// 32bit sigma = 4
static const uint32_t DistForSampler_2_KY[] = {
428361011u,
830363458u,
756054532u,
646687711u,
519628174u,
392235924u,
278136918u,
185279131u,
115944717u,
68160383u,
37641738u,
19528276u,
9517321u,
4357347u,
1874071u,
757193u,
287398u,
102474u,
34324u,
10800u,
3192u,
3192u,
};

// z, sigma = 2.5582018962155022023807759978808462619781494140625
static const uint32_t DistForBaseSampler_CDT[] = {
41861532u,
119426299u,
181103206u,
223197130u,
247854967u,
260252329u,
265602187u,
267583687u,
268213606u,
268385481u,
268425733u,
268433824u,
268435219u,
268435426u,
268435452u,
268435455u
};

//c = 0, sigma = 1.5, z+, 28bits
static const uint32_t DistForBSampler_3_CDT[] = {
71393501u,
185728396u,
244429888u,
263754008u,
267832790u,
268384793u,
268432693u,
268435358u,
268435453u,
268435455u,
};

//c = 0, sigma = 0.8, z+, 28bits
static const uint32_t DistForBSampler_4_CDT[] = {
133861942u,
256434869u,
268197835u,
268434457u,
268435455u,
};

// 64bits sigma = 0.75
static const uint64_t DistForSampler_1_LUT_64[] = {
0x882b0eea011c0800u,
0xf82108677a9c6000u,
0xffe890d4670d3800u,
0xfffff41cbe3a7800u,
0xfffffffefab90800u,
0xfffffffffffc3800u,
0xffffffffffffffffu
};

static const double DistForSampler_1_Reject[] = {
0.53190701687809405218132496884209103882312774658203125,
0.437347024091196912021217713117948733270168304443359375,
0.0303883806158780911399475144207826815545558929443359375,
0.0003568698513634223075212392739530287144589237868785858154296875,
0.0000007083258373455691967227575156496044428422464989125728607177734375,
0.000000000237616846807191745526457224069378539044183895612150081433355808258056640625,
0.00000000000001347231703693802467800338413244726765692106662530846961089991964399814605712890625,
0.00000000000000000012910060610088588833092745214697578012021821598580804957290268930591992102563381195068359375
};
/*******************Pre-Computation Table********************/

// Fixed sigma = 0.75 and center = 0
int sampler_1_CDT(void* ctx) {
	sampler_context* spc;
	spc = ctx;
	int z = 0;
	uint32_t r = prng_get_u32(&spc->p) >> 4;

	while ((DistForSampler_1_CDT[z] - r) >> 31) //未设置边界，在dist中需要保证最后一个数为max
	{
		z = z + 1;
	}

	int s = prng_get_u8(&spc->p) & 1;
	return z = s == 1 ? -z : z;
}

// Fixed sigma = 4 and center = 0
int sampler_1_KY(void* ctx) {
	sampler_context* spc;
	spc = ctx;

	int z = 0, d = 0, hit = 0, col = 0;

	while (hit == 0)
	{
		int r0 = prng_get_u8(&spc->p) & 1;
		d = 2 * d + r0;
		for (int i = (sizeof DistForSampler_2_KY) / sizeof(DistForSampler_2_KY[0]) - 1; i >= 0; i--)
		{
			int pij = (int)(DistForSampler_2_KY[i] << col >> 31);
			d = d - pij;
			if (d == -1)
			{
				z = i;
				hit = 1;
				break;
			}
		}
		col++;
		if (col == 32)
		{
			z = 0;
			break;
		}
	}

	int s = prng_get_u8(&spc->p) & 1;

	return z = s == 0 ? -z : z;
}

// Fixed sigma = 0.75 and center = 0
int sampler_1_LUT(void* ctx)
{
	sampler_context* spc;
	spc = ctx;
	int z = 0;
	uint64_t r = prng_get_u64(&spc->p);
	
	int length = sizeof(DistForSampler_1_LUT_64) / sizeof(DistForSampler_1_LUT_64[0]);

	for (int i = 0; i <= length; i++)
	{
		if (r < DistForSampler_1_LUT_64[i])
		{
			z = i;
		}
	}

	int s = prng_get_u8(&spc->p) & 1;
	return z = (s == 0) ? -z : z;
}


// Fixed sigma = 0.75 and center = 0
int sampler_1_Reject(void* ctx)
{
	int z = 0;
	sampler_context* spc;
	spc = ctx;

	while (1)
	{
		int t = prng_get_u8(&spc->p) & 7;
		double u = prng_get_rand(&spc->p);

		if (u < DistForSampler_1_Reject[t])
		{
			int s = prng_get_u8(&spc->p) & 1;
			return z = (s == 0) ? -t : t;
			break;
		}
		else
		{
			continue;
		}
	}
}

int BaseSampler2_CDT(prng* p)
{
	int z = 0;
	uint32_t r = prng_get_u32(p) >> 4; // 28 bit
	int temp = 0;

	while ((DistForBaseSampler_CDT[z] - r) >> 31) //未设置边界，在dist中需要保证最后一个数为max
	{
		z = z + 1;
	}

	int s = prng_get_u8(p) & 1;
	return z = s == 1 ? -z : z;
}

// Fixed sigma = 1024 and center = 0
int sampler_2(void* ctx) {
	sampler_context* spc;
	spc = ctx;
	int z = 0;

	int x1, x2, x3, x4, x5, x6, x7, x8;
	x1 = BaseSampler2_CDT(&spc->p);
	x2 = BaseSampler2_CDT(&spc->p);
	x3 = BaseSampler2_CDT(&spc->p);
	x4 = BaseSampler2_CDT(&spc->p);
	x5 = BaseSampler2_CDT(&spc->p);
	x6 = BaseSampler2_CDT(&spc->p);
	x7 = BaseSampler2_CDT(&spc->p);
	x8 = BaseSampler2_CDT(&spc->p);

	z = x1 + 4 * x2 + 8 * x3 + 32 * x4 + 12 * x5 + 48 * x6 + 96 * x7 + 384 * x8;

	return z;
}

// 接受采样，返回0或1
int AcceptSample(prng* pp, double sis, double x)
{
	double p = sis * expm_p63(-x);;
	
	int i = 1;
	uint32_t u, v;

	do {
		i = i * 0xFF;
		u = prng_get_u8(pp);
		v = (int)(p * i) & 0xFF; //强制类型转换，用于向下取整
	} while (u == v);

	return (int)((u - v) >> 31);
}

// Fixed sigma = 1.5 and center c is uniformly distributed over [0,1)
int sampler_3(void* ctx) {
	double sigma = 1.5;

	int z = 0;


	return z;

}

int BaseSampler4(prng* p)
{
	int z = 0;
	uint32_t r = prng_get_u32(p) >> 4; // 32bit
	int temp = 0;

	while ((DistForBSampler_4_CDT[z] - r) >> 31)
	{
		z = z + 1;
	}

	return z;
}

// sigma is uniformly distributed over (0.8,1.6) and center is uniformly distributed over [0,1)
int sampler_4(void* ctx) {
	sampler_context* spc;
	spc = ctx;
	int z = 0;

	double isigma_min = 1 / spc->sigma_min; // 1 / sigma_min
	double isigma = 1 / spc->sigma; // 1 / sigma
	double sis = spc->sigma_min * isigma; // sigma_min / sigma

	while (1)
	{
		int z0 = BaseSampler4(&spc->p);
		int b = prng_get_u8(&spc->p) & 1;
		z = (2 * b - 1) * z0 + b;
		double x = (z0 * z0 * isigma_min * isigma_min) / 8 - (z - spc->center) * (z - spc->center) * isigma * isigma / 2;
		if (AcceptSample(&spc->p, sis, x))
		{
			return z;
			break;
		}
	}
}
