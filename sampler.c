#include "sampler.h"
#include <immintrin.h>
#include <math.h> //快很多

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

 //c = 0, sigma = 0.75, 28bit
static const uint32_t DistForSampler_1_CDT[] = {
142782702u,
260182150u,
268339469u,
268435265u,
268435455u
};

// z, sigma = 2.5582018962155022023807759978808462619781494140625, 28bits
static const uint32_t DistForSampler_2_CDT[] = {
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

//c = 0, sigma = 1.5, z+, 30bits
static const uint32_t DistForBSampler_3_CDT[] = {
451157485u,
812416158u,
997892545u,
1058950071u,
1071837610u,
1073581749u,
1073733096u,
1073741516u,
1073741817u,
1073741823u
};

//c = 0, sigma = 0.9--1.6, z+, 30bits
static const uint32_t CDT4_16[] = {
428587674u,
781134279u,
977356015u,
1051253798u,
1070084626u,
1073331468u,
1073710265u,
1073740167u,
1073741764u,
1073741822u,
1073741823u
};
static const uint32_t CDT4_15[] = {
451157485u,
812416158u,
997892545u,
1058950071u,
1071837610u,
1073581749u,
1073733096u,
1073741516u,
1073741817u,
1073741823u
};
static const uint32_t CDT4_14[] = {
476236526u,
845242412u,
1016900815u,
1064843023u,
1072881840u,
1073691099u,
1073740010u,
1073741784u,
1073741823u
};
static const uint32_t CDT4_13[] = {
504267880u,
879389259u,
1033809185u,
1068985948u,
1073420293u,
1073729624u,
1073741565u,
1073741820u,
1073741823u
};
static const uint32_t CDT4_12[] = {
535805472u,
914431486u,
1048035764u,
1071577413u,
1073648795u,
1073739805u,
1073741802u,
1073741823u
};
static const uint32_t CDT4_11[] = {
571551042u,
949640434u,
1059089712u,
1072954405u,
1073722981u,
1073741625u,
1073741823u
};
static const uint32_t CDT4_10[] = {
612406982u,
983850593u,
1066730865u,
1073534092u,
1073739532u,
1073741814u,
1073741823u
};
static const uint32_t CDT4_09[] = {
659553470u,
1015321563u,
1071158033u,
1073707814u,
1073741692u,
1073741823u
};

static const double sigma_minT[] = {
0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5
};

// 1 / (2 * sigma_max^2)
static const double isigma_maxT[] = {
0.617283950617284027373443677789,
0.500000000000000000000000000000,
0.413223140495867724553136213217,
0.347222222222222154375259606240,
0.295857988165680318992656339105,
0.255102040816326369743194391049,
0.222222222222222126619683990612,
0.195312499999999861222121921855
};

// Fixed sigma = 0.75 and center = 0
int sampler_1(void* ctx) {
	sampler_context* sc = ctx; //统一变量名为sc
	int z = 0;
	uint32_t r = prng_get_u32(&sc->p);
	int s = (int)r & 1; //sign
	r = r >> 4;
	while ((DistForSampler_1_CDT[z] - r) >> 31) //边界为dist最后一个
	{
		z = z + 1;
	}
	return z = s == 1 ? -z : z;
}

// AVX2降频CPU
// Fixed sigma = 1024 and center = 0
static inline void BaseSampler_Vector(prng* p, __m256i* z_out)
{
	const __m256i v_zero = _mm256_setzero_si256();
	const __m256i v_one = _mm256_set1_epi32(1);

	__m256i v_r = _mm256_set_epi64x(prng_get_u64(p), prng_get_u64(p), prng_get_u64(p), prng_get_u64(p));
	__m256i v_r_shifted = _mm256_srli_epi32(v_r, 4);
	__m256i v_z = v_zero;

	for (size_t k = 0; k < sizeof(DistForSampler_2_CDT) / sizeof(DistForSampler_2_CDT[0]); k++) {
		__m256i v_dist = _mm256_set1_epi32(DistForSampler_2_CDT[k]);
		__m256i mask = _mm256_cmpgt_epi32(v_r_shifted, v_dist);
		v_z = _mm256_add_epi32(v_z, _mm256_and_si256(mask, v_one));
		if (_mm256_testz_si256(mask, _mm256_set1_epi32(-1))) { break; }
	}

	__m256i isHighestBitZero = _mm256_cmpeq_epi32(_mm256_and_si256(v_r, v_one), v_zero);
	v_z = _mm256_blendv_epi8(v_z, _mm256_sub_epi32(v_zero, v_z), isHighestBitZero);

	*z_out = v_z;
}

int sampler_2(void* ctx)
{
	sampler_context* sc = ctx;
	int z = 0;

	__m256i v_z;
	BaseSampler_Vector(&sc->p, &v_z);
	__m256i p_z = _mm256_set_epi32(384, 96, 48, 12, 32, 8, 4, 1);
	__m256i product = _mm256_mullo_epi32(p_z, v_z); // 两个向量的逐元素相乘

	int* i = (int*)&product;
	for (int j = 0; j < 8; j++)
	{
		z = z + i[j];
	}
	return z;
}

static inline int BaseSampler3(prng* p)
{
	int z = 0;
	uint32_t r = prng_get_u32(p) >> 2;
	while ((DistForBSampler_3_CDT[z] - r) >> 31)
	{
		z = z + 1;
	}
	return z;
}

// Fixed sigma = 1.5 and center c is uniformly distributed over [0,1)
int sampler_3(void* ctx) {
	sampler_context* sc = ctx;
	double sigma = 1.5;
	int z = 0;
	double isigma = 1 / (2 * sigma * sigma);

	while (1)
	{
		int z0 = BaseSampler3(&sc->p);
		int b = (int)prng_get_u8(&sc->p) & 1;
		z = (2 * b - 1) * z0 + b;
		double x = (z - sc->center) * (z - sc->center);
		x = x - z0 * z0;
		x = x * isigma;
		double p = expm_p63(x);
		int i = 1;
		uint8_t u, v;

		do {
			i = i * 0xff;
			u = prng_get_u8(&sc->p);
			v = (int)(p * i) & 0xff;
		} while (u == v);

		if (u < v)
		{
			return z;
		}
	}
}

static inline int BaseSampler4(prng* p, int mark)
{
	int z = 0;
	uint32_t r = prng_get_u32(p) >> 2; // 30bit
	const uint32_t* DistForBSampler4; // 指向DistForBSampler的指针

	switch (mark) {
	case 0: DistForBSampler4 = CDT4_09; break;
	case 1: DistForBSampler4 = CDT4_10; break;
	case 2: DistForBSampler4 = CDT4_11; break;
	case 3: DistForBSampler4 = CDT4_12; break;
	case 4: DistForBSampler4 = CDT4_13; break;
	case 5: DistForBSampler4 = CDT4_14; break;
	case 6: DistForBSampler4 = CDT4_15; break;
	case 7: DistForBSampler4 = CDT4_16; break;
	default: DistForBSampler4 = CDT4_16;
	}

	while ((DistForBSampler4[z] - r) >> 31) {
		z = z + 1;
	}
	return z;
}

// 接受采样，返回0或1
static int AcceptSample(prng* pp, double sis, double x)
{
	double p = sis * expm_p63(-x);

	int i = 1;
	uint16_t u, v;
	do {
		i = i * 0xff;
		u = prng_get_u8(pp);
		v = (int)(p * i) & 0xff;
	} while (u == v);
	return u < v;
}

// sigma is uniformly distributed over (0.8,1.6) and center is uniformly distributed over [0,1)
int sampler_4(void* ctx) {
	sampler_context* sc = ctx;
	int z = 0;

	int mark = 7;
	if (sc->sigma > 0.8 && sc->sigma <= 0.9) { mark = 0; }
	else if (sc->sigma > 0.9 && sc->sigma <= 1.0) { mark = 1; }
	else if (sc->sigma > 1.0 && sc->sigma <= 1.1) { mark = 2; }
	else if (sc->sigma > 1.1 && sc->sigma <= 1.2) { mark = 3; }
	else if (sc->sigma > 1.2 && sc->sigma <= 1.3) { mark = 4; }
	else if (sc->sigma > 1.3 && sc->sigma <= 1.4) { mark = 5; }
	else if (sc->sigma > 1.4 && sc->sigma <= 1.5) { mark = 6; }
	else { mark = 7; }

	double isigma = 1 / sc->sigma;
	double sis = sigma_minT[mark] * isigma;

	while (1)
	{
		int z0 = BaseSampler4(&sc->p, mark);
		int b = prng_get_u8(&sc->p) >> 7;
		z = (2 * b - 1) * z0 + b;
		double x = z0 * z0 * isigma_maxT[mark];
		x = x - (z - sc->center) * (z - sc->center) * isigma * isigma * 0.5;
		if (AcceptSample(&sc->p, sis, x))
		{
			return z;
		}
	}
}

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

// Fixed sigma = 4 and center = 0
int sampler_1_KY(void* ctx) {
	sampler_context* sc = ctx;

	int z = 0, d = 0, hit = 0, col = 0;

	while (hit == 0)
	{
		int r0 = prng_get_u8(&sc->p) & 1;
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

	int s = prng_get_u8(&sc->p) & 1;
	return z = s == 0 ? -z : z;
}

// Fixed sigma = 0.75 and center = 0
int sampler_1_LUT(void* ctx)
{
	sampler_context* sc = ctx;
	int z = 0;
	uint32_t r = prng_get_u32(&sc->p);
	int s = (int)r & 1;
	r = r >> 4;
	int length = sizeof(DistForSampler_1_CDT) / sizeof(DistForSampler_1_CDT[0]);

	for (int i = 0; i < length; i++)
	{
		if ((r - DistForSampler_1_CDT[i]) >> 31)
		{
			z = i;
			break;
		}
	}

	return z = (s == 0) ? -z : z;
}

// sigma = 0.75, center = 0, len = 7
static double DistForSampler_1_Reject[] = {
0.53192304053524364082505826445413,
0.43736019913598300767887394613354,
0.030389296063459844687848487865267,
0.00035688060203969431370685239102158,
0.00000070834717560584336899825097164096,
0.00000000023762400499605740608216507442384,
0.000000000000013472722889431712901342109797289,
0.00000000000000000012910449524651443744064140859293
};

// Fixed sigma = 0.75 and center = 0
int sampler_1_Reject(void* ctx)
{
	int z = 0;
	sampler_context* sc = ctx;
	uint8_t u, v;

	while (1)
	{
		int i = 1;
		int k = (int)prng_get_u8(&sc->p) & 7; //取0-7的均匀随机数
		do
		{
			i = i * 0xff; //8bit, 全1, 256
			u = (uint8_t)prng_get_u8(&sc->p); //0-256的均匀随机数
			v = (uint8_t)(DistForSampler_1_Reject[k] * i) & 0xff;
		} while (u == v);
		if (u < v) //不知道为什么不能使用位移比较大小，下同
		{
			int s = prng_get_u8(&sc->p) & 1;
			return z = s == 0 ? -k : k;
			break;
		}
	}
}

// 向量化采样器1，sigma = 0.75, center = 0，一次输出8个结果
int sampler_1_vector(void* ctx)  {
	sampler_context* sc = ctx;
}

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

static int BaseSampler2_CDT(prng* p)
{
	int z = 0;
	uint32_t r = prng_get_u32(p);
	int s = (int)r & 1; //符号位
	r = r >> 4;

	while ((DistForBaseSampler_CDT[z] - r) >> 31) //未设置边界，在dist中需要保证最后一个数为max
	{
		z = z + 1;
	}

	return z = s == 1 ? -z : z;
}

// Fixed sigma = 1024 and center = 0
int sampler_2_ori(void* ctx) {
	sampler_context* sc = ctx;
	int z = 0;

	int x1, x2, x3, x4, x5, x6, x7, x8;
	x1 = BaseSampler2_CDT(&sc->p);
	x2 = BaseSampler2_CDT(&sc->p);
	x3 = BaseSampler2_CDT(&sc->p);
	x4 = BaseSampler2_CDT(&sc->p);
	x5 = BaseSampler2_CDT(&sc->p);
	x6 = BaseSampler2_CDT(&sc->p);
	x7 = BaseSampler2_CDT(&sc->p);
	x8 = BaseSampler2_CDT(&sc->p);

	z = x1 + 4 * x2 + 8 * x3 + 32 * x4 + 12 * x5 + 48 * x6 + 96 * x7 + 384 * x8;

	return z;
}

// 改进Karney算法，中心固定
int sampler_2_karney(void* ctx) {
	sampler_context* sc = ctx;
}

//对任意标准差(> 1)，任意中心的Karney算法
int sampler_4_karney(void* ctx)
{
	sampler_context* sc = ctx;
	int z = 0;

	return z;
}

//对任意标准差，任意中心的实时计算拒绝采样法
int sampler_5(void* ctx)
{
	int z = 0;
	sampler_context* sc = ctx;

	uint8_t u, v;
	double x, p;
	double isigma = 1 / sc->sigma;
	int cut = (int)(3 * sc->sigma) + 1; //截断范围
	int i, k = 0;

	while (1)
	{
		do
		{
			k = (int)(prng_get_u8(&sc->p) & 15);
		} while (k > 2 * cut);
		k = k - cut;

		x = 0.5 * (k - sc->center) * (k - sc->center) * isigma * isigma;
		// p = exp(-x); // 需使用math.h库
		p = expm_p63(x);

		i = 1;
		do
		{
			i = i * 0xff;
			u = (uint8_t)prng_get_u8(&sc->p);
			v = (uint8_t)(p * i) & 0xff;
		} while (u == v);

		if (u < v)
		{
			return z = k;
		}
	}
}