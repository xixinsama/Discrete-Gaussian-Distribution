#include "sampler.h"
#include <immintrin.h>
#include <math.h>

// 测试用函数
// 8 / 16 / 24 / 32 / 64
// 递增数列以全1最大数为界，速度最优，因为概率以此递减
// 递减数列以0为界，减少了一个静态变量内存，速度降低

// Fixed sigma = 0.75 and center = 0, 22bits
static const uint32_t 
Dist1[] = {
2230979u,
4065346u,
4192804u,
4194301u,
4194303u
};
// 16bits
static const uint16_t 
Dist2[] = {
34859u,
63521u,
65512u,
65535u
};
int s1(void* ctx) {
	sampler_context* sc = ctx;
	int z = 0;
	uint16_t r = prng_get_u16(&sc->p);
	int s = (int)r & 1; //sign
	// r = r >> 1;
	while (Dist2[z] < r) //(Dist2[z] - r) >> 15
	{
		z++;
	}
	return z = s == 1 ? -z : z;
}

// 14bits
static const uint16_t
Dist3[] = {
17429u,
31760u,
32756u,
32767u
};
// 向量化版本，一次输出16个样本，必须留出一个符号位
int s1_v(void* ctx, int* samples) {
	sampler_context* sc = ctx;
	const __m256i v_zero = _mm256_setzero_si256();
	const __m256i v_one = _mm256_set1_epi16(1);

	__m256i v_r = _mm256_set_epi64x(prng_get_u64(&sc->p), prng_get_u64(&sc->p), prng_get_u64(&sc->p), prng_get_u64(&sc->p));
	__m256i v_r_shifted = _mm256_srli_epi16(v_r, 1);
	__m256i v_z = v_zero;

	for (size_t k = 0; k < sizeof(Dist3) / sizeof(Dist3[0]); k++) {
		__m256i v_dist = _mm256_set1_epi16(Dist3[k]);
		__m256i mask = _mm256_cmpgt_epi16(v_r_shifted, v_dist);
		v_z = _mm256_add_epi16(v_z, _mm256_and_si256(mask, v_one));
		if (_mm256_testz_si256(mask, _mm256_set1_epi16(-1))) { break; }
	}

	__m256i isHighestBitZero = _mm256_cmpeq_epi16(_mm256_and_si256(v_r, v_one), v_zero);
	v_z = _mm256_blendv_epi8(v_z, _mm256_sub_epi16(v_zero, v_z), isHighestBitZero);

	_mm256_storeu_si256((__m256i*)samples, v_z);
	return 0;
}

// sigma = 2.5582018962155022023807759978808462619781494140625, 15bits
static const uint16_t Dist4[] = {
5110u,
14578u,
22107u,
27245u,
30255u,
31769u,
32422u,
32664u,
32740u,
32761u,
32767u
};
static void bs2_v(prng* p, __m256i* z_out) {
	const __m256i v_zero = _mm256_setzero_si256();
	const __m256i v_one = _mm256_set1_epi16(1);

	__m256i v_r = _mm256_set_epi64x(prng_get_u64(p), prng_get_u64(p), prng_get_u64(p), prng_get_u64(p));
	__m256i v_r_shifted = _mm256_srli_epi16(v_r, 1);
	__m256i v_z = v_zero;

	for (size_t k = 0; k < sizeof(Dist4) / sizeof(Dist4[0]); k++) {
		__m256i v_dist = _mm256_set1_epi16(Dist4[k]);
		__m256i mask = _mm256_cmpgt_epi16(v_r_shifted, v_dist);
		v_z = _mm256_add_epi16(v_z, _mm256_and_si256(mask, v_one));
		if (_mm256_testz_si256(mask, _mm256_set1_epi16(-1))) { break; }
	}

	__m256i isHighestBitZero = _mm256_cmpeq_epi16(_mm256_and_si256(v_r, v_one), v_zero);
	v_z = _mm256_blendv_epi8(v_z, _mm256_sub_epi16(v_zero, v_z), isHighestBitZero);

	*z_out = v_z;
}
int s2_v(void* ctx, int16_t* samples) {
	sampler_context* sc = ctx;
	__m256i v_z;
	bs2_v(&sc->p, &v_z);
	__m256i p_z = _mm256_set_epi16(384, 96, 48, 12, 32, 8, 4, 1, 384, 96, 48, 12, 32, 8, 4, 1);
	__m256i product = _mm256_mullo_epi16(p_z, v_z); // 两个向量的逐元素相乘

	int16_t* i = (int16_t*)&product;
	samples[0] = samples[1] = 0;
	samples[0] += i[0];
	samples[0] += i[1];
	samples[0] += i[2];
	samples[0] += i[3];
	samples[0] += i[4];
	samples[0] += i[5];
	samples[0] += i[6];
	samples[0] += i[7];
	samples[1] += i[8];
	samples[1] += i[9];
	samples[1] += i[10];
	samples[1] += i[11];
	samples[1] += i[12];
	samples[1] += i[13];
	samples[1] += i[14];
	samples[1] += i[15];

	return 0;
}

// sigma = 1.5, 16bits
static const uint16_t
Dist5[] = {
27536u,
49585u,
60906u,
64633u,
65419u,
65526u,
65535u
};
static const uint16_t
Dist6[] = {
13768u,
24792u,
30453u,
32316u,
32709u,
32763u,
32767u
};
//c = 0, sigma = 1.5, z+, 30bits
static const uint32_t 
Dist7[] = {
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
static void bs3_v(prng* p, __m256i* z_out) {
	const __m256i v_zero = _mm256_setzero_si256();
	//const __m256i v_one = _mm256_set1_epi16(1);
	const __m256i v_one = _mm256_set1_epi16(1);

	__m256i v_r = _mm256_set_epi64x(prng_get_u64(p), prng_get_u64(p), prng_get_u64(p), prng_get_u64(p));
	__m256i v_r_shifted = _mm256_srli_epi16(v_r, 1);
	__m256i v_z = v_zero;

	for (size_t k = 0; k < sizeof(Dist6) / sizeof(Dist6[0]); k++) {
		__m256i v_dist = _mm256_set1_epi16(Dist6[k]);
		__m256i mask = _mm256_cmpgt_epi16(v_r_shifted, v_dist); // 生成遮罩
		v_z = _mm256_add_epi16(v_z, _mm256_and_si256(mask, v_one));
		if (_mm256_testz_si256(mask, _mm256_set1_epi16(-1))) { break; }
	}

	*z_out = v_z;
}
// Fixed sigma = 1.5 and center c is uniformly distributed over [0,1)
// 返回成功采样的个数，采样结果存储在samples中
int s3_v(void* ctx, int* samples) {
	sampler_context* sc = ctx;
	double isigma = 2.0 / 9;
	int z[16] = { 0 };

	int16_t sample[16] = { 0 };
	bs3_v(&sc->p, sample); // 生成16个样本
	int n = 0; // zs的索引
	int	m = 0; // sample的索引
	double x, p;
	int i, b;
	uint8_t u, v; // 惰性浮点
	
	while (1) {
		i = 0xff;
		u = prng_get_u8(&sc->p);
		b = (int)u & 1;// 0 or 1
		z[n] = (2 * b - 1) * sample[m] + b;
		x = (z[n] - sc->center) * (z[n] - sc->center);
		x = x - sample[m] * sample[m];
		x = x * isigma;
		p = expm_p63(x);
		v = (int)(p * i) & 0xff;

		while (u == v) {
			i = i * 0xff;
			u = prng_get_u8(&sc->p);
			v = (int)(p * i) & 0xff;
		}

		if (u < v) {
			samples[n] = z[n];
			n++;
			m++;
		}
		else { m++; }
		if (m == 16) { break; }
	}
	
	return n;
}

static const double 
sigma_minT[] = {
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
// sigma = 0.9, 15bits
static const uint16_t
CDT4_09[] = {
20127u,
30985u,
32689u,
32766u,
32767u
};
static const uint16_t
CDT4_10[] = {
18689u,
30024u,
32554u,
32761u,
32767u
};
static const uint16_t
CDT4_11[] = {
17442u,
28980u,
32320u,
32743u,
32767u
};
static const uint16_t
CDT4_12[] = {
16351u,
27906u,
31983u,
32701u,
32765u,
32767u
};
static const uint16_t
CDT4_13[] = {
15389u,
26836u,
31549u,
32622u,
32758u,
32767u
};
static const uint16_t
CDT4_14[] = {
14533u,
25794u,
31033u,
32496u,
32741u,
32766u,
32767u
};
static const uint16_t
CDT4_15[] = {
13768u,
24792u,
30453u,
32316u,
32709u,
32763u,
32767u
};
static const uint16_t
CDT4_16[] = {
13079u,
23838u,
29826u,
32081u,
32656u,
32755u,
32767u
};
static void 
bs4_v(prng* p, __m256i* z_out, int mark) {
	const uint16_t* Dist_4; // 指向DistForBSampler的指针

	switch (mark) {
	case 0: Dist_4 = CDT4_09; break;
	case 1: Dist_4 = CDT4_10; break;
	case 2: Dist_4 = CDT4_11; break;
	case 3: Dist_4 = CDT4_12; break;
	case 4: Dist_4 = CDT4_13; break;
	case 5: Dist_4 = CDT4_14; break;
	case 6: Dist_4 = CDT4_15; break;
	case 7: Dist_4 = CDT4_16; break;
	default: Dist_4 = CDT4_16;
	}

	const __m256i v_zero = _mm256_setzero_si256();
	const __m256i v_one = _mm256_set1_epi16(1);

	__m256i v_r = _mm256_set_epi64x(prng_get_u64(p), prng_get_u64(p), prng_get_u64(p), prng_get_u64(p));
	__m256i v_r_shifted = _mm256_srli_epi16(v_r, 1);
	__m256i v_z = v_zero;

	//for (size_t k = 0; k < sizeof(Dist_4) / sizeof(Dist_4[0]); k++) {
	//	__m256i v_dist = _mm256_set1_epi16(Dist_4[k]);
	//	__m256i mask = _mm256_cmpgt_epi16(v_r_shifted, v_dist); // 生成遮罩
	//	v_z = _mm256_add_epi16(v_z, _mm256_and_si256(mask, v_one));
	//	if (_mm256_testz_si256(mask, _mm256_set1_epi16(-1))) { break; }
	//}
	__m256i v_dist, mask;
	int k = 0;
	do {
		v_dist = _mm256_set1_epi16(Dist_4[k]);
		mask = _mm256_cmpgt_epi16(v_r_shifted, v_dist); // 生成遮罩
		v_z = _mm256_add_epi16(v_z, _mm256_and_si256(mask, v_one));
		k++;
	} while (!_mm256_testz_si256(mask, _mm256_set1_epi16(-1)));

	*z_out = v_z;
}
int s4_v(void* ctx, int* samples) {
	sampler_context* sc = ctx;
	int z[16] = { 0 };

	int16_t sample[16] = { 0 };
	int n = 0; // zs的索引
	int	m = 0; // sample的索引
	
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

	bs4_v(&sc->p, sample, mark); // 生成16个样本
	double x, p;
	int i, b;
	uint8_t u, v; // 惰性浮点

	while (1) {
		i = 0xff;
		u = prng_get_u8(&sc->p);
		b = (int)u & 1;// 0 or 1
		z[n] = (2 * b - 1) * sample[m] + b;
		x = sample[m] * sample[m] * isigma_maxT[mark];
		x = (z[n] - sc->center) * (z[n] - sc->center) * isigma * isigma * 0.5 - x;
		p = sis * expm_p63(x);
		v = (int)(p * i) & 0xff;

		while (u == v) {
			i = i * 0xff;
			u = prng_get_u8(&sc->p);
			v = (int)(p * i) & 0xff;
		}

		if (u < v) {
			samples[n] = z[n];
			n++;
			m++;
		}
		else { m++; }
		if (m == 16) { break; }
	}

	return n;
}

