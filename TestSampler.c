#include "sampler.h"
#include <immintrin.h>
#include <math.h>

// 测试用函数
// Fixed sigma = 0.75 and center = 0, 22bits
static const uint32_t Dist1[] = {
2230979u,
4065346u,
4192804u,
4194301u,
4194303u
};
// 14bits
static const uint16_t Dist2[] = {
8714u,
15880u,
16378u,
16383u
};
// 7bits
static const uint8_t Dist3[] = {
68u,
124u,
127u
};
int s1(void* ctx) {
	sampler_context* sc = ctx;
	int z = 0;
	uint8_t r = prng_get_u8(&sc->p);
	int s = (int)r & 1; //sign
	r = r >> 1; //24-22
	while ((Dist3[z] - r) >> 7) //16-1
	{
		z++;
	}
	return z = s == 1 ? -z : z;
	return;
}

// 向量化版本，一次输出16个样本
int s1_v(void* ctx, int* samples) {
	sampler_context* sc = ctx;
	const __m256i v_zero = _mm256_setzero_si256();
	const __m256i v_one = _mm256_set1_epi16(1);

	__m256i v_r = _mm256_set_epi64x(prng_get_u64(&sc->p), prng_get_u64(&sc->p), prng_get_u64(&sc->p), prng_get_u64(&sc->p));
	__m256i v_r_shifted = _mm256_srli_epi16(v_r, 2);
	__m256i v_z = v_zero;

	for (size_t k = 0; k < sizeof(Dist2) / sizeof(Dist2[0]); k++) {
		__m256i v_dist = _mm256_set1_epi16(Dist2[k]);
		__m256i mask = _mm256_cmpgt_epi16(v_r_shifted, v_dist);
		v_z = _mm256_add_epi16(v_z, _mm256_and_si256(mask, v_one));
		if (_mm256_testz_si256(mask, _mm256_set1_epi16(-1))) { break; }
	}

	__m256i isHighestBitZero = _mm256_cmpeq_epi16(_mm256_and_si256(v_r, v_one), v_zero);
	v_z = _mm256_blendv_epi8(v_z, _mm256_sub_epi16(v_zero, v_z), isHighestBitZero);

	_mm256_storeu_si256((__m256i*)samples, v_z);
	return 0;
}