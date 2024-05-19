#ifndef INTEGER_GAUSSIAN_SAMPLER__
#define INTEGER_GAUSSIAN_SAMPLER__

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include<immintrin.h>
#include<omp.h>

/**************************************************************************
 ************************     API     *************************************
 **************************************************************************
 */


int sampler_1_CDT(void* ctx);
int sampler_1_KY(void* ctx);
int sampler_1_LUT(void* ctx);
int sampler_1_Reject(void* ctx);

int sampler_2(void* ctx);

int sampler_3(void* ctx);

int sampler_4(void* ctx);

int sampler_karney(void* ctx);

/*
 * "Naming" macro used to apply a consistent prefix over all global
 * symbols.
 */
#define Zf(name)             Zf_(FALCON_PREFIX, name)
#define Zf_(prefix, name)    Zf__(prefix, name)
#define Zf__(prefix, name)   prefix ## _ ## name


 /* ==================================================================== */
 /*
  * SHAKE256 implementation (shake.c).
  *
  * API is defined to be easily replaced with the fips202.h API defined
  * as part of PQClean.
  */

typedef struct {
	union {
		uint64_t A[25];
		uint8_t dbuf[200];
	} st;
	uint64_t dptr;
} sampler_shake256_context;

#define sampler_shake256_init      Zf(i_shake256_init)
#define sampler_shake256_inject    Zf(i_shake256_inject)
#define sampler_shake256_flip      Zf(i_shake256_flip)
#define sampler_shake256_extract   Zf(i_shake256_extract)

void Zf(i_shake256_init)(
	sampler_shake256_context* sc);
void Zf(i_shake256_inject)(
	sampler_shake256_context* sc, const uint8_t* in, size_t len);
void Zf(i_shake256_flip)(
	sampler_shake256_context* sc);
void Zf(i_shake256_extract)(
	sampler_shake256_context* sc, uint8_t* out, size_t len);

/* ==================================================================== */
/*
 * RNG (rng.c).
 *
 * A PRNG based on ChaCha20 is implemented; it is seeded from a SHAKE256
 * context (flipped) and is used for bulk pseudorandom generation.
 * A system-dependent seed generator is also provided.
 * 实现了基于ChaCha20的PRNG； 它是从 SHAKE256 上下文（翻转）播种的，用于批量伪随机生成。
 * 还提供了依赖于系统的种子生成器。
 */

 /*
  * Structure for a PRNG. This includes a large buffer so that values
  * get generated in advance. The 'state' is used to keep the current
  * PRNG algorithm state (contents depend on the selected algorithm).
  * PRNG 的结构。 这包括一个大缓冲区，以便提前生成值。 ‘state’用于保存当前PRNG算法状态（内容取决于所选算法）。
  * 与“dummy_u64”的联合是为了确保正确对齐
  *
  * The unions with 'dummy_u64' are there to ensure proper alignment for
  * 64-bit direct access.
  */
typedef struct {
	union {
		uint8_t d[512]; /* MUST be 512, exactly */
		uint64_t dummy_u64;
	} buf;
	size_t ptr;
	union {
		uint8_t d[256];
		uint64_t dummy_u64;
	} state;
	int type;
} prng;

/*
 * Instantiate a PRNG. That PRNG will feed over the provided SHAKE256
 * context (in "flipped" state) to obtain its initial state.
 */
void Zf(prng_init)(prng* p, sampler_shake256_context* src);

/*
 * Refill the PRNG buffer. This is normally invoked automatically, and
 * is declared here only so that prng_get_u64() may be inlined.
 */
void Zf(prng_refill)(prng* p);

/*
 * Get some bytes from a PRNG.
 */
void Zf(prng_get_bytes)(prng* p, void* dst, size_t len);

/*
 * Get a 64-bit random value from a PRNG.
 */
static inline uint64_t
prng_get_u64(prng* p)
{
	size_t u;

	/*
	 * If there are less than 9 bytes in the buffer, we refill it.
	 * This means that we may drop the last few bytes, but this allows
	 * for faster extraction code. Also, it means that we never leave
	 * an empty buffer.
	 * 如果缓冲区中的字节数少于 9 个字节，我们将重新填充它。
	 * 这意味着我们可能会删除最后几个字节，但这可以更快地提取代码。
	 * 另外，这意味着我们永远不会留下空缓冲区。
	 */
	u = p->ptr;
	if (u >= (sizeof p->buf.d) - 9) {
		Zf(prng_refill)(p);
		u = 0;
	}
	p->ptr = u + 8;

	/*
	 * On systems that use little-endian encoding and allow
	 * unaligned accesses, we can simply read the data where it is.
	 */
	return (uint64_t)p->buf.d[u + 0]
		| ((uint64_t)p->buf.d[u + 1] << 8)
		| ((uint64_t)p->buf.d[u + 2] << 16)
		| ((uint64_t)p->buf.d[u + 3] << 24)
		| ((uint64_t)p->buf.d[u + 4] << 32)
		| ((uint64_t)p->buf.d[u + 5] << 40)
		| ((uint64_t)p->buf.d[u + 6] << 48)
		| ((uint64_t)p->buf.d[u + 7] << 56);
}

/*
 * Get a 32-bit random value from a PRNG.
 */
static inline uint32_t
prng_get_u32(prng* p)
{
	size_t u;

	u = p->ptr;
	if (u >= (sizeof p->buf.d) - 9) {
		Zf(prng_refill)(p);
		u = 0;
	}
	p->ptr = u + 4;

	return (uint32_t)p->buf.d[u + 0]
		| ((uint32_t)p->buf.d[u + 1] << 8)
		| ((uint32_t)p->buf.d[u + 2] << 16)
		| ((uint32_t)p->buf.d[u + 3] << 24);
}

/*
 * Get an 8-bit random value from a PRNG.
 */
static inline uint8_t
prng_get_u8(prng* p)
{
	size_t v;

	v = p->buf.d[p->ptr++];
	if (p->ptr == sizeof p->buf.d) {
		Zf(prng_refill)(p);
	}
	return v;
}


/*
 * Get a 256-bit random value (8x32-bit) from a PRNG.
 */
static inline __m256i
prng_get_u32_v(prng* p)
{
	size_t u;
	__m256i v = _mm256_setzero_si256();

	for (int i = 0; i < 8; ++i) {
		u = p->ptr;
		if (u >= (sizeof p->buf.d) - 9) {
			Zf(prng_refill)(p);
			u = 0;
		}
		p->ptr = u + 4;

		uint32_t rand_val = (uint32_t)p->buf.d[u + 0]
			| ((uint32_t)p->buf.d[u + 1] << 8)
			| ((uint32_t)p->buf.d[u + 2] << 16)
			| ((uint32_t)p->buf.d[u + 3] << 24);

		switch (i) {
		case 0: v = _mm256_insert_epi32(v, rand_val, 0); break;
		case 1: v = _mm256_insert_epi32(v, rand_val, 1); break;
		case 2: v = _mm256_insert_epi32(v, rand_val, 2); break;
		case 3: v = _mm256_insert_epi32(v, rand_val, 3); break;
		case 4: v = _mm256_insert_epi32(v, rand_val, 4); break;
		case 5: v = _mm256_insert_epi32(v, rand_val, 5); break;
		case 6: v = _mm256_insert_epi32(v, rand_val, 6); break;
		case 7: v = _mm256_insert_epi32(v, rand_val, 7); break;
		}
	}

	return v;
}

/*
 * Get uniformly distributed random numbers in the interval [0,1) from a PRNG.
 */
union uint64_double_union {
	uint64_t u;
	double d;
};

static inline double
prng_get_rand(prng* p)
{
	union uint64_double_union ud;

	size_t w;
	w = p->ptr;
	if (w >= (sizeof p->buf.d) - 9) {
		Zf(prng_refill)(p);
		w = 0;
	}
	p->ptr = w + 8;

	ud.u = (uint64_t)p->buf.d[w + 0]
		| ((uint64_t)p->buf.d[w + 1] << 8)
		| ((uint64_t)p->buf.d[w + 2] << 16)
		| ((uint64_t)p->buf.d[w + 3] << 24)
		| ((uint64_t)p->buf.d[w + 4] << 32)
		| ((uint64_t)p->buf.d[w + 5] << 40)
		| ((uint64_t)p->buf.d[w + 6] >> 4 << 48)
		| 0x3FF0000000000000ULL;

	return ud.d - 1.0; // Subtract 1 to get a value in [0,1)
}

typedef struct {
	prng p;
	double sigma_min;
	double center;
	double sigma;
} sampler_context;

// 计算exp(-x)的近似值
static inline double expm_p63(double x)
{
	double y;
	y = 0.000000002073772366009083061987;
	y = 0.000000025299506379442070029551 - y * x;
	y = 0.000000275607356160477811864927 - y * x;
	y = 0.000002755586350219122514855659 - y * x;
	y = 0.000024801566833585381209939524 - y * x;
	y = 0.000198412739277311890541063977 - y * x;
	y = 0.001388888894063186997887560103 - y * x;
	y = 0.008333333327800835146903501993 - y * x;
	y = 0.041666666666110491190622155955 - y * x;
	y = 0.166666666666984014666397229121 - y * x;
	y = 0.500000000000019206858326015208 - y * x;
	y = 0.999999999999994892974086724280 - y * x;
	y = 1.000000000000000000000000000000 - y * x;
	return y;
}

static inline uint64_t
fpr_expm_p63(double x, double ccs)
{
	/*
	 * Polynomial approximation of exp(-x) is taken from FACCT:
	 *   https://eprint.iacr.org/2018/1234
	 * Specifically, values are extracted from the implementation
	 * referenced from the FACCT article, and available at:
	 *   https://github.com/raykzhao/gaussian
	 * Tests over more than 24 billions of random inputs in the
	 * 0..log(2) range have never shown a deviation larger than
	 * 2^(-50) from the true mathematical value.
	 */


	 /*
	  * Normal implementation uses Horner's method, which minimizes
	  * the number of operations.
	  */

	double d, y;

	d = x;
	y = 0.000000002073772366009083061987;
	y = 0.000000025299506379442070029551 - y * d;
	y = 0.000000275607356160477811864927 - y * d;
	y = 0.000002755586350219122514855659 - y * d;
	y = 0.000024801566833585381209939524 - y * d;
	y = 0.000198412739277311890541063977 - y * d;
	y = 0.001388888894063186997887560103 - y * d;
	y = 0.008333333327800835146903501993 - y * d;
	y = 0.041666666666110491190622155955 - y * d;
	y = 0.166666666666984014666397229121 - y * d;
	y = 0.500000000000019206858326015208 - y * d;
	y = 0.999999999999994892974086724280 - y * d;
	y = 1.000000000000000000000000000000 - y * d;
	y *= ccs;
	return (uint64_t)(y * 9223372036854775808.0);

}

/* ==================================================================== */


#endif
