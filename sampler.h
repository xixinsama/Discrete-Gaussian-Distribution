#ifndef INTEGER_GAUSSIAN_SAMPLER__
#define INTEGER_GAUSSIAN_SAMPLER__

#include <stdint.h>

/**************************************************************************
 ************************     API     *************************************
 **************************************************************************
 */


 /* ==================================================================== */
 /*
  * SHAKE256 implementation (shake.c).
  *
  * API is defined to be easily replaced with the fips202.h API defined
  * as part of PQClean.
  */
  // 用于存储 SHAKE256 算法的上下文，包括状态数组 st 和数据指针 dptr。
typedef struct
{
	union
	{
		uint64_t A[25];
		uint8_t dbuf[200];
	} st;
	uint64_t dptr;
} sampler_shake256_context;

void i_shake256_init(sampler_shake256_context* sc);
void i_shake256_inject(sampler_shake256_context* sc, const uint8_t* in, size_t len);
void i_shake256_flip(sampler_shake256_context* sc);
void i_shake256_extract(sampler_shake256_context* sc, uint8_t* out, size_t len);

// 从 SHAKE256 上下文中获取一个 64 位的随机数。
static inline uint64_t
shake256_get_u64(sampler_shake256_context* p) {
	uint8_t random_bytes[8];

	i_shake256_extract(p, random_bytes, 8);
	uint64_t result;
	for (int i = 0; i < 8; i++) {
		result |= ((uint64_t)random_bytes[i] << (8 * i));
	}
	return result;
}

// 从 SHAKE256 上下文中获取一个 32 位的随机数。
static inline uint32_t
shake256_get_u32(sampler_shake256_context* p) {
	uint8_t random_bytes[4];

	i_shake256_extract(p, random_bytes, 4);
	uint32_t result;
	for (int i = 0; i < 4; i++) {
		result |= ((uint32_t)random_bytes[i] << (8 * i));
	}
	return result;
}

// 从 SHAKE256 上下文中获取一个 8 位的随机数。
static inline uint8_t
shake256_get_u8(sampler_shake256_context* p) {
	uint8_t random_bytes[1];
	i_shake256_extract(p, random_bytes, 1);
	return random_bytes[0];
}


/* ==================================================================== */
/*
 * RNG (rng.c).
 *
 * A PRNG based on ChaCha20 is implemented; it is seeded from a SHAKE256
 * context (flipped) and is used for bulk pseudorandom generation.
 * A system-dependent seed generator is also provided.
 */

 /*
  * Structure for a PRNG. This includes a large buffer so that values
  * get generated in advance. The 'state' is used to keep the current
  * PRNG algorithm state (contents depend on the selected algorithm).
  * The unions with 'dummy_u64' are there to ensure proper alignment for
  * 64-bit direct access.
  */
  // 用于存储伪随机数生成器（PRNG）的状态和缓冲区。
typedef struct
{
	union
	{
		uint8_t d[512]; /* MUST be 512, exactly */
		uint64_t dummy_u64;
	} buf;
	size_t ptr;
	union
	{
		uint8_t d[256];
		uint64_t dummy_u64;
	} state;
	int type;
} prng;

/*
 * Instantiate a PRNG. That PRNG will feed over the provided SHAKE256
 * context (in "flipped" state) to obtain its initial state.
 */
 // 初始化 PRNG p，使用 SHAKE256 上下文 src 作为种子。
void prng_init(prng* p, sampler_shake256_context* src);

/*
 * Refill the PRNG buffer. This is normally invoked automatically, and
 * is declared here only so that prng_get_u64() may be inlined.
 */
 // 重新填充 PRNG p 的缓冲区，确保有足够的随机数据可供使用。
void prng_refill(prng* p);

/*
 * Get some bytes from a PRNG.
 */
 // 从 PRNG p 中获取 len 字节的随机数据，存储到 dst 中。
void prng_get_bytes(prng* p, void* dst, size_t len);

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
	 */
	u = p->ptr;
	if (u >= (sizeof p->buf.d) - 9)
	{
		prng_refill(p);
		u = 0; //u = p->ptr;
	}
	p->ptr = u + 8;

	/*
	 * On systems that use little-endian encoding and allow
	 * unaligned accesses, we can simply read the data where it is.
	 */
	return (uint64_t)p->buf.d[u + 0] |
		((uint64_t)p->buf.d[u + 1] << 8) |
		((uint64_t)p->buf.d[u + 2] << 16) |
		((uint64_t)p->buf.d[u + 3] << 24) |
		((uint64_t)p->buf.d[u + 4] << 32) |
		((uint64_t)p->buf.d[u + 5] << 40) |
		((uint64_t)p->buf.d[u + 6] << 48) |
		((uint64_t)p->buf.d[u + 7] << 56);
}

/*
 * Get a 32-bit random value from a PRNG.
 */
static inline uint32_t
prng_get_u32(prng* p)
{
	size_t u;

	u = p->ptr;
	if (u >= (sizeof p->buf.d) - 5)
	{
		prng_refill(p);
		u = 0;
	}
	p->ptr = u + 4;

	return (uint32_t)p->buf.d[u + 0] |
		((uint32_t)p->buf.d[u + 1] << 8) |
		((uint32_t)p->buf.d[u + 2] << 16) |
		((uint32_t)p->buf.d[u + 3] << 24);
}

/*
 * Get an 8-bit random value from a PRNG.
 */
//static inline size_t
//prng_get_u8(prng* p)
//{
//	size_t v;
//
//	v = p->buf.d[p->ptr++];
//	if (p->ptr == sizeof p->buf.d)
//	{
//		prng_refill(p);
//	}
//	return v;
//}

static inline uint8_t 
prng_get_u8(prng* p) {
	size_t u = p->ptr;
	if (u >= (sizeof p->buf.d) - 2) {
		prng_refill(p);
		u = 0;
	}
	p->ptr = u + 1;
	return p->buf.d[u];
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
	if (w >= (sizeof p->buf.d) - 8) {
		prng_refill(p);
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

// compute exp(-x)
static inline double expm_p63(double x)
{
	double y;
	if (x > 1) {
		int n = (int)(x)+1;
		double e_x = expm_p63(x / n);
		y = e_x;
		for (int i = 0; i < n - 1; i++) {
			y = y * e_x;
		}
	}
	else {
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
	}
	return y;
}

/* ==================================================================== */

typedef struct {
	prng p;
	sampler_shake256_context p1;
	double center;
	double sigma;
} sampler_context;

int sampler_1(void* ctx); // 对应 sampler_1_CDT(void* ctx);
int sampler_2(void* ctx); // 对应 sampler_2_Vector(void* ctx);
int sampler_3(void* ctx);
int sampler_4(void* ctx);

int sampler_1_KY(void* ctx);
int sampler_1_LUT(void* ctx);
int sampler_1_Reject(void* ctx);
int sampler_1_vector(void* ctx, int* samplers);

int sampler_2_ori(void* ctx);
int sampler_2_karney(void* ctx);

int sampler_4_karney(void* ctx);
int sampler_5(void* ctx);

// 测试函数
int s1(void* ctx);
int s1_v(void* ctx, int* samples);

int s2_v(void* ctx, int16_t* samples);

int s3_v(void* ctx, int* samples);


/*
 * Get a 24-bit random value from a PRNG.
 */
static inline uint32_t
prng_get_u24(prng* p)
{
	size_t u;

	u = p->ptr;
	if (u >= (sizeof p->buf.d) - 4)
	{
		prng_refill(p);
		u = 0;
	}
	p->ptr = u + 3;

	return (uint32_t)p->buf.d[u + 0] |
		((uint32_t)p->buf.d[u + 1] << 8) |
		((uint32_t)p->buf.d[u + 2] << 16);
}

/*
 * Get a 16-bit random value from a PRNG.
 */
static inline uint16_t
prng_get_u16(prng* p)
{
	size_t u;

	u = p->ptr;
	if (u >= (sizeof p->buf.d) - 3)
	{
		prng_refill(p);
		u = 0;
	}
	p->ptr = u + 2;

	return (uint16_t)p->buf.d[u + 0] |
		((uint16_t)p->buf.d[u + 1] << 8);
}

#endif
