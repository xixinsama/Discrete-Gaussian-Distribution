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
112789371u,
203104039u,
249473136u,
264737517u,
267959402u,
268395437u,
268433274u,
268435379u,
268435454u,
268435455u
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
1073741823u,
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
1073741823u,
};


// 28bits sigma = 0.75
static const uint64_t DistForSampler_1_LUT[] = {
142782702u,
260182149u,
268339468u,
268435264u,
268435454u,
268435455u
};

static const double sigma_minT[] = {
0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5
};

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

	uint32_t r = prng_get_u32(&spc->p) >> 4;

	int length = sizeof(DistForSampler_1_LUT) / sizeof(DistForSampler_1_LUT[0]);

	for (int i = 0; i < length; i++)
	{
		if ((r - DistForSampler_1_LUT[i]) >> 31)
		{
			z = i;
			break;
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

static int BaseSampler2_CDT(prng* p)
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
static int AcceptSample(prng* pp, double sis, double x)
{
	double p = sis * expm_p63(-x);;

	int i = 1;
	uint16_t u, v;

	//惰性浮点伯努利采样
	do {
		i = i * 0xff;
		u = prng_get_u8(pp);
		v = (int)(p * i) & 0xff; //强制类型转换，用于向下取整
	} while (u == v);

	return (u - v) >> 31;
}

static int BaseSampler3(prng* p)
{
	int z = 0;
	uint32_t r = prng_get_u32(p) >> 4; // 32bit
	int temp = 0;

	while ((DistForBSampler_3_CDT[z] - r) >> 31)
	{
		z = z + 1;
	}

	return z;
}

/*
 * Sample a bit with probability exp(-x) for some x >= 0.
 */
static int
BerExp(prng* p, double x, double ccs)
{
	int s, i;
	double r;
	uint32_t sw, w;
	uint64_t z;

	/*
	 * Reduce x modulo log(2): x = s*log(2) + r, with s an integer,
	 * and 0 <= r < log(2). Since x >= 0, we can use fpr_trunc().
	 */
	s = (int)(uint64_t)(x * 1.4426950408889634073599246810);
	r = x - s * 0.69314718055994530941723212146;

	/*
	 * It may happen (quite rarely) that s >= 64; if sigma = 1.2
	 * (the minimum value for sigma), r = 0 and b = 1, then we get
	 * s >= 64 if the half-Gaussian produced a z >= 13, which happens
	 * with probability about 0.000000000230383991, which is
	 * approximatively equal to 2^(-32). In any case, if s >= 64,
	 * then BerExp will be non-zero with probability less than
	 * 2^(-64), so we can simply saturate s at 63.
	 */
	sw = (uint32_t)s;
	sw ^= (sw ^ 63) & -(int)((63 - sw) >> 31);  //这里跟原来的不同
	s = (int)sw;

	/*
	 * Compute exp(-r); we know that 0 <= r < log(2) at this point, so
	 * we can use fpr_expm_p63(), which yields a result scaled to 2^63.
	 * We scale it up to 2^64, then right-shift it by s bits because
	 * we really want exp(-x) = 2^(-s)*exp(-r).
	 *
	 * The "-1" operation makes sure that the value fits on 64 bits
	 * (i.e. if r = 0, we may get 2^64, and we prefer 2^64-1 in that
	 * case). The bias is negligible since fpr_expm_p63() only computes
	 * with 51 bits of precision or so.
	 */
	z = ((fpr_expm_p63(r, ccs) << 1) - 1) >> s;

	/*
	 * Sample a bit with probability exp(-x). Since x = s*log(2) + r,
	 * exp(-x) = 2^-s * exp(-r), we compare lazily exp(-x) with the
	 * PRNG output to limit its consumption, the sign of the difference
	 * yields the expected result.
	 */
	i = 64;
	do {
		i -= 8;
		w = prng_get_u8(p) - ((uint32_t)(z >> i) & 0xFF);
	} while (!w && i > 0);
	return (int)(w >> 31);
}

// Fixed sigma = 1.5 and center c is uniformly distributed over [0,1)
int sampler_3(void* ctx) {
	// double sigma = 1.5;

	sampler_context* spc;
	spc = ctx;
	int s;
	double r, dss, ccs, dss0;

	s = (int)(spc->center); //=0
	r = spc->center - s;

	double isigma = 1 / spc->sigma;
	dss = 0.5 * isigma * isigma; //=1 / 2*sigma^2

	ccs = isigma * spc->sigma_min;

	dss0 = 0.125 / (spc->sigma_min * spc->sigma_min); //=1 / 8*sigma_min^2

	for (;;) {
		int z0, z, b;
		double x;

		z0 = BaseSampler3(&spc->p);
		b = (int)prng_get_u8(&spc->p) & 1;
		z = b + ((b << 1) - 1) * z0;

		/*
		 * Rejection sampling. We want a Gaussian centered on r;
		 * but we sampled against a Gaussian centered on b (0 or
		 * 1). But we know that z is always in the range where
		 * our sampling distribution is greater than the Gaussian
		 * distribution, so rejection works.
		 *
		 * We got z with distribution:
		 *    G(z) = exp(-((z-b)^2)/(2*sigma0^2))
		 * We target distribution:
		 *    S(z) = exp(-((z-r)^2)/(2*sigma^2))
		 * Rejection sampling works by keeping the value z with
		 * probability S(z)/G(z), and starting again otherwise.
		 * This requires S(z) <= G(z), which is the case here.
		 * Thus, we simply need to keep our z with probability:
		 *    P = exp(-x)
		 * where:
		 *    x = ((z-r)^2)/(2*sigma^2) - ((z-b)^2)/(2*sigma0^2)
		 *
		 * Here, we scale up the Bernouilli distribution, which
		 * makes rejection more probable, but makes rejection
		 * rate sufficiently decorrelated from the Gaussian
		 * center and standard deviation that the whole sampler
		 * can be said to be constant-time.
		 */

		x = (z - r) * (z - r) * dss;
		x = x - z0 * z0 * dss0;
		if (BerExp(&spc->p, x, ccs)) {
			/*
			 * Rejection sampling was centered on r, but the
			 * actual center is mu = s + r.
			 */
			return s + z;
		}
	}

}

static int BaseSampler4(prng* p, int mark)
{
	int z = 0;
	uint32_t r = prng_get_u32(p) >> 2; // 30bit
	int temp = 0;
	uint32_t* DistForBSampler4; // 指向DistForBSampler的指针

	// 根据mark的值选择不同的DistForBSampler
	switch (mark) {
	case 0:
		DistForBSampler4 = CDT4_09;
		break;
	case 1:
		DistForBSampler4 = CDT4_10;
		break;
	case 2:
		DistForBSampler4 = CDT4_11;
		break;
	case 3:
		DistForBSampler4 = CDT4_12;
		break;
	case 4:
		DistForBSampler4 = CDT4_13;
		break;
	case 5:
		DistForBSampler4 = CDT4_14;
		break;
	case 6:
		DistForBSampler4 = CDT4_15;
		break;
	case 7:
		DistForBSampler4 = CDT4_16;
		break;
	default:
		DistForBSampler4 = CDT4_16;
	}

	while ((DistForBSampler4[z] - r) >> 31) {
		z = z + 1;
	}
	return z;
}

// sigma is uniformly distributed over (0.8,1.6) and center is uniformly distributed over [0,1)
int sampler_4(void* ctx) {
	sampler_context* spc;
	spc = ctx;
	int z = 0;

	int mark = 7;
	if (spc->sigma > 0.8 && spc->sigma <= 0.9) {
		mark = 0;
	}
	else if (spc->sigma > 0.9 && spc->sigma <= 1.0) {
		mark = 1;
	}
	else if (spc->sigma > 1.0 && spc->sigma <= 1.1) {
		mark = 2;
	}
	else if (spc->sigma > 1.1 && spc->sigma <= 1.2) {
		mark = 3;
	}
	else if (spc->sigma > 1.2 && spc->sigma <= 1.3) {
		mark = 4;
	}
	else if (spc->sigma > 1.3 && spc->sigma <= 1.4) {
		mark = 5;
	}
	else if (spc->sigma > 1.4 && spc->sigma <= 1.5) {
		mark = 6;

	}
	else if (spc->sigma > 1.5 && spc->sigma <= 1.6) {
		mark = 7;
	}

	double isigma = 1 / spc->sigma; // 1 / sigma
	double sis = sigma_minT[mark] * isigma; // sigma_min / sigma

	while (1)
	{
		int z0 = BaseSampler4(&spc->p, mark);
		int b = prng_get_u8(&spc->p) >> 7;
		z = (2 * b - 1) * z0 + b;
		double x = z0 * z0 * isigma_maxT[mark];
		x = x - (z - spc->center) * (z - spc->center) * isigma * isigma * 0.5;
		if (AcceptSample(&spc->p, sis, x))
		{
			return z;
			break;
		}
	}
}

int sampler_karney(void* ctx)
{
	sampler_context* spc;
	spc = ctx;
	int z = 0;

	return z;
}

