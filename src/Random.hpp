#ifndef RANDOM_H
#define RANDOM_H

#include <random>

/// Random takes care of random number generation

class Random
{
	private:
	std::uniform_real_distribution<> m_dist; ///< Uniform distribution in [0, 1)
	std::mt19937 m_re;                       ///< the random number generator
	std::random_device m_rd;                 ///< device for random seeding

	public:
	/// Always seed with the random device and seed if wanted later
	Random();

	/// Return a random number in the interval [0, 1);
	double rand();

	/// Seed the random number generator.
	/** \param seed <0 random; >=0 seed with seed */
	void seed(int seed);
};

#endif
