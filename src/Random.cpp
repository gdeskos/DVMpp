#include "Random.hpp"

Random::Random()
{
	m_re.seed(m_rd());
}

double Random::rand()
{
	return m_dist(m_re);
}

void Random::seed(int seed)
{
	if (seed >= 0) {
		m_re.seed(seed);
	}
}
