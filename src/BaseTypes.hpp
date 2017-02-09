#ifndef BASETYPES
#define BASETYPES

#include <armadillo>

using Matrix = arma::mat;
using Vector = arma::vec;
using Vector_un = arma::Col<unsigned>;

namespace math
{
const double pi = arma::datum::pi;
const double rpi2 = 1.0 / (2.0 * pi);
};

#endif
