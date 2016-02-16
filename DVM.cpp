#include "DVMBase.hpp"
#include <time.h>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <string>
#include <iomanip>
#include <iostream>
#include <vector>
#include <sstream>
#include <armadillo>

DVMBase::DVMBase()
{
    m_pi=4.0*atan(1.0);
}

