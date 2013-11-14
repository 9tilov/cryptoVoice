#ifndef FFT_H_
#define FFT_H_

#include <vector>

#define _USE_MATH_DEFINES
#include <math.h>

const double TwoPi = 2*M_PI;

void FFTAnalysis(const std::vector<double>& input, std::vector<double>& output, int Nvl, int Nft);

#endif /* FFT_H_ */