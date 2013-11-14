#include "fft.h"

void FFTAnalysis(const std::vector<double>& input, std::vector<double>& output, int Nvl, int Nft) {
  int i, j, n, m, Mmax, Istp;
  double Tmpr, Tmpi, Wtmp, Theta;
  double Wpr, Wpi, Wr, Wi;
  std::vector<double> temp(2 * Nvl, 0);
 
  n = Nvl * 2;
 
  for (int k = 0; k < Nvl; ++k) {
	temp.insert(temp.begin() + k * 2 , 0);
	temp.insert(temp.begin() + k * 2 + 1, input[k]);
  }
 
  i = 1; j = 1;
  while (i < n) {
    if (j > i) {
      Tmpr = temp[i];
	   temp[i] = temp[j];
	   temp[j] = Tmpr;
	   Tmpr = temp[i+1];
	   temp[i+1] = temp[j+1];
	   temp[j+1] = Tmpr;
    }
    i = i + 2;
	m = Nvl;
    while ((m >= 2) && (j > m)) {
      j = j - m; 
	  m = m >> 1;
    }
    j = j + m;
  }
 
  Mmax = 2;
  while (n > Mmax) {
    Theta = -TwoPi / Mmax; 
 Wpi = sin(Theta);
    Wtmp = sin(Theta / 2); 
 Wpr = Wtmp * Wtmp * 2;
    Istp = Mmax * 2;
 Wr = 1;
 Wi = 0;
 m = 1;
 
    while (m < Mmax) {
      i = m;
   m = m + 2; 
   Tmpr = Wr; 
   Tmpi = Wi;
      Wr = Wr - Tmpr * Wpr - Tmpi * Wpi;
      Wi = Wi + Tmpr * Wpi - Tmpi * Wpr;
 
      while (i < n) {
        j = i + Mmax;
        Tmpr = Wr * temp[j] - Wi * temp[j-1];
        Tmpi = Wi * temp[j] + Wr * temp[j-1];
  
  
		temp[j] = temp[i] - Tmpr;
		temp[j-1] = temp[i-1] - Tmpi;
		temp[i] = temp[i] + Tmpr;
		temp[i-1] = temp[i-1] + Tmpi;
        i = i + Istp;
      }
    }
 
    Mmax = Istp;
  }
 
  for (i = 0; i < Nft; i++) {
    j = i * 2;
	output[Nft - i - 1] = sqrt(pow(temp[j], 2) + pow(temp[j+1], 2));
  }
}