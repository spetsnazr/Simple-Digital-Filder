#ifndef _FUN_H_
#define _FUN_H_
#include <vector>
#include <complex>
#include <string>
#include "mbed.h"
#include "stdio.h"
#include "SPI_TFT_ILI9341.h"
#include "font_big.h"
#include <istream>
#include <stdlib.h>
#include <utility>
#include <ctime>
#include <deque>
#include <algorithm>
#define _CSTD   ::

void ImpulsniOdziv(std::vector<std::complex<float> > &y,  std::vector<float> &naz, std::vector<float> &br);
const float PI(4*atan(1.));
const unsigned short sz(31);
void Poruka(string *por);
void FFT(std::vector<std::complex<float> > &y, int i, int N, bool mod);
void NFFilter(std::vector<float> &num, std::vector<float> &den, double &fg, bool mod);
void VFFilter(std::vector<float> &num, std::vector<float> &den,double &fg, bool mod);
std::complex<float> Wn(int x, int m);
void BandPassFilter(std::vector<float> &num, std::vector<float> &den,double f1, double f2, bool mod);
void BandStopFilter(std::vector<float> &num, std::vector<float> &den,double f1, double f2, bool mod);
//pretvara string u double (uzeto iz biblioteke string sa C++11 standardom)
inline double stod(const string& _Str, size_t *_Idx = 0){   
    const char *_Ptr = _Str.c_str();
    char *_Eptr;
    double _Ans = _CSTD strtod(_Ptr, &_Eptr);
    return (_Ans);
}
std::pair<bool,std::complex<double> > Laguerre(std::vector<std::complex<double> > &koef, std::complex<double> x, int max,double eps,int red);
std::vector<std::complex<double> >PolyRoots(std::vector<std::complex<double> > &coefficients, double eps = 1e-10,int maxiters = 100, int maxtrials = 10);
#endif

