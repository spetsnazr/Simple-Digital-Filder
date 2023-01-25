#include "Fun.h"
void FFT(std::vector<std::complex<float> > &y, int i, int N, bool mod) {
    int m(2);
    if (mod) {
        y.resize(N);
        for (int i(N - 1), j(1); i > N / 2; i--, j++) y[i] = std::conj(y[j]);
    };
    int l(0);
    for (int i(0); i <= N - 2; i++) {
        if (i < l) {
            std::complex<float> temp(y[i]);
            y[i] = y[l], y[l] = temp;
        }
        int k(N / 2);
        while (k <= l) l -= k, k /= 2;
        l += k;
    }
    while (m <= y.size()) {
        std::complex<float> mi(1, 0), w(Wn(i, m));
        for (int j(0); j < m / 2; j++) {
            for (int k(j); k <= j + N - m; k += m) {
                std::complex<float>u(y[k]), v((mi*y[k + m / 2]));
                y[k] = u + v, y[k + m / 2] = u - v;
            }
            mi *= w;
        }
        m <<= 1;
    }
    if (i == 1) {
        for (int i(0); i < y.size(); i++) y[i] /= N;
    }
}
std::complex<float> Wn(int x, int m){
   return std::complex<float>(std::cos(2 * std::abs(x)*PI / m), std::sin(2 * x*PI / m));
}
void ImpulsniOdziv(std::vector<std::complex<float> > &y,  std::vector<float> &naz, std::vector<float> &br){
    unsigned short vel;
    if(br.size() > naz.size()) vel = br.size();
    else vel = naz.size(); 
    int k(4);
    bool nijestepen(false);
    if(vel&(vel-1)){
        while(k < vel) k<<=1;    
        nijestepen = true;
    }
    else k = vel;
    br.resize(vel);
    naz.resize(vel);
    y.resize(k);
   
    for(int i(0); i < vel; i++){
        std::complex<float> nazivnik(0), brojnik(0), temp;
        for(int j(0); j < vel; j++) nazivnik += naz[j]*(temp = Wn(i*j,k)),brojnik += br[j]*temp;
        y[i] = brojnik/nazivnik; 
    }
    if(nijestepen) for(int i(k-1), j(1); i >= vel; i--,j++) y[i] = std::conj(y[j]);
    FFT(y,1,k,false);   
}
void NFFilter(std::vector<float> &num, std::vector<float> &den, double &fg, bool mod){
    if(mod){
        for (int i(0); i < sz; i++)
            num[i] = sin(fg*(i-int(sz/2.))) /(PI * (i-int(sz/2.)));
    }else{ 
        double alfa((1-sin(fg))/cos(fg)); 
        num.push_back(1-alfa),num.push_back(1-alfa),den.push_back(2),den.push_back(-2*alfa);
    }
}
void VFFilter(std::vector<float> &num, std::vector<float> &den,double &fg, bool mod){
    if(mod){
           for (int i(0); i < sz; i++)
                num[i] = -sin(fg*i) / (PI * i); 
            num[0] = 1.-fg/PI;
    }else{
        double alfa((1-sin(fg))/cos(fg)); 
        num.push_back(1+alfa),num.push_back(-1-alfa),den.push_back(2),den.push_back(-2*alfa);
    }
}
void BandPassFilter(std::vector<float> &num, std::vector<float> &den,double f1, double f2, bool mod){
    if(mod){
        for (int i(0); i < sz; i++) num[i] = sin(f2*i)/(PI * i)-sin(f1*i) / (PI * i);
    }else{
        double beta(cos((f1+f2)/2)),a1((2+sqrt(4-4*cos(f2-f1)*cos(f2-f1)))/(2*cos(f2-f1))),a2((2-sqrt(4-4*cos(f2-f1)*cos(f2-f1)))/(2*cos(f2-f1))),pravi; 
        std::vector<std::complex<double> > koef;
        koef.push_back(a1*2),koef.push_back(-beta*2*(1+a1)),koef.push_back(2);
        koef = PolyRoots(koef);
        bool test(true);
        for(int i(0); i < koef.size(); i++) if(sqrt(koef[i].real()*koef[i].real() + koef[i].imag()*koef[i].imag()) > 1) pravi = a2,test = false;
        if(test) pravi = a1;
        num.push_back(1-pravi),num.push_back(0),num.push_back(pravi-1),den.push_back(2),den.push_back(-2*beta*(1+pravi)),den.push_back(2*pravi);
    }
}
void BandStopFilter(std::vector<float> &num, std::vector<float> &den,double f1, double f2, bool mod){
    if(mod){
        for (int i(1); i < sz; i++) num[i] = sin(f1*i)/(PI * i)-sin(f2*i) / (PI * i);
        num[0] = 1-(f2-f1)/PI;
    }else{
         double beta(cos((f1+f2)/2)),a1((2+sqrt(4-4*cos(f2-f1)*cos(f2-f1)))/(2*cos(f2-f1))),a2((2-sqrt(4-4*cos(f2-f1)*cos(f2-f1)))/(2*cos(f2-f1))),pravi; 
        std::vector<std::complex<double> > koef;
        koef.push_back(a1*2),koef.push_back(-beta*2*(1+a1)),koef.push_back(2);
        koef = PolyRoots(koef);
        bool test(true);
        for(int i(0); i < koef.size(); i++) if(sqrt(koef[i].real()*koef[i].real() + koef[i].imag()*koef[i].imag()) > 1) pravi = a2,test = false;
        if(test) pravi = a1;
        num.push_back(1-pravi),num.push_back(0),num.push_back(pravi-1),den.push_back(2),den.push_back(-2*beta*(1+pravi)),den.push_back(2*pravi);
    }
}
std::pair<bool,std::complex<double> > Laguerre(std::vector<std::complex<double> > &koef, std::complex<double> x, int max,double eps,int red){
    std::complex<double> dx(std::numeric_limits<double>::infinity());
    int k(1);
    int n(red);
    while (std::fabs(dx.real()) > eps && k < max) {
        std::complex<double> f(koef[red]);
        std::complex<double> d(0), s(0);
        for (int i(red-1); i > -1; i--) s = s*x + std::complex<double>(2, 0)*d, d = d*x + f, f = f*x + koef[i];
        if (std::fabs((f - std::complex<double>(0)).real()) < std::numeric_limits<double>::epsilon()) return std::pair<bool, std::complex<double> >(true, x);
        std::complex<double> r(std::sqrt(std::complex<double>(n-1, 0)*(std::complex<double>(n-1, 0)*d*d - std::complex<double>(n, 0)*f*s)));
        if (std::fabs((d + r).real()) > std::fabs(((d - r)).real())) dx = std::complex<double>(n, 0)*f / (d + r);
        else dx = (std::complex<double>(n, 0)*f / (d - r));
        x -= dx, k++;
    }
    if (dx.real() < eps || std::fabs(dx.real() - eps) < std::numeric_limits<double>::epsilon()) return std::pair<bool, std::complex<double> >(true, x);
    return std::pair<bool, std::complex<double> >(false, x);
}
std::vector<std::complex<double> >PolyRoots(std::vector<std::complex<double> > &coefficients, double eps,int maxiters, int maxtrials){
    std::pair<bool, std::complex<double> > x;
    std::vector<std::complex<double> > x2(coefficients.size()-1);
    std::complex<double> v, w;
    for (int i(coefficients.size() - 1); i > 0; i--) {
        double t(1);
        x.first = false;
        while (!x.first && t < maxtrials) {
            x.second = std::complex<double>(rand() % 20 - 10, rand() % 20 - 10);
            x = Laguerre(coefficients, x.second, maxiters, eps,i);
            t++;
        }
        if (std::fabs(x.second.imag() - eps) < std::numeric_limits<double>::epsilon() || std::fabs(x.second.imag()) < eps) x.second = x.second.real();
        x2[i-1] = x.second;
        v = coefficients[i];
        for (int j(i - 1); j > -1; j--) {
            w = coefficients[j];
            coefficients[j] = v;
            v = w + x.second*v;
        }
    }
    return x2;
}