#include "Fun.h"
Serial pc(USBTX,USBRX);

SPI_TFT_ILI9341 TFT(PTD2,PTD3,PTD1,PTE3,PTA20,PTD5,"TFT");
AnalogOut S(PTE30);
AnalogIn U(PTB0);
const unsigned short brdijelova(256);
double DajKoef();
void PosaljiNaIzlaz(std::vector<float> &den, std::vector<float> &num, int fz);
void PopuniDif(std::vector<float> &den, std::vector<float> &num);
void PopuniUzorke(std::vector<std::complex<float> > &signal, int fz);
void BrzaKonvolucija(std::vector<std::complex<float> > &signal, std::vector<std::complex<float> > *h, int fz, char &izbor);

int main() {
    pc.baud(9600);
    TFT.claim(stdout);
    TFT.set_orientation(2);
    TFT.background(Black);     
    TFT.foreground(White);    
    TFT.cls();  
    char t, izbor = 'k';
    while(!pc.readable());
    string  *por;
    int fz;
    char izb('u');
    while(1) {
        por = new string("Sa slovom E se unos moze potvrditi\r\n");
        Poruka(por);
        por = new string("Odaberite opciju (Unesite broj 0-9):\r\n1.Niskopropusni filter\r\n2.Visokopropusni filter\r\n3.Filter propusnik opsega frekvencija\r\n4.Filter nepropusnik opsega frekvencija\r\n5.Unosenje koeficijenata diferentne jednacine\r\n(Ispod ovog su brzi načini rada (preporučeno za visokofrekventne signale ili ako ima mnogo koeficijenata)\r\n6.Niskopropusni filter\r\n7.Visokopropusni filter\r\n8.Filter propusnik opsega frekvencija\r\n9.Filter nepropusnik opsega frekvencija\r\n0.Unosenje koeficijenata diferentne jednacine\r\n Ako se odabere brzi nacin rada morate uzeti u obzir slj. ograničenje\r\nFrekvencija uzorkovanja moze biti maksimalno 1025-(broj koeficijenata impulsnog odziva) jer je to ogranicenje memorije (algoritmi za ovaj rad zahtjevaju stepen dvojke)\r\n");
        Poruka(por);
        while(!pc.readable());
        while(true){
            if(pc.readable()) t = pc.getc();
            if(t >= '0' && t <= '9') izbor = t;
            if(t == 'E' || t == 'e') break;;
        }
        por = new string("\r\nUnesite frekvenciju uzorkovanja: ");
        Poruka(por);
        fz = DajKoef();
        switch(izbor){
            case '1':{
            por = new string("\r\nUnesite granicnu frekvenciju: ");
            Poruka(por);
            double fg(2*PI*DajKoef());
            std::vector<float> num(0),den(0);
            NFFilter(den,num,fg/=fz,false);
            PosaljiNaIzlaz(num,den,fz);
            break;
            }
            case '2':{
            por = new string("\r\nUnesite granicnu frekvenciju: ");
            Poruka(por); 
            double fg(2*PI*DajKoef());
            std::vector<float> num,den;
            VFFilter(den,num,fg/=fz,false);
            PosaljiNaIzlaz(num,den,fz);
            break;
            }
            case '3':{
            por = new string("\r\nUnesite gornju granicnu frekvenciju");
            Poruka(por);  
            double fg2(2*PI*DajKoef());
            por = new string("\r\nUnesite donju granicnu frekvenciju");
            Poruka(por);  
            double fg1(2*PI*DajKoef()); 
            std::vector<float> num,den;
            BandPassFilter(num,den,fg1/fz,fg2/fz,false);
            PosaljiNaIzlaz(num,den,fz);
            break;
            }
            case '4':{
            
            por = new string("\r\nUnesite gornju granicnu frekvenciju");
            Poruka(por);  
            double fg2(2*PI*DajKoef());
            por = new string("\r\nUnesite donju granicnu frekvenciju");
            Poruka(por);   
            double fg1(2*PI*DajKoef()); 
            std::vector<float> num,den;
            BandStopFilter(num,den,fg1/fz,fg2/fz,false);
            PosaljiNaIzlaz(num,den,fz);
            break;
            }
            case '5':{
            std::vector<float> num,den;
            PopuniDif(num,den);
            PosaljiNaIzlaz(num,den,fz);
            break;
            }
            case '6':{
            por = new string("\r\nUnesite granicnu frekvenciju: ");
            Poruka(por);
            std::vector<std::complex<float> > *h(new std::vector<std::complex<float> > (sz));
            double fg(2*PI*DajKoef());
            {
                std::vector<float> imps(sz);
                {
                    std::vector<float> den;
                    NFFilter(imps,den,fg/=fz,true);
                }
                for(int i(0); i < sz; i++) (*h)[i] = imps[i];
            }
            std::vector<std::complex<float> > signal(brdijelova-h->size()+1);
                BrzaKonvolucija(signal,h,fz,izb);
            break;
            }
            case '7':{
            por = new string("\r\nUnesite granicnu frekvenciju: ");
            Poruka(por);
            double fg(2*PI*DajKoef());
            std::vector<std::complex<float> > *h(new std::vector<std::complex<float> > (sz));
            {
                std::vector<float> imps(sz);
                {
                    std::vector<float> den;
                    VFFilter(imps,den,fg/=fz,true);
                }
                for(int i(0); i < sz; i++) (*h)[i] = imps[i];
            }
            std::vector<std::complex<float> > signal(brdijelova-h->size()+1);
                BrzaKonvolucija(signal,h,fz,izb);
            break;
            }
            case '8':{
                por = new string("\r\nUnesite gornju granicnu frekvenciju");
                Poruka(por);  
                double fg2(2*PI*DajKoef());
                por = new string("\r\nUnesite donju granicnu frekvenciju");
                Poruka(por);  
                double fg1(2*PI*DajKoef()); 
                std::vector<std::complex<float> > *h(new std::vector<std::complex<float> > (sz));
                {
                    std::vector<float> imps(sz);
                    {
                        std::vector<float> den;
                        BandPassFilter(imps,den,fg1/=fz,fg2/=fz,true);
                    }
                    for(int i(0); i < sz; i++) (*h)[i] = imps[i];
                }
                std::vector<std::complex<float> > signal(brdijelova-h->size()+1);
                BrzaKonvolucija(signal,h,fz,izb);
                break;
            }
            case '9':{
                por = new string("\r\nUnesite gornju granicnu frekvenciju");
                Poruka(por);  
                double fg2(2*PI*DajKoef());
                por = new string("\r\nUnesite donju granicnu frekvenciju");
                Poruka(por);  
                double fg1(2*PI*DajKoef()); 
                std::vector<std::complex<float> > *h(new std::vector<std::complex<float> > (sz));
                {
                    std::vector<float> imps(sz);
                    {
                        std::vector<float> den;
                        BandStopFilter(imps,den,fg1/=fz,fg2/=fz,true);
                    }
                    for(int i(0); i < sz; i++) (*h)[i] = imps[i];
                }
                std::vector<std::complex<float> > signal(brdijelova-h->size()+1);
                
                BrzaKonvolucija(signal,h,fz,izb);
                break;
            }
            case '0':{
            std::vector<float> num,den;
            std::vector<std::complex<float> > *h(new std::vector<std::complex<float> > (sz));
            PopuniDif(num,den);
            ImpulsniOdziv(*h,num,den);
            std::vector<std::complex<float> > signal(brdijelova-h->size()+1);
             while(true) {
                BrzaKonvolucija(signal,h,fz,izb);
                if(izb == 'N' || izb == 'n') break;
            }
            break;
            }
        }
    }
}

double DajKoef(){
    char t = '|';
    string broj;
    bool neg(false);
    while(true){
        if(pc.readable()){ 
            t = pc.getc();
            if(t>='0' && t <='9' || t == '.') broj += t;
            if(t == '-') neg = true;
            if(t == 'e' || t == 'E') break;
        }
    } 
    double ret(stod(broj));
    if(neg) return -ret;
    else return ret;
}
void PopuniDif(std::vector<float> &num, std::vector<float> &den){
    string broj;
    char t;
    string *por;
    por = new string("\r\nKoeficijenti se redaju na slj. način (prvo se unosi lijeva strana) a0*y[n] + a1*y[n-1] + ... ad*y[n-d]");
    Poruka(por);
    while(!pc.readable());
    bool uneseno(false);
    while(true){
         if(pc.readable()) t = pc.getc(),uneseno = false;
         if(t == ' ') den.push_back(stod(broj)),broj = "", t = 'k';
         if((t>='0' && t <='9' || t == '.') && !uneseno) broj += t,uneseno = true;
         if(t == 'e' || t == 'E') break;
     }
     broj = "";
     por = new string("\r\nSada koeficijenti sa desne strane a0*x[n] + a1*x[n-1] + ... ad*x[n-d]");
     Poruka(por);
     while(!pc.readable());
     while(true){
        if(pc.readable()) t = pc.getc(),uneseno = false;
        if(t == ' ') num.push_back(stod(broj)),broj = "", t = 'k';
        if((t>='0' && t <='9' || t == '.') && !uneseno) broj += t,uneseno = true;
        if(t == 'e' || t == 'E') break;
     }
}

void Poruka(string *por){
    pc.printf("%s",por->c_str());
    delete por;
}
void PosaljiNaIzlaz(std::vector<float> &den, std::vector<float> &num, int fz){
    std::deque<double> y,x;
    Timer t1;
    t1.start();
    int vel1(0),vel2(0),pix(5);
    float preth(0);
    x.push_front(U*3.3);
    string *por;
    por = new string("\r\nslovo N za prekid");
    Poruka(por);
    char izb;
    int k = num.size();
    while(true){
         
         if(pc.readable()) izb = pc.getc();
         if(izb == 'N' || izb == 'n') break;
         if(t1.read_us() > 1000/fz){
               if(vel1){
                   if(pc.readable()) izb = pc.getc(); //čita izbor
                   x.push_front(U*3.3); //u x unosi ulazni napon kada ga očita
                   y.push_front(preth);// u y unosi prethodnu vrijednost 
                   if(vel2 >= den.size()) y.erase(y.end()-1),vel2--; //ako smo popunili onoliko koliko ima vrijednosti sa lijeve strane diferentne jednačine obriši najstarije
                   if(vel1 >= num.size()) x.erase(x.end()-1),vel1--; //obriši najstariju desnu vrijednost diferentne jednačine (koliko će ih biti ovisi koliko koeficijenata unesemo)
                   double izl(0);
                   for(int i(0); i <= vel1; i++) izl += num[i]*x[i]; //prvo sabeir sve ulazne vrijednost
                   for(int i(1); i <= vel2; i++) izl -= den[i]*y[i]; //sada se od izlaza oduzimaju vrijednosit prethodnih signala
                   S = preth = (izl/den[0]+1)/2; //nađi prethodnu vrijednost signala
                   vel1++,vel2++;
                   //if(pix <= 320) TFT.line(y[1],pix-5,y[0],pix,Red),pix+=5;          
                }
                else preth = S = (x[0]/den[0]+1)/2,vel2++,vel1++;
                
            t1.reset();
        }
    }
}
void PopuniUzorke(std::vector<std::complex<float> > &signal, int fz){
    Timer t2;
    t2.start();
    //mala optimizacija brzine
    const unsigned short vel = signal.size();
    for(int i(0); i < vel;) if(t2.read_ms() > 1000./fz) signal[i] = U*3.3,i++,t2.reset();
}
void BrzaKonvolucija(std::vector<std::complex<float> > &signal, std::vector<std::complex<float> > *h, int fz, char &izbor){
    PopuniUzorke(signal,fz);
    signal.resize(brdijelova);
    FFT(signal,-1,brdijelova,false);
    signal.resize((brdijelova/2)+1);
    h->resize(brdijelova);
    FFT((*h),-1,brdijelova,false);
    h->resize((brdijelova/2)+1);
    for(int i(0); i < (brdijelova/2)+1; i++) signal[i] *= (*h)[i];
    delete h;
    char izb;
        FFT(signal,1,brdijelova,true);
        unsigned short i(0),j(signal.size());
        string *por;
        por = new string("\r\nslovo N za prekid");
        Poruka(por);
       /*for(int i(0); i < j; i+= ceil(j/320.)){
            if(i != 0) TFT.line(signal[i-ceil(j/320.)].real(),i-ceil(j/320.),signal[i].real(),i,Red);
            wait(1./fz);
        }*/
        while(true){
             if(signal[i].real() > 1e-2 && signal[i].real() < 3.3-1e-2) S = (signal[i].real()+1)/2;
             wait(1./fz);
             i++;
             if(i >= brdijelova) i = 0; 
             if(pc.readable()) izb = izbor = pc.getc();
             if(izb == 'N' || izb == 'n') break;
        }    
}