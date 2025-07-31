#include "fft_analyse.h"

using namespace std;
typedef std::complex<double> Complex;

fft_analyse::fft_analyse()
{
	//构造函数
}  

fft_analyse::~fft_analyse()
{
	//析构函数
} 

void fft_analyse::print_on()
{   
	//成员函数
    is_print = true;
}

void fft_analyse::print_off()
{
    is_print = false;
}

void fft_analyse::fft_setting(int sampling_frenquency_setting)
{
    sampling_frenquency = sampling_frenquency_setting;
}

void fft_analyse::fft(double *data, int data_len)
{
    int N, p=0, t;

    while(data_len>pow(2,p))p++;
    N = pow(2,p);

    for(int i=data_len;i<N;i++)data[i]=0;

    for(int i=0;i<N;i++)
    {
        t=reverse(i,p);
        fft_res[i].real(data[t]);
        fft_res[i].imag(0);
    }

    for(int k=1;k<=p;k++)
        for(int l=0;l<N/pow(2,k);l++)
            for(int i=0;i<pow(2,k-1);i++)
                f(l*pow(2,k)+i,l*pow(2,k)+i+pow(2,k-1),w(i,pow(2,k)));

    if (is_print)
    {
        cout << "FFT处理后结果：";
        for(int i=0; i<N; i++) cout << fft_res[i] << ' ';
        cout << endl;
    }

    double * amplitude_spectrum = new double(data_len/2+1);
    for (int i=0; i<=data_len/2; i++)
    {
        amplitude_spectrum[i] = 2* sqrt(fft_res[i].real()*fft_res[i].real() + fft_res[i].imag()*fft_res[i].imag())/data_len;
    }

    if (is_print)
    {
        cout << "amplitude_spectrum: ";
        for (int i=0; i<=data_len/2; i++)
        {
            cout << amplitude_spectrum[i] << " ";
        }
        cout << endl;
    }

    pair<double, size_t> max_pair = findMax(amplitude_spectrum, data_len/2+1);

    estimated_frequency = max_pair.second*((double)sampling_frenquency/data_len);
    estimated_amplitude = max_pair.first;

    estimated_phase_rad = arg(fft_res[max_pair.second]);
    estimated_phase_deg = estimated_phase_rad*180.0/M_PI;

    if (is_print)
    {
        cout << "estimated frequency: " << estimated_frequency << endl;
        cout << "estimated amplitude: " << estimated_amplitude << endl;
        cout << "estimated phase: " << estimated_phase_rad << "(" << estimated_phase_deg << ")" << endl;
    }
}

std::pair<double, size_t> fft_analyse::findMax(const double arr[], size_t size) {
    if (size == 0)
    {
        throw std::invalid_argument("数组大小不能为零");
    }

    double max_val = arr[0];
    size_t max_idx = 0;

    for (size_t i = 1; i < size; ++i)
    {
        if (arr[i] > max_val)
        {
            max_val = arr[i];
            max_idx = i;
        }
    }

    return {max_val, max_idx};
}

void fft_analyse::f(int i,int j,Complex w)
{
    Complex t  = fft_res[i];
    fft_res[i] = fft_res[i]+w*fft_res[j];
    fft_res[j] = t-w*fft_res[j];
}

Complex fft_analyse::w(int m,int N)
{
    Complex c;
    c.real(cos(2*acos(-1)*m/N));
    c.imag(sin(-2*acos(-1)*m/N));
    return c;
}

int fft_analyse::reverse(int n,int l)
{
    int r=0, t;
    for (int i=0; i<l; i++)
    {
        t = n&1;
        n = n>>1;
        r+= t*pow(2,l-1-i);
    }
    return r;
}
