/**
 * @file fft_analysis
 * @brief fft数据处理库，和matlab中fft的功能相似
 *
 * Copyright (c) 2025 [JUYI LI].
 * All rights reserved. Licensed under the MIT License.
 *
 * author: JUYI LI
 * created: 2025-07-29
 */

#include "fft_analysis.h"

using namespace std;
typedef std::complex<double> Complex;

fft_analysis::fft_analysis()
{
	//构造函数
}  

fft_analysis::~fft_analysis()
{
	//析构函数
} 

void fft_analysis::print_on()
{   
	//成员函数
    is_print = true;
}

void fft_analysis::print_off()
{
    is_print = false;
}

void fft_analysis::fft_setting(int sampling_frenquency_setting)
{
    sampling_frenquency = sampling_frenquency_setting;
}

void fft_analysis::fft(std::vector<double> data)
{
    int N, p=0, t;

    while(data.size()>pow(2,p))
        p++;
    N = pow(2,p);
//    cout << "N:" << N << endl;

    for(int i=data.size();i<N;i++)
        data.push_back(0);
//    cout << "N:" << N << endl;

    for(int i=0;i<N;i++)
    {
        t=reverse(i,p);
        fft_res[i].real(data[t]);
        fft_res[i].imag(0);
    }
//    cout << "fft_res:" << fft_res[0] << endl;

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

    vector<double> amplitude_spectrum;
    for (int i=0; i<=(int)data.size()/2; i++)
    {
        amplitude_spectrum.push_back(2* sqrt(fft_res[i].real()*fft_res[i].real() + fft_res[i].imag()*fft_res[i].imag())/data.size());
    }

    if (is_print)
    {
        cout << "amplitude_spectrum: ";
        for (int i=0; i<=(int)data.size()/2; i++)
        {
            cout << amplitude_spectrum[i] << " ";
        }
        cout << endl;
    }

    pair<double, size_t> max_pair = findMax(amplitude_spectrum, data.size()/2+1);

    estimated_frequency = max_pair.second*((double)sampling_frenquency/data.size());
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

double fft_analysis::phase_difference_deg(vector<double> data1, vector<double> data2)
{
    fft(data1);
    double pd_1 = estimated_phase_deg;
    fft(data2);
    double pd_2 = estimated_phase_deg;

    double pd = pd_1-pd_2;
    if (pd>180)
        pd = pd - 360;

    if (is_print)
        cout << "Phase difference in deg: " << pd << endl;

    return pd;
}

double fft_analysis::amplitude_dB(vector<double> data1, vector<double> data2)
{
    fft(data1);
    double ea_1 = estimated_amplitude;
    fft(data2);
    double ea_2 = estimated_amplitude;

    double a_dB = 20*log10(ea_1/ea_2);
    if (is_print)
        cout << "Amplitude in dB: " << a_dB << endl;
    return a_dB;
}

std::pair<double, size_t> fft_analysis::findMax(const vector<double> arr, size_t size) {
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

void fft_analysis::f(int i,int j,Complex w)
{
    Complex t  = fft_res[i];
    fft_res[i] = fft_res[i]+w*fft_res[j];
    fft_res[j] = t-w*fft_res[j];
}

Complex fft_analysis::w(int m,int N)
{
    Complex c;
    c.real(cos(2*acos(-1)*m/N));
    c.imag(sin(-2*acos(-1)*m/N));
    return c;
}

int fft_analysis::reverse(int n,int l)
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
