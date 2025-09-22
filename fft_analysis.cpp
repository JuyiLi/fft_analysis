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
    int N, p=0, t, L=data.size();

    if (is_print)
    {
        cout << "size of data: " << data.size() << endl;
    }

    while(data.size()>pow(2,p))
        p++;
    N = pow(2,p);

    for(int i=data.size();i<N;i++)
        data.push_back(0);

    for(int i=0;i<data.size();i++)
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
//        cout << "FFT处理后结果：";
//        for(int i=0; i<N; i++) cout << fft_res[i] << ' ';
//        cout << endl;
    }

    vector<double> amplitude_spectrum;
    for (int i=0; i<=(int)data.size()/2; i++)
    {
        amplitude_spectrum.push_back(2* sqrt(fft_res[i].real()*fft_res[i].real() + fft_res[i].imag()*fft_res[i].imag())/L);
    }

    if (is_print)
    {
//        cout << "amplitude_spectrum: ";
//        for (int i=0; i<=(int)data.size()/2; i++)
//        {
//            cout << amplitude_spectrum[i] << " ";
//        }
//        cout << endl;
    }

    pair<double, size_t> max_pair = findMax(amplitude_spectrum, L/2+1, 1);

    estimated_frequency = (max_pair.second)*((double)sampling_frenquency/L);
    estimated_amplitude = max_pair.first;
    if (is_print)
    {
        cout << sampling_frenquency << ", " << L << "," << (double)sampling_frenquency/L << endl;
        cout << "the chosen pair: " << max_pair.first << ", " << max_pair.second << endl;
    }

    estimated_phase_rad = arg(fft_res[max_pair.second]);
    estimated_phase_deg = estimated_phase_rad*180.0/3.14159265358979323846;

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
    while(pd>0)
        pd = pd - 360;
//    while(pd<-180)
//        pd = pd + 90;

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

double fft_analysis::amplitude_dB_with_peak(std::vector<double> data1, std::vector<double> data2)
{
    double ea_1 = find_max_in_vector(data1);
    double ea_2 = find_max_in_vector(data2);

    double a_dB = 20*log10(ea_1/ea_2);
    if (is_print)
        cout << "Amplitude in dB: " << a_dB << endl;
    return a_dB;
}

double fft_analysis::amplitude_dB_with_trough(std::vector<double> data1, std::vector<double> data2)
{
    double ea_1 = find_min_in_vector(data1);
    double ea_2 = find_min_in_vector(data2);

    double a_dB = 20*log10(ea_1/ea_2);
    if (is_print)
        cout << "Amplitude in dB: " << a_dB << endl;
    return a_dB;
}

double fft_analysis::amplitude_dB_with_trough(std::vector<double> data1_1, std::vector<double> data1_2, std::vector<double> data2)
{
    vector<double> data1_1_aftermax = getElementsAfterMax(data1_1);
    vector<double> data1_2_beforemax = getElementsBeforeMax(data1_2);

    data1_1_aftermax.insert(data1_1_aftermax.end(), data1_2_beforemax.begin(), data1_2_beforemax.end());

    return amplitude_dB_with_trough(data1_1_aftermax, data2);
}

double fft_analysis::amplitude_dB_with_peak_and_trough(std::vector<double> data1, std::vector<double> data2)
{
    double ea_1_max = find_max_in_vector(data1);
    double ea_2_max = find_max_in_vector(data2);

    double ea_1_min = find_min_in_vector(data1);
    double ea_2_min = find_min_in_vector(data2);

    double a_dB = 20*log10((ea_1_max-ea_1_min)/(ea_2_max-ea_2_min));
    if (is_print)
        cout << "Amplitude in dB: " << a_dB << endl;
    return a_dB;
}

double fft_analysis::amplitude_dB_with_peak_and_trough(std::vector<double> data1_1, std::vector<double> data1_2, std::vector<double> data2)
{
    vector<double> data1_1_aftermax = getElementsAfterMax(data1_1);
    vector<double> data1_2_beforemax = getElementsBeforeMax(data1_2);

    data1_1_aftermax.insert(data1_1_aftermax.end(), data1_2_beforemax.begin(), data1_2_beforemax.end());

    return amplitude_dB_with_peak_and_trough(data1_1_aftermax, data2);
}

std::pair<double, size_t> fft_analysis::findMax(const vector<double> arr, size_t size, int start_pos) {
    if (size <= start_pos)
    {
        throw std::invalid_argument("数组大小不能为零");
    }

    double max_val = arr[start_pos];
    size_t max_idx = start_pos;

    for (size_t i = start_pos+1; i < size; ++i)
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

double fft_analysis::frequency_sweeping(std::vector<double> data1, std::vector<double> data2, int start_freq, int end_freq, int freq_step, int amplitude_dB_method, int data_modify)
{

    modified_data1.clear();
    modified_data2.clear();

    bool is_valid_data = false;
    int freq_count = 0;
    int zero_count = 0;
    double frequency_tested = -1;
    int cycle_tested = -1;
    map<int, vector<double>> p_and_a;

    vector<double> data1_valid;
    vector<double> data2_valid;

    vector<double> data1_save = {0};
    vector<double> data2_save = {0};
    vector<double> data1_save_save = {0};
    vector<double> data2_save_save = {0};
    double pdd_save = 0;

    int data1_move_forward = 0;

    for (int i=0; i<(int)data2.size(); i++)
    {
        if (!data_modify)
        {
            modified_data1.push_back(data1.at(i+data1_move_forward));
            modified_data2.push_back(data2.at(i));
        }
        if (!is_valid_data)
        {
            if (data2.at(i) !=0)
            {
                for (int j=0; j<data1_move_forward; j++)
                {
                    modified_data1.at(data1.at(i+data1_move_forward-j)) = modified_data1.at(data1.at(i+data1_move_forward-j-1));
                }
                data1_move_forward = 0;
                data1_valid.push_back(data1.at(i+data1_move_forward));
                data2_valid.push_back(data2.at(i));
                is_valid_data = true;
                zero_count++;
                freq_count++;
                frequency_tested = start_freq+(freq_count-1)*freq_step;
                p_and_a.clear();
            }
            else
            {
                if (data_modify)
                {
                    modified_data1.push_back(data1.at(i+data1_move_forward));
                    modified_data2.push_back(data2.at(i));
                }

            }
        }
        else
        {
            data1_valid.push_back(data1.at(i+data1_move_forward));
            data2_valid.push_back(data2.at(i));
            if (data2.at(i) * data2.at(i+1) <0 or (data2.at(i) == 0 and data2.at(i-1) != 0))
            {
                zero_count++;

                if (zero_count%2 == 1)
                {
                    vector<double> data1_exp = resizeVector(data1_valid,256);
                    vector<double> data2_exp = resizeVector(data2_valid,256);

                    cycle_tested = zero_count/2;
                    double pdd = phase_difference_deg(data1_exp, data2_exp);
                    if (data_modify == 2)
                    {
                        while (pdd < -90 and frequency_tested <= frequency_limit)
                        {
                            data1_move_forward++;
                            data1_valid.push_back(data1.at(i+data1_move_forward));
                            data1_valid.erase(data1_valid.begin());

                            data1_exp = resizeVector(data1_valid,256);
                            pdd = phase_difference_deg(data1_exp, data2_exp);
                        }
                    }

                    double ad = 0;
                    if (amplitude_dB_method == 1)
                    {
                        ad = amplitude_dB(data1_exp, data2_exp);
                    }
                    else if (amplitude_dB_method == 2)
                    {
                        ad = amplitude_dB_with_peak(data1_valid, data2_valid);
                    }
                    else if (amplitude_dB_method == 3)
                    {
                        ad = amplitude_dB_with_trough(data1_save, data1_valid, data2_valid);
                    }
                    else if (amplitude_dB_method == 4)
                    {
                        ad = amplitude_dB_with_peak_and_trough(data1_save, data1_valid, data2_valid);
                    }
                    else if (amplitude_dB_method == 5)
                    {
                        vector<double> data1_save_exp = resizeVector(data1_save,256);
                        vector<double> data2_save_exp = resizeVector(data2_save,256);
                        data1_save_exp.insert(data1_save_exp.end(),data1_exp.begin(),data1_exp.end());
                        data2_save_exp.insert(data2_save_exp.end(),data2_exp.begin(),data2_exp.end());

                        ad = amplitude_dB(data1_save_exp, data2_save_exp);
                    }
                    else if (amplitude_dB_method == 6)
                    {
                        vector<double> data1_save_exp = resizeVector(data1_save,256);
                        vector<double> data2_save_exp = resizeVector(data2_save,256);
                        vector<double> data1_save_save_exp = resizeVector(data1_save_save,256);
                        vector<double> data2_save_save_exp = resizeVector(data2_save_save,256);
                        data1_save_exp.insert(data1_save_exp.end(),data1_exp.begin(),data1_exp.end());
                        data2_save_exp.insert(data2_save_exp.end(),data2_exp.begin(),data2_exp.end());
                        data1_save_save_exp.insert(data1_save_save_exp.end(),data1_save_exp.begin(),data1_save_exp.end());
                        data2_save_save_exp.insert(data2_save_save_exp.end(),data2_save_exp.begin(),data2_save_exp.end());

                        ad = amplitude_dB(data1_save_save_exp, data2_save_save_exp);
                    }

                    if (data_modify)
                    {
                        while (ad>0.4)
                        {
//                            cout << "start ad: " << ad;
                            if (amplitude_dB_method == 1)
                            {
                                for (int j=0; j<(int)data1_valid.size(); j++)
                                {
                                    data1_valid.at(j) = data1_valid.at(j)*0.98;
                                }

                                data1_exp = resizeVector(data1_valid,256);
                                ad = amplitude_dB(data1_exp, data2_exp);
                            }
//                            cout << ", end ad: " << ad << endl;
                        }
                        while (ad<-3 and frequency_tested <= frequency_limit)
                        {
//                            cout << "start ad: " << ad;
                            if (amplitude_dB_method == 1)
                            {
                                for (int j=0; j<(int)data1_valid.size(); j++)
                                {
                                    data1_valid.at(j) = data1_valid.at(j)*1.02;
                                }

                                data1_exp = resizeVector(data1_valid,256);
                                ad = amplitude_dB(data1_exp, data2_exp);
                            }
//                            cout << ", end ad: " << ad << endl;
                        }
                        modified_data1.insert(modified_data1.end(), data1_valid.begin(), data1_valid.end());
                        modified_data2.insert(modified_data2.end(), data2_valid.begin(), data2_valid.end());
                    }

                    if (is_print)
                    {
                        cout << "i:" << i << "*******************************************************" << endl;

                        cout << "data1_valid:" << endl;
                        for (auto data_1: data1_valid)
                        {
                            cout << data_1 << ", ";
                        }
                        cout << endl;

                        cout << "data2_valid:";
                        for (auto data_2: data2_valid)
                        {
                            cout << data_2 << ", ";
                        }
                        cout << endl;

                        cout << "frequency:" << frequency_tested << endl;
                        cout << "cycle:" << cycle_tested << endl;
                        cout << "phase_difference_deg: " << pdd << endl;
                        cout << "amplitude_dB: " << ad << endl << endl;

                    }

                    if (amplitude_dB_method == 3 or amplitude_dB_method == 4 or amplitude_dB_method == 5)
                    {
                        if (cycle_tested>0)
                        {
                            vector<double> res = {pdd_save, ad};
                            p_and_a[cycle_tested] = res;
                        }
                    }
                    else if (amplitude_dB_method == 6)
                    {
                        if (cycle_tested>1)
                        {
                            vector<double> res = {pdd_save, ad};
                            p_and_a[cycle_tested-1] = res;
                        }
                    }
                    else
                    {
                        vector<double> res = {pdd, ad};
                        p_and_a[cycle_tested] = res;
                    }

                    pdd_save = pdd;
                    data1_save_save=data1_save;
                    data2_save_save=data2_save;
                    data1_save=data1_valid;
                    data2_save=data2_valid;
                    data1_valid.clear();
                    data2_valid.clear();
                }
                else
                {

                }
            }

            if (zero_count >=11)
            {
                is_valid_data = false;
                zero_count = 0;
                frequency_sweeping_res[frequency_tested] = p_and_a;
                data1_save_save={0};
                data2_save_save={0};
                data1_save={0};
                data2_save={0};
            }
        }
    }

    return 0;
}

double fft_analysis::frequency_sweeping_phase_difference_deg(double frequency, int cycle)
{
    return frequency_sweeping_res[frequency][cycle][0];
}

double fft_analysis::frequency_sweeping_amplitude_dB(double frequency, int cycle)
{
    return frequency_sweeping_res[frequency][cycle][1];
}


std::vector<double> fft_analysis::expandVector(const std::vector<double>& input, int expansionFactor) {
    std::vector<double> result;

    if (input.empty() || expansionFactor <= 0) {
        return result;
    }

    if (expansionFactor == 1) {
        return input; // 如果扩展倍数为1，直接返回原vector
    }

    // 遍历原始数据
    for (size_t i = 0; i < input.size(); ++i) {
        // 添加当前原始数据点
        result.push_back(input[i]);

        // 如果不是最后一个元素，则插入扩展数据
        if (i < input.size() - 1) {
            // 计算当前点和下一个点之间的差值
            double diff = input[i+1] - input[i];

            // 插入扩展数据点
            for (int j = 1; j < expansionFactor; ++j) {
                double interpolatedValue = input[i] + (diff * j) / expansionFactor;
                result.push_back(interpolatedValue);
            }
        }
    }

    return result;
}

std::vector<double> fft_analysis::resizeVector(const std::vector<double>& input, size_t targetSize) {
    if (targetSize == 0) {
        return std::vector<double>();
    }

    if (input.empty()) {
        throw std::invalid_argument("Input vector cannot be empty");
    }

    if (targetSize == 1) {
        // 如果目标大小为1，返回输入向量的平均值
        double sum = 0.0;
        for (double val : input) {
            sum += val;
        }
        return std::vector<double>{sum / input.size()};
    }

    std::vector<double> result(targetSize);

    // 特殊情况：如果目标大小等于输入大小，直接返回输入
    if (targetSize == input.size()) {
        return input;
    }

    // 计算输入和输出索引之间的映射关系
    for (size_t i = 0; i < targetSize; ++i) {
        // 计算当前输出位置对应的输入位置
        double inputPos = static_cast<double>(i) * (input.size() - 1) / (targetSize - 1);

        // 找到输入位置所在的区间
        size_t lowerIndex = static_cast<size_t>(inputPos);
        size_t upperIndex = std::min(lowerIndex + 1, input.size() - 1);

        // 计算插值权重
        double weight = inputPos - lowerIndex;

        // 线性插值计算新值
        result[i] = input[lowerIndex] * (1 - weight) + input[upperIndex] * weight;
    }

    return result;
}

double fft_analysis::find_max_in_vector(std::vector<double> data1)
{
    if (data1.empty()) {
        std::cout << "向量为空" << std::endl;
        return -1;
    }

    double max = data1[0];
    for (size_t i = 1; i < data1.size(); ++i) {
        if (data1[i] > max) {
            max = data1[i];
        }
    }

    return max;
}

double fft_analysis::find_min_in_vector(std::vector<double> data1)
{
    if (data1.empty()) {
        std::cout << "向量为空" << std::endl;
        return -1;
    }

    double min = data1[0];
    for (size_t i = 1; i < data1.size(); ++i) {
        if (data1[i] < min) {
            min = data1[i];
        }
    }

    return min;
}

vector<double> fft_analysis::getElementsAfterMax(const std::vector<double> data1)
{
        std::vector<double> result;

        if (data1.empty()) {
            return result;
        }

        auto max_it = data1.begin();
        for (auto i = data1.begin(); i != data1.end(); i++) {
            if (*i >= *max_it) {
                max_it = i;
            }
        }

        if (max_it == data1.end() - 1) {
            return result;
        }

        result.assign(max_it + 1, data1.end());
        return result;
}

vector<double> fft_analysis::getElementsBeforeMax(const std::vector<double> data1)
{
        std::vector<double> result;

        if (data1.empty()) {
            return result;
        }

        auto max_it = data1.begin();
        for (auto i = data1.begin(); i != data1.end(); i++) {
            if (*i >= *max_it) {
                max_it = i;
            }
        }

        if (max_it == data1.begin()) {
            return result;
        }

        result.assign(data1.begin(), max_it-1);
        return result;
}

void fft_analysis::set_frequency_limit(double fl)
{
    frequency_limit = fl;
    return;
}
