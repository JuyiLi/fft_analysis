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

#ifndef FFT_ANALYSIS_H
#define FFT_ANALYSIS_H

#include <iostream>
#include<stdio.h>
#include<fstream>
#include<complex>
#include<cmath>
#include<vector>

class fft_analysis
{
public:	
    //构造和析构函数
    /**
     * @brief fft_analysis
     */
    fft_analysis();

    ~fft_analysis();

	//public members here;
	//including variables and functions you want to call in scripts;

    /**
     * @brief 打开打印
     */
    void print_on();

    /**
     * @brief 关闭打印
     */
    void print_off();

    void fft_setting(int sampling_frenquency_setting);

    /**
     * @brief fft
     * @param data 需要fft处理的数组
     */
    void fft(std::vector<double> data);

    /**
     * @brief data1和data2之间的相位差,单位deg
     * @param data1 第一组数据
     * @param data2 第二组数据
     * @return 相位差
     */
    double phase_difference_deg(std::vector<double> data1, std::vector<double> data2);

    /**
     * @brief data1和data2之间的振幅衰减单位dB
     * @param data1 第一组数据
     * @param data2 第二组数据
     * @return  振幅衰减，单位dB不
     */
    double amplitude_dB(std::vector<double> data1, std::vector<double> data2);
	
protected:
	//protected members here;
	//normally, most of variables and inside functions are here;

    bool is_print = true;

    std::complex<double> fft_res[100001];

    int sampling_frenquency=1000;

    double estimated_frequency = -1;
    double estimated_amplitude = -1;
    double estimated_phase_rad = -1;
    double estimated_phase_deg = -1;

    /**
     * @brief 查找double数组的最大值及其位置
     * @param arr 要查找的double数组
     * @param size double数组的大小
     * @return 由最大值和位置组成的pair
     */
    std::pair<double, size_t> findMax(const std::vector<double> arr, size_t size);

    void f(int i,int j,std::complex<double> w);

    std::complex<double> w(int m,int N);

    int reverse(int n,int l);


private:
	//private members here;
};

#endif
