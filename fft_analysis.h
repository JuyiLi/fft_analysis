/**
 * @file fft_analyse
 * @brief fft数据处理库，和matlab中fft的功能相似
 *
 * Copyright (c) 2025 [JUYI LI].
 * All rights reserved. Licensed under the MIT License.
 *
 * author: JUYI LI
 * created: 2025-07-29
 */

#ifndef FFT_ANALYSE_H
#define FFT_ANALYSE_H

#include <iostream>
#include<stdio.h>
#include<fstream>
#include<complex>
#include<cmath>
 
#ifdef VISION_STATE
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#endif 
//#include other .h 

class fft_analyse
{
public:	
    //构造和析构函数
    /**
     * @brief fft_analyse
     */
    fft_analyse();

    ~fft_analyse();

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
     * @param data_len 数组的长度
     */
    void fft(double * data, int data_len);
	
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
    std::pair<double, size_t> findMax(const double arr[], size_t size);

    void f(int i,int j,std::complex<double> w);

    std::complex<double> w(int m,int N);

    int reverse(int n,int l);


private:
	//private members here;
};

#endif
