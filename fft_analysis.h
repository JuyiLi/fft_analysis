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
#include<map>
 
#ifdef VISION_STATE
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#endif 
//#include other .h 

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

    /**
     * @brief 设置sampling_frenquency,1ms周期对应1000
     * @param sampling_frenquency_setting 设置sampling_frenquency,1ms周期对应1000
     */
    void fft_setting(int sampling_frenquency_setting);

    /**
     * @brief fft
     * @param data 需要fft处理的数组
     */
    void fft(std::vector<double> data);

    /**
     * @brief data1和data2之间的相位差,单位deg
     * @param data1 第一组数据，输出信号
     * @param data2 第二组数据，输入信号
     * @return 相位差
     */
    double phase_difference_deg(std::vector<double> data1, std::vector<double> data2);

    /**
     * @brief data1和data2之间的振幅衰减单位dB
     * @param data1 第一组数据，输出信号
     * @param data2 第二组数据，输入信号
     * @return  振幅衰减，单位dB
     */
    double amplitude_dB(std::vector<double> data1, std::vector<double> data2);

    /**
     * @brief 扫频分析，每个频率扫五个周期，每个周期之间输入信号均为0
     * @param data1，输出信号
     * @param data2，输入信号
     * @param start_freq 开始扫频的频域
     * @param end_freq 结束扫频的频域
     * @param freq_step 每次的频域增加量
     * @param amplitude_dB_method 幅值衰减算法，1：fft算法， 2：峰值算法， 3：谷值算法，4：波峰波谷振幅算法，5：双周期fft算法（2：1-2，3：2-3，4：3-4）,6:三周期fft算法（2：1-3，3：2-4，4：3-5）
     * @return
     */
    double frequency_sweeping(std::vector<double> data1, std::vector<double> data2, int start_freq, int end_freq, int freq_step, int amplitude_dB_method = 1, int data_modify = 0);

    /**
     * @brief 获取相位差。需要先执行frequency_sweeping
     * @param frequency 指定的频率
     * @param cycle 指定的周期，周期范围2-4
     * @return 特定频率和周期下的相位差
     */
    double frequency_sweeping_phase_difference_deg(double frequency, int cycle);

    /**
     * @brief 获取衰减。需要先执行frequency_sweeping
     * @param frequency 指定的频率
     * @param cycle 指定的周期，周期范围2-4
     * @return 特定频率和周期下的衰减
     */
    double frequency_sweeping_amplitude_dB(double frequency, int cycle);

    /**
     * @brief 通过波峰值判断data1和data2之间的振幅衰减单位dB
     * @param data1 输出信号
     * @param data2 输入信号
     * @return 振幅衰减，单位dB
     */
    double amplitude_dB_with_peak(std::vector<double> data1, std::vector<double> data2);

    /**
     * @brief 通过波谷值判断data1和data2之间的振幅衰减单位dB
     * @param data1 输出信号
     * @param data2 输入信号
     * @return 振幅衰减，单位dB
     */
    double amplitude_dB_with_trough(std::vector<double> data1, std::vector<double> data2);

    /**
     * @brief 通过波谷值判断data1和data2之间的振幅衰减单位dB，由于可能一个输出信号无法找到波谷，所以采用两个输出信号融合的信号寻找波谷
     * @param data1_1 输出信号前半部分
     * @param data1_2 输出信号后半部分
     * @param data2 输入信号
     * @return 振幅衰减，单位dB
     */
    double amplitude_dB_with_trough(std::vector<double> data1_1, std::vector<double> data1_2, std::vector<double> data2);

    double amplitude_dB_with_peak_and_trough(std::vector<double> data1, std::vector<double> data2);

    double amplitude_dB_with_peak_and_trough(std::vector<double> data1_1, std::vector<double> data1_2, std::vector<double> data2);


    std::vector<double> modified_data1;
    std::vector<double> modified_data2;

    /**
     * @brief set_frequency_limit
     * @param fl
     */
    void set_frequency_limit(double fl);
	
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

    std::map<double, std::map<int, std::vector<double>>> frequency_sweeping_res;

    /**
     * @brief 查找double数组的最大值及其位置
     * @param arr 要查找的double数组
     * @param size double数组的大小
     * @param start_pos 从start_pos开始寻找最大值，start_pos默认为0
     * @return 由最大值和位置组成的pair
     */
    std::pair<double, size_t> findMax(const std::vector<double> arr, size_t size, int start_pos = 0);

    void f(int i,int j,std::complex<double> w);

    std::complex<double> w(int m,int N);

    int reverse(int n,int l);

    /**
     * @brief 扩充vector
     * @param input 原始vector
     * @param expansionFactor 将原始vector的元素数扩大的倍数
     * @return 扩充后的vector
     */
    std::vector<double> expandVector(const std::vector<double>& input, int expansionFactor);

    /**
     * @brief 重新规划大小的vector
     * @param input 原始的vector
     * @param targetSize vector重新规划后的元素数
     * @return 重新规划大小后的vector
     */
    std::vector<double> resizeVector(const std::vector<double>& input, size_t targetSize);

    /**
     * @brief 找到vector中的最大值
     * @param data1 要处理的vector
     * @return vector中的最大值
     */
    double find_max_in_vector(std::vector<double> data1);

    /**
     * @brief 找到vector中的最小值
     * @param data1 要处理的vector
     * @return vector中的最小值
     */
    double find_min_in_vector(std::vector<double> data1);

    std::vector<double> getElementsAfterMax(const std::vector<double> data1);
    std::vector<double> getElementsBeforeMax(const std::vector<double> data1);

    double frequency_limit = 25;


private:
	//private members here;
};

#endif
