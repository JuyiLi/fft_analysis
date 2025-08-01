#include <iostream>
#include<stdio.h>
#include<fstream>
#include<complex>
#include<cmath>
#include "fft_analysis.h"
#include<vector>

using namespace std;

int main(int argc, char** argv)
{
    vector<double> data = {1, 7, 8, 2, 1, 2, 8, 2, 3, 4, 2, 8, 9, 8, 1, 2, 1};
    fft_analysis fft_a;

    fft_a.print_on();
    fft_a.fft(data);
}
