#include <iostream>
#include<stdio.h>
#include<fstream>
#include<complex>
#include<cmath>
#include "fft_analysis.h"

using namespace std;
//typedef std:: complex<double> Complex;

//Complex fft[100001];

//void f(int i,int j,Complex w)
//{
//    Complex t=fft[i];
//    fft[i]=fft[i]+w*fft[j];
//    fft[j]=t-w*fft[j];
//}

//Complex w(int m,int N)
//{
//    Complex c;
//    c.real(cos(2*acos(-1)*m/N));
//    c.imag(sin(-2*acos(-1)*m/N));
//    return c;
//}

//int reverse(int n,int l)
//{
//    int i=0,r=0,t;
//    while(i<l){
//        t=n&1;
//        n=n>>1;
//        r+=t*pow(2,l-1-i);
//        i++;
//    }
//    return r;
//}

//int main(int argc, char** argv)
//{
//    //读取文件
//    ifstream infile("../data.txt",ios::in);
//    double data[100001];
//    int len=0;
//    while(infile)infile>>data[len++];
//    len-=1;
//    cout<<"原始数据："<<endl;
//    for(int i=0;i<len;i++)cout<<data[i]<<' ';
//    cout<<endl<<endl;

//    int N,p=0,k,l,i,t;
//    while(len>pow(2,p))p++;
//    N=pow(2,p);
//    for(i=len;i<N;i++)data[i]=0;
//    for(i=0;i<N;i++){
//        t=reverse(i,p);
//        fft[i].real(data[t]);
//        fft[i].imag(0);
//    }
//    for(k=1;k<=p;k++)
//        for(l=0;l<N/pow(2,k);l++)
//            for(i=0;i<pow(2,k-1);i++)
//                f(l*pow(2,k)+i,l*pow(2,k)+i+pow(2,k-1),w(i,pow(2,k)));

//    //打印结果
//    i=0;
//    cout<<"FFT处理后结果："<<endl;
//    while(i<N)cout<<fft[i++]<<' ';
//    cout << endl;
//    return 0;
//}

int main(int argc, char** argv)
{
    vector<double> data = {1, 7, 8, 2, 1, 2, 8, 2, 3, 4, 2, 8, 9, 8, 1, 2, 1};
    fft_analysis fft_a;

    fft_a.print_on();
    fft_a.fft(data);
}
