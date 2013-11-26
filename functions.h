#pragma once

#include <vector>
#include <stdio.h>
#include <iomanip>
#include <cstdint>
#include <iostream>
#include <math.h>
#include <complex>
#include <conio.h>
//#include <time.h>
//#include <list>

const double PI = 3.14159265;
const double TwoPi = 6.2831853;
const int frame = 5646;
const int coeffs = 16;
const int freq_dis = 44100;

void getAmlitude(FILE *fp, std::vector<double> &amplitude);
void getFrames(const std::vector<double> &amplitude, std::vector<std::vector<double>> &frames);
void addZeroes(std::vector<double> &amplitude);
void get(const std::vector <double> &amplitude, std::vector<std::vector<double>> &frames);
void Hamming(std::vector<std::vector<double>>& frames);
void newHamming(std::vector<std::vector<double>>& frames);
void fourierTransform(std::vector<std::vector<double>>& fourierFrame, std::vector<std::vector<double>>& frames);
void FourierTransform(std::vector<std::vector<double>>& fourierFrame, const std::vector<std::vector<double>>& frames);
std::vector<double> ComparingAmplitudes(const std::vector<std::vector<double>>& first_file, const std::vector<std::vector<double>>& second_file);
double newComparingAmplitudes(const std::vector<double>& first_file, const std::vector<double>& second_file);
void newfourierTransformWithAmplitudes(const std::vector<double>& amplitude, std::vector<double>& fourier);
void FFTAnalysis(const std::vector<double>& input, std::vector<double>& output);
void cutAmplitude(std::vector<double>& standart_ampl, std::vector<double>& test_ampl);
void melCepstral(const std::vector<std::vector<double>>& fourier, std::vector<double>& coefficients);