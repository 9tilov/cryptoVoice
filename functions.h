#pragma once

#include <vector>
#include <stdio.h>
#include <iomanip>
#include <cstdint>
#include <iostream>
#include <math.h>
#include <complex>

#include <boost/circular_buffer.hpp>

const double PI = 3.14159265;
const int frame = 5646;

void getAmlitude(FILE *fp, std::vector<double> &amplitude);
void getFrames(const std::vector<double> &amplitude, std::vector<std::vector<double>> &frames);
void addZeroes(std::vector<double> &amplitude);
void get(const std::vector <double> &amplitude, std::vector<std::vector<double>> &frames);
void Hamming(std::vector<std::vector<double>>& frames);
void fourierTransform(std::vector<std::vector<double>>& fourierFrame, std::vector<std::vector<double>>& frames);
void newFourierTransform(std::vector<std::vector<double>>& fourierFrame, const std::vector<std::vector<double>>& frames);
std::vector<double> ComparingAmplitudes(const std::vector<std::vector<double>>& first_file, const std::vector<std::vector<double>>& second_file);
double newComparingAmplitudes(const std::vector<double>& first_file, const std::vector<double>& second_file);
void fourierTransformWithAmplitudes(const std::vector<double>& amplitude, std::vector<double>& fourier);
