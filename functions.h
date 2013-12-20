#pragma once

#include <fstream>
#include <vector>
#include <stdio.h>
#include <iomanip>
#include <cstdint>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <ctime>
#include <cstdlib>
#include <conio.h>
#include "bass.h"

const double PI = 3.14159265;
const double TwoPi = 6.2831853;
const int frame = 8192;
const int coeffs = 20;
const int freq_dis = 44100;
const double limit = 1.8;
const int CHANS = 2;
const int BUFSTEP = 170000;
const int times = 10;

static int device;
static int record_input;
static DWORD reclen;
static HRECORD rchan = 0;
static HSTREAM chan=0;
static char *recbuf=NULL;

void getAmplitude(char* ampl, std::vector<double>& amplitude);
void getFrames(const std::vector<double> &amplitude, std::vector<std::vector<double>> &frames);
void addZeroes(std::vector<double> &amplitude);
void get(const std::vector <double> &amplitude, std::vector<std::vector<double>> &frames);
void Hamming(std::vector<double>& frames);
void newHamming(std::vector<std::vector<double>>& frames);
void fourierTransform(std::vector<std::vector<double>>& fourierFrame, std::vector<std::vector<double>>& frames);
void FourierTransform(std::vector<std::vector<double>>& fourierFrame, const std::vector<std::vector<double>>& frames);
std::vector<double> ComparingAmplitudes(const std::vector<std::vector<double>>& first_file, const std::vector<std::vector<double>>& second_file);
double newComparingAmplitudes(const std::vector<double>& first_file, const std::vector<double>& second_file);
void newfourierTransformWithAmplitudes(const std::vector<double>& amplitude, std::vector<double>& fourier);
void FFTAnalysis(const std::vector<double>& input, std::vector<double>& output);
void cutAmplitude(std::vector<double>& standart_ampl, std::vector<double>& test_ampl);
void melCepstral(const std::vector<double>& fourier, std::vector<double>& coefficients);
double measureFrames(const std::vector<double>& standart_sample, const std::vector<double>& test_sample);
double townMeasure(const std::vector<double>& standart_sample, const std::vector<double>& test_sample);
double deltaMeasure(const std::vector<double>& standart_sample, const std::vector<double>& test_sample);
double summMeasure(const std::vector<double>& standart_sample, const std::vector<double>& test_sample);
void normalAmplitudes(std::vector<double>& amplitude);
double degreeMeasure(const std::vector<double>& standart_sample, const std::vector<double>& test_sample);
BOOL CALLBACK MyRecordProc(HRECORD handle, const void *buffer, DWORD length, void *user);
BOOL CALLBACK RecordingCallback(HRECORD handle, const void *buffer, DWORD length, void *user);
void StartRecording();
void StopRecording();
BOOL InitDevice(int device);
char* Record();
