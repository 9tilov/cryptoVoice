#include <iostream>
#include <vector>
#include <fstream>

#include "functions.h"

using namespace std;

const int frames = 5646;

int main(){

    FILE *fp;
    FILE *fp1;
    fp = fopen("/home/dmitry/projects/audio/test5.wav", "rb");
    fp1 = fopen("/home/dmitry/projects/audio/test2.wav", "rb");

    std::ofstream out("result1.txt");
    if (!out){
        std::cout << "Can't create file" << std::endl;
    }

    std::vector<double> amplitude;
    std::vector<double> amplitude1;

    getAmlitude(fp, amplitude);
    getAmlitude(fp1, amplitude1);
    std::ofstream out_amplitude("amp1.txt");
    for (std::size_t i = 0; i < amplitude.size(); ++i){
        out_amplitude << amplitude[i] << std::endl;
    }
    std::ofstream out_amplitude1("amp2.txt");
    for (std::size_t i = 0; i < amplitude.size(); ++i){
        out_amplitude1 << amplitude1[i] << std::endl;
    }

    std::cout << "first = " << amplitude.size() << std::endl;
    addZeroes(amplitude);
    addZeroes(amplitude1);
    for (const auto& it: amplitude){
        out << it << std::endl;
    }
    std::cout << "zeroes = " << amplitude.size() << std::endl;
    int size_vector_of_frame = (static_cast<int> (2 * amplitude.size())) / frames;
    std::vector<std::vector<double>> frames(size_vector_of_frame);
    std::vector<std::vector<double>> frames1(size_vector_of_frame);
    get(amplitude, frames);
    get(amplitude1, frames1);
    std::cout << "getFrame = " << frames.size() <<std::endl;
    std::ofstream out_frames1("frames1.txt");
    for (std::size_t i = 0; i < frames.size(); ++i){
        for (std::size_t j = 0; j < frames[i].size(); ++j){
            out_frames1 << frames[i][j] << std::endl;
        }
        out_frames1 << std::endl;
    }
    Hamming(frames);
    Hamming(frames1);
    std::vector<std::vector<double>> fourierFrame(size_vector_of_frame);
    std::vector<std::vector<double>> fourierFrame1(size_vector_of_frame);
    //    fourierTransform(fourierFrame, frames);
    //    fourierTransform(fourierFrame1, frames1);
    newFourierTransform(fourierFrame, frames);
    newFourierTransform(fourierFrame1, frames1);
    std::vector<double> result = newComparingAmplitudes(fourierFrame, fourierFrame1);
    double summa = 0;
    for (std::size_t i = 0; i < result.size(); ++i){
        summa += result[i];
    }
    summa /= frame;
    std::cout << "summa = " << summa << std::endl;
    std::ofstream out_result("result24.txt");
    for (std::size_t i = 0; i < fourierFrame.size(); ++i){
        out_result << result[i] << std::endl;
    }

    //    comparingAmplitudes(fourierFrame, fourierFrame1);
    std::cout << "framesHamming = " << frames.size() << std::endl;

    std::ofstream out_fourer("fourer.txt");
    for (std::size_t i = 0; i < fourierFrame.size(); ++i){
        for (std::size_t j = 0; j < fourierFrame[i].size(); ++j){
            out_fourer << fourierFrame[i][j] << std::endl;
        }
        out_fourer<< std::endl;
    }

    std::ofstream out_frames2("frames2.txt");
    for (std::size_t i = 0; i < frames.size(); ++i){
        for (std::size_t j = 0; j < frames[i].size(); ++j){
            out_frames2 << frames[i][j] << std::endl;
        }
        out_frames2 << std::endl;
    }

    std::ofstream out_frames3("frames3.txt");
    for (std::size_t i = 0; i < frames.size(); ++i){
        for (std::size_t j = 0; j < frames[i].size(); ++j){
            out_frames3 << frames[i][j] << std::endl;
        }
        out_frames3 << std::endl;
    }
    return 0;
}
