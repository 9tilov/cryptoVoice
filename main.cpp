#include <iostream>
#include <vector>
#include <fstream>

#include "functions.h"

using namespace std;

const int frames = 5646;

int main(){

    FILE *fp;
    FILE *fp1;

    fp = fopen("/home/dmitry/projects/audio/test2.wav", "rb");
    fp1 = fopen("/home/dmitry/projects/audio/parol.wav", "rb");


    std::ofstream out("result1.txt");
    if (!out){
        std::cout << "Can't create file" << std::endl;
    }

    std::vector<double> amplitude;
    std::vector<double> amplitude1;

    getAmlitude(fp, amplitude);
    getAmlitude(fp1, amplitude1);
    //    newfourierTransformWithAmplitudes(amplitude, amplitude1);


    //    double result_finish = newComparingAmplitudes(amplitude, amplitude1);
    //    std::cout << "result_finish = " << result_finish <<std::endl;


    std::cout << "first = " << amplitude.size() << std::endl;
    addZeroes(amplitude);
    addZeroes(amplitude1);

    std::cout << amplitude.size() << std::endl;

//    for (const auto& it: amplitude){
//        out << it << std::endl;
//    }
//    std::cout << "zeroes = " << amplitude.size() << std::endl;
//    int size_vector_of_frame = (static_cast<int> (2 * amplitude.size())) / frames;
//    std::vector<std::vector<double>> frames(size_vector_of_frame);
//    std::vector<std::vector<double>> frames1(size_vector_of_frame);
//    get(amplitude, frames);
//    get(amplitude1, frames1);
//    std::cout << "getFrame = " << frames.size() <<std::endl;


//    newHamming(frames);
//    newHamming(frames1);
////    Hamming(frames);
////    Hamming(frames1);
//    std::vector<std::vector<double>> fourierFrame(size_vector_of_frame);
//    std::vector<std::vector<double>> fourierFrame1(size_vector_of_frame);
//    FourierTransform(fourierFrame, frames);
//    FourierTransform(fourierFrame1, frames1);

//    std::ofstream out_frames1("frames1.txt");

//    std::vector<double> result = ComparingAmplitudes(fourierFrame, fourierFrame1);

//    double summa = 0;
//    for (std::size_t i = 0; i < result.size(); ++i){
//        summa += result[i];
//    }
//    summa /= result.size();
//    std::cout << "summa = " << summa << std::endl;
//    std::ofstream out_result("result24.txt");
//    for (std::size_t i = 0; i < fourierFrame.size(); ++i){
//        out_result << result[i] << std::endl;
//    }

    return 0;
}
