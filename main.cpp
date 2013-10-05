#include <iostream>
#include <vector>
#include <fstream>

#include <functions.h>

using namespace std;

const int frames = 5646;

int main(){

    FILE *fp;
    fp = fopen("/home/dmitry/projects/audio/speech.1.wav", "rb");
    std::ofstream out("result.txt");
    if (!out){
        std::cout << "Can't create file" << std::endl;
    }

    std::vector<double> amplitude;
    getAmlitude(fp, amplitude);
    std::ofstream out_amplitude("amp.txt");
    for (std::size_t i = 0; i < amplitude.size(); ++i){
        out_amplitude << amplitude[i] << std::endl;
    }
    std::cout << "first = " << amplitude.size() << std::endl;
    addZeroes(amplitude);
    for (const auto& it: amplitude){
        out << it << std::endl;
    }
    std::cout << "zeroes = " << amplitude.size() << std::endl;
    int size_vector_of_frame = (static_cast<int> (2 * amplitude.size())) / frames;
    std::vector<std::vector<double>> frames(size_vector_of_frame);
    get(amplitude, frames);
    std::cout << "getFrame = " << frames.size() <<std::endl;
    std::ofstream out_frames1("frames1.txt");
    for (std::size_t i = 0; i < frames.size(); ++i){
        for (std::size_t j = 0; j < frames[i].size(); ++j){
            out_frames1 << frames[i][j] << std::endl;
        }
        out_frames1 << std::endl;
    }
    Hamming(frames);
    std::vector<std::vector<double>> fourierFrame(size_vector_of_frame);
    fourierTransform(fourierFrame, frames);
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
