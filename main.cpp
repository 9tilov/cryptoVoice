#include <iostream>
#include <vector>
#include <fstream>

#include "functions.h"

using namespace std;

const int frames = 5646;

int main(int argc, char** argv){

    FILE *standart_fp;
    FILE *test_fp;

    standart_fp = fopen(argv[1], "rb");
    test_fp = fopen(argv[2], "rb");


    std::ofstream out("result1.txt");
    if (!out){
        std::cout << "Can't create file" << std::endl;
    }

    std::vector<double> standart_amplitude;
    std::vector<double> test_amplitude;

    getAmlitude(standart_fp, standart_amplitude);
    getAmlitude(test_fp, test_amplitude);

	cutAmplitude(standart_amplitude, test_amplitude);
	std::cout << "standart_size = " << standart_amplitude.size() << " test_size = " << test_amplitude.size() << std::endl;
	std::cout << "delta = " << test_amplitude.size() / frame << std::endl;



    //    newfourierTransformWithAmplitudes(amplitude, amplitude1);


    //    double result_finish = newComparingAmplitudes(amplitude, amplitude1);
    //    std::cout << "result_finish = " << result_finish <<std::endl;


   // addZeroes(amplitude);
   // addZeroes(amplitude1);




//    for (const auto& it: amplitude){
//        out << it << std::endl;
//    }
//    std::cout << "zeroes = " << amplitude.size() << std::endl;
    int size_vector_of_frame = (static_cast<int> (2 * standart_amplitude.size())) / frames - 1;
    std::vector<std::vector<double>> standart_frames(size_vector_of_frame);
    std::vector<std::vector<double>> test_frames(size_vector_of_frame);
    get(standart_amplitude, standart_frames);
    get(test_amplitude, test_frames);

	std::vector<std::vector<double>> output_standart_frames((int)standart_amplitude.size(), 0);
	std::vector<std::vector<double>> output_test_frames((int)standart_amplitude.size(), 0);

	for (std::size_t i = 0; i < output_standart_frames.size(); ++i){
		output_standart_frames[i].resize(frame);
		output_test_frames[i].resize(frame);
	}

	for (std::size_t i = 0; i < standart_frames.size(); ++i){
		FFTAnalysis(standart_frames[i], output_standart_frames[i], frame, frame);
		FFTAnalysis(test_frames[i], output_test_frames[i], frame, frame);
	}

	for(std::size_t i = 0; i < standart_frames.size(); ++i){
		//std::cout << standart_frames[i].size() << std::endl;
	}

	for(std::size_t i = 0; i < test_frames.size(); ++i){
		//std::cout << test_frames[i].size() << std::endl;
	}
	std::vector<double> coefficients;
	melCepstral(output_standart_frames, coefficients);
	
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
	std::ofstream standart_result("standart_result.txt");
	std::ofstream test_result("test_result.txt");
	for (std::size_t i = 0; i < standart_frames.size(); ++i){
		for (std::size_t j = 0; j < standart_frames[i].size(); ++j){
			standart_result << standart_frames[i][j] << std::endl;
			test_result << test_frames[i][j] << std::endl;
		}
		 standart_result << "======================================================================" << std::endl;
		 test_result << "======================================================================" << std::endl;
	}

	fclose(standart_fp);
	fclose(test_fp);
	getch();
    return 0;
}