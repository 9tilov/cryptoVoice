#include <iostream>
#include <vector>
#include <fstream>

#include "functions.h"

using namespace std;

int main(int argc, char** argv){

    FILE *standart_fp;
    FILE *test_fp;

    standart_fp = fopen(argv[1], "rb");
    test_fp = fopen(argv[2], "rb");

    std::vector<double> standart_amplitude;
    std::vector<double> test_amplitude;

	/*++++++++++++++++++++++++++++ Получаем амплитуды ++++++++++++++++++++++++++++*/
    getAmlitude(standart_fp, standart_amplitude);
    getAmlitude(test_fp, test_amplitude);
	/*++++++++++++++++++++++++++++ Получаем амплитуды ++++++++++++++++++++++++++++*/

	/*++++++++++++++++++++++++++++++++++++ Режем на фреймы +++++++++++++++++++++++++++++++++*/
	std::cout << "before : standart_size = " << standart_amplitude.size() << " test_size = " << test_amplitude.size() << std::endl;
	cutAmplitude(standart_amplitude, test_amplitude);
	std::cout << "standart_size = " << standart_amplitude.size() << " test_size = " << test_amplitude.size() << std::endl;
	std::cout << "delta = " << test_amplitude.size() / frame << std::endl;
	/*++++++++++++++++++++++++++++++++++++ Режем на фреймы +++++++++++++++++++++++++++++++++*/

    int size_vector_of_frame = (static_cast<int> (2 * standart_amplitude.size())) / frame - 1;
    
	std::vector<std::vector<double>> standart_frames;
	standart_frames.resize(size_vector_of_frame);
    
	std::vector<std::vector<double>> test_frames;
	test_frames.resize(size_vector_of_frame);
    
	get(standart_amplitude, standart_frames);
    get(test_amplitude, test_frames);

	std::vector<std::vector<double>> output_standart_frames;
	output_standart_frames.resize(size_vector_of_frame);

	std::vector<std::vector<double>> output_test_frames;
	output_test_frames.resize(size_vector_of_frame);

	for (std::size_t i = 0; i < output_standart_frames.size(); ++i){
		output_standart_frames[i].resize(frame);
		output_test_frames[i].resize(frame);
	}
	
	double t1 = clock();
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ФУРЬЁбанаротблясука ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	for (std::size_t i = 0; i < standart_frames.size(); ++i){
		FFTAnalysis(standart_frames[i], output_standart_frames[i]);
		FFTAnalysis(test_frames[i], output_test_frames[i]);
	}
	double t2 = clock();
	std::cout << "FFT done! " << t2 - t1 << " sizes : " << output_test_frames.size() << std::endl;
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ФУРЬЁ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	
	std::vector<double> standart_coefficients;
	std::vector<double> test_coefficients;

	melCepstral(output_standart_frames, standart_coefficients);
	melCepstral(output_test_frames, test_coefficients);
	 
	std::cout << "size 1 = " << standart_coefficients.size() << " size 2  = " << test_coefficients.size() << std::endl;

	double result = 0;
	result = measureFrames(standart_coefficients, test_coefficients);

	std::cout << "result = " << result << std::endl;

	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ВЫВОД ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	/*std::ofstream standart_mel("standart_mel.txt");
	std::ofstream test_mel("test_mel.txt");
	for (std::size_t i = 0; i < standart_coefficients.size(); ++i){
		standart_mel << standart_coefficients[i] << std::endl;
			test_mel << test_coefficients[i] << std::endl;
		}

	
	std::ofstream standart_result("standart_result.txt");
	std::ofstream test_result("test_result.txt");
	for (std::size_t i = 0; i < output_standart_frames.size(); ++i){
		for (std::size_t j = 0; j < output_standart_frames[i].size(); ++j){
			standart_result << output_standart_frames[i][j] << std::endl;
			test_result << output_test_frames[i][j] << std::endl;
		}
		 standart_result << "======================================================================" << std::endl;
		 test_result << "======================================================================" << std::endl;
	}

	std::ofstream standart_result1("result1.txt");
	std::ofstream test_result1("result2.txt");
	for (std::size_t i = 0; i < standart_frames.size(); ++i){
		for (std::size_t j = 0; j < standart_frames[i].size(); ++j){
			standart_result1 << standart_frames[i][j] << std::endl;
			test_result1 << test_frames[i][j] << std::endl;
		}
		 standart_result1 << "======================================================================" << std::endl;
		 test_result1 << "======================================================================" << std::endl;
	}
	std::cout << "write to file!" << std::endl;*/
	
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ВЫВОД ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


	fclose(standart_fp);
	fclose(test_fp);
	getch();
    return 0;
}