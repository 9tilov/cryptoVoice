#include "functions.h"

using namespace std;

int main(int argc, char** argv){

    FILE *standart_fp;
    FILE *test_fp;


    if (argc != 3) {
        //std::cout << "Enter correct path\n";
	return 1;
    }
    standart_fp = fopen(argv[1], "rb");
    test_fp = fopen(argv[2], "rb");
   if (standart_fp == NULL) {
        std::cout << "file 1 not found" << std::endl;
        exit (-1);
    } else if (test_fp == NULL)        {
        std::cout << "file 2 not found" << std::endl;
        exit (-1);
    }
    std::vector<double> standart_amplitude;
    std::vector<double> test_amplitude;
    /*++++++++++++++++++++++++++++ Ïîëó÷àåì àìïëèòóäû ++++++++++++++++++++++++++++*/
    getAmlitude(standart_fp, standart_amplitude);
    getAmlitude(test_fp, test_amplitude);
    //std::cout << standart_amplitude.size() << std::endl;
    /*++++++++++++++++++++++++++++ Ïîëó÷àåì àìïëèòóäû ++++++++++++++++++++++++++++*/
    //standart_amplitude.resize(65536);
    //test_amplitude.resize(65536);
    addZeroes(standart_amplitude);
    addZeroes(test_amplitude);

    std::vector<double> output_standart_frames;
    output_standart_frames.resize((int)standart_amplitude.size());

    std::vector<double> output_test_frames;
    output_test_frames.resize((int)standart_amplitude.size());


    double t1 = clock();
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ÔÓÐÜ¨áàíàðîòáëÿñóêà ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    FFTAnalysis(standart_amplitude, output_standart_frames);
    FFTAnalysis(test_amplitude, output_test_frames);
    double t2 = clock();
    //std::cout << "FFT done! " << t2 - t1 << " sizes : " << standart_amplitude.size() << " sdf = " << output_standart_frames.size() << std::endl;
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ÔÓÐÜ¨ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//    output_standart_frames.resize((int)standart_amplitude.size());
//    output_test_frames.resize((int)standart_amplitude.size());
    std::vector<double> standart_coefficients;
    std::vector<double> test_coefficients;

    melCepstral(output_standart_frames, standart_coefficients);
    melCepstral(output_test_frames, test_coefficients);
	
    double sum1 = 0, sum2 = 0;
    for (std::size_t i = 0; i < standart_coefficients.size(); ++i) {
    sum1 += fabs(standart_coefficients[i]);
    sum2 += fabs(test_coefficients[i]);
//        std::cout << std::setw(9) << standart_coefficients[i] << " % " << test_coefficients[i] << std::endl;
	
    }
    //std::cout << "s1 = " << sum1 << " s2 = " << sum2 << std::endl;

//    std::cout << "size 1 = " << standart_coefficients.size() << " size 2  = " << test_coefficients.size() << std::endl;

    double result = 0;
    //result = measureFrames(standart_coefficients, test_coefficients);
    //result = townMeasure(standart_coefficients, test_coefficients);
    result = deltaMeasure(standart_coefficients, test_coefficients);
   // result = summMeasure(standart_coefficients, test_coefficients);
    std::cout << std::left << std::setw(25)<< argv[1] << argv[2] << "    result = " << std::setw(8) << result << ((result < 2.5)? " -ok":"") << std::endl;

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ÂÛÂÎÄ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//    std::ofstream standart_mel("standart_mel.txt");
//    std::ofstream test_mel("test_mel.txt");
//    for (std::size_t i = 0; i < standart_coefficients.size(); ++i){
//        standart_mel << standart_coefficients[i] << std::endl;
//        test_mel << test_coefficients[i] << std::endl;
//    }


//    std::ofstream standart_result("standart_result.txt");
//    std::ofstream test_result("test_result.txt");
//    for (std::size_t i = 0; i < standart_amplitude.size(); ++i){
//        standart_result << standart_amplitude[i] << std::endl;
//        test_result << test_amplitude[i] << std::endl;
//    }

//    std::ofstream standart_result1("result1.txt");
//    std::ofstream test_result1("result2.txt");
//    for (std::size_t i = 0; i < standart_frames.size(); ++i){
//        for (std::size_t j = 0; j < standart_frames[i].size(); ++j){
//            standart_result1 << standart_frames[i][j] << std::endl;
//            test_result1 << test_frames[i][j] << std::endl;
//        }
//        standart_result1 << "======================================================================" << std::endl;
//        test_result1 << "======================================================================" << std::endl;
//    }
    //std::cout << "write to file!" << std::endl;

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ÂÛÂÎÄ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


    fclose(standart_fp);
    fclose(test_fp);
    return 0;
}
