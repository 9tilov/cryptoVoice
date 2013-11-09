#include "functions.h"

void getAmlitude(FILE *fp, std::vector<double> &amplitude){
    char byte_1, byte_2;
    int byteInt_1, byteInt_2, byteSum;
    fseek(fp, 0, SEEK_END);
    int sizeOfDataFile = ftell(fp) - 44;

    fseek(fp, 44, SEEK_SET);

    for (int i = 0; i < sizeOfDataFile / 2; i += 2){
        fscanf(fp, "%c%c", &byte_1, &byte_2);
        byteInt_1 = static_cast<int> (byte_1);
        byteInt_2 = static_cast<int> (byte_2);
        //        std::cout << ((int)byte_1) << "  " << ((int)byte_2) << std::endl;
        byteSum = (byteInt_2 << 8) | (byteInt_1 & 0x000000FF);
        double amp = byteSum / 32768.0;
        amplitude.push_back(amp);
        fseek (fp, 46 + 2 * i, SEEK_SET);
    }
}

//void addZeroes(std::vector<double> &amplitude){
//    int k = static_cast<int> (amplitude.size()) % frame;
//    for (int i_add = 0; i_add < frame - k; ++i_add){
//        amplitude.emplace_back(0);
//    }
//}

void addZeroes (std::vector<double> &amplitudes){
    int a = 1, degree = 0;
    for (degree = 0;;++degree){
        if (a > (int)amplitudes.size()){
            break;
        }
        a <<= 1;
    }
    int delta = pow(2, degree) - amplitudes.size();
    for (int i = 0; i < delta; ++i){
        amplitudes.push_back(0);
    }
    std::cout << degree << std::endl;
}

void getFrames(const std::vector<double> &amplitude, std::vector<std::vector<double>> &frames){
    int j = 0;
    int size_of_vector_frame = (static_cast<int> (amplitude.size())) / frame;
    while (j < size_of_vector_frame){
        for (std::size_t ampl = 1; ampl < amplitude.size() + 1; ++ampl){
            frames[j].push_back(amplitude[ampl - 1]);
            if (ampl % frame == 0){
                j++;
            }
        }
    }
}

//\ TODO для облегчения расчетов

void get(const std::vector <double> &amplitude, std::vector<std::vector<double>> &frames){
    std::size_t k = 0;
    int left = 0, right = frame;
    while (k < 2 * amplitude.size() / frame){
        for (int i = left + 1; i < right + 1; ++i){
            frames[k].push_back(amplitude[i - 1]);
        }
        left += frame / 2;
        right += frame / 2;
        k++;
    }
}

void Hamming(std::vector<std::vector<double>>& frames){
    for (std::size_t i = 0; i < frames.size(); ++i){
        for (std::size_t j = 0; j < frames[i].size(); ++j){
            frames[i][j] *= (0.53836 - 0.46164 * cos((2 * PI * j) / ((frame - 1) / 2)));
        }
    }
}

void newHamming(std::vector<std::vector<double>>& frames){
    for (std::size_t i = 0; i < frames.size(); ++i){
        for (std::size_t j = 1; j < frames[i].size(); ++j){
            frames[i][j] = (frames[i][j] - 0.9 * frames[i][j - 1]) * (0.53836 - 0.46164 * cos((2 * PI * j) / 180));
        }
    }
}

void fourierTransform(std::vector<std::vector<double>>& fourierFrame, std::vector<std::vector<double>>& frames){
    double sum = 0;
    for (std::size_t i_frames = 0; i_frames < frames.size(); ++i_frames){
        for (int k = 0; k < frame; ++k){
            for (int n = 0; n < frame; ++n){
                sum += (frames[i_frames][k]) * (cos(-2 * PI * k * n / frame));
            }
            fourierFrame[i_frames].push_back(sum);
            sum = 0;
        }
    }
}


void FourierTransform(std::vector<std::vector<double>>& fourierFrame, const std::vector<std::vector<double>>& frames){
    double sum = 0;
    double even_sum = 0;
    double odd_sum = 0;
    for (std::size_t i_frames = 0; i_frames < frames.size(); ++i_frames){
        for (int k = 0; k < frame; ++k){
            for (int n = 0; n < frame; ++n){
                if (n % 2 == 0){
                    even_sum += (frames[i_frames][n]) * exp(-4 * PI * k * n / frame);
                }else{
                    odd_sum += (frames[i_frames][n]) * exp(-4 * PI * k * n / frame);
                }
                sum += even_sum + (exp(-2 * PI * k * n / frame)) * odd_sum;
            }
            fourierFrame[i_frames].push_back(sum);
            sum = 0;
        }
    }
}



std::vector<double> ComparingAmplitudes(const std::vector<std::vector<double>>& first_file, const std::vector<std::vector<double>>& second_file){
    std::vector<double> first_sum, second_sum, result;
    double summa_first = 0, summa_second = 0, high_sum = 0, low_left_sum = 0, low_right_sum = 0;
    for (std::size_t i = 0; i < first_file.size(); ++i){
        for (std::size_t j = 0; j < second_file.size(); ++j){
            summa_first += first_file[i][j];
            summa_second += second_file[i][j];
        }
        summa_first /= frame;
        summa_second /= frame;
        first_sum.push_back(summa_first);
        second_sum.push_back(summa_second);
    }
    for (std::size_t i = 0; i < first_file.size(); ++i){
        for (std::size_t j = 0; j < second_file.size(); ++j){
            high_sum += (first_file[i][j] - first_sum[i]) * (second_file[i][j] - second_sum[i]);
            low_left_sum += pow((first_file[i][j] - first_sum[i]), 2);
            low_right_sum += pow((second_file[i][j] - second_sum[i]), 2);
        }
        low_left_sum = sqrt(low_left_sum);
        low_right_sum = sqrt(low_right_sum);
        result.push_back(fabs(high_sum / (low_left_sum * low_right_sum)));
        high_sum = 0;
        low_left_sum = 0;
        low_right_sum = 0;
    }
    return result;
}


//void newfourierTransformWithAmplitudes(const std::vector<double>& amplitude, std::vector<double>& fourier){
//    double sum = 0;
//    std::cout << "enter to fourier" << std::endl;
//    double even_sum = 0;
//    double odd_sum = 0;
//    for (std::size_t i = 0; i < amplitude.size(); ++i){
//        for (std::size_t k = 0; k < amplitude.size(); ++k){
//            if (k % 2 == 0){
//                even_sum += (amplitude[k]) * exp(-4 * PI * i * k / frame);
//            }else{
//                odd_sum += (amplitude[k]) * exp(-4 * PI * i * k / frame);
//            }
//            sum += even_sum + (exp(-2 * PI * k * i / frame)) * odd_sum;
//        }
//        fourier.push_back(sum);
//        sum = 0;
//        std::cout << "push" << " i = " << i << std::endl;
//    }
//    std::cout << "exit fourier" << std::endl;
//}

//double newComparingAmplitudes(const std::vector<double>& first_file, const std::vector<double>& second_file){
//    double result;
//    std::cout << "endter to comparing" << std::endl;
//    double summa_first = 0, summa_second = 0, high_sum = 0, low_left_sum = 0, low_right_sum = 0;
//    for (std::size_t i = 0; i < first_file.size(); ++i){
//        summa_first += first_file[i];
//        summa_second += second_file[i];
//    }
//    summa_first /= first_file.size();
//    summa_second /= second_file.size();
//    for (std::size_t i = 0; i < first_file.size(); ++i){
//        high_sum += (first_file[i] - summa_first) * (second_file[i] - summa_second);
//        low_left_sum += pow((first_file[i] - summa_first), 2);
//        low_right_sum += pow((second_file[i] - summa_second), 2);
//    }
//    low_left_sum = sqrt(low_left_sum);
//    low_right_sum = sqrt(low_right_sum);
//    result = fabs(high_sum / (low_left_sum * low_right_sum));
//    std::cout << "exit comparing" << std::endl;
//    return result;

//}
