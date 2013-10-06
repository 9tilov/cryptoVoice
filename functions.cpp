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

void addZeroes(std::vector<double> &amplitude){
    int k = static_cast<int> (amplitude.size()) % frame;
    for (int i_add = 0; i_add < frame - k; ++i_add){
        amplitude.emplace_back(0);
    }
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


void newFourierTransform(std::vector<std::vector<double>>& fourierFrame, std::vector<std::vector<double>>& frames){
    double sum = 0;
    double even_sum = 0;
    double odd_sum = 0;
    for (std::size_t i_frames = 0; i_frames < frames.size(); ++i_frames){
        for (int k = 0; k < frame; ++k){
            for (int n = 0; n < frame; ++n){
                if (k % 2 == 0){
                    even_sum += (frames[i_frames][k]) * exp(-4 * PI * k * n / frame);
                }else{
                    odd_sum += (frames[i_frames][k]) * exp(-4 * PI * k * n / frame);
                }
                sum += even_sum + (exp(-2 * PI * k * n / frame)) * odd_sum;
            }
            fourierFrame[i_frames].push_back(sum);
            sum = 0;
        }
    }
}

void comparingAmplitudes(const std::vector<std::vector<double>>& first_file, const std::vector<std::vector<double>>& second_file){
    int k = 0;
    for (std::size_t i = 0; i < first_file.size(); ++i){
        for (std::size_t j = 0; j < first_file[i].size(); ++j){
            if (abs(first_file[i][j] - second_file[i][j]) < 0.000005){
                k++;
            }
        }
        if (k > frame - 5){
            std::cout << "OK!" << std::endl;
        }else{
            std::cout << "NE OK!" << std::endl;
        }
    }

}

std::vector<double> newComparingAmplitudes(const std::vector<std::vector<double>>& first_file, const std::vector<std::vector<double>>& second_file){
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
