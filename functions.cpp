#include "functions.h"
//void getAmlitude(FILE *fp, std::vector<double> &amplitude){
//    char byte_1, byte_2;
//    int byteInt_1, byteInt_2, byteSum;
//    fseek(fp, 0, SEEK_END);
//    int sizeOfDataFile = ftell(fp) - 44;
//
//    fseek(fp, 44, SEEK_SET);
//    amplitude.reserve(140000);
//    for (int i = 0; i < sizeOfDataFile / 2; i += 2){
//        fscanf(fp, "%c%c", &byte_1, &byte_2);
//        byteInt_1 = static_cast<int> (byte_1);
//        byteInt_2 = static_cast<int> (byte_2);
//        byteSum = (byteInt_2 << 8) | (byteInt_1 & 0x000000FF);
//        double amp = byteSum / 32768.0;
//
//        amplitude.push_back(amp);
//        fseek (fp, 46 + 2 * i, SEEK_SET);
//    }
//    amplitude.erase(amplitude.begin(), amplitude.begin() + 9000);
//    amplitude.resize(57000);
//}

void getAmplitude(char* ampl, std::vector<double>& amplitude) {
	char byte_1, byte_2;
	int byteInt_1, byteInt_2, byteSum;
	amplitude.reserve(170004);
	std::cout << reclen << std::endl;
	for (int i = 0; i < (int)reclen; i += 4){
		byteInt_1 = (int) *(ampl + i);
		byteInt_2 = (int) *(ampl + i + 1);
		byteSum = (byteInt_2 << 8) | (byteInt_1 & 0x000000FF);
		double amp = byteSum / 32768.0;
		amplitude.push_back(amp);
		std::cout << "i = " << i << std::endl;
	}
	amplitude.erase(amplitude.begin(), amplitude.begin() + 9000);
	amplitude.resize(57000);
}

void StopRecording()
{
	std::cerr<<"stop";
	BASS_ChannelStop(rchan);
	rchan = 0;
	BASS_Free();
}

BOOL CALLBACK MyRecordProc(HRECORD handle, const void* buffer, DWORD length, void* user)
{
	std::cerr<<"me_record";
	if (reclen >= BUFSTEP){	
		std::cout << "100%" << std::endl;
		//recbuf = "\0";
		StopRecording();
		std::cerr << reclen; 
		return FALSE;
	}
	memcpy(recbuf + reclen, buffer, (reclen + length >= BUFSTEP)? (BUFSTEP - reclen) : length);
	std::cout << ((double)reclen / BUFSTEP) * 100 << "%" << std::endl;
	reclen += length;
	return TRUE;
}

void StartRecording()
{
	std::cerr<<"start";
	WAVEFORMATEX *wf;
	if (recbuf) { // free old recording
		BASS_StreamFree(chan);
		chan=0;
		free(recbuf);
		recbuf=NULL;
		BASS_Free();
	}
	recbuf=(char*)calloc(BUFSTEP, 1);

	DWORD retval = BASS_ChannelIsActive(rchan);
	if (retval == BASS_ACTIVE_STOPPED){
		std::cout << "BASS_ACTIVE_STOPPED" << std::endl;
	} else if (retval == BASS_ACTIVE_PLAYING){
		std::cout << "BASS_ACTIVE_PLAYING" << std::endl;
	}else if (retval == BASS_ACTIVE_PAUSED){
		std::cout << "BASS_ACTIVE_PAUSED" << std::endl;
	}else if (retval == BASS_ACTIVE_STALLED){
		std::cout << "BASS_ACTIVE_STALLED" << std::endl;
	}

	std::cout<< "rchan = " << rchan << std::endl;

	rchan=BASS_RecordStart(freq_dis,CHANS,0,MyRecordProc,0);
	std::cout<< "rchan = " << rchan << std::endl;

	float fft[128]; // fft data buffer
	for (int a=0; a<128; a++)
		retval = BASS_ChannelIsActive(rchan);
	if (retval == BASS_ACTIVE_STOPPED){
		std::cout << "BASS_ACTIVE_STOPPED" << std::endl;
	} else if (retval == BASS_ACTIVE_PLAYING){
		std::cout << "BASS_ACTIVE_PLAYING" << std::endl;
	}else if (retval == BASS_ACTIVE_PAUSED){
		std::cout << "BASS_ACTIVE_PAUSED" << std::endl;
	}else if (retval == BASS_ACTIVE_STALLED){
		std::cout << "BASS_ACTIVE_STALLED" << std::endl;
	}
	BASS_SetVolume(0.5);
	std::cout << "volume 1 = " << BASS_SetVolume(0.5) << std::endl;

	std::cout << "rchan = " << rchan << std::endl;
	if (!rchan) {
		std::cout << "Couldn't start recording" << std::endl;
		free(recbuf);
		recbuf=0;
		return;
	} else {
		std::cout << "In StartRecording()" << std::endl;
	}
}


BOOL InitDevice(int device)
{
	std::cerr<<"init";
	BASS_RecordFree(); // free current device (and recording channel) if there is one
	if (!BASS_RecordInit(device)) {
		std::cout<<"Can't initialize recording device"<<std::endl;
		return FALSE;
	}
	{ // get list of inputs
		int c;
		const char *i;
		for (c=0;i=BASS_RecordGetInputName(c);c++) {
			if (!(BASS_RecordGetInput(c,NULL)&BASS_INPUT_OFF)) { // this one is currently "on"
				record_input=c;
			}
		}
	}
	return TRUE;
}

char* Record() {
	std::cerr<< "record";
	BASS_DEVICEINFO di;
	int def;
	for (int c = 0; BASS_RecordGetDeviceInfo(c, &di); ++c) {
		if (di.flags&BASS_DEVICE_DEFAULT) { // got the default device
			def = c;
		}
	}
	std::cout << "Default device : " << def << std::endl;
	InitDevice(def);
	std::cout << "InitDevice : " << def <<" ...";
	std::cout<< ((InitDevice(def)==TRUE)? "ok":"not ok" )<< std::endl;
	for (int i = 0; BASS_RecordSetInput(i, BASS_INPUT_OFF, -1); ++i); // 1st disable all inputs, then...
	std::cout << "BASS_RecordSetInput : " << def <<" ..." << ((BASS_RecordSetInput(def,BASS_INPUT_ON,-1)==TRUE)? "ok":"not ok" )<< std::endl;
	WAVEFORMATEX *wf;
	if (recbuf) { // free old recording
		BASS_StreamFree(chan);
		chan=0;
		free(recbuf);
		recbuf=NULL;
		BASS_Free();
	}
	recbuf=(char*)malloc(BUFSTEP);
	reclen = 0;
	/*wf=(WAVEFORMATEX*)(recbuf+20);
	wf->wFormatTag=1;
	wf->nChannels=CHANS;
	wf->wBitsPerSample=16;
	wf->nSamplesPerSec=freq_dis;
	wf->nBlockAlign=wf->nChannels*wf->wBitsPerSample/8;
	wf->nAvgBytesPerSec=wf->nSamplesPerSec*wf->nBlockAlign;
	FILE *my_record = fopen("jabloko_denis2.wav", "wb");
	fwrite(recbuf, reclen, 1, my_record);*/
	std::cout << "Запись начнется через ";
	for (int i = 3; i > 0; --i) { 
		std::cout << i << " ";
		Sleep(1000);
	}
	std::cout << "Говорите!" << std::endl;

	HRECORD record=BASS_RecordStart(44100, 2, 0, MyRecordProc, recbuf); // start recording
	std::cerr<<"end";
	return recbuf;
}

void cutAmplitude(std::vector<double>& standart_ampl, std::vector<double>& test_ampl){
	if (standart_ampl.size() > test_ampl.size()) {
		int total = test_ampl.size() / frame;
		standart_ampl.resize(total * frame);
		test_ampl.resize(total * frame);
	} else {
		int total = standart_ampl.size() / frame;
		test_ampl.resize(total * frame);
		standart_ampl.resize(total * frame);
	}
}

void addZeroes (std::vector<double> &amplitudes){
	int a = 1, degree = 0;
	for (degree = 0;;++degree){
		if (a > (int)amplitudes.size()){
			break;
		}
		a <<= 1;
	}
	int delta = (int)(pow(2.0, (double)degree) - (double)amplitudes.size());
	for (int i = 0; i < delta; ++i){
		amplitudes.push_back(0);
	}
}

void getFrames(const std::vector<double> &amplitude, std::vector<std::vector<double>> &frames){
	int j = 0;
	int size_of_vector_frame = (static_cast<int>(amplitude.size())) / frame;
	frames.reserve(size_of_vector_frame);
	while (j < size_of_vector_frame){
		for (std::size_t ampl = 1; ampl < amplitude.size() + 1; ++ampl){
			frames[j].push_back(amplitude[ampl - 1]);
			if (ampl % frame == 0){
				j++;
			}
		}
	}
}

void normalAmplitudes(std::vector<double>& amplitude){
	double sum = 0;
	for (std::size_t i = 0; i < amplitude.size(); ++i) {
		sum += pow(amplitude[i], 2);
	}
	sum = sqrt(sum);
	for (std::size_t i = 0; i < amplitude.size(); ++i) {
		amplitude[i] /= sum;
	}
}

void get(const std::vector <double> &amplitude, std::vector<std::vector<double>> &frames){
	std::size_t left = 0, right = frame;

	for (std::size_t i = 0; i < frames.size(); ++i){
		frames[i].reserve(frame);
	}

	for (std::size_t k = 0; k < 2 * amplitude.size() / frame; ++k){
		if (right <= amplitude.size()){
			for (std::size_t i = left; i < right; ++i){
				frames[k].push_back(amplitude[i]);
			}
			left += frame / 2;
			right += frame / 2;
		}
	}
}

void Hamming(std::vector<double>& frames){
	for (std::size_t i = 0; i < frames.size(); ++i){
		frames[i] *= (0.53836 - 0.46164 * cos((2 * PI * i) / ((frame - 1) / 2)));
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

void FFTAnalysis(const std::vector<double>& input, std::vector<double>& output) {
	int i = 0, j = 0, n = 0, m = 0, Mmax = 0, Istp = 0, size = (int)input.size();
	double Tmpr = 0, Tmpi = 0, Wtmp = 0, Theta = 0;
	double Wpr = 0, Wpi = 0, Wr = 0, Wi = 0;

	double *temp;
	n = size * 2;

	temp = new double[2 * n];

	for (int g = 0; g < 2 * n; ++g){
		temp[g] = 0;
	}

	for (i = 0; i < size; ++i) {
		j = i * 2;
		temp[j] = 0;
		temp[j+1] = input[i];
	}
	i = 1;
	j = 1;

	while (i < n) {
		if (j > i) {
			Tmpr = temp[i];
			temp[i] = temp[j];
			temp[j] = Tmpr;
			Tmpr = temp[i+1];
			temp[i+1] = temp[j+1];
			temp[j+1] = Tmpr;
		}
		i = i + 2;
		m = size;
		while ((m >= 2) && (j > m)) {
			j = j - m;
			m = m >> 1;
		}
		j = j + m;
	}
	Mmax = 2;
	while (n > Mmax) {
		Theta = -TwoPi / Mmax;
		Wpi = sin(Theta);
		Wtmp = sin(Theta / 2);
		Wpr = Wtmp * Wtmp * 2;
		Istp = Mmax * 2;
		Wr = 1;
		Wi = 0;
		m = 1;
		while (m < Mmax) {
			i = m;
			m = m + 2;
			Tmpr = Wr;
			Tmpi = Wi;
			Wr = Wr - Tmpr * Wpr - Tmpi * Wpi;
			Wi = Wi + Tmpr * Wpi - Tmpi * Wpr;

			while (i < n) {
				j = i + Mmax;
				Tmpr = Wr * temp[j] - Wi * temp[j-1];
				Tmpi = Wi * temp[j] + Wr * temp[j-1];

				temp[j] = temp[i] - Tmpr; temp[j-1] = temp[i-1] - Tmpi;
				temp[i] = temp[i] + Tmpr; temp[i-1] = temp[i-1] + Tmpi;
				i = i + Istp;
			}
		}
		Mmax = Istp;
	}

	for (i = 0; i < size; ++i) {
		j = i * 2;
		output[size - i - 1] = sqrt(pow(temp[j], 2) + pow(temp[j + 1], 2));
	}
	//output.resize((int)input.size());
	delete []temp;
}

void melCepstral(const std::vector<double>& fourier, std::vector<double>& coefficients) {
	coefficients.reserve(coeffs);
	double freq_low = freq_dis / 300, freq_high = freq_dis / 42, mel_low = 1127 * log(1 + freq_low / 700), mel_high = 1127 * log(1 + freq_high / 700), mel_dist = (mel_high - mel_low) / (coeffs + 1);

	double mel_centries[coeffs], freq_centies[coeffs], x[coeffs], temp_coefficients = 0;
	int freq_samples[coeffs];
	for (int i = 0; i < coeffs; ++i) {
		x[i] = 0;
	}
	for (int i = 0; i < coeffs; ++i){
		mel_centries[i] = mel_low + mel_dist * (i + 1);
		freq_centies[i] = 700 * (exp(mel_centries[i] / 1127) - 1);
		freq_samples[i] = (int)((freq_centies[i]) * ((int)fourier.size()) / (2 * freq_dis));
	}

	for (int i = 0; i < coeffs; ++i) {
		if (i == 0) {
			for (int k = freq_samples[i]; k < freq_samples[i + 1]; ++k) {
				x[i] += (fourier[k]) * (freq_samples[i + 1] - k) / (freq_samples[i + 1] - freq_samples[i]);
			}
		} else if (i == coeffs - 1) {
			for (int k = freq_samples[i - 1]; k < freq_samples[i]; ++k) {
				x[i] += (fourier[k]) * (k - freq_samples[i - 1]) / (freq_samples[i] - freq_samples[i - 1]);
			}
		} else {
			for (int k = freq_samples[i - 1]; k < freq_samples[i + 1]; ++k) {
				if ((k < freq_samples[i]) && (k > freq_samples[i - 1]))
					x[i] += (fourier[k]) * (k - freq_samples[i - 1]) / (freq_samples[i] - freq_samples[i - 1]);
				else if ((k > freq_samples[i]) && (k < freq_samples[i + 1]))
					x[i] += (fourier[k]) * (freq_samples[i + 1] - k) / (freq_samples[i + 1] - freq_samples[i]);
			}
		}
	}
	for (int i = 0; i < coeffs; ++i) {
		x[i] = log(x[i]);
	}
	for (int j = 1; j <= coeffs; ++j) {
		temp_coefficients = 0;
		for (int k = 1; k <= coeffs; ++k){
			temp_coefficients += x[k - 1] * cos(j * (k - 1.0/2) * PI / coeffs);
		}
		coefficients.push_back(temp_coefficients);
	}
}

double measureFrames(const std::vector<double>& standart_sample, const std::vector<double>& test_sample) {
	double result = 0, up_sum = 0, left_sum = 0, right_sum = 0, standart_expected_value = 0, test_expected_value = 0;
	for (std::size_t i = 0; i < standart_sample.size(); ++i) {
		standart_expected_value += standart_sample[i];
		test_expected_value += test_sample[i];
	}
	standart_expected_value /= standart_sample.size();
	test_expected_value /= test_sample.size();

	for (std::size_t i = 0; i < standart_sample.size(); ++i) {
		up_sum += (standart_sample[i] - standart_expected_value) * (test_sample[i] - test_expected_value);
		left_sum += pow((standart_sample[i] - standart_expected_value), 2);
		right_sum += pow((test_sample[i] - test_expected_value), 2);
	}

	result = fabs(up_sum / (sqrt(left_sum) * sqrt(right_sum)));
	return result;
}

double townMeasure(const std::vector<double>& standart_sample, const std::vector<double>& test_sample) {
	double temp = 0;
	for (std::size_t i = 0; i < standart_sample.size(); ++i) {
		temp += fabs(standart_sample[i] - test_sample[i]);
	}

	return temp;
}

double deltaMeasure(const std::vector<double>& standart_sample, const std::vector<double>& test_sample) {
	std::vector<double> delta;
	double sum = 0, result = 0;
	delta.reserve((int) standart_sample.size());
	for (std::size_t i = 0; i < standart_sample.size(); ++i) {
		delta.push_back(fabs(fabs(standart_sample[i])-fabs(test_sample[i])));
		sum += delta[i];
	}
	result = sum / (double)(delta.size());
	return result;
}

double summMeasure(const std::vector<double>& standart_sample, const std::vector<double>& test_sample) {
	double sum1 = 0, sum2 = 0, result = 0;
	for (std::size_t i = 0; i < standart_sample.size(); ++i) {
		sum1 += fabs(standart_sample[i]);
		sum2 += fabs(test_sample[i]);
	}
	result = fabs(sum1 - sum2);
	return result;
}

double degreeMeasure(const std::vector<double>& standart_sample, const std::vector<double>& test_sample) {
	double sum = 0, p = 2.0, r = 5.0;
	for(std::size_t i = 0; i < standart_sample.size(); ++i) {
		sum += exp(p*log(fabs(standart_sample[i] - test_sample[i])));
	}
	sum = exp((1/r) * log(sum));
	return sum;
}