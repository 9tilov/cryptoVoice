clear all% Очистка памяти

ampl = fopen('vioce_records\den_o.txt', 'r');
key = fopen('vioce_records\den_o2.txt', 'r');
if( ampl == -1 )
  error( 'File ampl is not opened' );
end;
if( key == -1 )
  error( 'File key is not opened' );
end;


wav_signal = (0);
cnt = 1;
while ~feof(ampl)
     a = fgetl(ampl);
     wav_signal(cnt) = str2double(a);
     cnt=cnt+1;
end;
wav_signal(cnt) = 0;

wav_key = (0);
cnt = 1;
while ~feof(key)
     a = fgetl(key);
     wav_key(cnt) = str2double(a);
     cnt=cnt+1;
end;
wav_key(cnt) = 0;

%disp(length(wav_signal));
%disp(length(wav_key));

fclose(key);
fclose(ampl);

%% нормировка

% signal_mean=0;
% key_mean=0;
% ratio=0;
% for cnt = 1:1:length(wav_signal)
%     signal_mean = signal_mean + abs(wav_signal(cnt));
% end;
% signal_mean = signal_mean / length(wav_signal);
% 
% for cnt = 1:1:length(wav_key)
%     key_mean = key_mean + abs(wav_key(cnt));
% end;
% key_mean = key_mean / length(wav_key);
% 
% 
% if (signal_mean > key_mean)
%     ratio = signal_mean / key_mean;
%     disp(ratio);
%     for cnt = 1:1:length(wav_signal)
%         wav_signal(cnt) = wav_signal(cnt) / ratio;
%     end;
% else
%     ratio = key_mean / signal_mean;
%     disp(ratio);
%     for cnt = 1:1:length(wav_key)   
%         wav_signal(cnt) = wav_signal(cnt) / ratio;
%     end;
% end;


%% Hammings's windows
% for cnt = 2:1:length(wav_signal)
%     %wav_signal(cnt) = wav_signal(cnt)*(0.53836 - 0.46164*cos(cnt*2*pi/(length(wav_signal)-1)));
%     wav_signal(cnt) = (wav_signal(cnt)-0.9*wav_signal(cnt-1))*(0.54 - 0.46*cos((cnt-6)*2*pi/180));
% end;
% 
% for cnt = 2:1:length(wav_key)
%      %wav_key(cnt) = wav_key(cnt)*(0.53836 - 0.46164*cos(cnt*2*pi/(length(wav_key)-1)));
%      wav_key(cnt) = (wav_key(cnt)-0.9*wav_key(cnt-1))*(0.54 - 0.46*cos((cnt-6)*2*pi/180));     
% end;

%% Параметры
Tm=5644/44100;% Длина сигнала (с)
Fd=44100;% Частота дискретизации (Гц)
FftL=65536;% Количество линий Фурье спектра
% if (length(wav_signal) > length(wav_key))
%    Tm = length(wav_key)/Fd;
%    FftL = length(wav_key);
% else
%     Tm = length(wav_signal)/Fd;
%     FftL = length(wav_signal);
% end

%% Генерация рабочих массивов
T=0:1/Fd:Tm - 1/Fd;% Массив отсчетов времени
%disp(length(T));
%T=wav_signal;
%disp(length(T));
%disp(length(wav_signal));
Signal=wav_signal(1:length(T));
signal2=wav_key(1:length(T));
%disp(length(wav_key));

%% Спектральное представление сигнала
FftS=abs(fft(Signal, FftL));% Амплитуды преобразования Фурье сигнала
%FftS=2*FftS./FftL;% Нормировка спектра по амплитуде
%FftS(1)=FftS(1)/2;% Нормировка постоянной составляющей в спектре
FftSh=abs(fft(signal2, FftL));% Амплитуды преобразования Фурье смеси сигнал+шум
%FftSh=2*FftSh./FftL;% Нормировка спектра по амплитуде
%FftSh(1)=FftSh(1)/2;% Нормировка постоянной составляющей в спектре
%% Получаем спектральную мощность
P_ffts = FftS.^2;
%% фильтры:
P = 16;% число фильтров
Fl = Fd/215 - 1/FftL;% from 100 Hz
Fh = Fd/73 - 1/FftL;% to 1000 Hz
MEL_Fl = 1127*log(1+Fl/700);
MEL_Fh = 1127*log(1+Fh/700);
MEL_len = (MEL_Fh - MEL_Fl)/(P+1);

MEL_centres = (0);
for i=1:1:P
   MEL_centres(i) = MEL_Fl + MEL_len*i; 
end
%disp(MEL_centres);

% back to Hz
Fr_centres = 700*(exp(MEL_centres./1127)-1);
filt_smpls = Fr_centres.*((FftL/2)/Fd);
%% считаем кепстральные коэффициенты
x = (0);

for i=1:1:P
    x(i) = 0;
   for k=1:1:length(FftS)/2
       if ((i > 1) && (i < P))
           if ((k < filt_smpls(i)) && (k > filt_smpls(i-1)))          
               x(i) = x(i) + FftS(k)*((k - filt_smpls(i-1))/(filt_smpls(i) - filt_smpls(i-1)));
           elseif ((k < filt_smpls(i+1)) && (k > filt_smpls(i)))
               x(i) = x(i) + FftS(k)*((filt_smpls(i+1) - k)/(filt_smpls(i+1) - filt_smpls(i)));
           end           
       elseif (i == 1)
           if (k < filt_smpls(i))          
               x(i) = x(i) + FftS(k)*(k/filt_smpls(i));
           elseif ((k < filt_smpls(i+1)) && (k > filt_smpls(i)))
               x(i) = x(i) + FftS(k)*((filt_smpls(i+1) - k)/(filt_smpls(i+1) - filt_smpls(i)));
           end
       elseif (i == P)
           if ((k < filt_smpls(i)) && (k > filt_smpls(i-1)))          
               x(i) = x(i) + FftS(k)*((k - filt_smpls(i-1))/(filt_smpls(i) - filt_smpls(i-1)));
           end           
       end
   end
   x(i) = log(x(i));
end

c = (0);
for i=1:1:P
    c(i) = 0;
    for k=1:1:P
        c(i) = c(i) + x(k)*cos(i*(k-1/2)*pi/P);
    end
end

%disp(c);

%% то же самое для второго спектра

for i=1:1:P
    x(i) = 0;
   for k=1:1:length(FftSh)/2
       if ((i > 1) && (i < P))
           if ((k < filt_smpls(i)) && (k > filt_smpls(i-1)))          
               x(i) = x(i) + FftSh(k)*((k - filt_smpls(i-1))/(filt_smpls(i) - filt_smpls(i-1)));
           elseif ((k < filt_smpls(i+1)) && (k > filt_smpls(i)))
               x(i) = x(i) + FftSh(k)*((filt_smpls(i+1) - k)/(filt_smpls(i+1) - filt_smpls(i)));
           end           
       elseif (i == 1)
           if (k < filt_smpls(i))          
               x(i) = x(i) + FftSh(k)*(k/filt_smpls(i));
           elseif ((k < filt_smpls(i+1)) && (k > filt_smpls(i)))
               x(i) = x(i) + FftSh(k)*((filt_smpls(i+1) - k)/(filt_smpls(i+1) - filt_smpls(i)));
           end
       elseif (i == P)
           if ((k < filt_smpls(i)) && (k > filt_smpls(i-1)))          
               x(i) = x(i) + FftSh(k)*((k - filt_smpls(i-1))/(filt_smpls(i) - filt_smpls(i-1)));
           end           
       end
   end
   x(i) = log(x(i));
end

c1 = (0);
for i=1:1:P
    c1(i) = 0;
    for k=1:1:P
        c1(i) = c1(i) + x(k)*cos(i*(k-1/2)*pi/P);
    end
end

%diff_ = c-c1;
%disp(diff_);
%% kу
M_signal = 0;
M_key = 0;

for k = 1:1:P
    M_signal = M_signal + c(k);
    fprintf('%g, %g\n',c(k),c1(k));
end
M_signal = M_signal/P;
 
for k = 1:1:P
    M_key = M_key + c1(k);
end
M_key = M_key/P;

Up_sum = 0;
Left_sum = 0;
Right_sum = 0;
%disp(M_signal);
%disp(M_key);
for k = 1:1:P
    Up_sum = Up_sum + (c(k)-M_signal)*(c1(k)-M_key);
    %disp(Up_sum);
    Left_sum = Left_sum + (c(k)-M_signal)^2;
    Right_sum = Right_sum + (c1(k)-M_key)^2;
end

disp( abs(Up_sum/(sqrt(Left_sum)*sqrt(Right_sum))) );
%% Построение графиков
subplot(2,1,1);% Выбор области окна для построения
plot(T,Signal);% Построение сигнала
title('Сигнал 1');% Подпись графика
xlabel('Время (с)');% Подпись оси х графика
ylabel('Амплитуда (Попугаи)');% Подпись оси у графика
subplot(2,1,2);% Выбор области окна для построения
plot(T,signal2);% Построение смеси сигнал+шум
title('Сигнал 2');% Подпись графика
xlabel('Время (с)');% Подпись оси х графика
ylabel('Амплитуда (Попугаи)');% Подпись оси у графика

F=0:Fd/FftL:1000;%Fd/2-1/FftL;% Массив частот вычисляемого спектра Фурье
%F=0:Fd/FftL:Fd/2-1/FftL;
figure% Создаем новое окно
subplot(2,1,1);% Выбор области окна для построения
plot(F,FftS(1:length(F)));% Построение спектра Фурье сигнала
title('Спектр сигнала 1');% Подпись графика
xlabel('Частота (Гц)');% Подпись оси х графика
ylabel('Амплитуда (Попугаи)');% Подпись оси у графика
subplot(2,1,2);% Выбор области окна для построения
%F1 = 0:1:length(mels1)-1;
%disp(length(F1));
%disp(length(mels1));
plot(F,FftSh(1:length(F)));% Построение спектра Фурье сигнала
title('Спектр сигнала 2');% Подпись графика
xlabel('Частота (Гц)');% Подпись оси х графика
ylabel('Амплитуда (Попугаи)');% Подпись оси у графика

