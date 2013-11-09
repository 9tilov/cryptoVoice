clear all% Очистка памяти

ampl = fopen('artem3.txt', 'r');
key = fopen('dima3.txt', 'r');
if( ampl == -1 )
  error( 'File ampl is not opened' )
end;
if( key == -1 )
  error( 'File key is not opened' )
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

fclose(key);
fclose(ampl);
%% Hammings's windows
%for cnt = 1:1:length(wav_signal)
%    wav_signal(cnt) = wav_signal(cnt)*(0.53836 - 0.46164*cos(cnt*2*pi/(length(wav_signal)-1)));
%end;

%for cnt = 1:1:length(wav_key)
%    wav_key(cnt) = wav_key(cnt)*(0.53836 - 0.46164*cos(cnt*2*pi/(length(wav_key)-1)));
%end;

%% Параметры
Tm=3;% Длина сигнала (с)
Fd=44100;% Частота дискретизации (Гц)
FftL=65536;% Количество линий Фурье спектра

%% Генерация рабочих массивов
T=0:1/Fd:Tm;% Массив отсчетов времени
%T=wav_signal;
%disp(length(T));
%disp(length(wav_signal));
Signal=wav_signal;%Ak+A1*sind((F1*360).*T+Phi1)+A2*sind((F2*360).*T+Phi2);% Массив сигнала (смесь 2х синусоид и постоянной составляющей)
signal2=wav_key;
%disp(length(wav_key));

%% Спектральное представление сигнала
FftS=abs(fft(Signal,FftL));% Амплитуды преобразования Фурье сигнала
FftS=2*FftS./FftL;% Нормировка спектра по амплитуде
FftS(1)=FftS(1)/2;% Нормировка постоянной составляющей в спектре
FftSh=abs(fft(signal2,FftL));% Амплитуды преобразования Фурье смеси сигнал+шум
FftSh=2*FftSh./FftL;% Нормировка спектра по амплитуде
FftSh(1)=FftSh(1)/2;% Нормировка постоянной составляющей в спектре
%% Получаем спектральную мощность
P_ffts = FftS.^2;
%% фильтры:
P = 20;% число фильтров
Fl = 0;
Fh = Fd/2 - 1/FftL;
MEL_Fl = 0;
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
%% kуй
M_signal = 0;
M_key = 0;

for k = 1:1:P
    M_signal = M_signal + c(k);
end
M_signal = M_signal/P;
 
for k = 1:1:P
    M_key = M_key + c1(k);
end
M_key = M_key/P;

Up_sum = 0;
Left_sum = 0;
Right_sum = 0;
disp(M_signal);
disp(M_key);
for k = 1:1:P
    Up_sum = Up_sum + (c(k)-M_signal)*(c1(k)-M_key);
    %disp(Up_sum);
    Left_sum = Left_sum + (c(k)-M_signal)^2;
    Right_sum = Right_sum + (c1(k)-M_key)^2;
end

disp( abs(Up_sum/(sqrt(Left_sum)*sqrt(Right_sum))) );

%% Определяем похожесть
%M_signal = 0;
%M_key = 0;
%size = 1487;% 1000 Hz
%size = 1041;% 700 Hz
%size = 892;% 600 Hz

%for k = 1:1:size
%    M_signal = M_signal + FftS(k);
%end
%M_signal = M_signal/size;
 
%for k = 1:1:size
%    M_key = M_key + FftSh(k);
%end
%M_key = M_key/size;

% убираем частоты, амплитуда которых меньше средней в спетре (Зачем я это делаю?)
%for cnt = 1:1:size
%   if (FftS(cnt) < M_signal)
%       FftS(cnt) = 0;
%   end
%   if (FftSh(cnt) < M_signal)
%       FftSh(cnt) = 0;
%   end
%end

%Up_sum = 0;
%Left_sum = 0;
%Right_sum = 0;
%for k = 1:1:size
%    Up_sum = Up_sum + (FftS(k)-M_signal)*(FftSh(k)-M_key);
%    Left_sum = Left_sum + (FftS(k)-M_signal)^2;
%    Right_sum = Right_sum + (FftSh(k)-M_key)^2;
%end

%disp( (Up_sum/(sqrt(Left_sum)*sqrt(Right_sum))) );

%% Построение графиков
subplot(2,1,1);% Выбор области окна для построения
%plot(T,Signal);% Построение сигнала
title('Сигнал 1');% Подпись графика
xlabel('Время (с)');% Подпись оси х графика
ylabel('Амплитуда (Попугаи)');% Подпись оси у графика
subplot(2,1,2);% Выбор области окна для построения
%plot(T,signal2);% Построение смеси сигнал+шум
title('Сигнал 2');% Подпись графика
xlabel('Время (с)');% Подпись оси х графика
ylabel('Амплитуда (Попугаи)');% Подпись оси у графика

%F=0:Fd/FftL:Fd/2-1/FftL;% Массив частот вычисляемого спектра Фурье
F=0:Fd/FftL:Fd/2-1/FftL;
figure% Создаем новое окно
subplot(2,1,1);% Выбор области окна для построения
%plot(F,FftS(1:length(F)));% Построение спектра Фурье сигнала
title('Спектр сигнала 1');% Подпись графика
xlabel('Частота (Гц)');% Подпись оси х графика
ylabel('Амплитуда (Попугаи)');% Подпись оси у графика
subplot(2,1,2);% Выбор области окна для построения
F1 = 0:1:length(mels1)-1;
disp(length(F1));
disp(length(mels1));
%plot(F1,mels1);% Построение спектра Фурье сигнала
title('Спектр сигнала 2');% Подпись графика
xlabel('Частота (Гц)');% Подпись оси х графика
ylabel('Амплитуда (Попугаи)');% Подпись оси у графика

