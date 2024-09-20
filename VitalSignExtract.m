clc;
close all;

%基本参数设置
numOfADCSamples = 256;
% numOfFrames = 470;
numOfChirpsInFrame = 128;
adcSampleByte = 2;
IQ = 2;
numOfAllChirps = numOfChirpsInFrame*numOfFrames;
numOfADCBits = 16;
numOfTX = 1;
numOfRX = 4;
adcBinInfo = dir('BR20bpm-HR82bpm-1TX4RX.bin');
adcBinSize = adcBinInfo.bytes;
numOfFrames = adcBinSize/numOfADCSamples/numOfChirpsInFrame/numOfTX/numOfRX/adcSampleByte/IQ;
fs = 10000e3;%10000ksps
c = 3e8;
ts = numOfADCSamples/fs;
freqSlope = 30e12;%30MHz/us=30e3GHz/s
bandWidth = ts*freqSlope;
deltaR = c/(2*bandWidth);
startFrequency = 77e9;

%数组重组
% radarData = readDCA1000("testRange2M470F.bin");
radarData = readDCA1000("BR20bpm-HR82bpm-1TX4RX.bin");%  4RX*（256采样点*128chirps*525frames）
radarData1RX = reshape(radarData(1,:),numOfADCSamples,numOfAllChirps);%  256采样点*(128chirps*525frames）
radarData1RX1Chirp = zeros(numOfADCSamples,numOfAllChirps/numOfChirpsInFrame);%取每一帧中间的chirp 取第64个cirp
for chirpIndix = 64:128:numOfAllChirps
    radarData1RX1Chirp(:,(chirpIndix-64)/numOfChirpsInFrame+1) = radarData1RX(:,chirpIndix);% 256采样点*525chirps
end

[X,Y] = meshgrid((0:numOfADCSamples-1)*deltaR,(1:numOfFrames));%每帧中取一个chirp 帧数=chirp数
radarData1RX64thChirp_FFT = zeros(numOfFrames,numOfADCSamples);
for fftIndix = 1:numOfFrames
    radarData1RX64thChirp_FFT(fftIndix,:) = abs(fft(radarData1RX1Chirp(:,fftIndix)));
end
figure;
mesh(X,Y,radarData1RX64thChirp_FFT);%其频谱图中的频率分量可转为频率差，频率差小的频率分量代表与距离较近
xlabel('距离维度-m');
ylabel('脉冲数-n');
zlabel('幅度');
title('多脉冲距离维-FFT图');

%计算相位信息
phaseData = zeros(numOfADCSamples,numOfFrames);
for chirpIndix = 1:numOfFrames
    for samplesIndix = 1:numOfADCSamples
        phaseData(samplesIndix,chirpIndix) = angle(radarData1RX1Chirp(samplesIndix,chirpIndix));
    end
end

rangMaxIndex = 0;
rangMaxFFT = 0;
radarFFTtemp = zeros(1,numOfADCSamples);
for samplesIndix = 1:numOfADCSamples
    if(samplesIndix*deltaR < 1.5)
        for chirpIndix = 1:numOfFrames
            %对一个距离维的525个chirp幅度累加 非相关积累
            radarFFTtemp(samplesIndix) = radarFFTtemp(samplesIndix) + radarData1RX64thChirp_FFT(chirpIndix,samplesIndix);
        end
        if(radarFFTtemp(samplesIndix) > rangMaxFFT)
            rangMaxFFT = radarFFTtemp(samplesIndix);%保存最大的FFT幅度积累值
            rangMaxIndex = samplesIndix;%保存最大的索引
        end
    end
end

figure;
subplot(2, 2, 1);
plot(phaseData(rangMaxIndex,:));
xlabel('采样点数');
% ylabel('相位');
ylabel('归一化相位');
title('未展开相位信号');

%================================================处理相位BEGIN==========================================================%
%相位解缠
phaseTounwrap = phaseData(rangMaxIndex,:);
for chirpIndix = 1:numOfFrames-1
    phaseDiffTemp = phaseTounwrap(chirpIndix + 1) - phaseTounwrap(chirpIndix);
    while((phaseDiffTemp > pi) || (phaseDiffTemp < -pi))
        if phaseDiffTemp > pi
            phaseTounwrap(chirpIndix + 1) = phaseTounwrap(chirpIndix + 1) - 2*pi;
            phaseDiffTemp = phaseTounwrap(chirpIndix + 1) - phaseTounwrap(chirpIndix);
        elseif phaseDiffTemp < -pi
            phaseTounwrap(chirpIndix + 1) = phaseTounwrap(chirpIndix + 1) + 2*pi;
            phaseDiffTemp = phaseTounwrap(chirpIndix + 1) - phaseTounwrap(chirpIndix);
        end
    end
end

%HOC
phaseToHOC = phaseTounwrap; %待HOC变换向量
maxlag = fix((length(phaseToHOC)-1)/2);         % 最大滞后
% maxlag = 200;
nsamp = 256;          % 每个段的样本数
overlap = 50;         % 重叠百分比（50%）
flag = 'unbiased';    % 计算类型
k1 = -1;               % 固定滞后 k1
k2 = 1;               % 固定滞后 k2
% 调用函数计算四阶累积量
phaseTounwrap_HOC = cum4est(phaseToHOC, maxlag, nsamp, overlap, flag, k1, k2);
% phaseToDiff_HOC = cum4est(phaseToDiff, fix((length(phaseToDiff)-1)/2), 256, 50, 'unbiased', 0, 0);

subplot(2, 2, 2);
% plot(phaseTounwrap);
% hold on;
% plot(phaseTounwrap_HOC,'-r');
plot(rescale(phaseTounwrap,-1,1));
hold on;
plot(rescale(phaseTounwrap_HOC,-1,1),'-r');
xlabel('采样点数');
% ylabel('相位');
ylabel('归一化相位');
legend('正常方法','HOC方法');
title('展开相位信号');


%相位差分 得到相位差
phaseToDiff = zeros(1,numOfFrames);
phaseToDiff_HOC = zeros(1,numOfFrames);
for chirpIndix = 1:numOfFrames-1
    phaseToDiff(chirpIndix+1) = phaseTounwrap(chirpIndix+1) - phaseTounwrap(chirpIndix);
    phaseToDiff_HOC(chirpIndix+1) = phaseTounwrap_HOC(chirpIndix+1) - phaseTounwrap_HOC(chirpIndix);
end

subplot(2, 2, 3);
% plot(phaseToDiff);
% hold on;
% plot(phaseToDiff_HOC,'-r');
plot(rescale(phaseToDiff,-1,1));
hold on;
plot(rescale(phaseToDiff_HOC,-1,1),'-r');
xlabel('采样点数');
ylabel('相位');
legend('正常方法','HOC方法');
title('相位差信号');


%做平滑滤波
phaseSmooth = smoothdata(phaseToDiff,'movmean',5);
phaseSmooth_HOC = smoothdata(phaseToDiff_HOC,'movmean',5);
subplot(2, 2, 4);
% plot(phaseSmooth);
% hold on;
% plot(phaseSmooth_HOC,'-r');
plot(rescale(phaseSmooth,-1,1));
hold on;
plot(rescale(phaseSmooth_HOC,-1,1),'-r');
xlabel('采样点数');
% ylabel('相位');
ylabel('归一化相位');
legend('正常方法','HOC方法');
title('平滑相位差信号');
%================================================处理相位END==========================================================%

ff = (0:numOfFrames-1)/numOfFrames*25;%每帧周期为40ms   1s/40ms=25Hz
phaseSmooth_FFT = abs(fft(phaseSmooth));
phaseSmooth_HOC_FFT = abs(fft(phaseSmooth_HOC));
figure;
plot(ff(1:(numOfFrames-1)/2),rescale(phaseSmooth_FFT(1:(numOfFrames-1)/2),0,1));
hold on;
plot(ff(1:(numOfFrames-1)/2),rescale(phaseSmooth_HOC_FFT(1:(numOfFrames-1)/2),0,1),'-r');
xlabel('频率-Hz');
% ylabel('幅度');
ylabel('归一化幅度');
legend('正常方法','HOC方法');
title('相位信号-FFT');