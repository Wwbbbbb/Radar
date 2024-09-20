clc;
close all;

%������������
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

%��������
% radarData = readDCA1000("testRange2M470F.bin");
radarData = readDCA1000("BR20bpm-HR82bpm-1TX4RX.bin");%  4RX*��256������*128chirps*525frames��
radarData1RX = reshape(radarData(1,:),numOfADCSamples,numOfAllChirps);%  256������*(128chirps*525frames��
radarData1RX1Chirp = zeros(numOfADCSamples,numOfAllChirps/numOfChirpsInFrame);%ȡÿһ֡�м��chirp ȡ��64��cirp
for chirpIndix = 64:128:numOfAllChirps
    radarData1RX1Chirp(:,(chirpIndix-64)/numOfChirpsInFrame+1) = radarData1RX(:,chirpIndix);% 256������*525chirps
end

[X,Y] = meshgrid((0:numOfADCSamples-1)*deltaR,(1:numOfFrames));%ÿ֡��ȡһ��chirp ֡��=chirp��
radarData1RX64thChirp_FFT = zeros(numOfFrames,numOfADCSamples);
for fftIndix = 1:numOfFrames
    radarData1RX64thChirp_FFT(fftIndix,:) = abs(fft(radarData1RX1Chirp(:,fftIndix)));
end
figure;
mesh(X,Y,radarData1RX64thChirp_FFT);%��Ƶ��ͼ�е�Ƶ�ʷ�����תΪƵ�ʲƵ�ʲ�С��Ƶ�ʷ������������Ͻ�
xlabel('����ά��-m');
ylabel('������-n');
zlabel('����');
title('���������ά-FFTͼ');

%������λ��Ϣ
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
            %��һ������ά��525��chirp�����ۼ� ����ػ���
            radarFFTtemp(samplesIndix) = radarFFTtemp(samplesIndix) + radarData1RX64thChirp_FFT(chirpIndix,samplesIndix);
        end
        if(radarFFTtemp(samplesIndix) > rangMaxFFT)
            rangMaxFFT = radarFFTtemp(samplesIndix);%��������FFT���Ȼ���ֵ
            rangMaxIndex = samplesIndix;%������������
        end
    end
end

figure;
subplot(2, 2, 1);
plot(phaseData(rangMaxIndex,:));
xlabel('��������');
% ylabel('��λ');
ylabel('��һ����λ');
title('δչ����λ�ź�');

%================================================������λBEGIN==========================================================%
%��λ���
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
phaseToHOC = phaseTounwrap; %��HOC�任����
maxlag = fix((length(phaseToHOC)-1)/2);         % ����ͺ�
% maxlag = 200;
nsamp = 256;          % ÿ���ε�������
overlap = 50;         % �ص��ٷֱȣ�50%��
flag = 'unbiased';    % ��������
k1 = -1;               % �̶��ͺ� k1
k2 = 1;               % �̶��ͺ� k2
% ���ú��������Ľ��ۻ���
phaseTounwrap_HOC = cum4est(phaseToHOC, maxlag, nsamp, overlap, flag, k1, k2);
% phaseToDiff_HOC = cum4est(phaseToDiff, fix((length(phaseToDiff)-1)/2), 256, 50, 'unbiased', 0, 0);

subplot(2, 2, 2);
% plot(phaseTounwrap);
% hold on;
% plot(phaseTounwrap_HOC,'-r');
plot(rescale(phaseTounwrap,-1,1));
hold on;
plot(rescale(phaseTounwrap_HOC,-1,1),'-r');
xlabel('��������');
% ylabel('��λ');
ylabel('��һ����λ');
legend('��������','HOC����');
title('չ����λ�ź�');


%��λ��� �õ���λ��
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
xlabel('��������');
ylabel('��λ');
legend('��������','HOC����');
title('��λ���ź�');


%��ƽ���˲�
phaseSmooth = smoothdata(phaseToDiff,'movmean',5);
phaseSmooth_HOC = smoothdata(phaseToDiff_HOC,'movmean',5);
subplot(2, 2, 4);
% plot(phaseSmooth);
% hold on;
% plot(phaseSmooth_HOC,'-r');
plot(rescale(phaseSmooth,-1,1));
hold on;
plot(rescale(phaseSmooth_HOC,-1,1),'-r');
xlabel('��������');
% ylabel('��λ');
ylabel('��һ����λ');
legend('��������','HOC����');
title('ƽ����λ���ź�');
%================================================������λEND==========================================================%

ff = (0:numOfFrames-1)/numOfFrames*25;%ÿ֡����Ϊ40ms   1s/40ms=25Hz
phaseSmooth_FFT = abs(fft(phaseSmooth));
phaseSmooth_HOC_FFT = abs(fft(phaseSmooth_HOC));
figure;
plot(ff(1:(numOfFrames-1)/2),rescale(phaseSmooth_FFT(1:(numOfFrames-1)/2),0,1));
hold on;
plot(ff(1:(numOfFrames-1)/2),rescale(phaseSmooth_HOC_FFT(1:(numOfFrames-1)/2),0,1),'-r');
xlabel('Ƶ��-Hz');
% ylabel('����');
ylabel('��һ������');
legend('��������','HOC����');
title('��λ�ź�-FFT');