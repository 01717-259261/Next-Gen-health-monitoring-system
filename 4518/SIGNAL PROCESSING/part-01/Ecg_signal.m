clear all
close all
clc

load ('noisyecg.mat');
Ecgsignal = noisyECG_withTrend/200;
fs = 360;
t = (0:length(Ecgsignal)-1)/fs;
subplot(5,1,1)
plot(t,Ecgsignal)
xlabel('time (s)');
ylabel('Amplitude mv');
title('Raw ECG Signal')
f_low  = 5;
f_high = 15;
[b,a] = butter(2, [f_low,f_high]/(fs/2),'bandpass');
filtered_signal = filtfilt(b,a,Ecgsignal);
subplot(5,1,2)
plot(t,filtered_signal)
xlabel('time (s)');
ylabel('Amplitude mv');
title('Filtered Signal')

%Normalization of ECG signal in between range [-1,1]
normalized_signal = (filtered_signal - min(filtered_signal)) / (max(filtered_signal) - min(filtered_signal))*2-1;

%PLot of normalized signal
subplot(5,1,3)
plot(t,normalized_signal)
xlabel('time (s)');
ylabel('Amplitude mv');
title('Normalized Filtered Signal')

% Detection of R peaks for finding heart rate
% wt = modwt(Ecgsignal,4,'sym4');
% wtrec = zeros(size(wt));
% wtrec(3:4,:) = wt(3:4,:); %extractiong coefficients d1 d2
% y = imodwt(wtrec,'sym4'); %Inverse DWT with d1 d2
% y = abs(y).^2;
wt = modwt(Ecgsignal,5);
wtrec = zeros(size(wt));
wtrec(4:5,:) = wt(4:5,:);
y = imodwt(wtrec,'sym4');

% Finding peaks

timelimit = length(Ecgsignal)/fs;
 y = abs(y).^2;
avg = mean(y);
[qrspeaks,locs] = findpeaks(y,t,'MinPeakHeight',avg,'MinPeakDistance',0.150);

%Plotting

subplot(5,1,4)
plot(t,y)
hold on
plot(locs,qrspeaks,'ro')
xlabel('Seconds')
title('R Peaks Localized by Wavelet Transform with Automatic Annotations')
nobh = length(locs);
heart_beats = (nobh*60)/timelimit;
disp(strcat('Heart Rate = ',num2str(heart_beats)))
title(strcat('PQRST Detection and Heart Rate : ',num2str(heart_beats)))
subplot(5,1,5)
plot(t,y,'g')
hold on
plot(t,filtered_signal,'r')

for i = 1:nobh
    % Convert locs to integers using round or floor, depending on your needs.
    start_idx = (locs(i));
    end_idx = (locs(i) + 15);
    
    % Make sure the indices are within the valid range.
    if start_idx < 1
        start_idx = 1;
    end
    if end_idx > length(y)
        end_idx = length(y);
    end
    
    Speaks(i) = min(y(start_idx:end_idx));
    locs_s(i) = find(y == Speaks(i)) - 1;
end
for i = 1:nobh
    % Convert locs to integers using round or floor, depending on your needs.
    start_idx = (locs(i));
    end_idx = (locs(i) - 15);
    
    % Make sure the indices are within the valid range.
    if start_idx < 1
        start_idx = 1;
    end
    if end_idx < 1
        end_idx = 1;
    end
    
    Qpeaks(i) = min(y(end_idx:start_idx));
    locs_q(i) = find(y == Qpeaks(i));
end

for i = 1:nobh
    if locs_q(i) - 60 > 0
        start_idx = (locs_q(i) - 60);
        end_idx = (locs_q(i));
    else
        start_idx = (locs_q(i) - 10);
        end_idx = (locs_q(i));
    end
    
    % Make sure the indices are within the valid range.
    if start_idx < 1
        start_idx = 1;
    end
    if end_idx > length(y)
        end_idx = length(y);
    end
    
    Ppeaks(i) = max(y(start_idx:end_idx));
    locs_p(i) = find(y == Ppeaks(i));
end
for i = 1:nobh
    if locs_s(i) + 60 <= length(y)
        start_idx = (locs_s(i));
        end_idx = (locs_s(i) + 60);
    else
        start_idx = (locs_s(i));
        end_idx = length(y);
    end
    
    % Make sure the indices are within the valid range.
    if start_idx < 1
        start_idx = 1;
    end
    if end_idx > length(y)
        end_idx = length(y);
    end
    
    Tpeaks(i) = max(y(start_idx:end_idx));
    locs_t(i) = find(y == Tpeaks(i));
end

figure(2)
plot(t,normalized_signal )
grid on
xlim([0,nobh])
hold on
% xlim([0,length(Ecgsignal)])
% plot(locs_r,Rpeaks,'^r');
plot(locs_s,Speaks,'xm')
plot(locs_q,Qpeaks,'*g')
plot(locs_p,Ppeaks,'om')
plot(locs_t,Tpeaks,'pk')

xlabel('samples')
legend('','Speaks','Qpeaks','Ppeaks','Tpeaks')