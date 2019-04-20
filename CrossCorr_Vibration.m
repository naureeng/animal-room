% February 27, 2018 Data
% Automate process

load('avg_signal.mat');
time_stamp_vib = duration([09 34 57]);

% read in first video and parse time stamp
x = dir('*.txt'); 
fname = x(9).name; % first file name
hours_vid = str2double(fname(17:18));
mins_vid = str2double(fname(20:21));
secs_vid = str2double(fname(23:24));
time_stamp_vid = duration([hours_vid mins_vid secs_vid]);

% align by removing data in vibration
time_diff = abs(time_stamp_vid-time_stamp_vib);

% automate time data
str = split(char(time_diff));
newStr = split(str,':');
nums = cell2mat(newStr);
hours = str2double(nums(1,:));
mins = str2double(nums(2,:));
secs = str2double(nums(3,:));

% use time data to align
time_sec = (hours*60*60) + (mins*60) + (secs);
num_pts = round(time_sec./0.5178);
one_hour = round((60*60)./0.5178);
avg_vib_hour = avg_signal(num_pts:num_pts+round(one_hour));

% filter
signal = avg_vib_hour;
Fs = 800; % MMA8451 Default Rate = 800 Hz
Fn = Fs/2; % Nyquist freq
N = length(signal); % time

% remove 50 Hz noise
Wo = 50/(Fs/2); BW = Wo/35; 
[b,a] = iirnotch(Wo,BW);
Y = filter(b,a,signal); 

% design bandpass filter
Wp = [1.0 50]/Fn; % passband freq
Ws = [0.5 51]/Fn; % stopband freq
Rp = 1; % passband ripple (dB)
Rs = 150; % stopband ripple(dB)
[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs); % filter order
[z,p,k] = cheby2(n,Rs,Ws); % filter design
[sosbp,gbp] = zp2sos(z,p,k); % convert to second-order-section for stability

% freqz(sosbp, 2^16, Fs); % filter bode plot
Y1 = filtfilt(sosbp,gbp,Y);
vib_sec = resample(Y1, 30*60*60, length(Y1));
Max_vib = LocalMinima(-vib_sec,1000,-mean(vib_sec)-3*std(vib_sec));
Min_vib = LocalMinima(vib_sec,1000,mean(vib_sec)+3*std(vib_sec));
vib_peaks = unique([Max_vib;Min_vib]);

% read in video data
A = readmatrix(fname);
A = A(:,1); % take 1st col only
vid_sec = resample(A,10000,length(A));
vid_sec = resample(vid_sec,30*60*60,length(vid_sec));
Max_vid = LocalMinima(-vid_sec,1000,-mean(vid_sec)-3*std(vid_sec));
Min_vid = LocalMinima(vid_sec,1000,mean(vid_sec)+3*std(vid_sec));
vid_peaks = unique([Max_vid;Min_vid]);

[xcf,lags,bounds,h] = crosscorr(vid_peaks,vib_peaks);

% quantify positive signficant area
pos_sig_area = trapz(xcf)-( bounds(1)*40 );















