p = pwd;
addpath(genpath(d))

%% load demo data where the GC is from 1 to 2
load demo_data.mat;

% lfp: T x TRLS
% spk: structure for spike timing for TRLS

%% multitaper method parameter settings
mtm.Fs=1000; % sampling frequency
mtm.fpass=[0 120]; % band of frequencies to be kept
mtm.tapers=[3 5]; % taper parameters
mtm.movingwin=[0.15 0.01];

%% spike-lfp Granger causality
% spike 1 - LFP 2
[Ss1,Sp2,Cs1p2,t,f,GCs1p2,GCp2s1]= spklfp_granger(spike1,lfp2,mtm); 
figure;
plot(t, sum(GCs1p2,2))
hold on
plot(t, sum(GCp2s1,2),'r')


