function [S1,S2,C,tt,f,GCs2p,GCp2s]= spklfp_granger(dataspk,datalfp,mtm)
% 
% Evaluate spectrum, cross spectrum for a pair of spike trains and LFP from two different channels. 
% Input:
%   dataspk: structural array of spike times for each trial or 1d array of spike times
%   datalfp: continuous data in form of samples x trials
%   mtm: structure with fields tapers, Fs, fpass, movingwin
%           tapers   [TW K] where TW is the time-bandwidth product and K is the number of tapers to be used 
%           Fs    sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%           movingwin  in the form [length of movingwin, size of step] 

% Output:
%   C: coherence between spike - lfp
%   S1: trial averaged spike spectrogram in form of time x frequency x channels
%   S2: trial averaged LFP spectrogram in form of time x frequency x channels
%   GCs2p: Granger causality from spike to LFP
%   GCp2s: Granger causality from LFP to spike             
%   tt:    times
%   f:     frequencies
%   
% Xiajing Gong, Drexel University, 7/10
%
%% %%%%% Spectral denstiy matrix for 2 channels (channels x channels x frequency)
tapers=mtm.tapers;
Fs=mtm.Fs;
fpass=mtm.fpass;
movingwin=mtm.movingwin;

N=size(datalfp,1);Ntrl=size(datalfp,2); % size of data   
Nwin=round(Fs*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*Fs); % number of samples to step through
nfft=max(2^(nextpow2(Nwin)+2),Nwin);

% frequencies where spectrum is evaluated
f=0:Fs/nfft:Fs;
findx=find(f>=fpass(1) & f<=fpass(end));
f=f(findx);
Nf=length(f);

% obtain tapers
[tapers,eigs]=dpss(Nwin,tapers(1),tapers(2)); 
tapers = tapers*sqrt(Fs);

winstart=1:Nstep:N-Nwin+1;
nw=length(winstart);
% initialize 
   C=zeros(nw,Nf);
   S12=zeros(nw,Nf);
   S21=zeros(nw,Nf); 
   S1=zeros(nw,Nf); 
   S2=zeros(nw,Nf);

for n=1:nw
   indx=winstart(n):winstart(n)+Nwin-1;
   t=indx/Fs;
   data2=datalfp(indx,:);
   data1=extractspk(dataspk,[t(1) t(end)]);
   Ntap=size(tapers,2); % size of tapers
%%% fourier transform of discrete data
   H=fft(tapers,nfft,1);  % fft of tapers
   H=H(findx,:); % restrict fft of tapers to required frequencies
   w=2*pi*f; % angular frequencies at which ft is to be evaluated
   Nsp=zeros(1,Ntrl); Msp=zeros(1,Ntrl);
   for ch=1:Ntrl;
  if isstruct(data1);
     fnames=fieldnames(data1);
     eval(['dtmp=data1(ch).' fnames{1} ';'])
     indx=find(dtmp>=min(t)&dtmp<=max(t));
     if ~isempty(indx); dtmp=dtmp(indx);
     end;
  else
     dtmp=data1;
     indx=find(dtmp>=min(t)&dtmp<=max(t));
     if ~isempty(indx); dtmp=dtmp(indx);
     end;
  end;
  Nsp(ch)=length(dtmp);
  Msp(ch)=Nsp(ch)/length(t);
  if Msp(ch)~=0;
      data_proj=interp1(t',tapers,dtmp);
      exponential=exp(-i*w'*(dtmp-t(1))');
      J1(:,:,ch)=exponential*data_proj-H*Msp(ch);
  else
      J1(1:Nf,1:Ntap,ch)=0;
  end;
  end; 
%%% fourier transform of continuous data
   tapers2=tapers(:,:,ones(1,Ntrl)); % add channel indices to tapers
   data2=data2(:,:,ones(1,Ntap)); 
   data2=permute(data2,[1 3 2]); % reshape data to get dimensions to match those of tapers
   data2_tap=data2.*tapers2; % product of data with tapers
   J2=fft(data2_tap,nfft)/Fs; 
   J2=J2(findx,:,:); 
%%% evaluate spectral matrix of spike and LFP 
    s12=squeeze(mean(J1.*conj(J2),2)); % cross spectrum
    s21=squeeze(mean(J2.*conj(J1),2));
    s1=squeeze(mean(conj(J1).*J1,2)); % spectrum data 1
    s2=squeeze(mean(conj(J2).*J2,2)); % spectrum data 2
    s12=squeeze(mean(s12,2)); s21=squeeze(mean(s21,2));
    s1=squeeze(mean(s1,2)); s2=squeeze(mean(s2,2)); c=abs(s12./sqrt(s1.*s2));
    C(n,:)=c;S12(n,:)=s12;S21(n,:)=s21;S1(n,:)=s1;S2(n,:)=s2;  
end

%% %%%%% spectral matrix factorization to obtain transfer function and error covariance matrix using Wilsons algorithm: 
for kk=1:nw
S(1,2,:)=S12(kk,:);S(2,1,:)=S21(kk,:);
S(1,1,:)=S1(kk,:);S(2,2,:)=S2(kk,:); % spectral matrix for GC calculation
f_indx=0;
for jj=f
    f_indx=f_indx+1;
    Sarr(:,:,f_indx)=S(:,:,f_indx);
    if(f_indx>1)
%         Sarr(:,:,2*(Nf-1)+2-f_indx)=S(:,:,f_indx).';%'
      Sarr(:,:,2*(Nf-1)+2-f_indx)=S(:,:,f_indx).';%'
    end
end
% perform ifft to obtain gammas
for k1=1:2
    for k2=1:2
        gam(k1,k2,:)=real(ifft(squeeze(Sarr(k1,k2,:)))*Fs);
    end
end
gam0 = gam(:,:,1); h = chol(gam0); 

for ind =1:size(Sarr,3),
       psi(:,:,ind) = h; % initialization of  for the 1st iteration
end
for iter = 1:100
     for ind = 1:size(Sarr,3),
       g(:,:,ind)=inv(psi(:,:,ind))*Sarr(:,:,ind)*inv(psi(:,:,ind))'+eye(2);%'
     end
     gplus = PlusOperator(g,2,Fs,f); % Eq 3.1
   psiold=psi;
   for k = 1:size(Sarr,3),
          psi(:,:,k) = psi(:,:,k)*gplus(:,:,k);
          psierr(k)=norm(psi(:,:,k)-psiold(:,:,k),1);
   end
   iter;
   psierrf=mean(psierr);
   if(psierrf<1E-12),break;end;     
end 

for k = 1:Nf
      Snew(:,:,k) = squeeze(psi(:,:,k)*psi(:,:,k)'); %% Snew:improved spectral density matrix 
end

for k1=1:2
    for k2=1:2
        gamtmp(k1,k2,:)=real(ifft(squeeze(psi(k1,k2,:))));
    end
end
A0=squeeze(gamtmp(:,:,1)); % this is psi_1
A0inv=inv(A0);
Znew=A0*A0.'*Fs; % noise covariance(channels x channels x frequency),

for k = 1:Nf
      Hnew(:,:,k) = squeeze(psi(:,:,k))*A0inv; % transfer function
end

%%%%% estimate Granger causality from spectral density matrix, transfer function and noise covariance matrix.
clear f_indx jj;
f_indx = 0;
         for mm = f
             f_indx = f_indx + 1;
             Zp2s = Znew(2,2) - Znew(1,2)^2/Znew(1,1); %corrected noise covariance
             Zs2p = Znew(1,1) - Znew(2,1)^2/Znew(2,2);
             GCp2s(kk,f_indx) = log(abs(Snew(1,1,f_indx))/abs(Snew(1,1,f_indx)-(Hnew(1,2,f_indx)*Zp2s*conj(Hnew(1,2,f_indx)))/Fs)); %Geweke's original measure
             GCs2p(kk,f_indx) = log(abs(Snew(2,2,f_indx))/abs(Snew(2,2,f_indx)-(Hnew(2,1,f_indx)*Zs2p*conj(Hnew(2,1,f_indx)))/Fs));
%              GCp2s(kk,f_indx) = log(abs(Snew(1,1,f_indx))/abs(Snew(1,1,f_indx)-(abs(Hnew(1,2,f_indx))*Zp2s)/Fs)); %Geweke's original measure
%              GCs2p(kk,f_indx) = log(abs(Snew(2,2,f_indx))/abs(Snew(2,2,f_indx)-(abs(Hnew(2,1,f_indx))*Zs2p)/Fs));

         end
end

tt=(winstart+round(Nwin/2))/Fs;