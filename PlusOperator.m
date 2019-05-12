function gp = PlusOperator(g,m,Fs,freq)
%This function is for [ ]+operation
for k1=1:m
    for k2=1:m
          gam(k1,k2,:)= ifft(squeeze(g(k1,k2,:)));
    end
end
% take positive lags only and half of the zero lag
gamp = gam;beta0 = 0.5*gam(:,:,1); 
gamp(:,:,1) = triu(beta0);  %this is Stau
gamp(:,:,length(freq)+1:end) = 0;
% reconstitute
for k1=1:m
    for k2=1:m
         gp(k1,k2,:)= fft(squeeze(gamp(k1,k2,:)));
    end
end