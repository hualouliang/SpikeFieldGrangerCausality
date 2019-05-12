function data=extractspk(data,t)
% Extract segements of spike times between t(1) and t(2)
% Input:
% data: structural array of spike times for each channel/trial 
% t   : time as a 2d vector [t(1) t(2)]
% Output:
% data: spike times between t(1) and t(2)
% This code is adapted from point process extraction code in Chronux

if nargin < 2; error('Need data and times'); end;
if t(1) < 0 || t(2)<=t(1);
    error('times cannot be negative and t(2) has to greater than t(1)');
end;

if isstruct(data); 
    C=length(data);
elseif min(size(data))~=1; 
    error('Can only accept single vector data or a struct array'); 
else
    C=1;
    data=change_row_to_column(data);
end;
%fnames=fieldnames(data);
d2(1:C)=struct('times',[]);
for c=1:C,
    if isstruct(data)
       fnames=fieldnames(data);
       eval(['dtmp=data(c).' fnames{1} ';'])
    else
       dtmp=data(:);
    end
%     eval(['dtmp=data(c).' fnames{1} ';' ])
    sp=dtmp(dtmp>=t(1) & dtmp<t(2));
    d2(c).times=sp;
end
data=d2;