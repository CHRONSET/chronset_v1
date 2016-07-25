function [pc] = read_phono_coding_info(p2df,fn)
if nargin == 0
    p2df = '/bcbl/home/home_a-f/froux/chronset/data/SayWhen/';
    fn = 'file_and_phon.txt';
end;
%%
fid = fopen([p2df,fn],'r');
data = textscan(fid,'%s\t');
data = data{:};

data2 = cell(length(data)/2,2);
k = 0;
for it = 1:2:length(data)
    
    k= k+1;
    
    x = {data{it} data{it+1}};
    
    data2(k,:) = x;
    
end;

[pc] = data2;


return;