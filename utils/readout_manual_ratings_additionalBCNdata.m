%%
p2df = '/bcbl/home/home_a-f/froux/chronset/data/BCN/additional_data/';
%%
fn = [p2df,'additional_data_rater_coding.txt'];
fid = fopen(fn);
dat = textscan(fid,'%s');
dat = dat{:};
%%
dat = reshape(dat,[10 length(dat)/10 ]);
dat = dat';
%%
p2df = '/bcbl/home/home_a-f/froux/chronset/data/BCN/additional_data/raw_data/';

fID = cell(length(dat)-1,1);
man_RT = zeros(length(dat)-1,1);

k = 0;
for it = 2:length(dat)
    
    
    fn = [dat{it,4},'_*_',dat{it,3},'_*',dat{it,8},dat{it,2},'*.WAV'];
    
    chck = dir([p2df,fn]);
    if ~isempty(chck)
        k = k+1;
        fID{k} = chck.name;
        man_RT(k) = str2double(dat{it,7});
    end;
    
end;
man_RT(k+1:end) = [];
fID(k+1:end) = [];

return;