function [lds] = ds_comparison(ds1,ds2)

% algo: length x girth + angle of the tip


lds = ds1;

tip = abs(ds1-ds2);

girth = find(tip <= 55);

lds(girth) = ds2(girth);

