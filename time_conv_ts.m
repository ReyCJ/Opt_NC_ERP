% By Reza Mahini may 2017

function [v,w]=time_conv_ts(samples,start_ms,end_ms,tempStart,tempEnd)

step=samples/(end_ms-start_ms);

v=(tempStart-start_ms)*step; % time sample

w=(tempEnd-start_ms)*step; % time sample

end