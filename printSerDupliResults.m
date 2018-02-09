function [] = printSerDupliResults(thisfileName, Zs, Ss, elTimes, Z)


methods = 

[n,N] = size(Z)
if n >= N
    Z = Z';
end




fileID = fopen(thisfileName,'w');