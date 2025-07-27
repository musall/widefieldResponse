function binOut = fastDec2bin(decOut,n)
% usage: binOut = fastDec2bin(decOut,n)
% faster version of dec2bin but with binary output instead of string.
% decOut is a number of vector of values that should be translated into a
% binary secuence of 'n' bits.
% binOut is the resulting output matrix of size: length(decOut) x n 

if isvector(decOut)
    decOut = reshape(decOut, [], 1);
else
    error('Input has to be a number or 1-D vector');
end
    
[~,e] = log2(max(decOut)); % How many digits do we need to represent the numbers?
binOut = (rem(floor(decOut.*pow2(1-max(n,e):0)),2));