function loopReporter(cInd, totalInd, repFraction)
% function to report progress in a foor loop.
% cInd is the current loop index, totalInd is the total number of loops,
% and repFraction is the fraction of loops in which feedback should be
% give. e.g. repFraction of 10 will give feedback in 1/10 of loops.
%
% usage: loopReporter(cInd, totalInd, repFraction)

if rem(cInd, round(totalInd / repFraction)) == 0
    fprintf('%i/%i complete\n', cInd, totalInd);
end
    