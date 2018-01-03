function spikeTrain = generateSpikes(inRate, sampRate, meanRate, maxIter, prec)
% Generate a random spikeTrain (kind of Poisson) using an input to the
% process (inRate: range: 0-1) and generating a 'spikeTrain' with an
% overall mean firing rate equal to meanRate(dflt=1).
% inRate is (numCells x T). meanRate < sampRate.
%
% Usage: spikeTrain = generateSpikes(inRate, sampRate, meanRate, maxIter, prec)
% maxIter & prec are optional (dflts: 15,0.05) is how close (+-prec) to
% the meanRate we want to get to and how many (maxIter) we can use to do that
%
% Aman Saleem
% July 2013

if nargin < 2
    sampRate = 60;
end
if nargin < 3
    meanRate = 1;
end
if nargin < 4
    maxIter = 50;
end
if nargin < 5
    prec = 0.05;
end

if sampRate < meanRate
    error('Mean rate cannot be greater than sampling rate!');
end

spikeTrain = zeros(size(inRate));
numCells = size(inRate, 1);

if length(meanRate)<numCells
    meanRate(length(meanRate):numCells) = 1;
end

for icell = 1:numCells
    idx = 0;
    mRate = meanRate(icell)/sampRate;
    base = mRate;
    
    spikeTrain(icell,:) = 0;
    T = (1-base).*rand(1,size(spikeTrain,2));
    spikeTrain(icell,inRate(icell,:)>T) = 1;
    
    while (abs(mean(spikeTrain(icell,:)) - mRate)>prec/60) && (idx < maxIter)
        rateDiff = abs(mean(spikeTrain(icell,:)) - mRate);
        if mean(spikeTrain(icell,:)) < mRate
            base = base + rateDiff/2;
        else
            base = base - rateDiff/2;
        end
        spikeTrain(icell,:) = 0;
        T = (1-base).*rand(1,size(spikeTrain,2));
        spikeTrain(icell,inRate(icell,:)>T) = 1;
%         display(['Iteration index: ' num2str(idx) ', rate difference: ' num2str(60*(mean(spikeTrain(icell,:)) - mRate))]);
        idx = idx + 1;
    end
%     display(['Final iteration:' num2str(idx) ', Rate diff: ' num2str(abs(mean(spikeTrain(icell,:)) - mRate))]);
end