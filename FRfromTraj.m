function [firingRate, spikeTrain, nBins] = FRfromTraj(tuning, inputs, grid_size, meanRate, input_lims)

% firingRate = FRfromTraj(tuning, inputs, input_lims)
% tuning is the population tuning rates, possibly from makeFakePop.m
% inputs is the particular variable from es that is run through the tuning
% function and then used to get spikes by running through the spike
% generator (generateSpikes.m).
%
% inputs: NxT where there are N variables and T time bins
% tuning: Nx1 for each cell
% input_lims: are the limits to which the inputs are normalised

inputs = inputs';

if nargin<3 | isempty(grid_size)
    grid_size = 50*ones(1,size(inputs,1));
end

if nargin<4 | isempty(meanRate)
    meanRate = 1;
end

if nargin<5 | isempty(input_lims)
    input_lims(1,:) = min(inputs,[],2);
    input_lims(2,:) = max(inputs,[],2);
end

if sum(size(tuning)>1) ~= size(inputs,1)
    error('Tuning and inputs dimensions do not agree!')
end

% bins = size(tuning);
% bins = bins(1:sum(size(tuning)>1));

c = true(1,size(inputs,2));
for ivar = 1:size(inputs,1)
    c(isnan(inputs(ivar,:))) = false;
end

nor_inputs = zeros(size(inputs));

for ivar = 1:size(inputs,1)
    inputs(ivar,inputs(ivar,:)<input_lims(1,ivar)) = input_lims(1,ivar);
    inputs(ivar,inputs(ivar,:)>input_lims(2,ivar)) = input_lims(2,ivar);
    [nbins{ivar}, nor_inputs(ivar,:)]           = normalise_one_var(inputs(ivar,:), c,grid_size(ivar));
    nor_inputs(ivar,~c) = NaN;
end

firingRate = zeros(1,size(nor_inputs,2));
goodT = find(c);

% Change this line here
switch size(inputs,1)
    case str2num('1')
        firingRate(goodT) = tuning(nor_inputs(1,goodT));
    case str2num('2')
%         t1 = tic;
        for t = 1:length(goodT)
            firingRate(goodT(t)) = tuning(nor_inputs(1,goodT(t)),nor_inputs(2,goodT(t)));
        end
%         toc(t1)
    case str2num('3')
        firingRate(goodT) = tuning(nor_inputs(1,goodT),nor_inputs(2,goodT),nor_inputs(3,goodT));
    case str2num('4')
        firingRate(goodT) = tuning(nor_inputs(1,goodT),nor_inputs(2,goodT),nor_inputs(3,goodT),nor_inputs(4,goodT));
end

maxIter = 500;
prec    = 1;
spikeTrain = generateSpikes(firingRate,60,meanRate,maxIter, prec);

    function [bins, B]           = normalise_one_var(A, c,aGrid)
        maxA = max(A(c));
        minA = min(A(c));
        
        bins = minA:((maxA-minA)/(aGrid -1+eps)):maxA;
        A = (A-minA)/(maxA-minA);
        
        B = 1+floor(aGrid*A/(1+eps));
    end
% function [es, bins]         = normalise_es(es_input, c, pGrid, vrGrid)
%         es = es_input;
%
%         [bins.P_bins es.traj] = normalise_one_var(es.traj, c,pGrid);
%         [bins.V_bins es.trajspeed] = normalise_one_var(es.trajspeed, c,vrGrid);
%         [bins.R_bins es.ballspeed] = normalise_one_var(es.ballspeed, c,vrGrid);
%
%     end

%     function [es_out]           = runThroughPVR(es_in)
%         es_out = es_in;
%         es_out.spikeTrain = zeros(size(es_in.spikeTrain,1),numCells);
%         for icell = 1:numCells
%             c = (es_in.traj>0 & ~isnan(es_in.trajspeed) & ~isnan(es_in.traj) & ~isnan(es_in.ballspeed));
%             rate = zeros(size(es_in.traj));
%             rate(~c) = rand(1,sum(~c));
%             good = find(c);
%             for n = 1:length(good)
%                 rate(good(n)) = PVR(icell).tuning(es_in.traj(good(n)),es_in.trajspeed(good(n)),es_in.ballspeed(good(n)));
%             end
%             rand_gen = [1-rand(size(rate))];
%
%             es_out.spikeTrain(rate>rand_gen,icell) = 1;
% %             es_out.spikeTrain(~c,icell) = ;
%         end
end