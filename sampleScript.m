% Sample script for Asli

%% This is to generate a fake set of tuning curves for a population of cells
pop = makeFakePop(numCells, 'type', type, 'bins', bins,...
    'Pcen', Pcen, 'Pspread', Pspread,...
    'Rcen', rand(1,numCells)*30, 'Rspread',6);
% numCells: Number of cells is number of cells
% the other parameters are best explained in fakeCell.m which is an object,
% briefly it is to specific the type of model, which can be a position ('P'),
% V speed ('V'), or R speed ('R') only, different combinations 'PV',
% 'VR','PR', the full version 'PVR' and finally mixture ('M') of V and R,
% like in the paper (resp = aV + bR);
% The 'Xcen' is for the distribution of centres of the preferenses, and
% 'Xspread' the the distribution, distributed as a gaussian.

%% This is to generate a response from this population
for icell = 1:numCells
    [firingRate(:,icell), spikeTrain(:,icell)] ...
        = FRfromTraj(pop(icell).response, ...
        inputs + noise_common + noise_independent...
        , bin_sizes, meanRate(icell), []);
    % inputs are the input variables for the cell type, same as the
    % dimension of the cells, if you chose 'P' it would be a 1D array, 'PV'
    % 2D array, 'PVR' a 3D array. Look at lines 23-53 in genFakePopSpikes.m
    
    % noise_common and noise_independent are just if you want to add
    % something extra, this was to see if there were any weird errors that
    % I was checking them, again see genFakePopSpikes.m for more on that.
    
    % bin_sizes: bin_sizes = size(pop(1).response);
    % bin_sizes(bin_sizes==1) = []; <-- I don't remember why I had this
    % line.
    
    % for spikeTrain, the poission spike generator creates a spikeTrain
    % so that the mean rate is approzimately meanRate. You don't have to
    % use this function.
end

%% Plus
% You should add an extra kernel at the end, either to firingRate or
% spikeTrain