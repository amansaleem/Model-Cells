function [es pop] = genFakePopSpikes(es, numCells, type ...
    ,flag_common_noise, flag_independent_noise, meanRate,Pspread)

if nargin<5
    flag_common_noise = 0;
    flag_independent_noise = 0;
end

if nargin<2 | isempty(numCells)
    numCells = 20;
end
if nargin<3
    type = 'P';
end
if nargin<6
    meanRate = 5*ones(1,numCells); %exprnd(2,1,numCells);
    meanRate(meanRate>15) = 15;
    meanRate(meanRate<1)  = 1;
end
if nargin<7
    Pspread = 5;
end
switch type
    case 'P'
%         t = ~isnan(es.traj);
        inputs = es.traj; %(t);
    case 'V'
        t = ~isnan(es.trajspeed);
        inputs = es.trajspeed(t);
    case 'R'
        t = ~isnan(es.ballspeed);
        inputs = es.ballspeed(t);
    case 'PV'
        t = ~isnan(es.traj) & ~isnan(es.trajspeed);
        sV = smthInTime(es.trajspeed(t),60,150);
        inputs = [es.traj sV];
    case 'PR'
        t = ~isnan(es.traj) & ~isnan(es.ballspeed);
        sR = smthInTime(es.ballspeed(t),60,150);
        inputs = [es.traj(t) sR];
    case 'VR'
        t = ~isnan(es.trajspeed) & ~isnan(es.trajspeed);
        sV = smthInTime(es.trajspeed(t),60,150);
        sR = smthInTime(es.ballspeed(t),60,150);
        inputs = [sV sR];
    case 'PVR'
        t = ~isnan(es.traj) & ~isnan(es.trajspeed) & ~isnan(es.ballspeed);
        sV = smthInTime(es.trajspeed(t),60,150);
        sR = smthInTime(es.ballspeed(t),60,150);
        inputs = [es.traj(t) sV sR];
    otherwise
        error(message('Unrecognised type'))
end
bins = [50 15 15];
Pcen = (0:(bins(1)./(numCells-1)):bins(1))';

pop = makeFakePop(numCells, 'type', type, 'bins', bins,...
    'Pcen', Pcen, 'Pspread', Pspread,...
    'Rcen', rand(1,numCells)*30, 'Rspread',6);

bin_sizes = size(pop(1).response);
bin_sizes(bin_sizes==1) = [];
es.spikeTrain = zeros(length(inputs),numCells);
sigma_factor = zeros(size(es.traj));

sigma_low = 50;
sigma_med = 10;
sigma_high= 0;

sigma_factor(es.contrast==0.18) = sigma_low;
sigma_factor(es.contrast==0.6 ) = sigma_med;
sigma_factor(es.contrast>0.6) = sigma_high;

smthWin = 0;

if flag_common_noise
    noise_common = sigma_factor.*smthInTime(randn(size(es.traj)),60,smthWin);
else
    noise_common = 0.*randn(size(es.traj));
end

for icell = 1:numCells
    if flag_independent_noise
        noise_independent = sigma_factor.*smthInTime(randn(size(es.traj)),60,smthWin);
    else
        noise_independent = 0.*randn(size(es.traj));
    end
    [es.firingRate(:,icell), es.spikeTrain(:,icell)] ...
        = FRfromTraj(pop(icell).response, ...
        inputs + noise_common + noise_independent...
        , bin_sizes, meanRate(icell), []);
    pop(icell).meanRate = meanRate(icell);
end