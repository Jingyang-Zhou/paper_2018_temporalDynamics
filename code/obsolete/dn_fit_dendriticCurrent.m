function diff = dn_fit_dendriticCurrent(totalSpikes, data, param, dt, srate)


visualize = 0;

%% prep

tau   = param(1);
alpha = param(2);

bands = {[100, 170], 20}; 

normMax_range = @(x, range) (x - mean(x(range)))./max(x - mean(x(range)));

%% compute model prediction

% dendritic current I: is of dimension [number of trials x time course]
I = dn_computeDendriticCurrent(totalSpikes', tau, alpha, dt);

nTrials = size(I, 1);

bb   = [];
m_bb = [];


bb = extractBroadband(I', srate, 4, bands);


m_bb = median(bb, 2);

%% compute the difference between data and model fit

range = 1 : 300;

data = normMax_range(data, range);
m_bb = normMax_range(m_bb, range);

diff = sum((data(1 : 780) - m_bb(1 : 780)').^2);

%% visualization

if visualize == 1
    figure (100), clf
    plot(data(1 : 780), 'b-'), hold on
    plot(m_bb(1 : 780), 'r-'), drawnow
end

end