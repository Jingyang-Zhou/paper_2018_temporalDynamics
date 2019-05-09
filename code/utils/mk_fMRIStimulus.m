function [stimulus, time] = mk_fMRIStimulus(with0OrNot)

% function [stimulus, time] = importStimulus(with0OrNot)

% INPUT -------------------------------------------
% with0OrNot : 'with0' or 'without0'


% OUTPUT -------------------------------------------

% stimulus : 12 x: 4500
% time     : 1: 4500, unit: millisecond

% History :
% 03/07/2016, add zero condition to stimulus as the 13th condition

%% Example

% exampleOn = 0;
% exampleOn = checkExampleOn(exampleOn, mfilename);

% if exampleOn
%     with0OrNot = 'without0';
% end

%% Assess inputs

possibleInputs = {'with0', 'without0'};
assert(ismember(lower(with0OrNot), possibleInputs), 'Un-identifiable input.');

%% Pre-defined variables

trialTime = 4500; % in unit of millisecond
nStimulus = 12;
frameRate = 60; % Hz
tPattern  = [1, 2, 4, 8, 16, 32];
whichTisi = 4;

%% Secondary variables

frameTime = 1000 / frameRate;
dur       = round(tPattern * frameTime);

%% compute stimulus

switch with0OrNot
    case 'with0'
        stimulus  = zeros(nStimulus+1, trialTime);
    otherwise
        stimulus  = zeros(nStimulus, trialTime);
end

for k = 1: 12
    if k < 7
        stimulus(k, 1 : dur(k)) = 1;
    else
        stimulus(k, 1 : dur(whichTisi)) = 1;
        tStart = dur(whichTisi) + dur(k-6) + 1;
        tEnd   = dur(whichTisi) * 2 + dur(k-6);
        stimulus(k,  tStart : tEnd ) = 1;
    end
end

%% compute time

time = (1 : trialTime)./1000;

%% Plotting results

% if exampleOn
%     figure (1), clf
%     imagesc(stimulus)
%     title('12 stimuli')
%     xlabel('time(ms)')
%     ylabel('nth stimulus')
%     colorbar
% end

end