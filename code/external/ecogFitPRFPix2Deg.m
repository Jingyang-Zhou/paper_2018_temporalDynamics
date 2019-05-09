function [x, y, s] = ecogFitPRFPix2Deg(subj, x0, y0, s0)
%Convert pixels to degrees. 
%
% [x, y, s] = ecogFitPRFPix2Deg(subj, x0, y0, s0)
%
% We compute PRF models in pixels, but would like to report them in
% degrees.
%
% Inputs:   
%   subj: subject number (1-4)
%   x0:   x center of pRF in pixels
%   y0:   y center of pRF in pixels
%   s0:   pRF size in pixels (pRF size = sigma/sqrt(n))
%
% Outputs
%   x:  x center of pRF in degrees of visual angle (0 is fixation)
%   y:  y center of pRF in degrees of visual angle (0 is fixation)
%   s:  pRF size in degrees

% S1=RB=9:  61 cm
% S2=Dl=13: 53-55 cm
% S3=??=16: 55-60 cm
% S4=JT=17: 47 cm
% ecog23-kb 30 cm
% Make sure all input arguments are row vectors

nPoints = max([length(x0) length(y0) length(s0)]);
if isempty(x0), x0 = zeros(1,nPoints); end
if isempty(y0), y0 = zeros(1,nPoints); end
if isempty(s0), s0 = zeros(1,nPoints); end

if length(subj) == 1, subj = repmat(subj, nPoints, 1); end
x0 = x0(:)'; y0 = y0(:)'; s0 = s0(:)'; subj = subj(:)';

viewing_distance = zeros(size(subj));
screen_height    = zeros(size(subj));

%  viewing distance, screen height in cm (different for each subject)
for ii = 1:length(subj)
    switch subj(ii)
        case 9,  viewing_distance(ii) = 61;  screen_height(ii) = 20.7; % cm
        case 13, viewing_distance(ii) = 54;  screen_height(ii) = 20.7; % cm
        case 16, viewing_distance(ii) = 57.5;screen_height(ii) = 20.7; % cm
        case 17, viewing_distance(ii) = 47;  screen_height(ii) = 20.7; % cm
        case 19, viewing_distance(ii) = 50;  screen_height(ii) = 17.9; % cm
        case 23, viewing_distance(ii) = 30;  screen_height(ii) = 17.9; % cm% cm
        case 24, viewing_distance(ii) = 45;  screen_height(ii) = 17.9; % cm% cm
        otherwise
            error('Unknown viewing distance for subjec %d', subj);
    end
end

% height from center of screen to top of 15 inch apple macbook (gunjou or
% blue moon) used for ecog experiments
radius_in_cm  = screen_height/2; % radius in cm
radius_in_deg = rad2deg(atan(radius_in_cm./viewing_distance)); % radius in degrees

numpix  = 101; % number of pixels in stimulus description (one side) as input to fitprf
pix2deg = radius_in_deg*2/numpix;
xC      = (numpix+1)/2; % center location in pixels
yC      = xC;            
if subj(1) == 23, xC = 0; end

% convert pixels to degrees
x = pix2deg.*(x0-xC);
y = -pix2deg.*(y0-yC);
s = pix2deg.*s0;

end