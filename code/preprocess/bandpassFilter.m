function bp = bandpassFilter(x, srate, bands)

if ~exist('srate', 'var')  || isempty(srate),  srate = 1000; end
if ~exist('bands', 'var')  || isempty(bands),  bands = {[60 200], 20}; end

if isa(bands, 'cell')
    % Entire range for broadband
    band_rg  = bands{1};
    
    % Bin width
    band_w   = bands{2};
    
    % All bins
    lb       = band_rg(1):band_w:band_rg(2)-band_w;
    ub       = lb+band_w;
    bands   = [lb; ub]';
end

%% band pass filter each sub-band
bp  = zeros([size(x) size(bands,1)]);

for ii = 1:size(bands,1)
    bp(:,:, ii) = butterpass_eeglabdata(x,bands(ii,:),srate);
end


end