function data = dn_ecogNotch(data, srate, linenoise)

% This script is modified based on Dora's ecog_notch code

% data = ecog_notch(data,srate)
% notch filter data around 60Hz, 120Hz and 180Hz
% data is [time X electrodes]


%     Copyright (C) 2014  D Hermes
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% line noise and harmonics
noise = [linenoise, linenoise * 2, linenoise * 3];


% 5th order butterworth notch filter, [low/(samplerate/2) high/(samplerate/2)]
[n1_b, n1_a]=butter(5,2*[noise(1) - 1,  noise(1) + 1]/srate,'stop'); % low
[n2_b, n2_a]=butter(5,2*[noise(2) - 1,  noise(2) + 1]/srate,'stop'); % mid
[n3_b, n3_a]=butter(5,2*[noise(3) - 1,  noise(3) + 1]/srate,'stop'); % high
%disp(sprintf('notching out %d %d %d', noise(1), noise(2), noise(3)))
for elec=1:size(data,2)
    data(:,elec)=filtfilt(n1_b,n1_a,data(:,elec)); % low
    data(:,elec)=filtfilt(n2_b,n2_a,data(:,elec)); % mid
    data(:,elec)=filtfilt(n3_b,n3_a,data(:,elec)); % high
end