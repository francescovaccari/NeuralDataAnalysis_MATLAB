clear all;
close all;

dataFolder = 'C:\Users\frada\Desktop\WORK in PROGRESS\LAB\TouchScreen_project\Data\Datasets\';
filenameSpikes = [dataFolder 'ID2055_chs0_128_ave_bpass300_6.0k_Int16_r_4_unsorted.mat'];
filenameEvents = [dataFolder 'ID2055_data.csv'];
filenameCSMS = [dataFolder 'ID2055_CSMS_chs0_128_ave_bpass300_6.0k_Int16_r_4_unsorted.mat'];

eventType = "per_frame_events";

stateSeqs = [0 100];% 0 99];
matchOpt = "first&last";

otherMarkers = [1 4 5 6 7 8 9 10];

% stateSeqs = {[0 1 4 5 6 7 8 9 10 100]};
% matchOpt = "strict";


recordedData = RecordedData.makeCSMS(filenameSpikes, filenameEvents, eventType, stateSeqs, matchOpt, otherMarkers);

CS = recordedData.CS;
MS = recordedData.MS;
TC = recordedData.TC;
TCTS = recordedData.TCTS;
TRIAL_OUTCOME = recordedData.TRIAL_OUTCOME;

save (filenameCSMS,...
    "-v7.3", ...
    "CS", "MS", "TC", "TCTS", "TRIAL_OUTCOME")    ;