clear all;
close all;


filenameSpikes = "Z:\MEF29\Dati_neurali\DATI_PROCESSATI\DAT\DATI\FRA_TouchScreen\2026-05-05_13-29-17_MEF29_ID2020\Record Node 101\experiment1\recording1\continuous\Acquisition_Board-100.acquisition_board\continuous\DAT\chs128_192_ave_bpass300_6.0k_Int16_unsorted.mat";
filenameEvents = "Z:\MEF29\Dati_neurali\RAW\DATI\FRA_TouchScreen\CSV\Other\ID2020_data.csv";

eventType = "per_frame_events";

stateSeqs = [0 100];%; 0 99];
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

