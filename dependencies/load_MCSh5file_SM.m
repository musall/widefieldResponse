function [ephysData, analogData, digitalData, sRate] = load_MCSh5file_SM(targPath,fileName)
    % Script to read electrophysiological measurements from MCS system 
    % input: path and filename of an hdf5 file from MCS DataManager
    % output: ephysData is raw electrode data, analogData are analog inputs, digitalData are digital inputs
    % SM 13122023

    cFile = fullfile(targPath, fileName);   
    data = McsHDF5.McsData(cFile);
    
    for k = 1: length(data.Recording{1, 1}.AnalogStream) %find analog stream with electrode data
        if contains(data.Recording{1, 1}.AnalogStream{1, k}.Label,'Electrode Raw')
            rawStreamID = k;
        elseif contains(data.Recording{1, 1}.AnalogStream{1, k}.Label,'Digital Data')
            digitalStreamID = k;
        end
    end
    
    % get ephysData
    clear ephysData
    ephysData.streamObject = data.Recording{1}.AnalogStream{1,rawStreamID};
    ephysData.samples = ephysData.streamObject.ChannelData;    
    sRate = ephysData.streamObject.getSamplingRate; % sample frequency in Hz

    % get analogData - dont have a good example for this yet.. TBD
    analogData = [];
    
    % get ttlData
    clear digitalData
    digitalData.digitalObject = data.Recording{1}.AnalogStream{1,digitalStreamID};

    digitalLines = 16; %number of expected digital lines.
    bitRepresentation = fastDec2bin(single(digitalData.digitalObject.ChannelData), digitalLines);
    
    Cnt = 0;
    digitalData.events = cell(1, digitalLines);
    for b = digitalLines:-1:1
        Cnt = Cnt + 1;
        digitalData.events{Cnt} = digitalToTimestamp(bitRepresentation(:,b), sRate);
    end
end

