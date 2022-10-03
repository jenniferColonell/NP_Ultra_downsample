function callDownsampling

% for a list of input binary files, and a list of downsample patterns
% call poolChan with each pattern and a 'noiseModel', which is just the 
% channel map and fourier transform of noise data for a UHD probe.

% combining data across channels is a point by point calculation, so each
% block of data is analyzed independently (no overlap).

% Sites are combined by averaging the signal across sites. This reduces the
% amplifier and thermal noise by 1/sqrt(N sites). To simulate data with full
% amplifier noise, we add noise with the same frequency spectrum as a self 
% referencedd noise test in PBS. This data is stored in a 'noise model'
% structure generated by 'make_noise_model.m'.

% full path to noise model
nmPath = 'C:\Users\colonellj\Documents\matlab_ephys\NP_Ultra_downsample\UHD_extref_noiseModel.mat';
noiseModel = load(nmPath);       % loads structure nm
addNoise = 1;                    % set to 1 to add noise
noiseFrac = 0.8;    % UHD noise in saline ~ 6 uV; NP 2.0 ~ 8-8.5 uV; spec for NP1.0 = 5 uV.
%JIC note -- noiseFrac = 0.8 for the 'sortability/detectablity' vs downsampling 

% list of pattern types
patTypeArr = [2,3,4];
% need a list of start sites for each pattern type. For 2.0-like patterns,
% the pattern can be offset by one site in horizontal and vertical
% directions. (Note that the difference in signal is very small!)
stSiteList{1} = [0];
stSiteList{2} = [0];
stSiteList{3} = [0];


% SGLX data -- if true, will extract basename to create output and chanMap names
sglxData = 1;

% file locations for input binaries
% alm_list
% fileList{1} = 'G:\SC041_OUT\catgt_SC041_081120_g0\SC041_081120_g0_imec0\SC041_081120_g0_tcat.imec0.ap.bin';
% fileList{2} = 'G:\SC041_OUT\catgt_SC041_081120_g1\SC041_081120_g1_imec0\SC041_081120_g1_tcat.imec0.ap.bin';
% fileList{3} = 'G:\SC041_OUT\catgt_SC041_081120_g1\SC041_081120_g1_imec1\SC041_081120_g1_tcat.imec1.ap.bin';
% fileList{4} = 'G:\SC041_OUT\catgt_SC041_081120_g2\SC041_081120_g2_imec0\SC041_081120_g2_tcat.imec0.ap.bin';
% fileList{5} = 'G:\SC041_OUT\catgt_SC041_081120_g2\SC041_081120_g2_imec1\SC041_081120_g2_tcat.imec1.ap.bin';
% fileList{6} = 'G:\SC041_OUT\catgt_SC041_081120_g3\SC041_081120_g3_imec0\SC041_081120_g3_tcat.imec0.ap.bin';
% fileList{7} = 'G:\SC041_OUT\catgt_SC041_081220_g0\SC041_081220_g0_imec0\SC041_081220_g0_tcat.imec0.ap.bin';
% fileList{8} = 'G:\SC041_OUT\catgt_SC041_081220_g0\SC041_081220_g0_imec1\SC041_081220_g0_tcat.imec1.ap.bin';
% fileList{9} = 'G:\SC041_OUT\catgt_SC041_081220_g1\SC041_081220_g1_imec0\SC041_081220_g1_tcat.imec0.ap.bin';
% fileList{10} = 'G:\SC041_OUT\catgt_SC041_081220_g1\SC041_081220_g1_imec1\SC041_081220_g1_tcat.imec1.ap.bin';
% fileList{11} = 'G:\SC041_OUT\catgt_SC041_081220_g2\SC041_081220_g2_imec0\SC041_081220_g2_tcat.imec0.ap.bin';
% fileList{12} = 'G:\SC041_OUT\catgt_SC041_081220_g3\SC041_081220_g3_imec0\SC041_081220_g3_tcat.imec0.ap.bin';

% cerebellum list
fileList{1} = 'G:\SC041_OUT\catgt_SC041_081120_g0\SC041_081120_g0_imec2\SC041_081120_g0_tcat.imec2.ap.bin';
fileList{2} = 'G:\SC041_OUT\catgt_SC041_081120_g1\SC041_081120_g1_imec3\SC041_081120_g1_tcat.imec3.ap.bin';
fileList{3} = 'G:\SC041_OUT\catgt_SC041_081120_g3\SC041_081120_g3_imec3\SC041_081120_g3_tcat.imec3.ap.bin';
fileList{4} = 'G:\SC041_OUT\catgt_SC041_081220_g0\SC041_081220_g0_imec2\SC041_081220_g0_tcat.imec2.ap.bin';
fileList{5} = 'G:\SC041_OUT\catgt_SC041_081220_g0\SC041_081220_g0_imec3\SC041_081220_g0_tcat.imec3.ap.bin';
fileList{6} = 'G:\SC041_OUT\catgt_SC041_081220_g1\SC041_081220_g1_imec2\SC041_081220_g1_tcat.imec2.ap.bin';
fileList{7} = 'G:\SC041_OUT\catgt_SC041_081220_g2\SC041_081220_g2_imec2\SC041_081220_g2_tcat.imec2.ap.bin';
fileList{8} = 'G:\SC041_OUT\catgt_SC041_081220_g2\SC041_081220_g2_imec3\SC041_081220_g2_tcat.imec3.ap.bin';

% base for output; will go into director with base + probe name
outBase = 'D:\UHD_ds\SC041_curated';

nFile = numel(fileList);
NT = 100*30000;  % number of time points per batch

for i = 1:nFile 

    currFile = fileList{i};
    [currPath,currName,~] = fileparts(currFile);
    [~,prbFld] = fileparts(currPath);  % assumes standard catgt with -out_prb_fld
    outPath = fullfile(outBase,prbFld);

    for k = 1:numel(patTypeArr)
        
    newPatType = patTypeArr(k);
    stSiteArr = stSiteList{k};

    for j = 1:numel(stSiteArr)
        stSite = stSiteArr(j);
        if sglxData
            baseName = extractBefore(currName,'_g');
            gt_decor = sprintf('_g%s', extractAfter(currName,'_g'));
            outBinName = sprintf( '%s_pat%d_%d%s.bin', baseName,newPatType,stSite,gt_decor );
            outChanMapName = sprintf( '%s_pat%d_%d_chanMap.mat', baseName,newPatType,stSite );
            inputMetaName = sprintf('%s.meta', currName);
            outMetaName = sprintf( '%s_pat%d_%d%s.meta', baseName,newPatType,stSite,gt_decor );
    
            [meta] = ReadMeta(inputMetaName, currPath);
            rm_order = parseShankMap(meta);  % checking that channels are in row major order
            [AP,LF,SY] = ChannelCountsIM(meta);
            orig_NchanTOT = AP + SY;
    
        else
            outBinName = sprintf('pat%d_%d_out%d.bin', newPatType,stSite,i);
            outChanMapName = sprintf( '%pat%d_%d_out%d_chanMap.mat', newPatType,stSite,i );
            % assume default UHD configuration in binary
            orig_NchanTOT = 385;
            rm_order = 1:384;
        end
    
        need_reorder = issorted(rm_order);
    
        bytes       = get_file_size(fileList{i}); % size in bytes of raw binary
        nTimepoints = floor(bytes/orig_NchanTOT/2); % number of total timepoints
        Nbatch      = ceil(nTimepoints /NT); % number of data batches
    
        fid         = fopen(fileList{i}, 'r'); % open for reading raw data
        fidW        = fopen(fullfile(outPath,outBinName), 'w'); % open for writing processed data
    
        for ibatch = 1:Nbatch-1
            offset = 2*NT*(ibatch -1)*orig_NchanTOT;
            fseek(fid, offset, 'bof'); % fseek to batch start in raw file
            buff = fread(fid, [orig_NchanTOT NT], '*int16'); % read and reshape.
            if need_reorder
                buff = buff(rm_order,:);
            end
            if newPatType > 0
                [buff, chanMap, xCoord, yCoord] = poolChan(buff, addNoise, noiseModel, noiseFrac, newPatType, stSite);
            end
            fwrite(fidW, buff, 'int16'); % write this batch to binary file
        end
    
        % last batch
        offset = 2*NT*(Nbatch-1)*orig_NchanTOT;
        fseek(fid, offset, 'bof'); % fseek to batch start in raw file
        NTlast = nTimepoints - NT*(Nbatch-1);
        buff = fread(fid, [orig_NchanTOT NTlast], '*int16'); % read and reshape.
        if newPatType > 0
            [buff, chanMap, xCoord, yCoord, saveChanStr, saveChanArr] = poolChan(buff, addNoise, noiseModel, noiseFrac, newPatType, stSite);
        end
        fwrite(fidW, buff, 'int16'); % write this batch to binary file
    
        fclose(fid);
        fclose(fidW);
        if newPatType > 0
            % write a ks2 channel map for the new file
            NchanTOT = numel(chanMap);
            chanMap = 1:NchanTOT;
            chanMap0ind = chanMap -1;
            connected = ones([NchanTOT,1]);
            kcoords = ones([NchanTOT,1]);
            save(fullfile(outPath,outChanMapName), 'NchanTOT','chanMap','chanMap0ind', 'xCoord','yCoord','connected','kcoords')
        else
            NchanTOT = orig_NchanTOT;   
        end
    
        if sglxData
            % make a metadata file
            newShankMap = uhdShankMap(saveChanArr);
            fp = dir(fullfile(outPath,outBinName));        
            newTag = cell(5,1);
            newTag{1} = sprintf('%s%d', 'fileSizeBytes=', fp.bytes);
            newTag{2} = sprintf('%s%d', 'nSavedChans=', NchanTOT);
            newTag{3} = sprintf('%s%d%s', 'snsApLfSy=', NchanTOT, ',0,0');
            newTag{4} = sprintf('%s', saveChanStr);
            newTag{5} = sprintf('%s', newShankMap);
    
            repTags = cell(5,1);
            repTags{1} = 'fileSizeBytes';
            repTags{2} = 'nSavedChans';
            repTags{3} = 'snsApLfSy';
            repTags{4} = 'snsSaveChanSubset';
            repTags{5} = '~snsShankMap';
            
            fmodel = fopen( fullfile(currPath,inputMetaName), 'r');
            fmeta = fopen( fullfile(outPath,outMetaName), 'w');
            
            tline = fgetl(fmodel);
            while ischar(tline)
                currTag = extractBefore(tline,'=');
                tagFound = find(strcmp(repTags, currTag));
                if isempty(tagFound)
                    %copy over this line as is
                    fprintf(fmeta, '%s\n', tline );
                else
                    fprintf('found: %s\n', repTags{tagFound} );
                    fprintf(fmeta, '%s\n', newTag{tagFound} );
                end  
                tline = fgetl(fmodel);
            end
            fclose(fmeta);
            fclose(fmodel);
    
        end
    end  % end of loop over starting sites
    end  % end of loop over pattern type
end

end



function bytes = get_file_size(fname)
% gets file size in bytes, ensuring that symlinks are dereferenced
% MP: not sure who wrote this code, but they were careful to do it right on Linux
    bytes = NaN;
    if isunix
        cmd = sprintf('stat -Lc %%s %s', fname);
        [status, r] = system(cmd);
        if status == 0
            bytes = str2double(r);
        end
    end

    if isnan(bytes)
        o = dir(fname);
        bytes = o.bytes;
    end
end


function [meta] = ReadMeta(metaName, path)

    % Parse ini file into cell entries C{1}{i} = C{2}{i}    
    fid = fopen(fullfile(path, metaName), 'r');
% -------------------------------------------------------------
%    Need 'BufSize' adjustment for MATLAB earlier than 2014
%    C = textscan(fid, '%[^=] = %[^\r\n]', 'BufSize', 32768);
    C = textscan(fid, '%[^=] = %[^\r\n]');
% -------------------------------------------------------------
    fclose(fid);

    % New empty struct
    meta = struct();

    % Convert each cell entry into a struct entry
    for i = 1:length(C{1})
        tag = C{1}{i};
        if tag(1) == '~'
            % remake tag excluding first character
            tag = sprintf('%s', tag(2:end));
        end
        meta = setfield(meta, tag, C{2}{i});
    end
end % ReadMeta

function shankMapStr = uhdShankMap( saveChanArr )
    % saveChanArr is the channel numbers in a 0-based standard passive UHD
    % build a shank map for a passive UHD run saving just these sites
    % This is done from scratch to avoid issues with incorrect shank maps
    % in early UHD data. Note that to look at that data in the viewer, the 
    % metadata needs to be corrected.
    Nchan = numel(saveChanArr);
    shankMapStr = "~snsShankMap=(1,8,48)";
    for i = 1:Nchan
        row = floor(saveChanArr(i)/8);
        col = mod(saveChanArr(i),8);
        currEntry = sprintf('(0:%d:%d:1)', col, row);
        shankMapStr = sprintf('%s%s', shankMapStr, currEntry);
    end

end

function rm_order = parseShankMap( meta )
% parse shank map string to return array to map channels in row major order
% important for UHD2 downsampling; also downsampling simulated data.

    C = textscan(meta.snsShankMap, '(%d:%d:%d:%d', ...
                'EndOfLine', ')', 'HeaderLines', 1 );
    nElec = numel(C{1});
    rc_array = zeros([nElec,2]);
    rc_array(:,1) = double(cell2mat(C(3)));
    rc_array(:,2) = double(cell2mat(C(2)));
    [~,rm_order]= sortrows(rc_array);

end

function [AP,LF,SY] = ChannelCountsIM(meta)
    M = str2num(meta.snsApLfSy);
    AP = M(1);
    LF = M(2);
    SY = M(3);
end % ChannelCountsIM
