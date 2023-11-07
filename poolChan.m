function [datNew, chanMap, xc, yc, saveChanStr, saveChanArr] = poolChan(buff, addNoise, noiseModel, noiseFrac, newPatType, stSite)
% pool channels on an UHD probe to build a new pattern
% The patterns are hard coded here, pick which one with newPatType

% incoming datr is chan x time
[stChan, NT] = size(buff);
if stChan == 385
    nSite = 384;    % this is a real UHD recordign w/ SY channel
else
    nSite = stChan; % simulated data, no SY
end

% keep track of bottom left site 
saveChanStr = 'snsSaveChanSubset=';

switch newPatType
    
    case 1
        % pool and skip sets of 4 channels to make a checkerboard pattern
        % on the LHS of the probe. There will be 48 channels in the 
        % pooled data.
        Nchan = nSite/8;
        Ncomb = 4;
        datNew = zeros([Nchan,NT],'double');  % make double to match datr
        xc = zeros([Nchan,1],'double');
        yc = zeros([Nchan,1],'double');
        chanMap = zeros([Nchan,1],'double');
        chanMap(:,1) = (1:Nchan);
        connected = ones([Nchan,1],'double');
        for i = 0:2:Nchan-1       %step through rows in new pattern
            row = floor(i/2); %starts at zero
            bEven = ~(mod(row,2));
            if bEven 
                % rows 0, 2, 4 ...
                c0 = row*16 + 1;   % add 1 for matlab 
                xc(i+1) = 0;
                xc(i+2) = 24;
            else
                % rows 1, 3, 5 ...
                c0 = row*16 + 2 + 1;   % add 2 to offset from even rows, add 1 for matlab 
                xc(i+1) = 12;
                xc(i+2) = 36;
            end
            
            yc(i+1) = row*12;
            yc(i+2) = row*12;
            datNew(i+1,:) = (buff(c0,:) + buff(c0+1,:) + buff(c0+8,:) + buff(c0+9,:))/4;
            datNew(i+2,:) = (buff(c0+4,:) + buff(c0+5,:) + buff(c0+12,:) + buff(c0+13,:))/4;
            saveChanStr = sprintf('%s%d,%d,%d,%d,',saveChanStr,c0-1,c0+3);
        end
        nCol = 2;  % For spikeGLX shank map, staggered colums treated as two

    case 2
        % pool sets of 4 channels to make a NP 1.0-like pattern on the "RHS" of
        % the probe. There will be 24 channels in the final pattern
        % site 0 at x = 12 matches standard 1.0
        if ~ismember(stSite,[0,2,16,18])
            fprintf('stSite must be 0,2,16 or 18 for the NP 1.0 pattern')
            stSite = 0;
        end
        Nchan = nSite/16;
        Ncomb = 4;
        datNew = zeros([Nchan,NT],'double');  % make double to match datr
        xc = zeros([Nchan,1],'double');
        yc = zeros([Nchan,1],'double');
        chanMap = zeros([Nchan,1],'double');
        chanMap(:,1) = (1:Nchan);
        connected = ones([Nchan,1],'double');
        
        for i = 0:2:Nchan-1       %step through rows in new pattern

            switch stSite
                case 2
                    row = floor(i/2); %starts at zero
                    bEven = ~(mod(row,2));
                    if bEven 
                        % rows 0, 2, 4 ...
                        c0 = row*32 + 2 + 1;   % add 2 to offset to 3rd site in row, add 1 for matlab 
                        xc(i+1) = 12;
                        xc(i+2) = 36;
                    else
                        % rows 1, 3, 5 ...
                        c0 = row*32 + 1;   %  add 1 for matlab 
                        xc(i+1) = 0;
                        xc(i+2) = 24;
                    end
                    yc(i+1) = row*24;
                    yc(i+2) = row*24;
                case 18
                    row = floor(i/2); %starts at zero
                    bEven = ~(mod(row,2));
                    if bEven 
                        % rows 0, 2, 4 ...
                        c0 = row*32 + 18 + 1;   % add 16 + 2 to offset to 3rd site in 3rd uhd row, add 1 for matlab 
                        xc(i+1) = 12;
                        xc(i+2) = 36;
                    else
                        % rows 1, 3, 5 ...
                        c0 = row*32 + 16 + 1;   % add 16 to offset to 3rd uhd row, add 1 for matlab 
                        xc(i+1) = 0;
                        xc(i+2) = 24;
                    end
                    yc(i+1) = row*24 + 12;
                    yc(i+2) = row*24 + 12;
                case 0
                    row = floor(i/2); %starts at zero
                    bEven = ~(mod(row,2));
                    if bEven 
                        % rows 0, 2, 4 ...
                        c0 = row*32 + 1;   % add 1 for matlab 
                        xc(i+1) = 0;
                        xc(i+2) = 24;
                    else
                        % rows 1, 3, 5 ...
                        c0 = row*32 + 2 + 1;   % add 2 to offset to 3rd site in row, add 1 for matlab 
                        xc(i+1) = 12;
                        xc(i+2) = 36;
                    end
                    yc(i+1) = row*24;
                    yc(i+2) = row*24;
                case 16
                    row = floor(i/2); %starts at zero
                    bEven = ~(mod(row,2));
                    if bEven 
                        % rows 0, 2, 4 ...
                        c0 = row*32 + 16 + 1;   % add 16 to offset to 2nd row, add 1 for matlab 
                        xc(i+1) = 0;
                        xc(i+2) = 24;
                    else
                        % rows 1, 3, 5 ...
                        c0 = row*32 + 16 + 2 + 1;   % add 16 to offset to 2nd row, plus 2 for 3rd site, add 1 for matlab 
                        xc(i+1) = 12;
                        xc(i+2) = 36;
                    end
                    yc(i+1) = row*24 + 12;
                    yc(i+2) = row*24 + 12;
            end


            datNew(i+1,:) = (buff(c0,:) + buff(c0+1,:) + buff(c0+8,:) + buff(c0+9,:))/4;
            datNew(i+2,:) = (buff(c0+4,:) + buff(c0+5,:) + buff(c0+12,:) + buff(c0+13,:))/4;
            saveChanStr = sprintf('%s%d,%d,',saveChanStr,c0-1,c0+3);
            saveChanArr(i+1:i+2) = [c0-1,c0+3];
        end
        nCol = 2;  % For spikeGLX shank map, staggered colums treated as two
        
    case 3
        % pool sets of 4 channels to make a simulated two column probe with
        % separation of 1, 2, 3, or 4 sites 
        % (makes pitch of 18, 24, 30, or 36).
        % the probe. There will be 32 channels in the final pattern
        % starts from lower left hand site in stSite
        Nchan = floor(nSite/12);
        Ncomb = 4;
        datNew = zeros([Nchan,NT],'double');  % make double to match datr
        xc = zeros([Nchan,1],'double');
        yc = zeros([Nchan,1],'double');
        chanMap = zeros([Nchan,1],'double');
        chanMap(:,1) = (1:Nchan);
        connected = ones([Nchan,1],'double');
        saveChanArr = zeros([Nchan,1],'double');
        % which columns
        off1 = mod(stSite,8); % offset in sites to left hand column
        off2 = off1 + 5; % set to 5 for 30 um pitch, closest to 2.0
        row_off = floor(stSite/8); 
        
        for i = 0:2:Nchan-1       %step through rows in new pattern
            row = floor(i/2);      %starts at zero
            
            c0 = row_off*8 + row*24 + off1 + 1; % add off to get to site, plus 1 for matlab
            c1 = row_off*8 + row*24 + off2 + 1;
            
            xc(i+1) = off1*6;
            xc(i+2) = off2*6;
                 
            yc(i+1) = row*18;
            yc(i+2) = row*18;
            datNew(i+1,:) = (buff(c0,:) + buff(c0+1,:) + buff(c0+8,:) + buff(c0+9,:))/4;
            datNew(i+2,:) = (buff(c1,:) + buff(c1+1,:) + buff(c1+8,:) + buff(c1+9,:))/4;
            saveChanStr = sprintf('%s%d,%d,',saveChanStr,c0-1,c1-1);
            saveChanArr(i+1:i+2) = [c0-1,c1-1];            
        end
        nCol = 2;  % These patterns really are two columns
        
    case 4
        % pool sets of 4 channels to make a simulated two column probe that tiles
        % the same area as UHD. result has 24 rows x 4 columns
        Nchan = nSite/4;
        Ncomb = 4;
        datNew = zeros([Nchan,NT],'double');  % make double to match datr
        xc = zeros([Nchan,1],'double');
        yc = zeros([Nchan,1],'double');
        saveChanArr = zeros([Nchan,1],'double');
        chanMap = zeros([Nchan,1],'double');
        chanMap(:,1) = (1:Nchan);
        connected = ones([Nchan,1],'double');
        for i = 0:4:Nchan-1       %step through rows in new pattern
            row = floor(i/4); %starts at zero
           
            % indicies for the lower left hand site for each set of 4 to combine
            c0 = row*16 + 1; % add offset = [0,2,4,6] to get to site, plus 1 for matlab
            c1 = row*16 + 3;
            c2 = row*16 + 5;
            c3 = row*16 + 7;
            
            xc(i+1:i+4) = 12*[0,1,2,3];
            yc(i+1:i+4) = 12*row;

            datNew(i+1,:) = (buff(c0,:) + buff(c0+1,:) + buff(c0+8,:) + buff(c0+9,:))/4;
            datNew(i+2,:) = (buff(c1,:) + buff(c1+1,:) + buff(c1+8,:) + buff(c1+9,:))/4;
            datNew(i+3,:) = (buff(c2,:) + buff(c2+1,:) + buff(c2+8,:) + buff(c2+9,:))/4;
            datNew(i+4,:) = (buff(c3,:) + buff(c3+1,:) + buff(c3+8,:) + buff(c3+9,:))/4;
            saveChanArr(i+1:i+4) = [c0-1,c1-1,c2-1,c3-1];
            saveChanStr = sprintf('%s%d,%d,%d,%d,',saveChanStr,c0-1,c1-1,c2-1,c3-1);
        end
        nCol = 4;  % 8 sites in original, 2 in 4x4 combination

    otherwise
        fprintf('unrecoginzed pattern type: %d', newPatType)

end
    % Generate noise with the same power spectrum as a pre-analyzed
    % UHD dataset; spectrum stored in ops.noise_model
    
if addNoise
    eNoise = makeNoise( NT, noiseModel, chanMap, connected, Nchan );
    % fprintf( 'eNoise time: %.3f\n', toc );
    % adding extra noise to make the total variance ~0.8*UHD noise one one
    % site (to be roughly equivalent to a 1.0 probe, which has lower noise
    % than a UHD probe).
    % calculating how much noise we need to add:
    % desired variance ~ (noiseFrac * sigmaUHD)^2 = (sigma combined data)^2 + (scale*sigmaUHD)^2
    % (noiseFrac^2) * sigmaUHD^2 = sigmaUHD^2/Ncomb + scale^2*sigmaUHD^2
    % scale^2 = noiseFrac^2 - 1/Ncomb
    % If Ncomb = 1 and noiseFrac < 1, we don't want to add any noise, set
    % scale = 0.
    if ((noiseFrac^2) - (1/Ncomb)) > 0
        scale = sqrt( (noiseFrac^2) - (1/Ncomb) );
    else
        scale = 0;
    end
    % eNoise is returned in uV, scale to bits for the probe type of
    % interest (NP 1.0)
    scale = scale/2.438;    % NP 1.0, gain = 500, 2.438 uV/bit
    % fprintf( 'Channel 1 eNoise: %.3f mV\n', rms(eNoise(:,1)));
    datNew = datNew + scale*eNoise';
end

    % make new shank map string
    % replaced with indexing into original shank map
    % newShankMap = makeShankMap(xc, yc, nCol);
    
end

function eNoise = makeNoise( noiseSamp,noiseModel,chanMap,connected,NchanTOT )

    %if chanMap is a short version of a 3A probe, use the first
    %nChan good channels to generate noise, then copy that array 
    %into an NT X NChanTot array
    
    nChan = numel(chanMap);        %number of noise channels to generate
    goodChan = sum(connected);
    tempNoise = zeros( noiseSamp, goodChan, 'single' );
    nT_fft = noiseModel.nm.nt;         %number of time points in the original time series
    fftSamp = noiseModel.nm.fft;
    
    noiseBatch = ceil(noiseSamp/nT_fft);    
    lastWind = noiseSamp - (noiseBatch-1)*nT_fft; %in samples
    
    for j = 1:goodChan     
            for i = 1:noiseBatch-1
                tStart = (i-1)*nT_fft+1;
                tEnd = i * nT_fft;            
                tempNoise(tStart:tEnd,j) = fftnoise(fftSamp(:,j),1);
            end
            %for last batch, call one more time and truncate
            lastBatch = fftnoise(fftSamp(:,j),1);
            tStart = (noiseBatch-1)*nT_fft+1;
            tEnd = noiseSamp;
            tempNoise(tStart:tEnd,j) = lastBatch(1:lastWind);             
    end
    
    %unwhiten this array
    Wrot = noiseModel.nm.Wrot(1:goodChan,1:goodChan);
    tempNoise_unwh = tempNoise/Wrot;
    
    %scale to uV; will get scaled back to bits at the end
    tempNoise_unwh = tempNoise_unwh/noiseModel.nm.bitPerUV;
    
    %to get the final noise array, map to an array including all channels
    eNoise = zeros(noiseSamp, NchanTOT, 'single');
    %indicies of the good channels
    goodChanIndex = find(connected);
    eNoise(:,chanMap(goodChanIndex)) = tempNoise_unwh;
    
end

function noise=fftnoise(f,Nseries)
% Generate noise with a given power spectrum.
% Useful helper function for Monte Carlo null-hypothesis tests and confidence interval estimation.
%  
% noise=fftnoise(f[,Nseries])
%
% INPUTS:
% f: the fft of a time series (must be a column vector)
% Nseries: number of noise series to generate. (default=1)
% 
% OUTPUT:
% noise: surrogate series with same power spectrum as f. (each column is a surrogate).
%
%   --- Aslak Grinsted (2009)
%  
if nargin<2
    Nseries=1;
end
f=f(:);     %ensures f is a column vector
N=length(f); 
Np=floor((N-1)/2);
phases=rand(Np,Nseries)*2*pi;
phases=complex(cos(phases),sin(phases)); % this was the fastest alternative in my tests. 
f=repmat(f,1,Nseries);
f(2:Np+1,:)=f(2:Np+1,:).*phases;
f(end:-1:end-Np+1,:)=conj(f(2:Np+1,:));
noise=real(ifft(f,[],1)); 

end

function newMap = makeShankMap(xcoord, ycoord, nCol)
% for a set of x and y coordinate and number of columns, build a SGLX shank
% map string. Need to know the number of columns because staggered patterns
% are treated as 2 columns in SGLX

    yVal = unique(ycoord);
    nRow = numel(yVal);
    xVal = unique(xcoord);
    nElec = numel(xcoord);

    col = zeros([nElec,1]);
    row = zeros([nElec,1]);
    if nCol == 2
        col = xcoord > median(xVal);    %automatically zero based
        for i = 1:nElec 
            row(i) = find(yVal == ycoord(i))-1; % -1 for zero based
        end
    else
        for i = 1:nElec 
            col(i) = find(xVal == xcoord(i))-1;
            row(i) = find(yVal == ycoord(i))-1;
        end
    end

    %header line (shank, number of columns, number of rows)
    newMap = sprintf('(1,%d,%d)', nCol, nRow);
    for i = 1:nElec
        currEntry = sprintf('(0:%d:%d:1)', col(i), row(i)); 
        newMap = [newMap,currEntry];
    end
    

end

