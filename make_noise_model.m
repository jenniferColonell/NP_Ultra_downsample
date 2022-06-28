function make_noise_model
% Calculate the power spectrum and spatial correlation for all channels in 
% a "typical" dataset to use as a noise model. These noise models are used
% in two contexts:

% (1) In 'pure' simulations starting with model units, create noise with
% realistic frequency and spatial frequencies by creating a model from 
% in vivo recordings from a similar brain region. The spatial correlations
% of the background will be encoded in the whitening matrix (Wrot) and the
% frequency spectrum is encoded in the fft.
%
% (2) For simulating low density probes by averaging channels measured on
% high density probes (e.g. Neuropixels Ultra), create realistic single 
% channel noise levels by creating a model from an in vitro noise test
% recording -- the goal in this case is to capture the frequency spectrum
% of all of the non-biological noise in the recording.

% Data is filtered and whitened before fft, to make sure all data is 
% within the frequency range. This script uses the KS2 preprocessing code
% for filtering and whitening, copied to this file for convenience.
% Kilosort 2 repo: (https://github.com/MouseLand/Kilosort)

% start by building a rez structure -- will be used to call the whitening
% and filtering functions drawn from Kilosort2 

% params normally read in from a config file. These need to match the model
% data set from which the noise will be extracted.

%% USER PARAMS

%path to your data file; results will also be written to this folder
rootZ = 'E:\MEArec_for_localization\UHD512';
%name of the binary file
dataName = 'UHD512_120U_10uV_g0_t0.imec0.ap.bin';
% name of kilsort style chanMap file, located in the rootZ directory
mapName = 'UHD512_120U_10uV_kilosortChanMap.mat';


% sample rate
rez.ops.fs = 30000;  
rez.ops.bitPerUV = 0.4267; %(1/uVPerBit, 1/2.3437 for NP 1.0)
rez.ops.trange = [0 8]; % time range to use when extracting the whitening matrix
rez.ops.NchanTOT    = 512; % total number of channels in your recording, including digital

cm = load(fullfile(rootZ, mapName));
rez.ops.chanMap = cm.chanMap; % 1-based indicies of channels in file. If channels need to be excluded, they should be omitted from chanMap
rez.xc = cm.xcoords;
rez.yc = cm.ycoords;
rez.ops.whiteningRange      = 32; % number of channels to use for whitening each channel

% frequency for high pass filtering (150)
rez.ops.fshigh = 150;  
%Time to sample in sec, need to sample at least enough to capture freq down to the high pass limit
clipDur = 0.25; 
%directory for temporary whitened data file
rootH = 'D:\kilosort_datatemp\';
%%


%processing params usually read in from config. These should not be changed
rez.ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry
rez.ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
rez.ops.NT                  = 64*1024+ rez.ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory). 
rez.ops.nSkipCov            = 25; % compute whitening matrix from every N-th batch
rez.ops.scaleproc           = 100;   % int16 scaling of whitened data
rez.ops.CAR                 = 1; % do CAR

rez.ops.fproc       = fullfile(rootH, 'temp_wh.imec.ap.bin'); % proc file on a fast SSD

rez.ops.fbinary = fullfile(rootZ, dataName);

% call preprocessDataSub to filter data and create the whitening matrix.
[rez] = preprocessDataForNoise(rez);

%calculate fft for each channel in the whitened data (stored in ops.fproc)
tic;
[fftSamp, nT_fft] = sampFFT( rez, clipDur );
fprintf( 'time to calc fft: %f\n', toc );
%use the fft to generate noise for each channel for a time period 
%Can only make batches of nT points at a time

    noiseT = 2.2;           %in seconds, chosen for ~66000 points
    noiseSamp = noiseT*rez.ops.fs;
    nChan = numel(rez.ops.chanMap);
    eNoise = zeros( noiseSamp, nChan, 'single' );
    
    noiseBatch = ceil(noiseSamp/nT_fft);    
    lastWind = noiseSamp - (noiseBatch-1)*nT_fft; %in samples
    fprintf( 'noiseSamp, batchSamp, lastWind: %d, %d, %d\n', noiseSamp, nT_fft, lastWind);
    
    for j = 1:nChan     
            for i = 1:noiseBatch-1
                tStart = (i-1)*nT_fft+1;
                tEnd = i * nT_fft;            
                eNoise(tStart:tEnd,j) = fftnoise(fftSamp(:,j),1);
            end
            %for last batch, call one more time and truncate
            tempNoise = fftnoise(fftSamp(:,j),1);
            tStart = (noiseBatch-1)*nT_fft+1;
            tEnd = noiseSamp;
            eNoise(tStart:tEnd,j) = tempNoise(1:lastWind);             
    end
    
    selectChan = [10,11,12,13,14];  %channels to plot, no effect on model output
        
    tR = [4200:5200]; %just looking at 1000 samples
    h = figure(1);  
    for k = 1:numel(selectChan)  
        currNoise = eNoise(tR,selectChan(k));
        plot( tR, currNoise + k*500 );
        hold on;
    end  
    title("1000 samples of generated noise, before unwhitening (scaled by fscale)");
    hold off

    
    %unwhiten this noise array
    eNoise_unwh = eNoise/rez.Wrot;
    
    %to get the final noise array, map to an array including all channels
    nAllChan = max(rez.ops.chanMap);
    eNoise_final = zeros(noiseSamp, nAllChan, 'single');
    eNoise_final(:,rez.ops.chanMap) = eNoise_unwh;
    

    tR = [4200:5200]; %plot 1000 samples post unwhitening
    h = figure(2);  
    for k = 1:numel(selectChan)  
        currNoise = eNoise_final(tR,selectChan(k));
        plot( tR, currNoise + k*50 );
        hold on;
    end   
    title('1000 samples of generated noise, after unwhitening');
    hold off
   
    nm.fft = fftSamp;
    nm.nt = nT_fft;
    nm.Wrot = rez.Wrot;
    nm.chanMap = rez.ops.chanMap;
    nm.bitPerUV =  rez.ops.bitPerUV;
    
%save rez file
fprintf('Saving rez file  \n');
fname = fullfile(rootZ, 'rez.mat');
save(fname, 'rez', '-v7.3');

%save noise model file
fprintf('Saving nm file  \n');
fname = fullfile(rootZ, 'noiseModel.mat');
save(fname, 'nm', '-v7.3');

end



function [rez] = preprocessDataForNoise(rez)

%build rez structure -- code taken from preprocessDataSub in Kilsort 2
%Only changes are to skip checking for low spike rate samples and
%write out the data as rows = time -- just a convenience so it 
%it can be read in using the same code

rez.ops.nt0 	= 61;
rez.ops.nt0min  = 20;

NT       = rez.ops.NT ;
NchanTOT = rez.ops.NchanTOT;

bytes = get_file_size(rez.ops.fbinary);
nTimepoints = floor(bytes/NchanTOT/2);
rez.ops.tstart = ceil(rez.ops.trange(1) * rez.ops.fs);
rez.ops.tend   = min(nTimepoints, ceil(rez.ops.trange(2) * rez.ops.fs));
rez.ops.sampsToRead = rez.ops.tend-rez.ops.tstart;
rez.ops.twind = rez.ops.tstart * NchanTOT*2;

Nbatch      = ceil(rez.ops.sampsToRead /(NT-rez.ops.ntbuff));
rez.ops.Nbatch = Nbatch;


rez.ops.igood = true(size(rez.ops.chanMap));

rez.ops.Nchan = numel(rez.ops.chanMap);


NTbuff      = NT + 4*rez.ops.ntbuff;

% by how many bytes to offset all the batches
rez.ops.NTbuff = NTbuff;


%build whitening matrix. This function filters before calculating the
%cross correlation

tic;

fprintf('Time %3.0fs. Computing whitening matrix.. \n', toc);

Wrot = get_whitening_matrix(rez);

%Apply the whitening matrix to a subset of the original data

rez.ops.tstart = ceil(rez.ops.trange(1) * rez.ops.fs);
rez.ops.tend   = min(nTimepoints, ceil(rez.ops.trange(2) * rez.ops.fs));
rez.ops.sampsToRead = rez.ops.tend-rez.ops.tstart;
rez.ops.twind = rez.ops.tstart * NchanTOT*2;
Nbatch      = ceil(rez.ops.sampsToRead /(NT-rez.ops.ntbuff));
rez.ops.Nbatch = Nbatch;

fprintf('Time %3.0fs. Loading raw data and applying filters... \n', toc);

fid         = fopen(rez.ops.fbinary, 'r');

fidW        = fopen(rez.ops.fproc,   'w');

for ibatch = 1:Nbatch
    offset = max(0, rez.ops.twind + 2*NchanTOT*((NT - rez.ops.ntbuff) * (ibatch-1) - 2*rez.ops.ntbuff));
    if offset==0
        ioffset = 0;
    else
        ioffset = rez.ops.ntbuff;
    end
    fseek(fid, offset, 'bof');
    
    buff = fread(fid, [NchanTOT NTbuff], '*int16');
    if isempty(buff)
        break;
    end
    nsampcurr = size(buff,2);
    if nsampcurr<NTbuff
        buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr);
    end

    datr    = gpufilter(buff, rez.ops, rez.ops.chanMap); % apply filters and median subtraction

    datr    = datr(ioffset + (1:NT),:); % remove timepoints used as buffers

    datr    = datr * Wrot; % whiten the data and scale by 200 for int16 range

    datcpu  = gather(int16(datr)); % convert to int16, and gather on the CPU side
  
    %next step will be reading in this file to get FTs
    %Most convenient for that process to have the data in the original
    %format (row = time, column = channel index), so write the transpose
    % (differs from the standard preprocessing code)
    fwrite(fidW, datcpu', 'int16');
end

Wrot        = gather_try(Wrot);
rez.Wrot    = Wrot;

fclose(fidW);
fclose(fid);


end


function Wrot = get_whitening_matrix(rez)
% based on a subset of the data, compute a channel whitening matrix
% this requires temporal filtering first (gpufilter)
% drawn directly from Kilosort2 by Marius Pachitariu (https://github.com/MouseLand/Kilosort)
% copied here for convenience

ops = rez.ops;
Nbatch = ops.Nbatch;
twind = ops.twind;
NchanTOT = ops.NchanTOT;
NT = ops.NT;
NTbuff = ops.NTbuff;
chanMap = ops.chanMap;
Nchan = rez.ops.Nchan;
xc = rez.xc;
yc = rez.yc;

% % load data into patches, filter, compute covariance
% if isfield(ops,'fslow')&&ops.fslow<ops.fs/2
%     [b1, a1] = butter(3, [ops.fshigh/ops.fs,ops.fslow/ops.fs]*2, 'bandpass');
% else
%     [b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high');
% end

fprintf('Getting channel whitening matrix... \n');
fid = fopen(ops.fbinary, 'r');
CC = gpuArray.zeros( Nchan,  Nchan, 'single'); % we'll estimate the covariance from data batches, then add to this variable


ibatch = 1;
while ibatch<=Nbatch
    offset = max(0, twind + 2*NchanTOT*((NT - ops.ntbuff) * (ibatch-1) - 2*ops.ntbuff));
    fseek(fid, offset, 'bof');
    buff = fread(fid, [NchanTOT NTbuff], '*int16');

    if isempty(buff)
        break;
    end
    nsampcurr = size(buff,2);
    if nsampcurr<NTbuff
        buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr);
    end

    datr    = gpufilter(buff, ops, rez.ops.chanMap); % apply filters and median subtraction

    CC        = CC + (datr' * datr)/NT; % sample covariance

    ibatch = ibatch + ops.nSkipCov; % skip this many batches
end
CC = CC / ceil((Nbatch-1)/ops.nSkipCov); % normalize by number of batches

fclose(fid);

if ops.whiteningRange<Inf
    % if there are too many channels, a finite whiteningRange is more robust to noise in the estimation of the covariance
    ops.whiteningRange = min(ops.whiteningRange, Nchan);
    Wrot = whiteningLocal(gather(CC), yc, xc, ops.whiteningRange); % this function performs the same matrix inversions as below, just on subsets of channels around each channel
else
    Wrot = whiteningFromCovariance(CC);
end
Wrot    = ops.scaleproc * Wrot; % scale this from unit variance to int 16 range. The default value of 200 should be fine in most (all?) situations.

fprintf('Channel-whitening matrix computed. \n');

end


function [fftSamp, nT] = sampFFT( rez, clipDur )

    ops = rez.ops;
    
    %get file of whitened data
    datFile = ops.fproc;
    
    %read in some data
    nChansInFile = numel(ops.chanMap);  % channels in whitenened data, excludes ref chans in original 
    d = dir(ops.fproc); 
    nSamps = d.bytes/2/nChansInFile;
    tSkip = 1.0; % skip first tSkip seconds
    
    sampStart = round(ops.fs*tSkip); 
    nClipSamps = round(ops.fs*clipDur);
    mmf = memmapfile(datFile, 'Format', {'int16', [nChansInFile nSamps], 'x'});
    thisDat = (double(mmf.Data.x(:, (1:nClipSamps)+sampStart)));
    %subtract the DC
    thisDat = bsxfun(@minus, thisDat, mean(thisDat,2));
    [~,nT] = size(thisDat);
    %data is formatted as Nchan rows by nt columns. fft returns the
    %transform of the columns, so transpose
    fftSamp = fft(thisDat');
    
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

function datr = gpufilter(buff, ops, chanMap)
% filter this batch of data after common average referencing with the
% median
% buff is timepoints by channels
% chanMap are indices of the channels to be kep
% ops.fs and ops.fshigh are sampling and high-pass frequencies respectively
% if ops.fslow is present, it is used as low-pass frequency (discouraged)

% set up the parameters of the filter
if isfield(ops,'fslow')&&ops.fslow<ops.fs/2
    [b1, a1] = butter(3, [ops.fshigh/ops.fs,ops.fslow/ops.fs]*2, 'bandpass'); % butterworth filter with only 3 nodes (otherwise it's unstable for float32)
else
    [b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high'); % the default is to only do high-pass filtering at 150Hz
end

dataRAW = gpuArray(buff); % move int16 data to GPU
dataRAW = dataRAW';
dataRAW = single(dataRAW); % convert to float32 so GPU operations are fast
dataRAW = dataRAW(:, chanMap); % subsample only good channels

% subtract the mean from each channel
dataRAW = dataRAW - mean(dataRAW, 1); % subtract mean of each channel

% CAR, common average referencing by median
if ops.CAR
    dataRAW = dataRAW - median(dataRAW, 2); % subtract median across channels
end

%next four lines should be equivalent to filtfilt (which cannot be used because it requires float64)
datr = filter(b1, a1, dataRAW); % causal forward filter
datr = flipud(datr); % reverse time
datr = filter(b1, a1, datr); % causal forward filter again
datr = flipud(datr); % reverse time back

end

function Wrot = whiteningLocal(CC, yc, xc, nRange)
% function to perform local whitening of channels
% CC is a matrix of Nchan by Nchan correlations
% yc and xc are vector of Y and X positions of each channel
% nRange is the number of nearest channels to consider
    Wrot = zeros(size(CC,1), size(CC,1));
    for j = 1:size(CC,1)
        ds          = (xc - xc(j)).^2 + (yc - yc(j)).^2;
        [~, ilocal] = sort(ds, 'ascend');
        ilocal      = ilocal(1:nRange); % take the closest channels to the primary channel. First channel in this list will always be the primary channel.
    
        wrot0 = whiteningFromCovariance(CC(ilocal, ilocal));
        Wrot(ilocal, j)  = wrot0(:,1); % the first column of wrot0 is the whitening filter for the primary channel
    end
end


function Wrot = whiteningFromCovariance(CC)
% takes as input the matrix CC of channel pairwise correlations
% outputs a symmetric rotation matrix (also Nchan by Nchan) that rotates
% the data onto uncorrelated, unit-norm axes
    [E, D] 	= svd(CC); % covariance eigendecomposition (same as svd for positive-definite matrix)
    D       = diag(D); % take the non-zero values from the diagonal
    eps 	= 1e-6;
    Wrot 	= E * diag(1./(D + eps).^.5) * E'; % this is the symmetric whitening matrix (ZCA transform)
end


function x = gather_try(x)
    try
        x = gather(x);
    catch
    end
end