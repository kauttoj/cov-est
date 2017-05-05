
function result = klab_get_data(dataOut,class_labels,THRESHOLD,BINNING)
% Berkson's paradox

warning('on','all');

rng(666);

DT = 1/dataOut.AquisitionParams.framerate;

N_bin_frames = round(BINNING/DT);

N_labels = length(dataOut.blockLabels);

frames = [];
for trial = 1:length(dataOut.trial)
   for repeat = 1:dataOut.trial(trial).stimuliRepeats       
       for stim = 1:N_labels
            frames(end+1) = sum((dataOut.trial(trial).stimulus_vector==stim).*(dataOut.trial(trial).repetition_vector==repeat));
       end
   end
end
N_blocks = median(frames);

bin_indices = 1:N_bin_frames:N_blocks;
bin_indices(end+1) = round(N_blocks);
bin_indices = unique(bin_indices);

N_bin_indices = diff(bin_indices);

%bin_window = mean(N_bin_indices)*DT;

if N_bin_indices(end)<mean(N_bin_indices(1:end-1)*0.5)
    N_bin_indices(end) = [];
    bin_indices(end) = [];
    warning('Removed last bin with less than 50% of frames!!')
end

N_bins = length(bin_indices)-1;

N_max_bad_frames = floor(sum(N_bin_indices)*0.25);

result.N_max_bad_frames = N_max_bad_frames;
result.N_bin_indices=N_bin_indices;

label = nan(1,N_labels);

if not(length(class_labels)==N_labels)
    error('Number of stimulus blocks does not match!')
end

k=0;
for i=1:N_labels
    for j=1:length(class_labels)
        if strcmp(class_labels{j},dataOut.blockLabels{i})
            k=k+1;
            label(i) = k;
            label_str{i} = dataOut.blockLabels{i};
        end
    end
end
if k<N_labels
    error('not all labels were present in dataOut!');
end
valid = label>0;
label_ind = find(valid);
result.labels = dataOut.blockLabels(label_ind);

if ~isempty(THRESHOLD)
    if THRESHOLD==1
        CELLS = dataOut.stats.global.responsive_cells_p005_fdr_average(:)';
    else
        CELLS = THRESHOLD(:)';
    end
else
    CELLS = 1:dataOut.totalNumCells;
end
N_cells = length(CELLS);

result.CELL_IDs = CELLS;

DATA = nan(N_cells,length(dataOut.trial),length(label_ind),dataOut.trial(1).stimuliRepeats,N_bins);
DATA_isbad = zeros(size(DATA));
DATA_neuropil = DATA;
DATA_baseline = DATA;

for file = 1:length(dataOut.trial)           
    for stim = 1:length(label_ind)
        for repeat = 1:dataOut.trial(file).stimuliRepeats       
            vec = find((dataOut.trial(file).stimulus_vector==label_ind(stim)).*(dataOut.trial(file).repetition_vector==repeat));
            signal = dataOut.trial(file).signal_deconv(CELLS,vec);            
            
            isbad = ~dataOut.trial(file).isValidFrame(vec);
            for bin = 1:N_bins   
                binind = bin_indices(bin):(bin_indices(bin+1)-1);
                DATA(:,file,stim,repeat,bin) = mean(signal(:,binind),2);                
            end
                        
            vec1 = max(1,vec(1)-10):min(length(dataOut.trial(file).stimulus_vector),vec(end)+10);
            assert(length(vec1)>length(vec));
            signal_neuropil = dataOut.trial(file).signal_neuropil(CELLS,vec1);
            m = median(signal_neuropil,2);
            assert(all(m>0));
            DATA_neuropil(:,file,stim,repeat,:) = repmat(m,[1,N_bins]);
            
            if sum(isbad)>N_max_bad_frames
                DATA_isbad(:,file,stim,repeat,:) = 1;
            end
            
            signal_baseline = dataOut.trial(file).signal_raw(CELLS,vec1);
            m = arrayfun(@(i) median(signal_baseline(i,signal_baseline(i,:)<prctile(signal_baseline(i,:),50))),1:size(signal_baseline,1))';
            assert(all(m>0));
            DATA_baseline(:,file,stim,repeat,:) = repmat(m,[1,N_bins]);
            
            if sum(isbad)>N_max_bad_frames
                DATA_isbad(:,file,stim,repeat,:) = 1;
            end            
       end
   end
end

assert(nnz(isnan(DATA))==0);

result.DATA_dimensions = {'cell','file','stimulus','repeat','bin'};
result.DATA_isbad = DATA_isbad;
result.DATA_neuropil = DATA_neuropil;
result.DATA_baseline = DATA_baseline;

median_response = median(vectorize(DATA(DATA>0)));

result.median_response = median_response;

coeff = dataOut.trial(file).neuropil_coefficient(CELLS);

bad_cells = [];
for cell=1:size(DATA,1)
    a = mean(DATA(cell,:,:,:,:),5);
    a = vectorize(permute(a,[2,4,3,1]));
    bad_ratio = sum(a<median_response/1000)/length(a);
    if bad_ratio > 0.75 || coeff(cell)>1.2 || coeff(cell)<0.25
        bad_cells(end+1) = cell;
    end
end
bad_cells = unique(bad_cells);

result.bad_cells = CELLS(bad_cells);

result.CELL_IDs = CELLS;

result.max_allowed_bad_ratio_limit = 0.75;

DATA(bad_cells,:,:,:,:)=[];
CELLS(bad_cells)=[];

result.DATA = DATA;

result.CELL_IDs = CELLS;

if length(bad_cells)>0
   fprintf('\n\n!!!!! dropping %i bad cells (not responding or bad neuropil fit) !!!\n\n',length(bad_cells)); 
end

end

function mat = vectorize(mat)
    mat = mat(:);
end

function mat = nodiag(mat)
    mat = mat - diag(diag(mat));
end
