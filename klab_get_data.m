
function [DATA,CELLS,label_str] = klab_get_data(dataOut,class_labels,THRESHOLD,BINNING)
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

N_max_bad_frames = floor(N_bin_indices*0.25);

label = nan(1,N_labels);

if not(size(class_labels,1)==N_labels)
    error('Number of stimulus blocks does not match!')
end

k=0;
kk=0;
for i=1:N_labels
    for j=1:size(class_labels,1)
        if strcmp(class_labels{j,1},dataOut.blockLabels{i})
            k=k+1;
            if class_labels{j,2}>0
                kk=kk+1;
                label(i) = kk;
                label_str{i} = dataOut.blockLabels{i};
            else
                label(i) = 0;
            end
        end
    end
end
if k<N_labels
    error('not all labels were present in dataOut!');
end
valid = label>0;
label_ind = find(valid);
label_str = dataOut.blockLabels(label_ind);

N_labels = length(label_ind);

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

DATA = zeros(N_cells,N_labels,N_bins,dataOut.totalNumTrials);
DATA_isbad = DATA;

SESSION = zeros(N_cells,N_labels,N_bins,dataOut.totalNumTrials);

k=0;
for trial = 1:length(dataOut.trial)           
   for repeat = 1:dataOut.trial(trial).stimuliRepeats
       k=k+1;
       k1=0;
       for stim = label_ind
            k1=k1+1;
            vec = find((dataOut.trial(trial).stimulus_vector==stim).*(dataOut.trial(trial).repetition_vector==repeat));
            signal = dataOut.trial(trial).signal_deconv(CELLS,vec);            
            isbad = ~dataOut.trial(trial).isValidFrame(vec);
            for bin = 1:N_bins   
                binind = bin_indices(bin):(bin_indices(bin+1)-1);
                DATA(:,k1,bin,k) = mean(signal(:,binind),2);
                
                SESSION(:,k1,bin,k) = trial;
                
                if sum(isbad(binind))>N_max_bad_frames(bin)                  
                    DATA_isbad(:,k1,bin,k) = 1;
                end
            end
       end
   end
end

file_mean = nan(1,length(dataOut.trial));
for trial = 1:length(dataOut.trial)
    ind = SESSION==trial;
    file_mean(trial) = mean(vectorize(DATA(ind)));
    DATA(ind) = DATA(ind)-file_mean(trial);
end

fprintf('file-wise means: %s\n',num2str(file_mean,'%0.3f, '));

bad_cells = [];
for c=1:size(DATA,1)
    for stim=1:size(DATA,2)        
        sig = squeeze(DATA(c,stim,:,:));        
        if sum(sum(sig>1e-6)>0)<size(DATA,4)*0.15
            bad_cells(end+1)=c;
            break;
        end        
    end
end
bad_cells = unique(bad_cells);

DATA(bad_cells,:,:,:)=[];
CELLS(bad_cells)=[];

if length(bad_cells)>0
   fprintf('\n\n!!!!! dropping %i cells not responding enough !!!\n\n',length(bad_cells)); 
end

end

function mat = flatten(mat)
    mat = mat(:);
end

function mat = nodiag(mat)
    mat = mat - diag(diag(mat));
end
