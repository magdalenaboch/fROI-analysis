
% ---------------------------------------------------------------------
% Script to define individual functional regions of interest (fROIs)
% based on top % active voxels within defined search spaces
%
% Magdalena Boch, magdalena.boch@univie.ac.at
% ---------------------------------------------------------------------

clc
clear

disp('creating fROIs');

%% ======================================================

% This section needs your input

% subject IDs
% see also loop subj = 1:numel(subjects) to add letters
subjects = [2,5:7,9,11:14,16,18,21,23,27,28];

% percentage of active voxels to define fROI?
% can be one number or a range of percentages
percentages = [0.01,0.05:0.05:1];

% take all anatomical search spaces in search space folder?
allsearchspaces = false;

% or select only specific search spaces
if allsearchspaces == false
    clear sspaces; sspaces = {'caudal_ssg_L', 'caudal_ssg_R', 'medial_ssg_L', 'medial_ssg_R'};
end

% indicate the first level contrastname to define fROIs
% (assuming you already created folders for each contrast, containing the contrast maps from all
% subject after your first level analysis) 

% this contrast will mask the search space
% i.e. only voxels that have positive activation levels for faces > bodies
% within the given anatomical mask define the search space

sspace_masking = 'faces_bodies_ss1'; % within search space

% this is the contrast to define the fROI within the search space
contrast = 'faces_objects_ss1'; % top % voxel define fROI

catname = 'face_fROIs'; %  e.g. 'body_fROIs' or 'face_fROIs'
% this will determin to which category the fROI belongs to, it is important
% because for the dog data we use the same anatomical search spaces
% but the contrast name will also be saved in fROI name just to be sure

% ---------------------------------------------------------------------
% Directories

dirs.home = ('YOUR PROJECT PATH');

% directory containing scripts, e.g.
dirs.scripts = fullfile(dirs.home,'scripts'); 

% Directory containing contrast images 
% (has subfolders for each contrast
% containing  contrast maps from all subjects)

% e.g.
dirs.glm = fullfile(dirs.home, 'data_glm', 'contrasts');

% where I want to save the fROIs, e.g.
dirs.analysis = fullfile(dirs.home, 'fROI_analysis');

% create subfolder for fROI type (face/body)
dirs.froi = fullfile(dirs.analysis,'individual_fROIs', catname);

% directory containing search spaces (i.e., anatomical masks)
dirs.sspaces = fullfile(dirs.analysis,'search_spaces');

addpath(dirs.scripts); addpath(dirs.sspaces); addpath(dirs.glm); addpath(dirs.glm);

%% ===========================================================================
% let's get started

cd(dirs.sspaces);

% select all search spaces if you didn't specify specific ones
if allsearchspaces ==  true
    selectsspaces = dir('*.nii');
    sspaces = [];
    for r = 1:numel(selectsspaces)
        sspaces{r} = selectsspaces(r).name(1:end-4); %delete '.nii'
    end
else
end

% create empty table for subject level information (size for each
% percentage etc.
headers = {'subjname','contrast', 'searchspace', 'activevoxels', 'percentage','nvoxel'};
nlines = numel(sspaces)*numel(subjects)*numel(percentages);
infosubjlvl = cell2table(cell(nlines,numel(headers)),'VariableNames',headers);
cnttable = 1; % index for table

% create empty table for overview of overall active voxels per subject
% is also in more detailed table
header = {'subjname','contrast', 'searchspace', 'activevoxels'};
nlines = numel(sspaces)*numel(subjects);
infoactivevoxels = cell2table(cell(nlines,numel(header)),'VariableNames',header);
idxtable = 1;


%% START LOOPS
% we have 3 loops
% 1: search spaces, 2: subjects, 3: percentage

% ------------------------------------------------
% search space loop

for nsspaces = 1:numel(sspaces)
    
    % get anatomical mask or parcel
    cd(dirs.sspaces);
    clear filename; filename = dir([sspaces{nsspaces},'.nii']);
    
    disp(sspaces{nsspaces});
    
    % load it
    clear vol; vol = spm_vol(filename.name);
    clear mask; mask_sspace = spm_read_vols(vol);
    
    % make it binary
    clear ssbinary
    ssbinary = zeros(vol.dim(1), vol.dim(2), vol.dim(3));
    ssbinary(mask_sspace > 0.25) = 1;
    
    % ------------------------------------------------
    % subject-level loop: load contrasts per subject
    
    for subj = 1:numel(subjects)
        
        if subjects(subj) < 10; subjname = ['F0',num2str(subjects(subj))];
        else subjname = ['F',num2str(subjects(subj))]; end
        
        disp (['subject: ', subjname]);
        
        
        % ----------
        
        % first restrain search space to voxels with stronger signal for
        % for faces or bodies (i.e, sspace_masking)
        
        % get masking contrast 
        clear maskfolder; maskfolder = fullfile(dirs.glm, ['cons_',sspace_masking]);
        cd(maskfolder);
        
        %get contrast name
        clear con; con_mask = dir([subjname,'_',sspace_masking,'*.nii']);
        
        %read contrast
        clear vol; vol = spm_vol(con_mask.name);
        clear mask_condata; mask_condata = spm_read_vols(vol);
        
        % label negative values as NaNs
        mask_condata(mask_condata <= 0) = NaN;
        
        % make positive values to binary mask
        mask_condata(mask_condata > 0) = 1;
        
        % 1: create empty mask
        clear mbinary; mbinary = nan(vol.dim(1), vol.dim(2), vol.dim(3));
        
        % 2: get indices of overlapping voxel between contrast mask and search
        % space
        idx_m = find((mask_condata ==1)&(ssbinary == 1));
        
        
        % 3: "fill" overlap with activation data
        mbinary(idx_m) = mask_condata(idx_m);
    
    
        % ---------
        
        % get fROI defining contrast 
        clear confolder; confolder = fullfile(dirs.glm, ['cons_',contrast]);
        cd(confolder);
        
        %get contrast name
        clear con; con = dir([subjname,'_',contrast,'*.nii']);
        
        %read contrast
        clear vol; vol = spm_vol(con.name);
        clear condata; condata = spm_read_vols(vol);
        
        
        % get overlap between searchspace and contrast image
        
        % 1: create empty vol
        clear overlap; overlap = nan(vol.dim(1), vol.dim(2), vol.dim(3));
        
        % 2: get indices of overlapping voxel between contrasts and search
        % space
        idx = find(~isnan(condata)&(mbinary == 1));
        
        % 3: "fill" overlap with activation data
        overlap(idx) = condata(idx);
        
        % remove voxels with negative values from contrast
        % we only want voxel more sensitive for condition (i.e. faces/bodies)
        % e.g bodies > objects, i only want voxel with more activation for bodies
        
        % 1: get correct dimensions
        clear posoverlap; posoverlap = overlap;
        % 2: negative values are labelled as NaNs
        posoverlap(posoverlap < 0) = NaN;
        
        % 3: get all positiv activation values (v)
        v = posoverlap(~isnan(posoverlap));
        
        % 4: sort them in descending order
        v = sort(v, 'descend');
 
        % -------------------------------------------------
        % percentage loop: define fROI based on percentages
        
        for nperc = 1:numel(percentages)
            
            
            if numel(v) >= 10 % doesn't make sense if active voxels are less than 10
                
                %calculate number of voxels = xx % of all active voxels
                clear nvoxel; nvoxel = ceil(numel(v)* percentages(nperc));
                % we always round up (ceil()) because with low percentages
                % nvoxel could be < 0.05 and then it would round to 0
                
                %get critical value / cutoff (c)
                c = v(nvoxel); %nvoxel = number of voxels i want in ROI (set on top of script)
                
                %create empty functional ROI
                % this time zeros because I want a binary mask containing top
                % nVoxel
                clear froi; froi = zeros(vol.dim(1), vol.dim(2), vol.dim(3));
                froi(posoverlap >= c) = 1;
                
                %create empty based on contrast dimensions
                clear newvol; newvol = vol;
                
                newvol.fname = [subjname, '_fROI_', contrast,'_', ...
                    num2str(nvoxel),'voxel_',num2str(percentages(nperc)*100),'perc','.nii'];
                
                newvol.descrip = ['fROI based on ', newvol.descrip];
                
                % get fROI subfolder
                fROIsubfolder = fullfile(dirs.froi,contrast,sspaces{nsspaces}, subjname);
                if ~exist(fROIsubfolder, 'dir'); mkdir(fROIsubfolder); end
                cd(fROIsubfolder);
                
                % write fROI
                newvol = spm_write_vol(newvol,froi);
                
            elseif numel(v) < 10 % don't write fROIs but save data
                
                disp ([subjname, ' not enough active voxel']);
                
                 %calculate number of voxels = xx % of all active voxels
                clear nvoxel; nvoxel = ceil(numel(v)* percentages(nperc));
                % we always round up (ceil()) because with low percentages
                % nvoxel could be < 0.05 and then it would round to 0
                
                
            end % if nvoxel > 10 loop
            
            % write down subject level information            
            infosubjlvl.subjname{cnttable, 1} = subjname;
            infosubjlvl.activevoxels{cnttable, 1} = numel(v);
            infosubjlvl.contrast{cnttable, 1} = contrast;
            infosubjlvl.searchspace{cnttable, 1} = sspaces{nsspaces};
            infosubjlvl.percentage{cnttable, 1} = percentages(nperc)*100; 
            infosubjlvl.nvoxel{cnttable, 1} = nvoxel;
            
            cnttable = cnttable + 1;
            
        end % percentage loop
        
        infoactivevoxels.subjname{idxtable,1} = subjname;
        infoactivevoxels.activevoxels{idxtable,1} = numel (v);
        infoactivevoxels.contrast{idxtable,1} = contrast;
        infoactivevoxels.searchspace{idxtable,1} = sspaces{nsspaces};
        
        idxtable = idxtable +1;
    
    end % subject loop
    
end % search space loop


% got to fROI directory & save info tables
cd(fullfile(dirs.froi, contrast));

writetable(infosubjlvl,['con_', contrast,'_sizes_subjectlevel_data.txt']);
writetable(infoactivevoxels,['con_', contrast,'_active_voxels.txt']);


disp('done.')
