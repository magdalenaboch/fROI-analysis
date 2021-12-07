
% ---------------------------------------------------------------------
% Script to extract activations levels from individual
% functional regions of interest (fROIs) using REX toolbox
%
% uses 'create_fROIs.m' to create fROIs
%
% Magdalena Boch, magdalena.boch@univie.ac.at
% ---------------------------------------------------------------------

clear all
clc

%% ======================================================
% This section needs your input

% subject IDs
% see also in first loop: subj = 1:numel(subjects) to add letters
subjects = [2,5:7,9,11:14,16,18,21,23,27,28];

% contrast images of test dataset to extract activation levels from
% in this case it is faces / bodies / objects > scrambled images (run 2)
condition = {'faces_scr_ss2', 'bodies_scr_ss2', 'objects_scr_ss2'};

% for name of table.txt containing results
suffix = '-scrambled'; % or -species;

% Which fROIs do you want to use?
% 1: category name (face_fROI or body_fROI)
catname = 'face_fROIs';
% 2: fROI contrast (how did you define the fROI?)
contrast = 'faces_objects_ss1';

% List of all percentages of active voxels you have from each subject
percentages = [1,5:5:100]; % to create enough lines for data table


% ---------------------------------------------------------------------
% Directories


dirs.home = ('YOURPATH');

dirs.project = fullfile(dirs.home, 'PROJECTNAME');

% directory containing the rex toolbox
dirs.rex = fullfile(dirs.home, 'tools_general', 'Rex');

% directory containing this script
dirs.scripts = fullfile(dirs.project,'scripts');

% Directory containing contrast images 
% (contains subfolders for each contrast containing contrast maps from all subjects)
dirs.glm = fullfile(dirs.project, 'data_glm', 'contrasts');

% directory to save data
dirs.analysis = fullfile(dirs.project, 'fROI_analysis');

% create subfolder for fROI type (face/body)
dirs.activation = fullfile(dirs.analysis,'activation_levels', catname);
if ~exist(dirs.activation, 'dir'); mkdir(dirs.activation); end

% directory containing individual fROIs
dirs.frois = fullfile(dirs.analysis,'individual_fROIs', catname, contrast);

% add folders to path (otherwise already included in start up file)
addpath(dirs.scripts); addpath(dirs.glm); addpath(dirs.frois); addpath(dirs.activation);


%% ===========================================================================
% Let's get started

% get search spaces from category folder (e.g., body_fROIs/bodies_objects_ss1/)
clear searchspaces; searchspaces = dir(dirs.frois);

% get rid of non-directories and non-files
searchspaces(ismember( {searchspaces.name}, {'.', '..'})) = [];

% get rid of files in list
% 1: What is a directory in this list?
dirFlags = [searchspaces.isdir];
% 2: Extract only those that are directories.
searchspaces = searchspaces(dirFlags);

% create empty table to save all activation levels
header = {'subjname','fROIcategory', 'fROIcontrast', 'searchspace', 'hemisphere', 'condition', 'percentage', 'activation'};
nlines = numel(subjects)*numel(searchspaces)*numel(condition)*numel(percentages); 
activationlvls = cell2table(cell(nlines,numel(header)),'VariableNames',header);
idxtable = 1;


% Rex needs two inputs:
% SOURCES: character array with one row being the path to an image to
% extract signal from = "condition" in this script

% ROIS:  character array withone row being the path to the respective ROI .nii file
% in this script this are the fROIs (with different percentages for each
% subj and search space)

% We have 4 loops
% 1: subjects
% 2: search spaces
% 3: conditions
% 4: outputdata X containts activation levels for each fROI
%    in this loop we write it into overall table activationlvls

disp(catname);

% ------------------------------------------------
% subject-level loop

for nsubj = 1:numel(subjects)
    
    if subjects(nsubj) < 10; subjname = ['F0',num2str(subjects(nsubj))];
    elseif subjects(nsubj) > 9; subjname = ['F',num2str(subjects(nsubj))];
    end
    
    disp(['subjname: ',subjname]);
    
    % ------------------------------------------------
    % search space loop ("ROIS")
    
    for idxs = 1:numel(searchspaces)
        
        clear sspace; sspace = searchspaces(idxs).name;
        clear fROIs; fROIs = spm_select('FPList', fullfile(dirs.frois,sspace,subjname), [subjname, '.*nii']);
        
        disp(['search space: ', sspace]);
        
        % ------------------------------------------------
        % condition loop ("SOURCES")
        
        for ncond = 1:numel(condition)
            clear cond; cond = condition{ncond};
            
            % path where condition contrast image is save)
            conpath = fullfile(dirs.glm, ['cons_', cond]);
            % get con image
            subj_con = dir(fullfile(conpath,[subjname,'_', cond, '*']));
            
            % add path
            sources = fullfile(conpath,subj_con.name);
            
            if strcmp(suffix, '-species')
                
                % short name for condition column in datatable (face, body, object...)
                clear underscore_indx; underscore_indx = strfind(cond,'_');
                cond_short = cond(1:(underscore_indx(2)-1));
                
                disp(['condition: ', cond_short]);
            else
                % short name for condition column in datatable (face, body, object...)
                clear underscore_indx; underscore_indx = strfind(cond,'_');
                cond_short = cond(1:(underscore_indx(1)-1));
                
                disp(['condition: ', cond_short]);
            end
            
            
            X = rex(sources,fROIs,'output_type','none', 'gui',0,'select_clusters',0);
            
            %{
  ---Help for rex---
  MEANS=rex(SOURCES, ROIS);
  where SOURCES is a list of M source volume files (image files to extract from)
  and ROIS is list of N roi files (image or .tal files)
  returns the mean values of each of the source volume files at the voxels
  identified by each ROI in the matrix MEANS (with size M x N).
 
  rex(SOURCES, ROIS, 'paramname1',paramvalue1,'paramname2',paramvalue2,...);
    permits the specification of additional parameters:
        'summary_measure' :     choice of summary measure (across voxels) [{'mean'},'eigenvariate','median','weighted mean','count']
        'level' :               summarize across [{'rois'},'clusters','peaks','voxels']
        'scaling' :             type of scaling (for timeseries extraction) [{'none'},'global','roi']
        'conjunction_mask':     filename of conjunction mask volume(s)
        'output_type' :         choice of saving output ['none',{'save'},'saverex']
        'gui' :                 starts the gui [{0},1]
        'select_clusters'       asks user to select one or more clusters if multiple clusters exists within each ROI [0,{1}]
        'dims' :                (for 'eigenvariate' summary measure only): number of eigenvariates to extract from the data
        'mindist' :             (for 'peak' level only): minimum distance (mm) between peaks
        'maxpeak' :             (for 'peak' level only): maximum number of peaks per cluster
            %}
 

            % ------------------------------------------------
            % activation levels - outputdata loop
            % (X = output containing mean activation for each fROI 
            
            for noutput = 1:size(fROIs,1) % could also be numel(X)
                
                 % we already have all info we need for each data column, 
                 % except for percentages, we extract info from filename
                 
                clear froipath; froipath = fROIs(noutput,:);    
                clear underscore_indices; underscore_indices = strfind(froipath,'_');
                clear perc_indices; perc_indices = strfind(froipath,'perc');
                clear percentage; percentage = str2double(froipath(underscore_indices(end)+1:perc_indices(end)-1));

                % write 'down' all information and activation levels
              
                activationlvls.fROIcontrast{idxtable,1} = contrast;
                activationlvls.searchspace{idxtable,1} = sspace(1:end-2); % rm hemisphere info
                activationlvls.hemisphere{idxtable,1} = sspace(end); 
                activationlvls.fROIcategory{idxtable,1} = catname;
                activationlvls.subjname{idxtable,1} = subjname;
                activationlvls.condition{idxtable,1} = cond_short;
                
                %this changes within this loop
                activationlvls.percentage{idxtable,1} = percentage;
                activationlvls.activation{idxtable,1} = X(:,noutput);
                % X = output from 
                                
                idxtable = idxtable + 1;   
            
            end % output data loop
        end % condition loop
    end % search space loop
end % subject loop

cd(fullfile(dirs.activation));

tablename = [catname, '_', contrast, '_activationlevels', suffix, '.txt'];
writetable(activationlvls,tablename);


disp('done.')