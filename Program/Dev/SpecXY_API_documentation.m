
%% SpecXY API
%% SpecDB
%% public functions
%% insert_data(SpecDB,data,update_set)
data.DB.x= X_signal; % [n_x,1]
data.DB.y=y_signals; % [n_x,n_y]

data.Meta  %any information for assignment to Meta and MetaPrivat for example:

data.Meta.Sample='Sample'; % Meta
data.Meta.SampleName='Sample_Spectra_01'; % Meta
data.Meta.Mineral='CPX'; % Meta
data.Meta.Thickness=100; % Meta - thickness in micron for no thickness normalisation set to 10000 micron
data.Meta.MetaTable; % % MetaPrivat MetaTable containing column with x, y coordinates, Group, Fileindex, Type etc.

% etc...
update_set=true; % updates the UI and table
insert_data(SpecDB,data,update_set)

%% [x, y, SampleName, Thickness, Id, eID] = get_specs(SpecDB,selections,mask_specs)
% selections= Unique Ids from data table e.g.: [1  12 51 210 ...]
% mask_specs (default=false)
% Id-unique Ids from data table  e.g.: [1  12 51 210 ...]
% eId position of entrie in app.data.Meta(eId), app.MetaPrivat(eId) and
% app.data.DB(eId)
% x, y, SampleName, Thickness are returned as cell arrays with {n_selections}

%% answer=CalcIntAbs_cm(SpecDB,settings_return)
% this funciton calculates and extracts values from spectra and is used as
% return function from the IntegrationSettings Module
% settings_return needs to contain the following information:

settings_return.CalcAllSpecs=true;
settings_return.xlimits=[xlimMin xlimMax];
settings_return.fn_str='CalculationResultColName';
settings_return.nlines=xpos; % x positions used for calculation
settings_return.CallibrationFactor=1; %
settings_return.ABC='B';
% 'A' - integration between 2 xpos ;
% 'B' - linear baseline xpos-> integration between 2 xpos;
% 'C' - linear baseline xlim-> integration between 2 xpos;
% 'D' - value of xpos;
% 'E' - linear baseline between xlimits -> Value of xpos;
% 'F' - linear baseline between xlimits -> normalisation [0 1] -> Value of xpos;


%%  MatchWithDB(SpecDB,x_signal,y_signal,n_match,mmv_input,mmv_DB,baseline_cor)
% this funciton is experimental and needs more work
% at the monent it uses:
% testfit=100-(mean((y_DB-y_unknown).^2));
% [~,idx]=sort(testfit,"descend");
% results_id=selections(idx(1:n_match));
% if you need this feature or want to improve this funciton please contact me.

%% update_public(SpecDB,what2update)
what2update=0;% update_UI
what2update=1;% update_DB_table
what2update=2;% update_DB_list

%% OH calcultation function see code for more information
% OH_table = CalculateH2O(SpecDB,x_signal,y_signal,eID,sample,int_range,individual_peak)
% OH_table=Calc_OH_from_selection(SpecDB,selection)

%% PlotSelectionFun(SpecDB,selections,axis2plot)
% selections= Unique Ids from data table e.g.: [1  12 51 210 ...]
%% Public SpecDB Propoerties
   Data % data structure including all data
   edit_state % default false, if change in project since last save true


%% SpecMaps 
sampleinfo= [selected_sample selected_subsample selected_map selected_submap];
%% ThicknessMapReturn(app,Thickness_map,sampleinfo)
%function to insert Thickness_map
sampleinfo= [selected_sample];

%% ClassificationReturn(app,Classification,sampleinfo)
sampleinfo= [selected_sample];
Classification.Classes % matrix size of map with integers 1:n_classes
Classification.ClassesNames % {'classname_01','classname_...','classname_n'}

%% MapMatchReturn(app,matrix_to_match,matrix_to_matchNames,sampleinfo,type)
% function to insert referenced and interpolated numerical matrices (maps)
% with same size and extend as selected_sample
% type=150% Quantdata
sampleinfo= [selected_sample selected_subsample selected_map];
matrix_to_match(n).MapPlotMatrix =matrix;
matrix_to_matchNames{n};

%% matrix=ClassModPlot(app,Axis2Plot,sampleinfo,thickness_state,plotimage,matrix,alpha) %ATH rename
% plots and returns matrix of sampleinfo
sampleinfo= [selected_sample selected_subsample selected_map selected_submap];
thickness_state=true; % thickness correction before plotting
% optional inputs:
% plotimage=true; % plot selected optical image if available
matrix % UI matrix instead selected_map
alpha % transparency of matrix

%% MatrixLimits=update_MatixLimits(app,sampleinfo);
sampleinfo= [selected_sample selected_subsample selected_map selected_submap];
% returns the matrix xy min max coordinates


%% [x_signal,y_signal,range] = get_specs(app,sampleinfo,xlimits)
% returns the raw x_signal and y_signal
sampleinfo= [selected_sample selected_subsample];
% optional inputs:
% xlimits=[min_x max_x]

%% [x_signal,y_signal,range] = get_specs_cor(app,sampleinfo,xlimits,signal,thickness_state)
% returns corrected x_signal and y_signal
sampleinfo= [selected_sample selected_subsample];
% optional inputs:
xlimits=[min_x max_x];
signal =[x_signal,y_signal];% optional input to correct input signal 
thickness_state=true; % thickness correction before return

%% [ColorMap] = PlotMap_UpdateColormap(app,Resolution,~)
% returns the colormap as RGB matrix - XMapTools code

%% [app, sampleinfo]=import_Map(app,SampleName,SubsampleName,DataIn,mapsize_raw,SpecType,Mask,Thickness_map,MetaTable,Classification,etc)
% insert map from external source into the SpecMaps structure. 
% see function code for details
% (Image import not supported yet) 

%% ResetTree(app,Tree2Update)
% this function resets a tree

%% update_class_tree(app,Tree2Update)
% this function updates a tree with the classifiction structure of the
% selected_sample;

%% UpdateTree(app,Tree2Update)
% this function updates the sample and map structure of SpecMaps into the tree input 


%% SpecMaps_plot_fun(app,hFig,settings)
% the development of this funciton is in the beginning and more features
% will be added

settings(n).Axis2Plot = Axis2Plot; % nexxtile or subplot(x,x,x);
settings(n).sampleinfo = [selected_sample selected_subsample selected_map selected_submap];
settings(n).Type = 'Map'; % 'Profile', 'optical image'

%optional inputs:
settings(n).Options.ROIs=true; % plot ROIs;
plt_settings.cb_label='text'; % colorbar label
settings(n).Options.Scalebar=true;
settings(n).Options.FontSize=18;
settings(n).Options.Title=true; %adds SubSampleName  uses title
settings(n).Options.Title='title string'; %uses 
% 'title string' as title

% not supported yet but will be added later:
% submaps 
% deconvolution results
% peaks 
% classmap 
% thicknessmap
% maskmap

%% Public SpecMaps Propoerties
   SpecMap % data structure including all data
   ColorMaps
   ColorMapValues
   ColorScale
