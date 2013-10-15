runtype=input('\nWhat are you working on? (0 for 77933, 1 for ITER hybrid, 2 for ramp-up, 3 for JET hybrid, 4 for TS qs-scan, 5 for ASDEX hybrid, 6 for newJET hybrid archive, 7 for JET_ILW): ');
if runtype==0, workname='JET/77933/'; end
if runtype==1, workname='ITER_hybrid/'; end
if runtype==2, workname='JET_rampup/';    end
if runtype==3, workname='JET7962/'; end
if runtype==4, workname='TS43191/'; end
if runtype==5, workname='ASDEXhybnew/'; end
if runtype==6, workname='JEThybnew/'; end
if runtype==7, workname='JET_ILW/'; end

if isempty(workname), error('Runtype not recognized'); end


num=input('\nWhich run shall I load? ');
fprintf('\nYour wish is my command...\n\n');

fIN=sprintf('run%ddata',num); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load 'runXdata.mat' file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirs=dir([getenv('CRONOS_RUNS_FOLDER'), workname]);
dirname = 'none';

%load(['~/IntegratedModelling/cronos_abbrev/' workname fIN]) ; 
%load([getenv('CRONOS_RUNS_FOLDER'), workname, fIN]) 

testname = ['run' num2str(num) '_'];
for i=1:length(dirs)
    if strncmp(testname, dirs(i).name, length(testname))
        dirname = dirs(i).name
        fprintf(['Now loading ', [getenv('CRONOS_RUNS_FOLDER'), workname, dirname, '/', fIN], '\n'])
        
        load([getenv('CRONOS_RUNS_FOLDER'), workname, dirname, '/', fIN])
    end
end

if (strcmp(dirname,'none'))
    error(['Run ', num2str(num), ' not found!']);
end

%run.number = num;
%run = setfield(run, sprintf('n%d', num), evalin('caller', sprintf('ldr%d', num)))

%flag=input('Do you want the raw data as well? (1 for yes)');
%if flag
%    fIN=sprintf('run%draw',num); 
%    load(['~/CRONOS_result_files/' fIN]) ; 
%    fprintf('\nThere you go...\n\n');
%else
%    fprintf('\nNever mind then...\n\n');
%end
