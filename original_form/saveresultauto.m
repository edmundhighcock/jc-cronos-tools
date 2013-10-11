%latest version J.Citrin 23.7.2009

%FROM A CRONOS RESULT FILE CHOSEN BY THE USER, THIS FUNCTION AND ITS
%ACCOMPANYING FUNCTION DATSAV SAVES NUMEROUS IMPORTANT VARIABLES INTO A
%FILE RUN#DATA.MAT, WHEREAS # IS AN INTEGER (RUN INDEX) INPUTTED BY THE USER. IN
%RUN#DATA.MAT THE VARIABLE NAMES ALL INCLUDE A SUFFIX #

%THE RUN#DATA.MAT FILE IS SAVED IN A DIRECTORY GIVEN BY DATAPATH 

function saveresultauto

global RUNNUM
global DATATEMP


runtype=input('\nWhat are you working on? (0 for 77914, 1 for ITER hybrid, 2 for ramp-up, 3 for JET hybrid, 4 for TS qs-scan, 5 for ASDEX hybrid, 6 for newJET hybrid archive, 7 for JET_ILW): ');
if runtype==0, workname='JET/77914/'; end
if runtype==1, workname='ITER_hybrid/'; end
if runtype==2, workname='JET_rampup/';    end
if runtype==3, workname='JET7962/'; end
if runtype==4, workname='TS43191/'; end
if runtype==5, workname='ASDEXhybnew/'; end
if runtype==6, workname='JEThybnew/'; end
if runtype==7, workname='JET_ILW/'; end


if isempty(workname), error('Runtype not recognized'); end
    
%%%%
    datapath=[getenv('CRONOS_RUNS_FOLDER'),workname];
    fprintf(['\nDefault data directory is ',datapath]);
    runnums=input('\nChoose your favourite CRONOS runs for saving: ');
for j=1:length(runnums)
    %datapath=['~/IntegratedModelling/cronos/',workname];
    datapath=[getenv('CRONOS_RUNS_FOLDER'),workname];
    dirs=dir(datapath)
    testname=['run' num2str(runnums(j)) '_']
    %if length(testname) == 4
    %    testname=[testname '_'];
    %end
    len=length(testname);
    dirname='none';
	
    for i=1:length(dirs)
        compname=dirs(i).name;
        if strcmp(testname,compname(1:min(len,length(compname))))
            dirname=compname;
            dirname
            break
        end
    end
    if strcmp(dirname,'none')
        error(['Run ',num2str(runnums(j)),' not found!']);
    end

    %pathname=[datapath dirname '/']
    pathname=[datapath  '/' compname '/']
    filename=[compname,'_resultat.mat']
    %%%%

    %[filename,pathname]=uigetfile('*.gz','Choose your favourite CRONOS result file');
    %gunzip([pathname filename]);
    %newfilename=strrep(filename,'.gz','');
    %load(newfilename);
    %delete(newfilename);
	
	load([pathname filename]);	


    DATATEMP=data;

    structconv('DATATEMP');
    data=DATATEMP;

    %datapath=['~/IntegratedModelling/cronos_abbrev/',workname];
    datapath=[getenv('CRONOS_SAVEJETAUTO_FOLDER'),workname];
    fprintf(['\nDefault data save directory is ' datapath,'\n']);
    %datapath=input('Type in data directory (hit enter for default): ');
    %if isempty(datapath) datapath=datapathdef; end

    datsav(data,param,datapath,runnums(j))
    
end

function structconv(name)
%searches a structure for vectors which are in string form and converts them 
%to a numeric form via the eval function

global DATATEMP
nextlevel=fieldnames(eval(name));

for i=1:length(nextlevel)
    newname=[name '.' char(nextlevel(i))];
    if strcmp(class(eval(newname)),'struct')
        structconv(newname);
    else
        if strcmp(class(eval(newname)),'char')
            a=eval(newname);
            ws='caller';
            eval([newname,'=evalin(ws,a);']);    
        end
    end
end
