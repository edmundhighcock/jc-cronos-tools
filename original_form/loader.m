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


num=input('\nWhich run shall I load? ');
fprintf('\nYour wish is my command...\n\n');

fIN=sprintf('run%ddata',num); 

%load(['~/IntegratedModelling/cronos_abbrev/' workname fIN]) ; 
load([getenv('CRONOS_SAVEJETAUTO_FOLDER') workname fIN]) ; 

%flag=input('Do you want the raw data as well? (1 for yes)');
%if flag
%    fIN=sprintf('run%draw',num); 
%    load(['~/CRONOS_result_files/' fIN]) ; 
%    fprintf('\nThere you go...\n\n');
%else
%    fprintf('\nNever mind then...\n\n');
%end

%what is below is redundant JC 11.3.10

if num == 119
    otherq=q119(181,:);
    others=s119(181,:);
    otherpsi=psi119/(2*pi);
end

if num == 170
    otherq=q170(181,:);
    others=s170(181,:);
    otherpsi=psi170/(2*pi);
end
