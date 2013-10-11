%clear xth xsh teth tesh xth2 xsh2 teth2 tesh2
%clear tetemp titemp netemp vtortemp

%J.CITRIN 21.3.12
%THIS FUNCTION CARRIES OUT A STAND ALONE FIT OF JET DATA, BASED ON THE SPLINEFIT.M ITERATIVE WEIGHTED REGRESSION SCHEME WITH BREAKPOINTS AND A CONSTRAINT
%FOR VANISHING FIRST DERIVATIVES ON AXIS
%
%THIS VERSION EXPECTS TO HAVE A CRONOS DATA STRUCTURE IN THE WORKSPACE, SAY BASED FROM A NON-OPTIMIZED FIT BASED ON ZJET. THIS IS SINCE THE CXFM AND CXSE 
%MAPPINGS ARE CARRIED OUT HERE MANUALLY USING HELENA EQUILIBRIUM
%
%THIS VERSION EXPECTS TO FIND THE FOLLOWING DATA IN THE WORKSPACE: ZJET#DATA.MAT, WHEREAS # CORRESPONDS TO THE SUFFIX FROM THE SPECIFIC JET DISCHARGE
% 								   TIDATA.MAT, IN WHICH THE CXFM AND CXSE PPF DATA FOR THE VARIOUS SHOTS ARE STORED
%
%THE ZJET#DATA.MAT FILE CAN BE AUTOMATICALLY MADE BY THE M-FILE SAVEZJET ZJET


%IDEA FOR IMPROVEMENT: IF CX DATA CORRECTED BY CHEAP EXISTS, WE SHOULD READ THOSE DATAFILES FROM THE PPF. THIS IS TYPICALLY ONLY IMPORTANT FOR ROTATION


%MODIFICATIONS: 11/11/11 Allowance for gaps in the CX data made (copies last valid data to missing data timepoints) J.Citrin
%MODIFICATIONS: 21/3/12 Constraint for vanishing first derivates on axis added 
%MODIFICATIONS: 29/3/12 ECE x>0.8 filter is now a flagged option. In addition, any negative Te points are now taken care of (set equal to 1/2 the previous point in a loop)


%SUGGESTIONS FOR FUTURE IMPROVEMENT...
%
%TIME DEPENDENT INPUT FILE
%CHOOSE FROM WHICH TIME TO START THE FIT 
%SAVE CONFIGURATION OF THE FITS
%IMPROVE INCLUSION OF CXSE DATA WHEN EXISTS
%BYPASS ZJET 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%BELOW WE HAVE THE 'INPUT FILE' FOR THE FIT, i.e. THE VARIOUS BREAKPOINTS AND POLYNOMIAL ORDERS, DIAGNOSTICS TO BE USED IN THE FITS, AND VARIOUS
%ANNOYING LITTLE SETTINGS WHICH ARE SOMETIMES NEEDING FOR PARTICULAR RUNS (e.g. IF THE ECE DATA IS NOT VIABLE BELOW A CERTAIN TIME, WE CAN SET THIS TIME)

twait=41; %jump to this time 
tearly=45; %sets last timepoint of early time

cxb=[0 0.8 1]; %breakpoint for cx fit spline 
cxtin=7; %spline polynomial order for the Ti fit spline
cxvtorn=7; %spline polynomial order for the vtor fit spline
cxfilfac=1.5; %filtering parameter: ratio above average of first 3 points, above which points are ignored
cxaxfac=1.1; %ratio above first data point, at which the ti value at x=0 is set
cxedgefac=2; %ratio below the last 4 CXSE edge points (or last CXFM edge point), at which the x=1 ti value is set

cxnp=0; %3; %filtering parameter: the last cxnp+1 CXSE data points are filtered out in the fit (due to high resolution)
rmlowxti=1; %for removing the points in the CX data corresponding to the high-field size of the flux surfaces, which are not always consistent with the low field side
tetirat=1.2; %assumed ratio for Te/Ti if we are outside the range of CX measurements (1.2 roughly corresponds to an Ohmic phase)

teearly=200; %in ev, early boundary
teb=[0 0.75 1]; %breakpoint for te fit spline 
ten=8; %spline polynomial order for the te fit spline
tefilfac=3; %filtering parameter: ratio above average of first 5 te points, above which points are ignored
teaxfac=1.1; %ratio above first data point, at which the te value at x=0 is set
teedgefac=1.3; %ratio below the last 4 te edge points, at which the x=1 te value is set
tthtemponly=40.5; %time under which only the HRTS data is used
tshonly=1; %set 1 for ECE only 
tthonly=0; %set 1 for HRTS only
ecex08filter=0; %if set 1, then all ECE points beyond x=0.8 are ignored (useful when combining ECE and HRTS data)

neearly=1e18;
neb=[0 0.65 1]; %breakpoint for ne fit spline 
nen=5; %spline polynomial order for the ne fit spline
nefilfac=1.5; %filtering parameter: ratio above average of first 5 ne points, above which points are ignored
neaxfac=1.02; neaxfac=1.02; %ratio above first few data points, at which the ti value at x=0 is set
needgefac=5; %ratio below the last 4 ne edge points, at which the x=1 ne value is set
tLH=[55 57]; %period of LH transition typically with hollow profiles
nebLH=[0 0.6 0.8 0.9  1];
nenLH=4;

newti=0; %not all runs have CXSE data. The way we work is first fit all channels, and then if there is CXSE data, we rerun the script with newti=1
docxse=0; %if we want to include CXSE data in the actual fit of the second iteration, then set to 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

con.xc=0; con.yc=0; con.cc=zeros(2,1); con.cc(2)=1; %SETS THE DERIVATIVE CONSTRAINT FOR EACH PROFILE. SHOULD NOT BE CHANGED

doit=input('Enter 1 to input fits from previous JETnewfit run into CRONOS data file: '); %if satisfied with the fit, then doit
if doit==1
    data.prof.ti=titemp;
    data.prof.vtor_exp=vtortemp;
    param.gene.ti_invar=1; %should always be 1, this just ensures it
    data.mode.ae(:)=2; data.mode.brem(:)=2; data.mode.prad(:)=2;	 %they are not calculated by default for some reason. This ensures it
    if newti~=1
	    data.prof.te=tetemp;
	    data.prof.pe=tetemp.*netemp*1.6e-19;
	    data.prof.ne=netemp;
    end
    error('Not an error, just exiting');
end



if newti==1 %should carried out in a second iteration if CXSE data exists

    load('Tidata.mat'); %CXSE and CXFM data. JETnewfit expects to find the data for the specific discharge you're fitting with a 2 digit discharge number suffix

    titemp=data.prof.ti;
    vtortemp=data.prof.vtor_exp;
    runnum=input('Which run suffix? '); %needed in order to load the correct zjet data, needed for jetprof
    cxtime=input('When does the beam turn [on off]? '); %set time window of ti measurements
    
    
    eval(['load(',char(39),'zjet',num2str(runnum),'data.mat',char(39),');']); %load the jetprof data
    eval(['jetprof=jetprof',num2str(runnum),';']);
    
  %  cxtime=[jetprof.tcx(1)+0.1 jetprof.tcx(length(jetprof.tcx))-0.1];
    
    tte=eval(['tte',num2str(runnum)]);  %times for edge data from CXSE/TI   
    tie=eval(['tie',num2str(runnum)]);  %edge temperatures from CXSE/TI
    vtore=eval(['vtore',num2str(runnum)]); %vtor for edge data from CXSE/ANGF. Has same times as TI

    rc=eval(['r',num2str(runnum),'c']);    %radii from CXFM/RCOR
    rec=eval(['re',num2str(runnum),'c']);  %processed (using EFIT) edge radii on midplane from CXSE/RT
    re=eval(['re',num2str(runnum)]);  %raw edge radii from CXSE/TI
 
    iibeg=round(interp1(data.gene.temps,1:length(data.gene.temps),cxtime(1))); %set relevant time indexes from cronos data
    iiend=round(interp1(data.gene.temps,1:length(data.gene.temps),cxtime(2))); %set relevant time indexes from cronos data        
  
    figure;
    for ij=iibeg:iiend
        
	t=data.gene.temps(ij);
	i=round(interp1(jetprof.tcx,1:length(jetprof.tcx),t)); %find relevant time index in raw core data
	ie=round(interp1(tte,1:length(tte),t)); %find relevant index in raw edge data
	isim=round(interp1(data.gene.temps,1:length(data.gene.temps),t)); %find relevant index in cronos data
    
	
	
	if sum(jetprof.tirhoncx(i,:))==0 %in order to catch gaps in the CX measurements
		ij
	    	titemp(ij,:)=lastti;
		vtortemp(ij,:)=lastvtor;
	    	subplot(1,2,1); plot(xp,titemp(ij,:),'r');
        	title(['GAP IN CX, assumed ti \newline from last measured timepoint']);
	        subplot(1,2,2); plot(xp,vtortemp(ij,:),'r');
        	title(['GAP IN CX, assumed vtor \newline from last measured timepoint']);
		input(['Now ij=',num2str(ij),', hit enter for next']);
        	clf
    	else
	        ticx=jetprof.tirhoncx(i,:); %core ti from cxfm
	        rotcx=jetprof.rotrhoncx(i,:); %core rotation from cxfm

	        xcx=interp1(data.equi.R(isim,:,1),data.equi.rhoRZ(isim,:)./data.equi.rhoRZ(isim,101),rc(:,i)); %interpolation of core cx radial data to toroidal flux coordinate grid
    
	        ticxe=tie(:,ie); %rename for convenience edge ti and vtor data
	        rotcxe=vtore(:,ie);
        
	        xcxe=interp1(data.equi.R(isim,:,1),data.equi.rhoRZ(isim,:)./data.equi.rhoRZ(isim,101),rec(:,ie)); %interpolation of edge cx radial data to toroidal flux coordinate grid
  
	        xcx=double(xcx'); xcxe=double(xcxe'); ticxe=double(ticxe'); rotcxe=double(rotcxe'); %reframe the data into doubles and arrays of the right size
    
    		if docxse==1 
		        xcx2=[xcx xcxe]; ticx2=[ticx ticxe]; %combine the core and edge data into a new name for convenience
		        rotcx2=[rotcx rotcxe]; 
			xrot=[xcx xcxe];
                else
		 	xcx2=[xcx]; ticx2=[ticx]; %combine the core and edge data into a new name for convenience
		        rotcx2=[rotcx]; 
			xrot=[xcx];
                end
		
				  
	        %filter out np+1 last data points in the combined set      
	        xcx2(length(xcx2)-cxnp:length(xcx2))=[];
	        ticx2(length(ticx2)-cxnp:length(ticx2))=[];

		xrot(length(xrot)-cxnp:length(xrot))=[];
	        rotcx2(length(rotcx2)-cxnp:length(rotcx2))=[];

	      	%filter out NANs
		nanscx=isnan(xcx2);
	        nansrotcx=isnan(xrot);
             
	        xcx2(nanscx)=[];
	        ticx2(nanscx)=[];
	
	        xrot(nansrotcx)=[];
        	rotcx2(nansrotcx)=[];
	
	        %filter out spurious zeros
	        xcx2(ticx2<0.1)=[];
	        ticx2(ticx2<0.1)=[];

	        xrot(rotcx2<0.1)=[];
	        rotcx2(rotcx2<0.1)=[];

        
	        %filter out spurious high points
	        toohighcx=mean(ticx2(1:3))*cxfilfac;
        	toohighrotcx=mean(rotcx2(1:3))*cxfilfac;

        
	        xcx2(ticx2>toohighcx)=[]; 
	        ticx2(ticx2>toohighcx)=[]; 
    
        	xrot(rotcx2>toohighrotcx)=[]; 
	        rotcx2(rotcx2>toohighrotcx)=[]; 
        
        
	        %filter out spurious low points on axis (especially important for ECE)
        	ticx2(xcx2<0.001)=[];
	        xcx2(xcx2<0.001)=[];
        
	        rotcx2(xrot<0.001)=[];
        	xrot(xrot<0.001)=[];
        
	        %arrange CX data by removing points corresponding to the high field side of the flux surfaces, which are not always consistent with the low field side, leading
		%to multiple values at the same x-coorerrordinate. This point can perhaps be refined.
        	ind=0;
	        i=0;
        	while ind==0
	            i=i+1;
        	    if xcx2(i+1)>xcx2(i)
                	ind=1;
	            end
	        end
        
        	if rmlowxti == 1%actually remove the points found above, if desired
	           xcx2(1:i)=[]; ticx2(1:i)=[];
        	   xrot(1:i)=[]; rotcx2(1:i)=[];  
	        end
       
        	%set value on axis
	        ticxaxis=mean(ticx2(1))*cxaxfac;
	        rotcxaxis=mean(rotcx2(1))*cxaxfac;

        
        	%place points on 0 and 1 boundary of x-axis   
	        lencx=length(xcx2); lenrotcx=length(xrot);
        	minticx=mean(ticx2(lencx-4:lencx))/cxedgefac; 
	       % minrotcx=mean(rotcx2(lenrotcx-4:lenrotcx))*1.1; 
	       minrotcx=0; %simple boundary condition for rotation: 0-rotation at edge
        
	        xcx2=[0 xcx2 1]; 
	        xrot=[0 xrot 1];
    
	        ticx2=[ticxaxis ticx2 minticx];
	        rotcx2=[rotcxaxis rotcx2 minrotcx];
             
    
	        xp=0:0.01:1;  %make toroidal flux grid for fit

	        breaks=[0 cxb 1]; %set breakpoints
	        ppti = splinefit(xcx2,ticx2,breaks,cxtin,con); %fit ti
	        pprot = splinefit(xrot,rotcx2,breaks,cxvtorn,con); %fit rotation

	        ypti=ppval(ppti,xp);    %get fit on toroidal flux grid
	        yprot=ppval(pprot,xp);  
        
	        %plot fit
		subplot(1,2,1); plot(xcx,ticx,'.',xcx2,ticx2,'*',xcxe,ticxe,'*',xp,ypti);
	        if docxse == 1
			title(['Improved Ti for t=',num2str(t)]);
		else
			title(['Improved Ti for t=',num2str(t),' (fit without CXSE)']);
		end
        
	        subplot(1,2,2); plot(xcx,rotcx,'.',xrot,rotcx2,'*',xcxe,rotcxe,'*',xp,yprot);
	        if docxse == 1
			title(['Improved vtor for t=',num2str(t)]);
		else
			title(['Improved vtor for t=',num2str(t),' (fit without CXSE)']);
		end
        
	        titemp(ij,:)=ypti; vtortemp(ij,:)=yprot;
	        lastti=ypti; lastvtor=yprot; %in case there are gaps in the measurement
	
	        input(['Now ij=',num2str(ij),', hit enter for next']);
        	clf

    	end
    end
    
  nothing=error('Not an error, just moving along'); %exit fitting program here
  
end


clear tetemp titemp netemp vtortemp

runnum=input('Which run suffix? ');
cxtime=input('When does the beam turn [on off]? ');

eval(['load(',char(39),'zjet',num2str(runnum),'data.mat',char(39),');']); %load jetprof from zjet#data.mat file
eval(['jetprof=jetprof',num2str(runnum),';']);
%cxtime=[jetprof.tcx(1)+0.1 jetprof.tcx(length(jetprof.tcx))-0.1];
t=41;
if t>twait 
	figure;
end

%for k=168:201 %loop over CX times
for k=1:length(data.gene.temps) %loop over all times
       
    %FIT FOR TE   
    t=data.gene.temps(k);
	
    i=round(interp1(jetprof.tth,1:length(jetprof.tth),t));
    
    if tthonly==0
	    j=round(interp1(jetprof.tteshf,1:length(jetprof.tteshf),t));
    end
    
    teth=jetprof.terhonth(i,:); %HRTS data;
    xth=jetprof.rhonth(i,:);
    
    if tthonly==0
	    tesh=jetprof.terhonshf(j,:); %ECE data;
	    xsh=jetprof.rhonshf(j,:);
       
   	 %arrange ECE data (remove second set of higher points, and fix vector sizes
	    if length(xsh)<length(tesh)
	        tesh=tesh(1:length(xsh));
	    end
   
   	 ind=0;
	    i=0;
	    while ind==0
	        i=i+1;
	        if xsh(i+1)>xsh(i)
	               ind=1;
	        end
	    end
   	 xsh=xsh(i:length(xsh));
	 tesh=tesh(i:length(tesh));
	 xsh2=xsh; tesh2=tesh; 
	nanssh=isnan(xsh2); 
        xsh2(nanssh)=[];
        tesh2(nanssh)=[];
	 xsh2(tesh2<0.1)=[];
        tesh2(tesh2<0.1)=[];
        toohighsh=mean(tesh2(1:5))*tefilfac;
	xsh2(tesh2>toohighsh)=[]; 
        tesh2(tesh2>toohighsh)=[]; 
	 tesh2(xsh2<0.001)=[];
        xsh2(xsh2<0.001)=[];
	 teshaxis=mean(tesh2(1:5))*teaxfac;
	  
	if ecex08filter == 1 %filtering out edge points if desired
	  tesh2(xsh2>0.8)=[];
   	  xsh2(xsh2>0.8)=[];
	end  
   	
	lensh=length(xsh2);
	mintesh=mean(tesh2(lensh-4:lensh))/teedgefac;
	
	if t>tearly
		xsh2=[0 xsh2 1];
		tesh2=[teshaxis tesh2 mintesh];
	else
		xsh2=[0 xsh2(1:end-1) linspace(xsh2(end),1,5)];		
		tesh2=[teshaxis tesh2(1:end-1) linspace(tesh2(end),teearly,5)];
	
	end
	
    end
    xth2=xth; teth2=teth; 


    %filter out NANs
    nansth=isnan(xth2);
    xth2(nansth)=[];
    teth2(nansth)=[];
    
 
    %filter out spurious zeros
    xth2(teth2<0.1)=[];
    teth2(teth2<0.1)=[];
  
    %filter out spurious high points
    toohighth=mean(teth2(1:5))*tefilfac;
   

    xth2(teth2>toohighth)=[]; 
    teth2(teth2>toohighth)=[]; 
    
    
    %filter out spurious low points on axis (especially for ECE)
    teth2(xth2<0.001)=[];
    xth2(xth2<0.001)=[]; 
    
    tethaxis=mean(teth2(1:5))*teaxfac;
       
    %place points on 0 and 1 boundary of x-axis   
    lenth=length(xth2); 
    minteth=mean(teth2(lenth-4:lenth))/teedgefac; 

     if t> tearly
	    xth2=[0 xth2 1]; 
	    teth2=[tethaxis teth2 minteth];  %add boundary points
	else
		xth2=[0 xth2(1:end-1) linspace(xth2(end),1,5)];
		teth2=[tethaxis teth2(1:end-1) linspace(teth2(end),teearly,5)];
	end
    
    
    
    if xth2(length(xth2)-1) < 0.85 && tthonly==1 %sets breakpoint lower, to 0.7, if there are only 2 datapoints above 0.85 (could be the case for LIDAR)
        teb=0.7;
    end
        
    xp=0:0.01:1;  

    breaks=[0 teb 1];

	
   
    if t<tthtemponly || tthonly==1
        pp = splinefit(xth2,teth2,breaks,ten,con); %redo fit with Thompson data only, if in times defined in 'input file'
    else
        pp = splinefit([xth2 xsh2],[teth2 tesh2],breaks,ten,con);
        if tshonly
           pp = splinefit(xsh2,tesh2,breaks,ten,con);
        end
    end

    yp=ppval(pp,xp); 
	
	for ll=1:101 %cleans up negative points if they exist
		if yp(ll)<0
			yp(ll)=yp(ll-1)/2;
		 end
	end
	
   if t>twait
	    subplot(2,2,1);
    
    
	    if t<tthtemponly || tthonly==1
	        plot(xth,teth,'.',xp,yp,'r');
	    else
	        plot(xth,teth,'.',xsh,tesh,'*',xp,yp,'r');
	        if tshonly
	            plot(xsh,tesh,'*',xp,yp,'r');
	        end
	    end  
    
	    title(['te for t=',num2str(t)]);
   end
    
    %figure; plot(xth2,teth2,'.',xsh2,tesh2,'*',xp,yp);
    tetemp(k,:)=yp;

  
%CARRY OUT FIT FOR NE
    
    i=round(interp1(jetprof.tth,1:length(jetprof.neth),t));

    neth=jetprof.nerhonth(i,:);
    xth=jetprof.rhonth(i,:);
    
    xth2=xth; neth2=neth; 

    %filter out NANs
    nansth=isnan(xth2);

    xth2(nansth)=[];
    neth2(nansth)=[];
      
    %filter out spurious zeros
    xth2(neth2<0.1)=[];
    neth2(neth2<0.1)=[];

    %filter out spurious high points
    toohighth=mean(neth2(1:5))*nefilfac;

    xth2(neth2>toohighth)=[]; 
    neth2(neth2>toohighth)=[]; 
    
    %filter out spurious low points on axis (especially for ECE)
    neth2(xth2<0.001)=[];
    xth2(xth2<0.001)=[];
          
    nethaxis=mean(neth2(1:5))*neaxfac;
        
    %place points on 0 and 1 boundary of x-axis   
    lenth=length(xth2); 
    minneth=mean(neth2(lenth-4:lenth))/needgefac; 
	if t> tearly
	    xth2=[0 xth2 1]; 
	    neth2=[nethaxis neth2 minneth];
	else		
	    xth2=[0 xth2(1:end-1) linspace(xth2(end),1,5)];
	    neth2=[nethaxis neth2(1:end-1) linspace(neth2(end),neearly,5)];
	end

    xp=0:0.01:1;  

    if t > tLH(1) && t < tLH(2)
	    breaks=[0 nebLH 1];
	    pp = splinefit(xth2,neth2,breaks,nenLH,con);
    else
	    breaks=[0 neb 1];
	    pp = splinefit(xth2,neth2,breaks,nen,con);
    end
    
    yp=ppval(pp,xp); 
    if t>twait
	  subplot(2,2,2); plot(xth,neth,'.',xp,yp,'r');
	  title(['ne for t=',num2str(t)]);
    end
    netemp(k,:)=yp;

   
    %CARRY OUT TI FIT IF WE ARE IN THE TIME CORRESPONDING TO TI DATA MEASUREMENTS
    if t > cxtime(1) && t < cxtime(2)  
                  
        i=round(interp1(jetprof.tcx,1:length(jetprof.tcx),t));
	
	if sum(jetprof.tirhoncx(i,:))==0
	    	titemp(k,:)=lastti;
		vtortemp(k,:)=lastvtor;
		if t>twait
		    	subplot(2,2,3); plot(xp,titemp(k,:),'r');
	        	title(['GAP IN CX, assumed ti \newline from last measured timepoint']);
		        subplot(2,2,4); plot(xp,vtortemp(k,:),'r');
	        	title(['GAP IN CX, assumed vtor \newline from last measured timepoint']);
		end
    	else

	      
	      
	      ticx=jetprof.tirhoncx(i,:);
	      
	      %BEGINNING OF CHANGE
	      xcx=jetprof.rhoncx(i,:);
    	      %warning('NO CXSE but doing rcore interpolation');		
	      %isim=round(interp1(data.gene.temps,1:length(data.gene.temps),t)); %find relevant index in cronos data
   	      %xcx=interp1(data.equi.R(isim,:,1),data.equi.rhoRZ(isim,:)./data.equi.rhoRZ(isim,101),jetprof.rcxcor(i,:)); %interpolation of core cx radial data to toroidal flux coordinate grid
	      %xcx=double(xcx);
     	      %END OF CHANGE
	      
	      xcx2=xcx; ticx2=ticx; 

      	      %filter out NANs
       	      nanscx=isnan(xcx2);

              xcx2(nanscx)=[];
      	      ticx2(nanscx)=[];
      
      	      %filter out spurious zeros
	      xcx2(ticx2<0.1)=[];
	      ticx2(ticx2<0.1)=[];

	      %filter out spurious high points
	      toohighcx=mean(ticx2(1:3))*cxfilfac;

	      xcx2(ticx2>toohighcx)=[]; 
	      ticx2(ticx2>toohighcx)=[]; 
    
	      %filter out spurious low points on axis (especially for ECE)
	      ticx2(xcx2<0.001)=[];
	      xcx2(xcx2<0.001)=[];
        
	      %arrange CX data (remove second set of higher points)
	      ind=0;
     	      i=0;
       	      while ind==0
	       	 i=i+1;
               	 if xcx2(i+1)>xcx2(i)
	       	     ind=1;
	         end
	      end
      
              if rmlowxti == 1 %actually remove the points found above, if desired
           	 xcx2(1:i)=[]; ticx2(1:i)=[];
              end
  
   	      ticxaxis=mean(ticx2(1))*cxaxfac;
        
	      %place points on 0 and 1 boundary of x-axis   
	      lencx=length(xcx2); 
	      minticx=mean(ticx2(lencx))/cxedgefac; 

	      xcx2=[0 xcx2 1]; 
    
	      ticx2=[ticxaxis ticx2 minticx];
	      xp=0:0.01:1;  

	      breaks=[0 cxb 1];
	      pp = splinefit(xcx2,ticx2,breaks,cxtin,con);


	      yp=ppval(pp,xp); 

	      subplot(2,2,3); plot(xcx,ticx,'.',xcx2,ticx2,'*',xp,yp,'r');
	      title(['ti for t=',num2str(t)]);
	      titemp(k,:)=yp;
	      lastti=yp; %in case there are gaps in the measurement

      
	      %ROTATION FITTING
	      i=round(interp1(jetprof.tcx,1:length(jetprof.tcx),t));

   	      rotcx=jetprof.rotrhoncx(i,:);
  	      xcx=jetprof.rhoncx(i,:);
    
   	      xcx2=xcx; rotcx2=rotcx; 

  	      %filter out NANs
   	      nanscx=isnan(xcx2);

	      xcx2(nanscx)=[];
    	      rotcx2(nanscx)=[];
      
    	      %filter out spurious zeros
    	      xcx2(rotcx2<0.1)=[];
    	      rotcx2(rotcx2<0.1)=[];

    	      %filter out spurious high points
    	      toohighcx=mean(rotcx2(1:3))*cxfilfac;

    	      xcx2(rotcx2>toohighcx)=[]; 
    	      rotcx2(rotcx2>toohighcx)=[]; 
    
    	      %filter out spurious low points on axis (especially for ECE)
              rotcx2(xcx2<0.001)=[];
              xcx2(xcx2<0.001)=[];
        
    	      %arrange CX data (remove second set of higher points)
    	      ind=0;
    	      i=0;
     	      while ind==0
       	          i=i+1;
       	          if xcx2(i+1)>xcx2(i)
                     ind=1;
	          end
              end
     
      	      if rmlowxti == 1 %actually remove the points found above, if desired
        	   xcx2(1:i)=[]; rotcx2(1:i)=[];     
	      end
  
              rotcxaxis=mean(rotcx2(1))*cxaxfac;
        
     	      %place points on 0 and 1 boundary of x-axis   
     	      lencx=length(xcx2); 
     	      minrotcx=mean(rotcx2(lencx))/cxedgefac; 

     	      xcx2=[0 xcx2 1]; 
    
     	      rotcx2=[rotcxaxis rotcx2 minrotcx];
     	      xp=0:0.01:1;  

     	      breaks=[0 cxb 1];
	      
	      if runnum == 12 && t>46.5 && t<48.2
		  pp = splinefit(xcx2,rotcx2,[0 0.6 1],5,cxvtorn,con);
	      else
	         pp = splinefit(xcx2,rotcx2,breaks,cxvtorn,con);
	      end
    

              yp=ppval(pp,xp); 

              subplot(2,2,4); plot(xcx,rotcx,'*',xcx2,rotcx2,'.',xp,yp,'r');
              title(['vtor for t=',num2str(t)]);
        
              vtortemp(k,:)=yp;
	      lastvtor=yp; %in case there are gaps in the measurement
	
    	      lastk=k; %(see where it is used below for an explanation)

        end
    else %SET ASSUMPTIONS FOR TI AND VTOR IF WE ARE OUTSIDE THE MEASUREMENT TIME
        titemp(k,:)=tetemp(k,:)/tetirat;
        vtortemp(k,:)=zeros(1,101); %simple assumption for zero

        if t>cxtime(2)
            titemp(k,:)=titemp(lastk,:);       %IF WE ARE ABOVE THE MEASUREMENT TIME, THEN THE PROFILES ARE SIMPLY COPIED FROM THE LAST DATA TIME SLICE OF MEASUREMENTS
            vtortemp(k,:)=vtortemp(lastk,:);
        end

	if t>twait
	        subplot(2,2,3); plot(xp,titemp(k,:),'r');
	        title(['NO CX, assumed ti for t=',num2str(t)]);
	        subplot(2,2,4); plot(xp,vtortemp(k,:),'r');
	        title(['NO CX, assumed vtor for t=',num2str(t)]);
	end
    end
    if t>twait
	    nothing=input('Hit enter for next');
	    clf;
    end

end





