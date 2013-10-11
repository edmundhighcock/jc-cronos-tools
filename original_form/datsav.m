%latest version J.Citrin 20.4.2009
%see saveresult.m header comment... as and ev functions can be seen at 
%the bottom of the file
%DATSAV EXTRACTS CERTAIN PROFILES AND TIME-TRACES FROM THE CRONOS RESULT
%DATA STRUCTURE, AND SAVES THEM IN THE CALLER WORKSPACE (SAVERESULT
%WORKSPACE) WITH AN ADDED RUN SUFFIX (RUNNUM) AFTER EACH VARIABLE NAME

function datsav(data,param,datapath,runind);

global RUNNUM

%RUNNUM=input('Type in run number: ');

RUNNUM=runind;

%GENERAL 
as('tokamak','ITER');
as('M',2); %average mass of plasma in AMU for ASDEX and JET
as('S',680); %plasma surface area in m^2
as('Bout',4.3); %outer B in T
as('a',data.geo.a);
as('R',data.geo.r0); 
as('eps',ev('a')./ev('R')); 
as('k',data.geo.e1);
as('B',data.geo.b0); 
as('t',data.gene.temps);
as('r',param.gene.x); 
as('equirhoRZ',data.equi.rhoRZ);
as('equiR',data.equi.R);
as('equiZ',data.equi.Z);
as('equia',data.equi.a);
as('polflux',data.equi.psiRZ);
as('eta',data.neo.eta);
as('rhomax',data.equi.rhomax);
as('spr',data.equi.spr);
as('vpr',data.equi.vpr);
as('BZ',data.equi.BZ);
as('BR',data.equi.BR);

%calculate outboard R. NOTE THAT CSOU IS DEFINED FOR DEUTERIUM AND WITHOUT A SQRT(2)
[tmp,rhoend]=meshgrid(1:length(param.gene.x),data.equi.rhoRZ(:,end));
Rout=squeeze(data.equi.R(:,:,1));
rhonorm=data.equi.rhoRZ./rhoend;
for ii=1:length(data.gene.temps)
	%Routx(ii,:)=interp1(rhonorm(ii,:),Rout(ii,:),param.gene.x);
	
	%gammaE(ii,:)=-gradient(data.prof.vtor_exp(ii,:),Routx(ii,:)).*data.equi.a(ii,:)./data.prof.q(ii,:);
	gammaE(ii,:)=-gradient(data.prof.vtor_exp(ii,:),param.gene.x).*param.gene.x./data.prof.q(ii,:);
	csou(ii,:)=sqrt(data.prof.te(ii,:).*1.6e-19/(2*1.67e-27));
	
end
%as('Rx',Routx);
%as('gamE',gammaE./csou);

%0D parameters
as('beta',data.gene.beta);
as('betap',data.gene.betap);
as('tauth',data.gene.tauth);
as('tauei',data.gene.tauei);
as('tauj',data.gene.tauj);
as('wth',data.gene.wth);
as('wdia',data.gene.wdia);

as('Tim',data.gene.timoy/1000); %moy for average
as('Nem',data.gene.nemoy/1e19);
as('nbar',data.gene.nbar);
as('Tem',data.gene.temoy/1000);
as('zeffm',data.gene.zeffm); %average Zeff


%0D POWERS and CURRENTS

as('P_tot',data.gene.paddtot);  
as('P_loss',data.gene.ploss); %loss power (I think this means the plasma power 'exiting' via heat transport)  
as('P_fus',data.gene.paddfus);  
as('P_nbi',data.gene.paddidn);  
as('P_ic',data.gene.paddfci);  
as('P_ec',data.gene.paddfce);  
as('P_lh',data.gene.paddhyb);  
as('P_ohm',data.gene.paddohm);
as('dWdiadt',data.gene.dwdiadt); %diamagnetic energy derivative (total energy including suprathermal) 
as('P_aux',ev('P_tot')-ev('P_fus'));  
as('Prad',data.gene.prad);
as('Pbrem',data.gene.pbrem);
as('Pcyclo',data.gene.pcyclo);
as('Pradtot',ev('Prad')+ev('Pbrem')+ev('Pcyclo'));
as('Q',5*ev('P_fus')./(ev('P_tot')-ev('P_fus')));

as('ptotel',data.gene.pel);
as('ptotion',data.gene.pion);




as('Ip',data.gene.ip); %currents
as('Ip2',data.equi.ip);
as('Ip3',data.gene.ipjmoy);

as('Ioh',data.gene.ipohm);
as('Ilh',data.gene.ihyb); 
as('Ifce',data.gene.ifce);
as('Ibs',data.gene.iboot);
as('Ini',data.gene.ini); %non-inductive
as('Inbi',data.gene.iidn);

as('Ipmes',data.cons.ip); %other current related 'stuff'
as('li',data.equi.li);
as('li2',data.gene.li);

as('equi',data.equi);
as('pinis',data.cons.idn); %input NBI power
as('idncons',sum(data.cons.idn')'); %input NBI power
as('fcicons',data.cons.fci); %input ICRH power

as('Vmes',data.cons.vloop); %surface loop voltage measured
as('Vloop',data.gene.vres); %average loop voltage
as('fluxmes',data.cons.flux); %edge poloidal flux
as('vsurf',data.gene.vsurf);

%PROFILES

as('q',data.prof.q);
as('alpha',data.prof.alpha);
as('BPHI',data.equi.BPHI);
as('epar',data.prof.epar);
as('psid1',data.prof.psid1);
as('psid2',data.prof.psid2);
as('s',data.prof.shear);
as('sq',data.prof.shear./data.prof.q);
as('Te',data.prof.te/1000); %
as('Ti',data.prof.ti/1000); %
as('ne',data.prof.ne);
as('ni',data.prof.ni);
as('psi',data.prof.psi);
as('ptot',data.prof.ptot); %total pressure
as('pe',data.prof.pe); %electron pressure
as('pion',data.prof.pion); %ion pressure
as('impur',data.impur.impur); %impurity profile
as('ptot',data.prof.ptot); %total pressure
as('pth',data.prof.pe+data.prof.pion); %thermal pressure
as('F',data.equi.F); %diamagnetic function
as('equie',data.equi.e);
as('equitrh',data.equi.trh);
as('equitrl',data.equi.trl);
as('grho2r2',data.equi.grho2r2);
as('r2i',data.equi.r2i);
as('dpsidt',data.prof.dpsidt);
as('vtheta',double(squeeze(data.neo.vtheta(:,:,1))));
as('vtor',double(squeeze(data.neo.vtor(:,:,1))));
as('lti',data.prof.lti); %Ti gradient length
as('lte',data.prof.lte); %Te gradient length
as('lne',data.prof.lne); %ne gradient length
as('lni',data.prof.lni); %ni gradient length

as('vtorexp',data.prof.vtor_exp); %Toroidal velocity

as('nhe',data.impur.impur(:,:,3)); %helium density
as('zeff',data.impur.zeff); %zeff profile
as('xdur',data.prof.xdur); %suprathermal particle density


as('Pfce',data.source.fce.el); %source and sink deposition profiles
as('Plh',data.source.hyb.el);
as('Pnbi_i',data.source.idn.ion);
as('Pnbi_e',data.source.idn.el);
as('Pnbi_ne',data.source.idn.ne);
as('Pnbi_w',data.source.idn.w);
as('Pfci_i',data.source.fci.ion);
as('Pfci_e',data.source.fci.el);
as('Pfus_i',data.source.fus.ion);
as('Pfus_e',data.source.fus.el);
as('Pohm',data.source.ohm);
as('P_i',data.source.totale.ion);
as('P_e',data.source.totale.el);
as('q_ei',data.source.qei);  %sinks
as('Psyn',data.source.cyclo);
as('P_rad',data.source.prad);
as('P_brem',data.source.brem);

as('Psupra_nbi',data.source.idn.psupra);
as('Psupra_fci',data.source.fci.psupra);

if exist('data.source.idn.nfast');
	as('nfast_nbi',data.source.idn.nfast);
end

as('Pnbi',ev('Pnbi_i')+ev('Pnbi_e'));
as('Pfci',ev('Pfci_i')+ev('Pfci_e'));
as('Pfus',ev('Pfus_i')+ev('Pfus_e'));
as('Ptot',ev('Pnbi')+ev('Pfci')+ev('Pfus')+ev('Pohm')+ev('Pfce')+ev('Plh'));

as('jfce',data.source.fce.j); %current deposition profiles
as('jlh',data.source.hyb.j);
as('jnbi',data.source.idn.j);
as('jfci',data.source.fci.j);
as('j',data.prof.jmoy);
as('jbs',data.neo.jboot);
as('jphi',data.prof.jphi);
as('jpol',data.prof.jpol);
as('psi',data.prof.psi*2*pi); %for calculating flux consumption

jni=data.neo.jboot+data.source.fce.j+data.source.hyb.j+data.source.idn.j;

as('jni',jni);

%TRANSPORT COEFFICIENTS
%chie=chiecal; chii=chiical;
as('chie_int',data.prof.flux.kean./data.prof.ne); %effect el heat trans dif coeff
as('chii_int',data.prof.flux.kian./data.prof.ni);
as('chie_pre',data.coef.ee./data.prof.ne); 
as('chii_pre',data.coef.ii./data.prof.ni);
as('chie_neo',data.neo.coef.ee./data.prof.ne);  %neoclassical heat diffusitivity electrons
as('chii_neo',data.neo.coef.ii./data.prof.ni);  %neoclassical heat diffusitivity ions
as('qe',data.prof.flux.qe); %electron heat flux 
as('qi',data.prof.flux.qi); %ion heat flux 


B=ev('B'); B=B(1);
R=ev('R'); R=R(1);
as('rawqiGB',( ((ev('Ti').*1.6e-19*1e3).^2.5.*ev('ne').*sqrt(2.*1.67e-27) )./ ((1.6e-19*B*R).^2) ) ); %GB units

as('qiGB',data.prof.flux.qi./ ( ((ev('Ti').*1.6e-19*1e3).^2.5.*ev('ne').*sqrt(2.*1.67e-27) )./ ((1.6e-19*B*R).^2) ) ); %ion heat flux in GB units

as('qeiGB',data.source.qei./ ( ((ev('Ti').*1.6e-19*1e3).^2.5.*ev('ne').*sqrt(2.*1.67e-27) )./ ((1.6e-19*B*R).^2) ) ); %ion heat flux in GB units

%CALCULATED VALUES
as('betaN',100*ev('beta').*ev('a').*ev('B')./ev('Ip')*1e6); %
as('P_loss2_',ev('P_tot')-ev('Pbrem')-ev('Pcyclo')-1./3*ev('Prad'));

as('P_loss3_',ev('P_tot')-ev('Pbrem')-ev('Pcyclo')-1./3*ev('Prad')-ev('P_nbi')+ev('idncons')); %WITH INPUT NBI


as('f_b_s',ev('Ibs')./ev('Ip')); %
as('f_n_i',ev('Ini')./ev('Ip')); %
as('nG',(ev('Ip')./1e6)/pi./ev('a').^2); %%
as('f_G',(ev('nbar')/1e20)./ev('nG')); %

as('roman',4./3*(1+2.*ev('s')./ev('q')).*(1+ev('Ti')./ev('Te')));
as('romansq',4./3*(1+2.*ev('s')./ev('q')).*(1+ev('Ti')./ev('Ti')));

as('jenko',(4./3+1.91.*ev('s')./ev('q')).*(1+ev('zeff').*ev('Ti')./ev('Te')));

as('tauIPB',0.0562*(ev('Ip')/1e6).^0.93.*ev('B').^0.15.*(ev('P_loss')/1e6).^(-0.69).*(ev('nbar')/1e19).^0.41 ...
             .*ev('M')^0.19.*ev('R').^1.97.*ev('eps').^0.58.*ev('k').^0.78);

  as('tauIPB2_',0.0562*(ev('Ip')/1e6).^0.93.*ev('B').^0.15.*(ev('P_loss2_')/1e6).^(-0.69).*(ev('nbar')/1e19).^0.41 ...
             .*ev('M')^0.19.*ev('R').^1.97.*ev('eps').^0.58.*ev('k').^0.78);
       
     as('tauIPB3_',0.0562*(ev('Ip')/1e6).^0.93.*ev('B').^0.15.*(ev('P_loss3_')/1e6).^(-0.69).*(ev('nbar')/1e19).^0.41 ...
             .*ev('M')^0.19.*ev('R').^1.97.*ev('eps').^0.58.*ev('k').^0.78);

  as('tauth2_',(ev('wth')./(ev('P_tot')-ev('Pbrem')-ev('Pcyclo'))));
    as('tauth3_',(ev('wth')./(ev('P_tot')-ev('Pbrem')-ev('Pcyclo')-ev('P_nbi')+ev('idncons'))));



   P=ev('P_tot')-ev('P_nbi')+ev('idncons').*(1-exp(3.35-2./3.*ev('Ip')/1e6-0.2*ev('nbar')/1e19)/100)-0.472*1e6;
   P2=ev('P_tot')-ev('P_nbi')+ev('idncons').*(1-exp(3.35-2./3.*ev('Ip')/1e6-0.2*ev('nbar')/1e19)/100)-0.72*1e6;
   
   P3=ev('P_tot')-ev('P_nbi')+ev('idncons').*(1-exp(3.35-2./3.*ev('Ip')/1e6-0.2*ev('nbar')/1e19)/100)-ev('idncons')+ev('P_nbi');
   P4=ev('P_tot');
  
   
   as('tauth4_',(ev('wth')./P));
   as('tauIPB4_',0.0562*(ev('Ip')/1e6).^0.93.*ev('B').^0.15.*(P/1e6).^(-0.69).*(ev('nbar')/1e19).^0.41 ...
             .*ev('M')^0.19.*ev('R').^1.97.*ev('eps').^0.58.*ev('k').^0.78);

  as('tauth5_',(ev('wth')./P2));
   as('tauIPB5_',0.0562*(ev('Ip')/1e6).^0.93.*ev('B').^0.15.*(P2/1e6).^(-0.69).*(ev('nbar')/1e19).^0.41 ...
             .*ev('M')^0.19.*ev('R').^1.97.*ev('eps').^0.58.*ev('k').^0.78);

  as('tauth6_',(ev('wth')./P3));
   as('tauIPB6_',0.0562*(ev('Ip')/1e6).^0.93.*ev('B').^0.15.*(P3/1e6).^(-0.69).*(ev('nbar')/1e19).^0.41 ...
             .*ev('M')^0.19.*ev('R').^1.97.*ev('eps').^0.58.*ev('k').^0.78);
  
    as('tauth7_',(ev('wth')./P4));
   as('tauIPB7_',0.0562*(ev('Ip')/1e6).^0.93.*ev('B').^0.15.*(P4/1e6).^(-0.69).*(ev('nbar')/1e19).^0.41 ...
             .*ev('M')^0.19.*ev('R').^1.97.*ev('eps').^0.58.*ev('k').^0.78);
  
   
         
  as('H2_',ev('tauth2_')./ev('tauIPB2_')); %
  as('H3_',ev('tauth3_')./ev('tauIPB3_')); 
  
  as('H4_',ev('tauth4_')./ev('tauIPB4_')); 

   as('H5_',ev('tauth5_')./ev('tauIPB5_')); 

   as('H6_',ev('tauth6_')./ev('tauIPB6_'));  %FOR JET when using SPOT
 as('H7_',ev('tauth7_')./ev('tauIPB7_'));  %FOR ASDEX when using SPOT: no ion orbit loss taken into account
      
as('H',ev('tauth')./ev('tauIPB')); %
as('Plhtrans',2.76*ev('M').^(-1).*ev('B').^0.96.*(ev('nbar')/1e20).^0.77.*ev('R').^1.23.*ev('a').^0.76); %ITER physics basics
as('Plhtrans2_',0.072*ev('Bout')^0.7.*(ev('nbar')/1e20).^0.7*ev('S').^0.9.*(ev('zeffm')/2).^0.7.*2./ev('M'));  %ITPA H-mode power threshold database working group
%psi=ev('psi'); len=length(psi(:,1));endtime=length(psi(1,:));as('fluxcons',psi(len,endtime)-psi(len,1));
psi=ev('psi'); len=length(psi(1,:)); endtime=length(psi(:,1)); 

%as('totfluxcons',psi(endtime,len)-psi(1,len));
%as('rampfluxcons',psi(t100,len)-psi(1,len));
%as('flatfluxcons',psi(endtime,len)-psi(t100,len));
%as('thoufluxcons',psi(

fOUT=sprintf('run%ddata',RUNNUM);
clear data param
save([datapath fOUT]);

function as(str,val)
global RUNNUM
ws='caller';
assignin(ws,[str sprintf('%d',RUNNUM)],val);

function res=ev(str)
global RUNNUM
ws='caller';
res=evalin(ws,[str sprintf('%d',RUNNUM)]);


