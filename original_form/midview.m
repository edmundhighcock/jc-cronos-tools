%latest version: J.Citrin 11.6.09

%takes a look at the intermediate results of a run in progress at a given
%timestep.. 

%total currents (i.e. I_p, I_boot, various I_non-inductive)
%total powers (P_alpha, P_rad, various P_non_inductive)
%dimensionless quantities (<Z_eff>, l_i, beta_N, etc etc)
%miscellaneous time quantities (Q, H, V_loop)

%PROFILES (within a spatial range defined by user inputted xbord variable, at time snapshots defined by user inputed tsnaps vector)
%q-profile (also finds and displays time of q=1 surface appearance)
%T_e, T_i
%non-inductive current profiles
%total non-inductive current profile, ohmic current profile, total current profile
%power deposition profiles
%anomalous transport coefficients chi_i and chi_e

%function midview

runnum=input('\nRun number please: ');
indnum=input('\nIndex number please: ');

datapath=[getenv('CRONOS_TEMP')];

unix(['copyrun.sh ' num2str(runnum) ' ' num2str(indnum)]);

dirs=dir(datapath);
testname=['run' num2str(runnum)];
len=length(testname);
dirname='none';
for i=1:length(dirs)
    compname=dirs(i).name;
    if strcmp(testname,compname(1:min(len,length(compname))))
        dirname=compname;
        break
    end
end
if strcmp(dirname,'none')
    error(['Run ',num2str(runnum),' not found!']);
end
    
%filename=[compname,'_resultat_',num2str(indnum),'.mat'];
filename=[compname];
load([datapath filename]);

M=2.5; %average mass of plasma in AMU
S=680; %plasma surface area in m^2
Bout=4.3; %outer B in T

a=datakp1.geo.a;
R=datakp1.geo.r0; 
eps=a./R; 
k=datakp1.geo.e1;
B=datakp1.geo.b0; 
t=datakp1.gene.temps;
r=param.gene.x; 

%0D parameters
beta=datakp1.gene.beta;
tauth=datakp1.gene.tauth;
tauei=datakp1.gene.tauei;
tauj=datakp1.gene.tauj;
wth=datakp1.gene.wth;

Tim=datakp1.gene.timoy/1000; %moy for average
Nem=datakp1.gene.nemoy/1e19;
nbar=datakp1.gene.nbar;
Tem=datakp1.gene.temoy/1000;
Zeff=datakp1.gene.zeffm; %average Zeff

%0D POWERS and CURRENTS

P_tot=datakp1.gene.paddtot;  
P_loss=datakp1.gene.ploss; %loss power (I think this means the plasma power 'exiting' via heat transport)  
P_fus=datakp1.gene.paddfus;  
P_nbi=datakp1.gene.paddidn;  
P_ic=datakp1.gene.paddfci;  
P_ec=datakp1.gene.paddfce;  
P_lh=datakp1.gene.paddhyb;  
dWdiadt=datakp1.gene.dwdiadt; %diamagnetic energy derivative (total energy including suprathermal) 
P_aux=P_tot-P_fus;  
Prad=datakp1.gene.prad;
Pbrem=datakp1.gene.pbrem;
Pcyclo=datakp1.gene.pcyclo;
Pradtot=Prad+Pbrem+Pcyclo;
Q=5*P_fus./(P_tot-P_fus);

Ip=datakp1.gene.ip; %currents
Ioh=datakp1.gene.ipohm;
Ilh=datakp1.gene.ihyb; 
Ifce=datakp1.gene.ifce;
Ibs=datakp1.gene.iboot;
Ini=datakp1.gene.ini; %non-inductive
Inbi=datakp1.gene.iidn;

Ipmes=datakp1.cons.ip; %other current related 'stuff'
li=datakp1.equi.li;
Vmes=datakp1.cons.vloop; %surface loop voltage
Vloop=datakp1.gene.vres; %average loop voltage
fluxmes=datakp1.cons.flux; %edge poloidal flux

%PROFILES

q=datakp1.prof.q;
Bphi=datakp1.prof.bphi;
epar=datakp1.prof.epar;
s=datakp1.prof.shear;
Te=datakp1.prof.te/1000; %
Ti=datakp1.prof.ti/1000; %
ne=datakp1.prof.ne;
psi=datakp1.prof.psi;
ptot=datakp1.prof.ptot; %total pressure

Pfce=datakp1.source.fce.el; %source and sink deposition profiles
Plh=datakp1.source.hyb.el;
Pnbi_i=datakp1.source.idn.ion;
Pnbi_e=datakp1.source.idn.el;
Pfci_i=datakp1.source.fci.ion;
Pfci_e=datakp1.source.fci.el;
Pfus_i=datakp1.source.fus.ion;
Pfus_e=datakp1.source.fus.el;
Pohm=datakp1.source.ohm;
P_i=datakp1.source.totale.ion;
P_e=datakp1.source.totale.el;
q_ei=datakp1.source.qei;  %sinks
Psyn=datakp1.source.cyclo;
P_rad=datakp1.source.prad;
P_brem=datakp1.source.brem;

Pnbi=Pnbi_i+Pnbi_e;
Pfci=Pfci_i+Pfci_e;
Pfus=Pfus_i+Pfus_e;
Ptot=Pnbi+Pfci+Pfus+Pohm+Pfce+Plh;

jfce=datakp1.source.fce.j; %current deposition profiles
jlh=datakp1.source.hyb.j;
jnbi=datakp1.source.idn.j;
jfci=datakp1.source.fci.j;
j=datakp1.prof.jmoy;
jbs=datakp1.neo.jboot;
jphi=datakp1.prof.jphi;
psi=datakp1.prof.psi*2*pi; %for calculating flux consumption

%TRANSPORT COEFFICIENTS
%chie=chiecal; chii=chiical;
chie_int=datakp1.prof.flux.kean./datakp1.prof.ne; %effect el heat trans dif coeff
chii_int=datakp1.prof.flux.kian./datakp1.prof.ni;
chie_pre=datakp1.coef.ee./datakp1.prof.ne; 
chii_pre=datakp1.coef.ii./datakp1.prof.ni;
chie_neo=datakp1.neo.coef.ee./datakp1.prof.ne;  %neoclassical heat diffusitivity electrons
chii_neo=datakp1.neo.coef.ii./datakp1.prof.ni;  %neoclassical heat diffusitivity ions
qe=datakp1.prof.flux.qe; %electron heat flux 
qi=datakp1.prof.flux.qe; %ion heat flux 

%CALCULATED VALUES
betaN=100*beta.*a.*B./Ip*1e6; %
f_b_s=Ibs./Ip; %
f_n_i=Ini./Ip; %
nG=Ip./1e6/pi./a.^2; %%
f_G=(nbar/1e20)./nG; %
tauIPB=0.0562*(Ip/1e6).^0.93.*B.^0.15.*(P_loss/1e6).^(-0.69).*(nbar/1e19).^0.41 ...
             .*M^0.19.*R.^1.97.*eps.^0.58.*k.^0.78;
H=tauth./tauIPB; %
Plhtrans=2.76*M.^(-1).*B.^0.96.*(nbar/1e20).^0.77.*R.^1.23.*a.^0.76; %ITER physics basics
Plhtrans2_=0.072*Bout^0.7.*(nbar/1e20).^0.7*S.^0.9.*(Zeff/2).^0.7.*2./M;  %ITPA H-mode power threshold database working group


vstep=0.22;
hstep=0.3;
fsize=12;
vtop=1;

%TOTAL CURRENTS
figure; set(gcf,'PaperType','A4','PaperUnits','centimeters','PaperPosition',[0.63 0.63 19.72 28.41]);
subplot(4,2,1);
title(['run',num2str(runnum),' at t = ' num2str(t)],'FontSize',18)

text(0.1,vtop-vstep*1,['Ip=',num2str(round(100*Ip/1e6)/100),'MA'],'FontSize',fsize);
text(0.1,vtop-vstep*2,['I_{NBI}=',num2str(round(100*Inbi/1e6)/100),'MA'],'FontSize',fsize);
text(0.1,vtop-vstep*3,['I_{ECCD}=',num2str(round(100*Ifce/1e6)/100),'MA'],'FontSize',fsize);
text(0.1,vtop-vstep*4,['I_{LH}=',num2str(round(100*Ilh/1e6)/100),'MA'],'FontSize',fsize);

text(0.4,vtop-vstep*1,['I_{BS}=',num2str(round(100*Ibs/1e6)/100),'MA'],'FontSize',fsize);
text(0.4,vtop-vstep*2,['<n>=',num2str(round(100*nbar/1e20)/100),'\cdot10^{20}cm^{-3}'],'FontSize',fsize)
text(0.4,vtop-vstep*3,['P_{\alpha}=',num2str(round(100*P_fus/1e6)/100),'MW'],'FontSize',fsize);
text(0.4,vtop-vstep*4,['P_{rad}=',num2str(round(100*Pradtot/1e6)/100),'MW'],'FontSize',fsize);

text(0.7,vtop-vstep*1,['P_{NBI}=',num2str(round(100*P_nbi/1e6)/100),'MW'],'FontSize',fsize);
text(0.7,vtop-vstep*2,['P_{EC}=',num2str(round(100*P_ec/1e6)/100),'MW'],'FontSize',fsize);
text(0.7,vtop-vstep*3,['P_{LH}=',num2str(round(100*P_lh/1e6)/100),'MW'],'FontSize',fsize);
text(0.7,vtop-vstep*4,['P_{IC}=',num2str(round(100*P_ic/1e6)/100),'MW'],'FontSize',fsize)

subplot(4,2,2);

text(0.1,vtop-vstep*1,['Z_{eff}=',num2str(round(100*Zeff)/100)],'FontSize',fsize);
text(0.1,vtop-vstep*2,['l_i=',num2str(round(100*li)/100)],'FontSize',fsize);
text(0.1,vtop-vstep*3,['\beta_N=',num2str(round(100*betaN)/100)],'FontSize',fsize);
text(0.1,vtop-vstep*4,['f_{BS}=',num2str(round(100*f_b_s)/100)],'FontSize',fsize);

text(0.4,vtop-vstep*1,['f_{NI}=',num2str(round(100*f_n_i)/100)],'FontSize',fsize);
text(0.4,vtop-vstep*2,['f_G=',num2str(round(100*f_G)/100)],'FontSize',fsize)
text(0.4,vtop-vstep*3,['Q=',num2str(round(100*Q)/100)],'FontSize',fsize);
text(0.4,vtop-vstep*4,['H=',num2str(round(100*H)/100)],'FontSize',fsize);

text(0.7,vtop-vstep*1,['V_{loop}=',num2str(round(100*Vloop)/100)],'FontSize',fsize);
text(0.7,vtop-vstep*3,['P_{LH}=',num2str(round(100*Plhtrans)/100)],'FontSize',fsize);
text(0.7,vtop-vstep*4,['P_{loss}=',num2str(round(100*P_loss/1e6)/100)],'FontSize',fsize);

subplot(4,2,3);
xbord=[0,1]; dx=0.1;
plot(r,q,'b','linewidth',1);
xlabel('x','FontSize',fsize);
    ymin=floor(min(q));
    ymax=ceil(max(q));
    set(gca,'XLim',[0 1],'YLim',[ymin ymax],'XTick',0:0.1:1,'YTick',ymin:1:ymax,...
        'XGrid','on','YGrid','on','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.05]);  
    if max(q)>10
        set(gca,'XLim',[0 1],'YLim',[ymin ymax],'XTick',xbord(1):dx:xbord(2),'YTick',ymin:2:ymax,...
            'XGrid','on','YGrid','on','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.05]);
    end
%TEMPERATURES
subplot(4,2,4);    
plot(r,Ti,'b',r,Te,'r','linewidth',1);
ymax=ceil(max([max(Te) max(Ti)]));

set(gca,'XLim',[xbord(1) xbord(2)],'YLim',[0 ymax],'XTick',xbord(1):dx:xbord(2),'XGrid','on','YGrid','on',...
    'TickLength',[0.02 0.05]);  
xlabel('x','FontSize',fsize);
ylabel('T [KeV]','FontSize',fsize);  
legtemp=legend('T_i','T_e');
set(legtemp,'FontSize',12); 

%CURRENT PROFILES 1

subplot(4,2,5);
plot(r,jfce,r,jlh,r,jbs,r,jnbi,r,jfci,'linewidth',1); 

set(gca,'XLim',[xbord(1) xbord(2)],'YLim',[0 inf],'XTick',xbord(1):dx:xbord(2),...
        'XGrid','on','YGrid','on','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.05]); 
xlabel('x','FontSize',fsize);
ylabel('j [A/m^2]','FontSize',fsize);        

legtemp=legend('j_{ECCD}','j_{LH}','j_{BS}','j_{NBI}','j_{ICCD}');
set(legtemp,'FontSize',12,'Location','WestOutside'); 

%CURRENT PROFILES 2

subplot(4,2,6);
jni=jfce+jlh+jnbi+jfci+jbs;
johm=j-jni;
plot(r,jni,r,johm,r,j,'linewidth',1); 
set(gca,'XLim',[xbord(1) xbord(2)],'YLim',[0 inf],'XTick',xbord(1):dx:xbord(2),'XGrid','on','YGrid','on'...
        ,'XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.05]); 
xlabel('x','FontSize',fsize);
ylabel('j [A/m^2]','FontSize',fsize);        
legtemp=legend('j_{NI}','j_{OHM}','j_{TOT}');
set(legtemp,'FontSize',12,'Location','EastOutside'); 

%POWER DEPOSITION PROFILES
subplot(4,2,7)
plot(r,Pfce,r,Plh,r,Pfus,r,Pnbi,r,Pfci,r,Ptot,'linewidth',1); 
set(gca,'XLim',[xbord(1) xbord(2)],'YLim',[0 inf],'XTick',xbord(1):dx:xbord(2),'XGrid','on','YGrid',...
        'on','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.05]); 
xlabel('x','FontSize',fsize);
ylabel('P [W/m^3]','FontSize',fsize);        
legtemp=legend('P_{ECRH}','P_{LH}','P_{fus}','P_{NBI}','P_{ICRH}','P_{tot}');
set(legtemp,'FontSize',12,'Location','WestOutside'); 
set(gca,'FontSize',10);

%transport coefficients
subplot(4,2,8);
chie=chie_int;
chii=chii_int;
plot(r,chie,r,chii,'linewidth',1);
ymin=floor(min([min(chie) min(chii)])); ymax=ceil(max([max(chie) max(chii)]));
set(gca,'XLim',[xbord(1) xbord(2)],'YLim',[ymin ymax],'XTick',xbord(1):dx:xbord(2),'YTick'...
        ,ymin:1:ymax,'XGrid','on','YGrid','on','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.05]); 
xlabel('x','FontSize',fsize);
ylabel('\chi [m^2s^{-1}]','FontSize',fsize);
legtemp=legend('\chi_e','\chi_i');
set(legtemp,'FontSize',12,'Location','EastOutside'); 

unix(['rm ' datapath '*.*']);


