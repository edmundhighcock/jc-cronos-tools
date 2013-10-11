  clear rltitest ti dti dRtemp

%ind=111; %JET
%ind=161;
%ind=161; %JET2
ind=111;
%loc=34;
%ind=384; %JET
%ind=113; %AUG

%kk=110:112; %JET
%kk=172:189; %JET
%kk=110:112; %JET
%kk=159:161;
kk=101:121;
%kk=186:189; %JET GLF23 sims
%kk=103:113; %AUG
%kk=110:113; %AUG GLF23 sims

%ind=161; kk=159:161; 
loc=34;

%a=6.2; %ITER
a=2.86; %JET
%a=1.64; %AUG

%prof=datak.prof.ne;
x=linspace(0,1,101);
%
%prof=mean(ptot612(kk,:));

%prof=ti95mean;

%pion=mean(data.prof.pion(172:189,:));

%pionrot=mean(data.prof.ne(172:189,:).*2*1.67e-27.*data.prof.vtor_exp(172:189,:).^2.*9./2) ;

%prof=data.prof.vtor_exp(ind,:);
%prof=mean(Ti2001(kk,:));  %mean(data.prof.ti(kk,:));
prof=Ti4032(ind,:);

%prof=Ti55(ind,:);
%prof=vtorexp2004(ind,:);
rhomax=data.equi.rhomax(ind);
%prof=mean(data.prof.pion(172:189,:)+pionrot);
test=diff(prof); test1=[test(1) test]; test2=[test test(end)];
test=test1/2+test2/2; 
rlprofx=-a./(prof./(test./(rhomax*0.01)));



%prof=ti95mean/1e3;
%prof=mean(data.impur.impur(kk,:,1));
%prof=interp1(rho,angrotp_exp,x); 
%prof=vtorexp614(189,:); 
%prof(60:77)=linspace(2.308,1.89,18); %20995
%prof(60:76)=linspace(2.1,1.72,17); %20993
%prof(62:73)=linspace(3.1,2.76,12); %20993
%prof=mean(impur72(kk,:,1));
jetto=0;
if jetto==1
 xjet=r33_jet_ti;
 profjet=ti33_jet_eb0(:,6); prof=prof';
 prof=interp1(xjet,profjet,x);
end



%prof=linspace(1e19,0,101);

xrho=data.equi.rhoRZ(ind,:)./data.equi.rhoRZ(ind,101);

Rout=interp1(xrho,data.equi.R(ind,:,1),x);

%prof=prof2;

clear xR dR
xR=zeros(101,65);
for i=1:65
    xR(:,i)=interp1(xrho,data.equi.R(ind,:,i),x);
    xZ(:,i)=interp1(xrho,data.equi.Z(ind,:,i),x);
end

%xR=interp1(xrho,data.equi.R(ind,:,1),x);


rlti=a./data.prof.lti(ind,:);



ltitest=zeros(1,101); 
%ti=zeros(1,101);
%dti=zeros(1,101);
%dRtemp=zeros(1,101);

for i=2:100;
    dprof=(prof(i-1)-prof(i+1));
    dR(i)=mean(sqrt(((xR(i+1,:)-xR(i-1,:)).^2)+(xZ(i+1,:)-xZ(i-1,:)).^2));
    da(i)=data.equi.a(ind,i+1)-data.equi.a(ind,i-1);
    dRout(i)=Rout(i+1)-Rout(i-1);
    lprof(i)=prof(i)/(dprof/dR(i));
    lprofa(i)=prof(i)/(dprof/da(i));
    lprofRout(i)=prof(i)/(dprof/dRout(i));
end

for i=2:101;
    dprof=(prof(i-1)-prof(i));
    dR2(i)=mean(sqrt(((xR(i,:)-xR(i-1,:)).^2)+(xZ(i,:)-xZ(i-1,:)).^2));
     R2(i)=mean(sqrt(((xR(i,:)-xR(1,:)).^2)+(xZ(i,:)-xZ(1,:)).^2));
   % R2(1)=R2(2);
   % dR3(i)=R2(i)-R2(i-1);
    da2(i)=data.equi.a(ind,i)-data.equi.a(ind,i-1);    
    a2(i)=data.equi.a(ind,i);
    dxrho(i)=xrho(i)-xrho(i-1);
end



lprof(1)=lprof(2); lprof(101)=lprof(100);
lprofa(1)=lprofa(2); lprofa(101)=lprofa(100);
lprofRout(1)=lprofRout(2); lprofRout(101)=lprofRout(100);

dprof=prof./lprof;
dprofa=prof./lprofa;
dprofRout=prof./lprofRout;



%kk=input('Input radial position: ');
%lprof=a./lprof(kk*100+1)

%lprofa=a./lprofa(kk*100+1)

%lprofout=a./lprofRout(kk*100+1)

%lprofout=a./lprofa([36:10:76])



dR(1)=dR(2); dR(101)=dR(100);
da(1)=da(2); da(101)=da(100);
dR2(1)=dR2(2); dxrho(1)=dxrho(2);
da2(1)=da2(2); 
rhomax=1.2;
facR=rhomax.^2.*x.*0.01./R2./dR2;
faca=rhomax.^2.*x.*0.01./a2./da2;

rlprof2 = -a./prof.*rpdederive(x,prof,0,2,2,1) ./ data.equi.rhomax(ind,:) .* data.equi.grho(ind,:);

rlprof=a./lprof;
rlprofa=a./lprofa;
rlprofRout=a./lprofRout;

csou=mean(sqrt(2*data.prof.te(kk,:)*1.6e-19./2./1.67e-27));
gamE=mean(data.equi.a(kk,:))./mean(data.prof.q(kk,:)).*rlprofRout.*mean(data.prof.vtor_exp(kk,:))./a.*a./csou;

rlprofx(loc)
rlprof(loc)
rlprofa(loc)
rlprofRout(loc)

error('stopping')

%rlprof2=a.*gprof./prof;

%plot(x,rlti,x,rltitest); set(gca,'XLim',[0 0.8]);

%figure; plot(x,rlprof2); set(gca,'XLim',[0 0.85]);
%figure; plot(x,prof/1000); 


lprof2 = -prof.*rpdederive(x,prof,0,2,2,1) ./ data.equi.rhomax(ind,:) .* data.equi.grho(ind,:);

if flag==1
fsize=12;
plot(x,rlprof2,x,rlprofa,x,rlprofRout,x,roman97(189,:),'--'); set(gca,'XLim',[0.3 0.8]);
t1=xlabel('x'); t2=ylabel('R/L_{Ti}'); 
l=legend('With gradient from flux surface averagedr','With gradient from d(R_{max}-R_{min})/2','With gradient from dR_{max}','Romanelli formula','Location','NorthWest');
t3=title('79626 averaged between 46.1-46.6 s');
set(t1,'FontSize',fsize); set(t2,'FontSize',fsize); set(t3,'FontSize',fsize); set(l,'FontSize',fsize-2);  legend('boxoff');
print(gcf,'-dpsc','-append',['79626_variousRLti_methods.ps']);
end

if flag==2
figure;
if cronos==1
	rlprofa_cro=rlprofa;
	rlprofRout_cro=rlprofRout;
end
if jetto==1
	rlprofa_jet=rlprofa;
	rlprofRout_jet=rlprofRout;
end
plot(x,rlprofa_cro,x,rlprofRout_cro,x,rlprofa_jet,x,rlprofRout_jet); set(gca,'XLim',[0.3 0.75]);
t1=xlabel('x'); t2=ylabel('R/L_{Ti}'); 
l=legend('CRONOS: gradient from d(R_{max}-R_{min})/2','CRONOS: gradient from dR_{max}','JETTO: gradient from d(R_{max}-R_{min})/2','JETTO: gradient from dR_{max}','Location','NorthWest');
t3=title('R/L_{Ti} of CRONOS/JETTO predictions for 77933 at 49.6s');
set(t1,'FontSize',fsize); set(t2,'FontSize',fsize); set(t3,'FontSize',fsize); set(l,'FontSize',fsize-2);  legend('boxoff');

end
