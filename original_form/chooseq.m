%J.Citrin 10.3.2011

%GIVEN A Q-PROFILE, DEFINED AS q BELOW, THIS SCRIPT CALCULATES THE NECESSARY PSI PROFILE FOR ACHIEVING THAT Q-PROFILE. THE SCRIPT EXPECTS
%THE DATA STRUCTURE OF THE RELEVANT RUN IN THE WORKSPACE. BEST RESULTS ARE ACHIEVED IF THE EQUILIBRIUM QUANTITITES IN THE RELEVANT TIMEFRAME
%ARE CLOSE TO THE FINAL EQUILIRBIUM QUANTITIES ACHIEVED AFTER THE RUN. AN ITERATION CAN BE CARRIED OUT IF NECESSARY

%THE CALCULATION IS BASED ON AN INVERSION OF EQUATION B36 FROM THE CRONOS USERGUIDE
clear psi

mu0=4*pi*1e-7;
dx=0.01;
x=linspace(0,1,101);

clear jtpar1 jtpar2 jppar1 jppar2 psi djtpar2 djppar2 psi q
 
%Set the desired q-profile here below. 
 

%k is the time index range in which geometrical terms are averaged over (needed for calculating psi and/or j)
k=1:3;

%set the desired q profile below

%%%%% INSERT NEW Q%%%%%
newq=mean(jetqpol77914.qPOLx(123:126, :));
%%%%%%%%%%%%%%%%%%%%%%%

qpsi=newq; 
qj=newq';
for i=1:length(k)-1
	qj=[qj,newq'];
end

if length(k) == 1
	qj=qj';
end
qj=qj';
%calculate psi from q

if length(k)>1
	psipar=mean(-1./(4.*pi^2).*data.equi.r2i(k,:).*data.equi.vpr(k,:).*mean(data.equi.rhomax(k))...
	    .*data.equi.F(k,:));
else
	psipar=-1./(4.*pi^2).*data.equi.r2i(k,:).*data.equi.vpr(k,:).*mean(data.equi.rhomax(k))...
	    .*data.equi.F(k,:);
end

psi(1)=0; %sets arbitrary constant

for i=2:101
    psi(i)=trapz(0:dx:(dx*(i-1)),psipar(1:i)./qpsi(1:i));
end

vpr=data.equi.vpr;%(k,:);
rat=linspace(0,0.29,30);
%rat(1)=1e-1;
for i=k
%	vpr(i,1:30)=vpr(i,31).*(rat./0.3).^2;
%	vpr(i,1)=2e-1;
end
 
 
%calculate toroidal current from q

if length(k)>1
	jtpar1=mean(data.equi.F(k,:)./4./mu0./vpr(k,:)./pi^2);
	jtpar2=mean(data.equi.F(k,:).*(vpr(k,:).^2).*data.equi.r2i(k,:).*data.equi.grho2r2(k,:)./qj);
	djtpar2=diff(jtpar2)./(dx.*mean(data.equi.rhomax(k))); %djtpar2=[djtpar2 djtpar2(end)]; 
else
	jtpar1=data.equi.F(k,:)./4./mu0./vpr(k,:)./pi^2;
	jtpar2=data.equi.F(k,:).*(vpr(k,:).^2).*data.equi.r2i(k,:).*data.equi.grho2r2(k,:)./qj;
	djtpar2=diff(jtpar2)./(dx.*data.equi.rhomax(k)); %djtpar2=[djtpar2 djtpar2(end)]; 
end

djtpar2_v1=[0 djtpar2];
djtpar2_v2=[djtpar2 0];
djtpar2=(djtpar2_v1+djtpar2_v2)/2;

jt=jtpar1.*djtpar2; 
jt(101)=0;
jt(1)=jt(2); 
pp=splinefit(x,jt,[0 1],8);
yp=ppval(pp,x);

jt=yp;
 
%calculate poloidal current from q
if length(k)>1
	jppar1=-mean(1./pi^2./4./mu0./qj.*data.equi.F(k,:).*vpr(k,:).*data.equi.r2i(k,:).*data.equi.grho2r2(k,:));
	jppar2=mean(data.equi.F(k,:));
	djppar2=diff(jppar2)./(dx.*mean(data.equi.rhomax(k))); %djppar2=[djppar2 djppar2(end)]; 
else
	jppar1=-1./pi^2./4./mu0./qj.*data.equi.F(k,:).*vpr(k,:).*data.equi.r2i(k,:).*data.equi.grho2r2(k,:);
	jppar2=data.equi.F(k,:);
	djppar2=diff(jppar2)./(dx.*data.equi.rhomax(k)); %djppar2=[djppar2 djppar2(end)]; 
end
djppar2=[djppar2 0];

%djppar2=[djppar2(1) djppar2];
jp=jppar1.*djppar2; 

%calculate total current
if length(k)>1
	jmoy=sqrt(jt.^2+jp.^2)/mean(data.geo.b0(k));
else
	jmoy=sqrt(jt.^2+jp.^2)/data.geo.b0(k);
end
%dirty trick
%dirty=[linspace(1.4,1,15) linspace(1,1,86)];   
%jmoy=jmoy./dirty;

doit=input('Input 1 for changing psi in data structure, input 2 for changing jmoy in data structure: ');
if doit==1
    for i=k
        data.prof.psi(i,:)=psi;
        data.mode.psi(:)=1; %of course PSI must now be prescribed
    end
end

if doit==2
    for i=k
        data.prof.jmoy(i,:)=jmoy;
	data.mode.psi(:)=2;
    end
end



if doit==3

jtor_in=grho2r2339(231,:).*psid1339(231,:).*vpr339(231,:)./rhomax339(231);
jtor=diff(jtor_in)/0.01./rhomax339(231);
jtor_v1=[0 jtor]; jtor_v2=[jtor 0]; jtor=(jtor_v1+jtor_v2)/2;
jtor=-F339(231,:)./mu./vpr339(231,:).*jtor/5.3;

end

if doit==4

psi1co=log(grho2r2339(231,:).*vpr339(231,:)./F339(231,:));

psi1co=diff(psi1co)./0.01; 
psi1co_v1=[0 psi1co]; psi1co_v2=[psi1co 0];
psi1co=(psi1co_v1+psi1co_v2)/2; 

jnitest=-grho2r2339(231,:)./mu./5.3./rhomax339(231).^2.*F339(231,:).*(psid2339(231,:)+psi1co.*psid1339(231,:));

end

