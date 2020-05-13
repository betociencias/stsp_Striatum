%Graphs
clear all
clf 


file_path = 'C:\Users\Janet Barroso\Desktop\Work\Modelo\Matrices\';
%;'%Tf500Pinf05q02\'

strFac='tauFaci10000.0';
strPinf='pinf0.3';
strq='q0.7';
strext='.mat';


a=50:10:110;
[~, col]=size(a);
scale =1;
strTauRec='tauRec';
strdigits='.0';
rec2=[];
RecList={};
rec=scale.*a;
for m = 1:col 
    p=floor(rec(m))./rec(m);
    if p==1
        strRec=num2str(rec(m));
        strRec2=[strTauRec,strRec,strdigits];
    else
        strRec=num2str(rec(m));
        strRec2=[strTauRec, strRec];
    end
    RecList{m}=strRec2;
end

file_names={};
for k=1:col
    filename=[strFac,RecList{k}, strPinf, strq, strext];
    file_names{k}= filename;
end
files=[];
for i = 1:length(file_names)
    str_fullname=[file_path,file_names{i}];
    file = load(str_fullname);
    %file=load([file_path,file_names(i)])
    files=[files; file];
end
 


FirstSpike = 10;
Period = 1000/20;
TimeStart = 0;
TimeMax = 5000;
TimeStep = 0.1;
SampTimes = [TimeStart: TimeStep: TimeMax];
PreSpikes = [FirstSpike: Period: TimeMax];
BellTrain = TrainGauss(SampTimes, PreSpikes, 1);

subplot(1,2,1)
%plot(BellTrain(:)); hold on
%plot(files(2,1).orbit(1,:), 'b-') %Ocupancy of the release pool
%hold on
%plot(files(2,1).orbit(2,:), 'r-') % release probability
plot(files(2,1).orbit(1,:).*files(2,1).orbit(2,:), 'k-')
%ylim([0 1])
%xlim([0 5e5])
%hold off
IPSC = []; 
 
for j = 1:length(file_names)
    amp=[];
    
    for i= 1.18e004: 5e4: length(files(1,1).orbit)
    %for i= 8000: 5e4: (length(files(1,1).orbit))
        
         %a = min(files(j,1).orbit(1,i:(i+ (4.2e4)-1)).*files(j,1).orbit(2,i:(i+ (4.2e4)-1)));
         a = min(files(j,1).orbit(1,i:(i+ (3.82e4)-1)).*files(j,1).orbit(2,i:(i+ (3.82e4)-1)));
         amp = [amp, a];
     
    end
    
    IPSC = [IPSC; amp];
end
[row, col]= size(IPSC);


IpscnIpsc1=[];
for i = 1: row
    IpscnIpsc1(i,:) = IPSC(i,:)./IPSC(i,1);
end
subplot(1,2,2)
plot(IpscnIpsc1(:,:)')
hold on
%legend(file_names(:))



%Datos experimentales 

stim = 1: 1 : 10;
AjusteSTDc =(0.8036-1.36408*exp(-stim/0)).*(1.36408*exp(-stim/1.6197)+0.4968);%Doublecheck
AjusteSTFc = (1.3329-1.0202*exp(-stim/2.1527)).*(1.0202*exp(-stim/0)+ 1.3882);%Doublecheck
AjustePBc = (2.1348-3.1090*exp(-stim/1.6353)).*(3.1090*exp(-stim/1.8831)+0.3986);%Doublecheck
%AjusteSPNc = (1.2192-1.5183*exp(-stim/0)).*(1.5183*exp(-stim/0.8124)+ 0.3759);%Doublecheck
%AjusteSPNL =(0.7214-2.5426*exp(-stim/0)).*(2.5426*exp(-stim/1.0732)+ 0.3731);%Doublecheck
%AjusteSTDL = (0.8080-1.3902*exp(-stim/0)).*(1.3902*exp(-stim/1.5636)+0.4937);%Doublecheck
AjusteSTFL = (1.8867-1.7358*exp(-stim/8.4923)).*(1.7358*exp(-stim/0)+ 2.4098);%Doublecheck
AjustePBL = (0.1672-6.7600*exp(-stim/0.1969)).*(6.7600*exp(-stim/5.1625)+2.4272);%Doublecheck
TeoricoSTD =[1,0.740157187979938,0.608025509055085,0.529166744756153,0.477496788589588,0.441534359510047,0.415438196638262,0.395922340959858,0.380996207935543,0.369383045721080];%New model: logaritmic x(t)
TeoricoSTFc =[1,1.33573489527345,1.51587407118270,1.62686838009015,1.70079437938734,1.75151757520505,1.78708125280952,1.81244284502182,1.83077557753338,1.84416990148608]; %New model:logaritmic x(t)
TeoricoSTFL=[1,1.38254590708514,1.71114566273384,1.99346647322901,2.23607055393364,2.44457686828141,2.62380199278875,2.77787573571024,2.91034022600142,3.02423598346546]; %New model:logaritmic x(t)
%TeoricoPBc = [1,1.50117414726664,1.55169139262031,1.32112091103811,1.40385325034201,1.48906056392916,1.02302831954817,1.15051985071896,1.27975631327949,0.804787880489809];
TeoricoPBc=[1,1.50574095030732,1.26607123717361,1.13836161903652,1.08059513888682,1.05425369177588,1.04149319171292,1.03493408569692,1.03142134224737,1.02949396052672];
%TeoricoPBL=[1,1.23071458637905,1.03413667091095,0.859474535245435,0.744162925637215,0.676825557869262,0.639126070562074,0.617837429095097,0.605324350891986,0.597566943961941];
TeoricoPBL=[1,1.21165008937181,1.01058645715053,0.870306257521906,0.781374868420163,0.724259734387429,0.686482765489745,0.660797033965318,0.642942871630930,0.630323386405744]; %New model:logaritmic x(t)
plot(stim,AjusteSTDc,'k-','linewidth',2)
%hold on
plot(stim,AjusteSTFc,'k-','linewidth',2)
plot(stim,AjustePBc,'k-','linewidth',2)
plot(stim,AjustePBL,'k-','linewidth',2)
%plot(stim,AjusteSTFL,'k-','linewidth',2)
%plot(stim, TeoricoSTD, 'r-', 'linewidth',3)%**New Tf=2,Tr=64,pinf=0.45,q=0.9
%plot(stim, TeoricoSTFc, 'b-', 'linewidth',3)%**New Tf=650,Tr=1.7, pinf=0.2, q=0.2
%plot(stim, TeoricoSTFL, 'b-', 'linewidth', 3)%**New Tf =10000,Tr=0.01,pinf=0.5,%q=0.15
plot(stim, TeoricoPBL, 'g-', 'linewidth', 3)%**New Tf =50,Tr=50,pinf=0.6,q=0.3
%plot(stim, TeoricoPBc, 'g-', 'linewidth', 3)%Tf =700,Tr=90,pinf=0.8, q=0.5
%ylim([0.37 4])
hold off

Errores = [];
ErrMin = [];
for i=1:row
    
    Resta = [];
    for j=1:length(file_names)
        err =  AjusteSTDc-IpscnIpsc1(j,:);
        Resta =[Resta; err];
    end

    
    Err =Resta(i,:).*Resta(i,:);
    Errores = [Errores; Err];
    Min = min(sum(Errores(i,:)));
    
    ErrMin =[ErrMin; Min];
    
end
BestFitted = find(min(ErrMin));
file_names(BestFitted)


    