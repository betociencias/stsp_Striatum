%new
tic

Stim = 20; %hz
FirstSpike = 10;
Period = 1000/Stim;
TimeStart = 0;
TimeMax = 1000;
TimeStep = 1;
SampTimes = [TimeStart: TimeStep: TimeMax];
PreSpikes = [FirstSpike: Period: TimeMax];
BellTrain = TrainGauss(SampTimes, PreSpikes, 1);
TauFac = 100; 
Pinf = 0.2;
Xinf = 1;
q = 0.2;
n_inf=1;
p_inf = 0.2;
tau_n = 10;
tau_p = 100;
q = 0.2;
U0 = sparse(1,2002);
U0(1,1:2) = [1, 0.2]; %Initial conditions x(0) =1, p(0) = 0.2
tCrit = PreSpikes;
Tau_Rec =5;

% def TauREcSweep
%TauRecs = 75.*[0.2: 0.2: 3.001];
%nSimulations = length(TauRecs);
xi=[];
for m = 1: 21
 %   Tau_Rec = TauRecs(m)
    %pars = [Tau_Rec, TauFac, Xinf, Pinf, q, PreSpikes];
    [ai,b] = ode45(@(t,x)SynapticDynamics(t,x, BellTrain(m,:)'), SampTimes' ,U0); %Persistir hasta que salga. Ánimo !!,:
    xi(:,:,m) = b;
end

subplot(2,1,1)
plot(BellTrain(:))


% for j=1:21
%     plot(xi(:,:,j)) % ocupancy release pool
%     hold on
% 
% end
%     plot(xi(:,2), '-b') % probability of release
%     hold off
  %  hold on


bla2 = xi(:,:,1);
[y w z] = size(bla2);

for i = 1: 20
    [j k g]= size(bla2);
    bla2(j+1:(i+1)*y,:,1) = xi(:,:,i+1)
end
subplot(2,1,2)
plot(bla2(:,end))

    
toc    
    
    
    
  