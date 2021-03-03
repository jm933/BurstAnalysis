% % AIM: Analyse spike trains


clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bin size for the spike train
bs = 10;   % ms

% time-scale for smoothing the spike train
Ts = 200;   % ms

% burst onset is defined as when the smoothed firing rate r(t) rises above theta
% burst offset is defined as, after burst onset, when r(t) first dips below theta
theta = 3;    % spikes/s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The first set of data
fn = '2021-03-01-CS-traces.xlsx';
datafn = 'data.mat';
genotype = 1; %adjust based on number of genotypes being analysed 

if ~exist(datafn,'file')
    range{1} = {'B:B',"C:C",'D:D','E:E','F:F','G:G','H:H','I:I','J:J','K:K', 'L:L','M:M','N:N','O:O','P:P','Q:Q'};
%add ranges for other genotypes     
    Ntypes = numel(range);
    
    clear spiketime type
    for type = 1:Ntypes
        Ncells = numel(range{type});
        for i = 1:Ncells
            % add fake spike at 0
            spiketime{type}{i} = [xlsread(fn,type,range{type}{i})];
        end
    end
    clear type i
    save(datafn)
    
else
    load(datafn)
end
clear fn datafn range
%spiketime1 = spiketime;
%clear spiketime

%% analyse...
close all
clear h
% f1 = figure;
% f2 = figure;
% f3 = figure;
% f4 = figure;
% Ntypes = 1;
for type = 1:Ntypes
    Ncells = numel(spiketime{type});
    figure
    for i = 1:Ncells
        h(i) = subplot(Ncells,1,i);
        st = spiketime{type}{i};
        isi = diff(st);
        stb = zeros(1,ceil(max(st)/bs));
        tb = bs*(1:numel(stb))/1000;    % s
        stidx = ceil(st/bs);
        for idx = 1:numel(stidx)
            stb(stidx(idx)) = stb(stidx(idx))+1;
        end
        clear stidx idx idx1
        Lkernel = round(Ts/bs);
        kernel = ones(1,Lkernel)/Lkernel;
        rate = conv(stb,kernel,'same')*(1000/bs);   % spikes per second
        plot(tb,rate,'r')
        axis tight
        hold on
        if ~isempty(st)
            plot(st/1000,0,'k.')
        end
%         set(gca,'xticklabel','')
        z = zeros(size(rate));
        z(rate>=theta) = 1;
        burstonsetidx = find(diff([0 z 0])==1);
        burstoffsetidx = find(diff([0 z 0])==-1);
        if burstoffsetidx(end)>numel(tb)
            burstonsetidx(end) = [];
            burstoffsetidx(end) = [];
        end
        clear z
%         plot([tb(burstonsetidx)' tb(burstoffsetidx)'],[1 1],'b-')
%         plot(tb(burstonsetidx),1,'go')
%         plot(tb(burstoffsetidx),1,'ro')
        tburststart = zeros(size(burstonsetidx));
        tburststop = zeros(size(burstoffsetidx));
        check4spikes = true(size(burstonsetidx));
        spikesperburst = zeros(size(burstoffsetidx));
        for bidx = 1:numel(burstonsetidx)
            tstart = 1000*tb(burstonsetidx(bidx));
            tstop = 1000*tb(burstoffsetidx(bidx));
            tburststart(bidx) = st(find(st>=tstart,1));
            tburststop(bidx) = st(find(st<=tstop,1,'last'));
            if sum((st>=tstart) & (st<=tstop)) ==0
                check4spikes(bidx) = false;
            end
            spikesperburst(bidx) = sum((st>=tstart)&(st<=tstop));
        end
        tburststart = tburststart(check4spikes);
        tburststop = tburststop(check4spikes);
        spikesperburst = spikesperburst(check4spikes);
        clear check4spikes
%         if type==2 & i==1
%             tstart21 = tb(burstonsetidx);
%             tstop21 = tb(burstoffsetidx);
%             tburststart21 = tburststart;
%             tburststop21 = tburststop;
%         end
        burstdurations{type}{i} = tburststop-tburststart;   % ms
        burstdurmedian(type,i) = median(burstdurations{type}{i});
        burstdurmean(type,i) = mean(burstdurations{type}{i});
        burstdursd(type,i) = std(burstdurations{type}{i});
        burstduriqr(type,i) = iqr(burstdurations{type}{i});
        T = (tb(end)-tb(1));    % s
        firingrate(type,i) = sum(stb)/T;   % spikes/s
        burstrate(type,i) = numel(burstdurations{type}{i})/T;   % bursts/s
        spikesperburstmean(type,i) = mean(spikesperburst);

        plot([tburststart' tburststop']/1000,[1 1],'b-')
        plot(tburststart/1000,1,'go')
        plot(tburststop/1000,1,'ro')

        %         figure(f4)
%         isi = diff(st);
%         subplot(Ntypes,Ncells,(type-1)*Ncells+i);
%         loglog(isi(1:end-1),isi(2:end),'k.')
%         axis([1 3000 1 3000])
%         plot(isi(1:end-1),isi(2:end),'k.')
%         axis([0 3000 0 3000])
        
    end
    linkaxes(h,'x')
end
% linkaxes(h,'x')
clear type Ncells i st isi

%% summary statistics
% burstrate
% burstdurmean/1000
% burstdurmedian/1000
% burstdursd/1000
% burstduriqr/1000


figure
Nsp = 4;

subplot(2,genotype,1)
Nb = 1000;
tmp = [burstdurations{:}];
tmp = [tmp{:}]/1000;
[~,bc] = hist(tmp,Nb);
n = zeros(genotype,Nb);
for type = 1:genotype
    n(type,:) = hist([burstdurations{type}{:}]/1000,bc);
    n(type,:) = n(type,:)/sum(n(type,:));
end
xlabel('burst dur (s)')
%subplot(2,3,4)
plot(bc,n,'.-')
axis tight
axis([0 5 0 inf])
legend('1','2','3')
xlabel('burst dur (s)')
ylabel('freq (normalised)')
title('Burst duration histogram')

subplot(2,genotype,2), hold on
plot(1:genotype,burstdurmean/1000,'k.')
plot(1:genotype,mean(burstdurmean')/1000,'ko')
ylabel('mean burst duration (s)')
xlabel('cell type')
title('mean burst duration')

subplot(2,genotype,3), hold on
plot(1:genotype,burstdursd/1000,'k.')
plot(1:genotype,mean(burstdursd')/1000,'ko')
ylabel('SD burst duration (s)')
xlabel('cell type')
title('SD burst duration')

subplot(2,genotype,4), hold on
plot(1:genotype,burstrate,'k.')
plot(1:genotype,mean(burstrate'),'ko')
ylabel('bursts/s')
xlabel('cell type')
title('burst rate')

subplot(2,genotype,5), hold on
plot(1:genotype,firingrate,'k.')
plot(1:genotype,mean(firingrate'),'ko')
ylabel('spikes/s')
xlabel('cell type')
title('firing rate')

subplot(2,genotype,6), hold on
plot(1:genotype,spikesperburstmean,'k.')
plot(1:genotype,mean(spikesperburstmean'),'ko')
ylabel('spikes/burst')
xlabel('cell type')
title('Spikes per burst')
