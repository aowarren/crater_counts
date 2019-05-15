%resample_craters_for_sasha
%Edwin Kite (edwin.kite@gmail.com)
%Purpose: Construct shot-noise error bars by resampling for observed crater
%size-frequency distributions.
load data_option

picturesque=1
%load quickcheckpressureeffectsdump
if strcmp(data_option,'mawrth_2_2018')==1
    load mawrth2_pass1
x(1,1,1:length(cr)) = cr;
x(1,2,1:(length(crp)+length(cr))) = [cr,crp];
xu(2,1,1:length(cr)) = cl(1,:); %x uncertainties on measurement %SW:12/14/18 - cl and cpl have TWO ROWS - max and min measurements to give uncertainties around Rfit (cp/crp)
xu(2,3,1:(length(crp)+length(cr))) = [cl(1,:),cpl(1,:)];
xu(2,2,1:length(cr)) = cl(2,:); %x uncertainties on measurement
xu(2,4,1:(length(crp)+length(cr))) = [cl(2,:),cpl(2,:)];


elseif strcmp(data_option,'mawrth_phyllo_2018')==1
    load mawrth_phyllo
x(1,1,1:length(cr)) = cr;
x(1,2,1:(length(crp)+length(cr))) = [cr,crp];
xu(2,1,1:length(cr)) = cl(1,:); %x uncertainties on measurement %SW:12/14/18 - cl and cpl have TWO ROWS - max and min measurements to give uncertainties around Rfit (cp/crp)
xu(2,3,1:(length(crp)+length(cr))) = [cl(1,:),cpl(1,:)];
xu(2,2,1:length(cr)) = cl(2,:); %x uncertainties on measurement
xu(2,4,1:(length(crp)+length(cr))) = [cl(2,:),cpl(2,:)];

elseif strcmp(data_option,'meridiani_2_2018')==1
    load meridiani_2
x(1,1,1:length(cr)) = cr;
x(1,2,1:(length(crp)+length(cr))) = [cr,crp];
xu(2,1,1:length(cr)) = cl(1,:); %x uncertainties on measurement %SW:12/14/18 - cl and cpl have TWO ROWS - max and min measurements to give uncertainties around Rfit (cp/crp)
xu(2,3,1:(length(crp)+length(cr))) = [cl(1,:),cpl(1,:)];
xu(2,2,1:length(cr)) = cl(2,:); %x uncertainties on measurement
xu(2,4,1:(length(crp)+length(cr))) = [cl(2,:),cpl(2,:)];

end



%Resample:

nboots= 1000;%number of bootstraps

%EK 12/5/18 - There is likely a much more elegant way to do the below,
%using the MATLAB function https://www.mathworks.com/help/stats/bootstrp.html

%FOR EACH DTM:
add_measurement_error = 1;
for l = 1:nboots
    disp(l)
    for dtm = 1
        for ct = 1:size(x,2)
            nn = max(find(x(dtm,ct,:)>0));%number of craters
            for bscr = 1:nn %bootstrapped crater
               rnd1 = rand;
               whichcr(bscr) = ceil(rnd1*nn);
               rnd2 = rand;
               %returns raw or error-sampled crater diameter depending on
               %'add_measurement_error' switch
               withmeasureerror(bscr) = xu(dtm,ct*2-1,whichcr(bscr)) + rnd2* ...
                   (xu(dtm,ct*2,whichcr(bscr)) - xu(dtm,ct*2-1,whichcr(bscr))); %Taking account of measurement error
               bs(dtm,ct,bscr,l) = (1-add_measurement_error).*x(dtm,ct,whichcr(bscr))+ ...
                                   add_measurement_error.*(withmeasureerror(bscr));%bootstrap
               if rnd2==rnd1; error('Random number generator is misfiring'); end
                               %bootstrap histogram
            end
              bsh(dtm,ct,:,l) = cumsum(histc(squeeze(bs(dtm,ct,:,l)),[0:0.5:1500]));
              bsh(dtm,ct,:,l) = bsh(dtm,ct,:,l) - bsh(dtm,ct,1,l);
        end
    end
end

% %COMBINE DTMS:
% ndtm=2;
% for l = 1:nboots %Combine the DTMs bootstrap
%     disp(l)
%         for ct = 1:size(x,2)
%             nn = max(find(x(1,ct,:)>0)) + max(find(x(1,ct,:)>0));%number of craters
%                bootstrapthis = squeeze([x(1,ct,:),x(2,ct,:)]); %note change from '3' to 'ct'
%             for bscr = 1:nn %bootstrapped crater
%                bs(ndtm+1,ct,bscr,l) = bootstrapthis(ceil(rand*nn));%bootstrap
%                %bootstrap histogram
%             end
%               bsh(ndtm+1,ct,:,l) = cumsum(histc(squeeze(bs(ndtm+1,ct,:,l)),[0:0.5:1500]));
%               bsh(ndtm+1,ct,:,l) = bsh(ndtm+1,ct,:,l) - bsh(ndtm+1,ct,1,l);
%         end
% 
% end
% 

for ct = 1:size(x,2)
    for dtm = 1
      datahist(ct,dtm,:) = (cumsum(histc(squeeze(x(dtm,ct,:)),[0:0.5:1500])));
    end
end

bshsort = sort(bsh,4,'ascend');
errb = [0.159,0.841];%[0.025,0.975]; %error bars


if picturesque == 1;
%Plot to confirm:
figure;

hold on

myclrs = ['b','r','k';'g','m','k']; %colors for lines for 2 dtms

sbss = size(bshsort);
bshsort = cat(2,zeros(sbss(1),1,sbss(3),sbss(4)),bshsort); %This is a kludge to get around trouble with 'squeeze'

for dtm = 1
    for ct = 2:size(bshsort,2) %converting radii to diameters:
            toshow1=sum(squeeze(bshsort(dtm,1:ct,:,round(nboots*errb(1))))); %first bound
        plot([0:1:3000],toshow1./ ...
         max(toshow1),myclrs(dtm,1),'LineWidth',3-dtm);
        grid minor
        hold on
            toshow2=sum(squeeze(bshsort(dtm,1:ct,:,round(nboots*errb(2))))); %second bound
        plot([0:1:3000],toshow2./ ...
         max(toshow2),myclrs(dtm,2),'LineWidth',3-dtm);
        %set(gca,'xscale','log')
        if (ct-1) >= 2; datacr = squeeze(sum(datahist(1:(ct-1),dtm,:)));else
                    datacr = squeeze((datahist(1:(ct-1),dtm,:)));end
        datacr = datacr - datacr(1);
        plot([0:1:3000],datacr./max(datacr),myclrs(dtm,3),'LineWidth',3-dtm)
    end
end



grid minor
xlim([10 500])
xlabel('diameter (m)') %EK 12/5/18
ylabel('cumulative probability')
title({'Black lines are measurements,','colored lines are the 1-sigma Poisson-error envelope around the measurements'}) 
end