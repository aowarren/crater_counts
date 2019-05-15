%quickcheckpressureeffects_for_sasha
%Edwin Kite (edwin.kite@gmail.com)
%EK 12/5/18: Purpose: Fit crater size-frequency distribution measured in ArcGIS to 
%forward-models of crater size-frequency distributions (filtered by atmospheres
%of different thicknesses) that were calculated in
%Fortran using Jean-Pierre Williams' code
%EK 12/5/18: Script originally written: 2012 (for use in Kite et al. Nature Geoscience
%2014)
%EK 12/5/18: Comments added: Dec 2018 (for use in Sasha Warren's first-year project)
%EK 12/5/18: Note that for basically everything in the workflow, the crater
%size-frequency distributions are 'self-normalized' (we are looking at the
%ratio of big ancient craters to small ancient craters, not the absolute number of small
%anicent craters, or the size of the smallest ancient crater). 
%It *is* possible in principle to use the minimum diameter of
%ancient crater to set an upper bound on the paleo-atmospheric pressure
%(may be a useful back-of-envelope exercise),
%but this (a) throws out most of the data, (b) throws away the chance for
%some internal consistency checks, and (c) means that one has to be 
%very, very sure that the percentage of erosion-era ("modern") craters that
%are misclassified as ancient is tiny (so is less safe).
clear all
close all

data_option = 'mawrth_phyllo_2018' %'September2012' 'mawrth_phyllo_2018' 'mawrth_2_2018' 'meridiani_2_2018'

candidates = 0;

save data_option

load_previous = 0; %Set to zero to rebuild (takes about 5 min on EK laptop)

binsize = 0.5; %histogram width %EK 12/5/18: in meters, for crater diameter

show_not_fractal = 1; %EK 12/5/18: Set this to '1' to overplot dash-dot lines on the main (last-to-plot) figure showing the pressure
%predictions for the case where all the craters are on the same
%stratigraphic level.

fractal_correction = 'on'
% if strcmp(data_option,'mawrth_2_2018') == 1
%     fractal_correction = 'on'
% end

if strcmp(fractal_correction,'on')==1
    fractal_correction_vec = [10*binsize:binsize:3000];
    corrected_for_holsapple_K1_issue_flag = 0;
elseif strcmp(fractal_correction,'off')==1
    fractal_correction_vec = ones(1,5991);
    corrected_for_holsapple_K1_issue_flag = 0;
end 

%Load data:

disp('Loading files.')
if load_previous == 1;
        if exist('quickcheckpressureeffectsdump_December_2018_MAWRTH.mat','file') ~= 2; error;end
        disp('Loading previously assembled model output ...')
        load quickcheckpressureeffectsdump_October_2013_for_sasha.mat
elseif load_previous == 0; 
        for l = 1:8 %EK 12/5/18: These files are forward model results for impactor-atmosphere interactions for atmospheres of
            %differing at-the-surface densities 'rho'. 'rho' = 0.02
            %corresponds roughly to the present atmospheric pressure.

        % if l == 1; aa = importdata('results_ek_rho_0.02.txt');end
        % if l == 2; aa = importdata('results_ek_rho_0.83.txt');end
        % if l == 3; aa = importdata('results_ek_rho_1.500.txt');end
        % if l == 4; aa = importdata('results_ek_rho_3.3.txt');end
        % if l == 5; aa = importdata('results_ek_rho_6.6.txt');end
        % if l == 6; aa = importdata('results_ek_rho_9.9.txt');end
        % if l == 7; aa = importdata('results_ek_rho_13.2.txt');end
        % if l == 8; aa = importdata('results_ek_rho_16.5.txt');end

        % 
        % if l == 1; aa = importdata('results_inclsmall_ek_rho_0.02.txt');end
        % if l == 2; aa = importdata('results_inclsmall_ek_rho_0.415.txt');end
        % if l == 3; aa = importdata('results_inclsmall_ek_rho_0.83.txt');end
        % if l == 4; aa = importdata('results_inclsmall_ek_rho_1.500.txt');end
        % if l == 5; aa = importdata('results_inclsmall_ek_rho_3.3.txt');end
        % if l == 6; aa = importdata('results_inclsmall_ek_rho_6.6.txt');end
        % if l == 7; aa = importdata('results_inclsmall_ek_rho_9.9.txt');end
        % if l == 8; aa = importdata('results_inclsmall_ek_rho_16.5.txt');end
        % 
        % rhoa = ['0.02','0.415','0.83','1.500','3.3','6.6','9.9','13.2'];
        % if l == 1; aa = importdata('results_ek_large_rho_0.02.txt');end
        % if l == 2; aa = importdata('results_ek_large_rho_0.415.txt');end
        % if l == 3; aa = importdata('results_ek_large_rho_0.83.txt');end
        % if l == 4; aa = importdata('results_ek_large_rho_1.500.txt');end
        % if l == 5; aa = importdata('results_ek_large_rho_3.3.txt');end
        % if l == 6; aa = importdata('results_ek_large_rho_6.6.txt');end
        % if l == 7; aa = importdata('results_ek_large_rho_9.9.txt');end
        % if l == 8; aa = importdata('results_ek_large_rho_16.5.txt');end

        % rhoa = ['0.02','0.415','0.83','1.500','3.3','6.6'];
        % if l == 1; aa = importdata('results_ek_10cm_min_rho_0.02.txt');end
        % if l == 2; aa = importdata('results_ek_10cm_min_rho_0.415.txt');end
        % if l == 3; aa = importdata('results_ek_10cm_min_rho_0.83.txt');end
        % if l == 4; aa = importdata('results_ek_10cm_min_rho_1.500.txt');end
        % if l == 5; aa = importdata('results_ek_10cm_min_rho_3.3.txt');end
        % if l == 6; aa = importdata('results_ek_10cm_min_rho_6.6.txt');end
        % %if l == 7; aa = importdata('results_ek_10cm_min_rho_9.9.txt');end
        % %if l == 8; aa = importdata('results_ek_10cm_min_rho_16.5.txt');end

        % rhoa = ['0.02','0.415','0.83','1.500','3.3','6.6','9.9','16.5'];
        % if l == 1; aa = importdata('results_ek_33cm_min_rho_0.02.txt');end
        % if l == 2; aa = importdata('results_ek_33cm_min_rho_0.415.txt');end
        % if l == 3; aa = importdata('results_ek_33cm_min_rho_0.83.txt');end
        % if l == 4; aa = importdata('results_ek_33cm_min_rho_1.500.txt');end
        % if l == 5; aa = importdata('results_ek_33cm_min_rho_3.3.txt');end
        % if l == 6; aa = importdata('results_ek_33cm_min_rho_6.6.txt');end
        % if l == 7; aa = importdata('results_ek_33cm_min_rho_9.9.txt');end
        % if l == 8; aa = importdata('results_ek_33cm_min_rho_16.5.txt');end

%         rhoa = ['0.02','0.415','0.83','1.500','3.3','6.6','9.9','16.5'];
%         rhos = [0.02,0.415,0.83,1.500,3.3,6.6,9.9,16.5];
%         if l == 1; aa = importdata('results2ek1002den0.02.txt');end
%         if l == 2; aa = importdata('results2ek1002den0.415.txt');end
%         if l == 3; aa = importdata('results2ek1002den0.83.txt');end
%         if l == 4; aa = importdata('results2ek1002den1.500.txt');end
%         if l == 5; aa = importdata('results2ek1002den3.3.txt');end
%         if l == 6; aa = importdata('results2ek1002den6.6.txt');end
%         if l == 7; aa = importdata('results2ek1002den9.9.txt');end
%         if l == 8; aa = importdata('results2ek1002den16.5.txt');end

%EK 12/5/18: These files have a minimum cutoff size of impactor given
%in meters (e.g. results2ek1004den16.5min3.0.txt has no impactors smaller
%than 3.0 m). This is to speed up the forward model: I found in the earlier 
%forward-model runs that
%essentially zero impactors smaller than 3.0 m were making it to the
%surface at hypervelocity speeds for a den of 16.5 kg/m^3. As a result,
%almost all of my forward model "shots" were being wasted for thick atmospheres because the
%impactor size frequency distribution heavily weights small impactors, but I
%needed a large number of big impactors to build up a smooth crater size-frequency
%distribution at the crater diameters (hundreds of meters) that *do* make
%it through the 16.6 g/m^3 atmosphere.
        rhoa = ['0.02','0.415','0.83','1.500','3.3','6.6','9.9','16.5'];
        rhos = [0.02,0.415,0.83,1.500,3.3,6.6,9.9,16.5];
        if l == 1; aa = importdata('results2ek1004den0.02min0.1.txt');end
        if l == 2; aa = importdata('results2ek1004den0.415min0.25.txt');end
        if l == 3; aa = importdata('results2ek1004den0.83min0.5.txt');end
        if l == 4; aa = importdata('results2ek1004den1.65min1.0.txt');end
        if l == 5; aa = importdata('results2ek1004den3.3min1.5.txt');end
        if l == 6; aa = importdata('results2ek1004den6.6min2.5.txt');end
        if l == 7; aa = importdata('results2ek1004den9.9min3.0.txt');end
        if l == 8; aa = importdata('results2ek1004den16.5min3.0.txt');end
  
        
%WORKAROUND FOR HOLSAPPLE K1 ISSUE (SEE EMAILS FROM KEITH HOLSAPPLE (U. WASHINGTON)
%EK 12/5/18: The Holsapple K1 parameter is related to target strength
%effects on crater size. For craters >~1 km diameter, even those lacking a
%central peak, we are in the 'gravity regime' (see either of Melosh's
%textbooks) and target strength can be neglected. But all our craters are
%in the smaller-diameter 'strength regime' where hypervelocity craters formed in weak
%material can be a factor of 2 bigger than hypervelocity craters formed in
%strong, intact rock.
%EK 12/5/18: Email Edwin to ask for the Holsapple emails. The basic issue
%is that Holsappple published a paper with some important errors, that are fixed in an
%unpublished and not-obviously-related pdf on Holsapple's website. We added
%the correction to the workflow.
disp('Warning! Do not double-count the Holsapple K1 correction!') 
     
%DEPRECATED: diameter_correction_factor = (0.24/0.132)^(1/3); %cube root of the ratio of the crater-volume K1 factors

%if exist('corrected_for_holsapple_K1_issue_flag') == 0
%    c2 = c; d2 = d; e2 = e;
    load dcorr dcorr %Correction interpolated from new script per Pierre's email of 30 Oct 2013
    dcorr = nanmin(dcorr);
%      diameter_correction_factor = interp1([1:600],dcorr,[1:6001],'nearest','extrap')  
      %note that diameter_correction_factor can be greater than or less
      %than 1
%    for dind = 1:1:max(size(c))        
%      c2(:,dind) = c(:,round(dind*diameter_correction_factor(dind)));
%      d2(:,dind) = d(:,round(dind*diameter_correction_factor(dind)));
%      e2(:,dind) = e(:,round(dind*diameter_correction_factor(dind)));
%    end
%    
%    disp('Done')
%    c = c2; clear c2
%    d = d2; clear d2
%    e = e2; clear e2
%end
        c(l,:) = histc(aa.data(:,11).*interp1([1:600],dcorr,aa.data(:,11),'nearest','extrap'),[0:binsize:3000]);
        d(l,:) = histc(aa.data(:,12).*interp1([1:600],dcorr,aa.data(:,11),'nearest','extrap'),[0:binsize:3000]);
        e(l,:) = histc(aa.data(:,2) .*interp1([1:600],dcorr,aa.data(:,11),'nearest','extrap'),[0:binsize:3000]);

        cold(l,:) = histc(aa.data(:,11)  ,[0:binsize:3000]);
        dold(l,:) = histc(aa.data(:,12)  ,[0:binsize:3000]);
        eold(l,:) = histc(aa.data(:,2)   ,[0:binsize:3000]);
%END OF HOLSAPPLE K1 WORKAROUND

        
        m(l,1:size(aa.data,1),1:size(aa.data,2)) = aa.data;

        cols = aa.textdata;

        clear aa

        disp(l)

        end
end
disp('Loading of model output complete')


    
%Display 'model summary' plots:


%EK 12/5/18: This figure shows the minimum cutoff size of impact given
%the minimum size cutoff of impactor in meters (e.g. results2ek1004den16.5min3.0.txt has no impactors smaller
%than 3.0 m). This is to speed up the forward model: I found in the earlier 
%forward-model runs that
%essentially zero impactors smaller than 3.0 m were making it to the
%surface at hypervelocity speeds for a den of 16.5 kg/m^3. As a result,
%almost all of my forward model "shots" were being wasted for thick atmospheres because the
%impactor size frequency ditribution heavily weights small impactors, but I
%needed a large number of big impactors to build up a smooth crater size-frequency
%distribution at the crater diameters (hundreds of meters) that *do* make
%it through the 16.6 g/m^3 atmosphere.
figure
for l=1:size(c,1)
plot([0:binsize:3e3],c(l,:))
hold on
end
grid on
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'FontSize',16)
xlabel('Zero-atmosphere crater diameter (m)')
ylabel('Incremental frequency per 0.5m-diameter bin')
legend(num2str(rhos)) %EK 12/5/18 addition
title({'Forward model input','(no atmosphere filtering yet, see script comments)'}) %EK 12/5/18 addition

%EK 12/5/18: This figure shows the incremental distribution of craters
%making it to the surface.
figure
for l=1:size(d,1)
plot([0:binsize:3e3],d(l,:)/max(d(l,2:end)))
hold on
end
grid on
set(gca,'yscale','log')
set(gca,'xscale','log') %EK 12/5/18 addition
xlabel('Crater diameter (m)') %EK 12/5/18 addition
ylabel('Incremental frequency per 0.5m-diameter bin')
legend(num2str(rhos)) %EK 12/5/18 addition


figure
for l=1:size(e,1)
plot([0:binsize:3e3],e(l,:)/max(e(l,2:end)))
hold on
end
grid on
set(gca,'yscale','log')
legend(num2str(rhos)) %EK 12/5/18 addition
legend


figure
subplot(2,1,1)
contourf(log10(d))
xlim([0 100])
xlabel('Crater diameter (m)') %EK 12/5/18 addition
ylabel('Forward model number (small number = small atm density)') %EK 12/5/18 addition
title({'log10(number of craters at surface','forward model)'})
colorbar

subplot(2,1,2)
contourf(log10(c))
xlim([0 100])
xlabel('Zero-atmosphere crater diameter (m)') %EK 12/5/18 addition
ylabel('Forward model number (small number = small atm density)') %EK 12/5/18 addition
title({'log10(number of craters that would';'form with no atmosphere, forward model)'})
colorbar

%MAIN CDF: (Main Cumulative Distribution Function for data)

picturesque=0;
add_measurement_error = 0; %0 = shot-noise error only. 1 = include a boxcar prior for the true radius of each crater with bounds of 'crl(1,:)' and 'crl(2,:)'
resample_craters_f %find shot-noise error bar


%Fractal correction: 
%tll = tripletlist. tll(1,:) = density. tll(2,:) = cratersize. tll(3,:) =
%csdf.
tll = [rhos(1)*ones(1,1+3000/binsize); ...
      [0:binsize:3000]; ...
      [zeros(1,10),(cumsum((fractal_correction_vec.*d(1,11:end))')./sum((fractal_correction_vec.*d(1,11:end))))'] ];
for l = 2:size(rhos,2)
    tll = cat(2,tll,[rhos(l)*ones(1,1+3000/binsize); ...
          [0:binsize:3000]; ...
          [zeros(1,10),(cumsum((fractal_correction_vec.*d(l,11:end))')./sum((fractal_correction_vec.*d(l,11:end))))'] ]);
end

tll = tll'; %3 x n


%Interp onto a regular grid

        for l = 1:size(rhos,2)
            model(l,:) = tll(1+(6001*(l-1)):(6001*l),3);
            if min(find(model(l,:)>1e-10)) >1
                esr = min(find(model(l,:)>1e-2)); %end-slope root
                es = (log10(model(l,esr+5))-log10(model(l,esr))) ...
                    /(log10(esr+5) - log10(esr));%end slope
                fnz = min(find(model(l,:)>1e-10));%first non zerox(1,1,1:length(cr)) = cr;
x(1,2,1:(length(crp)+length(cr))) = [cr,crp];
xu(2,1,1:length(cr)) = cl(1,:); %x uncertainties on measurement %SW:12/14/18 - cl and cpl have TWO ROWS - max and min measurements to give uncertainties around Rfit (cp/crp)
xu(2,3,1:(length(crp)+length(cr))) = [cl(1,:),cpl(1,:)];
xu(2,2,1:length(cr)) = cl(2,:); %x uncertainties on measurement
xu(2,4,1:(length(crp)+length(cr))) = [cl(2,:),cpl(2,:)];
                dr = log10(fnz) - log10([1:(-1+fnz)]);%meters
                df = dr.*es;
                nv = 10.^(log10(model(l,fnz)) - df);%new values
                model(l,1:(-1+fnz)) = nv;
            end
        end
 
for l = 1:size(rhos,2)
    tll(1+(6001*(l-1)):(6001*l),3) = model(l,:);
end
    
%clear ZI %TEMPORARY

%if exist('ZI','var') ~= 1
% if corrected_for_holsapple_K1_issue_flag == 0
    disp('Griddata...')
    atmospoi = [10:10:5e3]%atmosp of interest
    %convert to density - take current rho, divide by P in mbars, multiply
    %by P - assumes no change in near surface atmosphere T (pv = nRT)
    atmospt = 0.017/6.6*atmospoi; 
    %Given an atmos and a diameter, tell me the cdf fraction.
    rhodex = 0;
    for atmosr=atmospt 
        atmosr
        rhodex = rhodex + 1; %note that the atmospoi are uncorrected for the elevation of the DTMs
     ZI(rhodex,:) = griddata(tll(:,1),tll(:,2),log10(tll(:,3)),atmosr,[0:binsize:3000]);  %
    end
    ZI = 10.^ZI;

    disp('Exited griddatan.')
% end
%end

%ZI = griddata(tll(:,2),tll(:,3),tll(:,1),XI(:,1),XI(:,2))


%figure
%scatter(XI(:,1),XI(:,2),7,ZI);



disp('Saving ...')


save quickcheckpressureeffectsdump_December_2018_MAWRTH

disp('Save complete.')

%DISPLAY RESULTS

mint=3;

%colscal = 1/-min(isfinite(log10(sum(m(mint:end,:,12)'>0)./length(m(mint:end,:,12)'))).*log10(sum(m(mint:end,:,12)'>0)./length(m(mint:end,:,12)')));


if strcmp(data_option,'mawrth_2_2018')==1
    elevation = 3100;
    dih = xlsread('dch_mawrth2_pass1.xlsx'); %the number of individual craters in each 0.5m-diameter bin
    show_individual_dtms = 0;
    show_this_dtm = [1 1];
    
    dicombined = dih; 
elseif strcmp(data_option,'mawrth_phyllo_2018')==1
    elevation = 3100;
    dih = xlsread('dch_mawrth_phyllo.xlsx'); %the number of individual craters in each 0.5m-diameter bin
    show_individual_dtms = 0;
    show_this_dtm = [1 1];
    
    dicombined = dih; 
elseif strcmp(data_option,'meridiani_2_2018')==1
    elevation = 670;
    dih = xlsread('dch_meridiani_2.xlsx');
    show_individual_dtms = 0;
    show_this_dtm = [1 1];
    
    dicombined = dih; 
end

%******OLD ANALYSIS METHOD, DEPRECATED ********
%FIND ATMOSPHERIC PRESSURE:
% data_opt = 1;%ancient craters only
% find_early_Mars_atmospheric_pressure
% %Plot error bar of bootstrapping on fit to data
% 
% do(1).bsu = bsu; do(1).bsl = bsl; %Save the results
% do(1).bss = bss; do(1).bsI = bsI; 
% do(1).bsb = bsb; do(1).bsC = bsC; 
% 
% 
% data_opt = 2; %Change to include rimmed circular mesas
% find_early_Mars_atmospheric_pressure    
%******OLD ANALYSIS METHOD, DEPRECATED ********



%******NEW ANALYSIS METHOD, LEGITIMATE ********
bayesian_palaeopressure_f
%non-combined "conservative"
do(1).bsu = round(errorhigh(1,1)/10)*10;
do(1).bsl = round(errorlow(1,1)/10)*10;
do(1).bsb = round(bestfit(1,1)/10)*10;
%combined "preferred"
bsl = round(errorhigh(1,3)/10)*10;
bsu = round(errorlow(1,3)/10)*10;
bsb = round(bestfit(1,3)/10)*10;
%******NEW ANALYSIS METHOD, LEGITIMATE ********

save('wegotthisfar')
load('wegotthisfar')


figure('Position', [100 100 1000 700]) %EK 12/5/18. This is the main figure (pretty similar to Fig 2 from the 2014 Nature Geoscience paper)
xlim([10 200])
%SHOWING MODEL OUTPUT:
colorscheme = 2;
mycolors = [255,0,0;255,215,0;... %;%orangered255,69,0; ...
           0,191,255;0,255,255; ...
          127,255,0;0,128,0];
mylinewidth = [2 2 4 4 6 6 8 8];
mypressures = [250,500,1000,2000,3000,4000];   

colz = parula(size(d-4,1));
colz = fliplr(colz)

newmap = brighten(parula(size(d,1)),.5);
colza = fliplr(newmap);

colz = [1,1,1;1,1,1;1,1,1;... %;%orangered255,69,0; ...
           0.6,0.8,1;0.3,0.5,0.9; ...
          0.3,0.3,0.9;0.1,0.01,0.7;1,1,1;1,1,1];

cust = customcolormap([0 0.4 0.55 0.7 1],[1 0 0;1 0 0;0.5 0 1;0 0.8 1;0 0.5 1])


if show_not_fractal == 1
    if strcmp(fractal_correction,'on')
for l =mint:size(d,1)
    if l == 1|l == 2 | l == 3 | l == 8| l == 9
        thecolor = 'none';
        bg = 'none';
    else
        thecolor = [0.7 0.7 0.7];
        bg = 'w';
    end
plot([0:binsize:2999.5],cumsum(d(l,2:end)')/sum(d(l,2:end)'), ...
    'Color',thecolor','LineWidth',2,'LineStyle','-.')
hold on

end

    end
end


for l =mint:size(d,1)
    if l == 1|l == 2 | l == 3| l == 8| l == 9
        thecolor = 'none';
        bg = 'none';
    else
        thecolor = 0.8*cust(l*5,:);
        bg = 'w';
    end
    if colorscheme == 1; colornow = [1-(l-mint+1)*1.25/(size(d,1)) (l-mint+1)*1.25/(size(d,1)) 0];
elseif colorscheme == 2; colornow = 1./255.*mycolors(l-mint+1,:);end
%Note elevation correction
[Cc,Ic] = (min(abs(atmospoi-mypressures(l-mint+1)./exp(-elevation/10700))));


if strcmp(fractal_correction,'on')
    
plot([0:binsize:3000],ZI(Ic,:), ...
    'Color',thecolor, ...
    'LineWidth',mylinewidth(l),'LineStyle','-')
            mystr = strcat(num2str(mypressures(l-mint+1)/1000),' bar')
        if l == 3; mystr = '250 mbar';end
        if l == 4; mystr = '0.5 bar'; end
        if l == 5; mystr = '1 bar'; end
        rtate = [0 0 60 45 45 45 45 40];
        ht(l-mint+1) = text(binsize*max(find(ZI(Ic,:)<0.5)),0.5,mystr,'FontSize',30', ... %binsize*find((cumsum(d(l,2:end)')/sum(d(l,2:end)'))<0.5,1,'last')
            'FontWeight','bold','Color',thecolor,'BackgroundColor',bg)
        set(ht(l-mint+1),'rotation',rtate(l))

elseif strcmp(fractal_correction,'off')
    plot([0:binsize:3000],ZI(Ic,:), ...
    'Color',thecolor, ...
    'LineWidth',mylinewidth(l),'LineStyle','-')
    mystr = strcat(num2str(mypressures(l-mint+1)/1000),' bar')
    if l == 3; mystr = '250 mbar';end
    if l == 4; mystr = '0.5 bar'; end
    if l == 5; mystr = '1 bar'; end
    rtate = [0 0 70 50 50 50 50 50];
    
    ht(l-mint+1) = text(binsize*max(find(ZI(Ic,:)<0.5)),0.5,mystr,'FontSize',30', ...
        'FontWeight','bold','Color',thecolor,'BackgroundColor',bg)
    set(ht(l-mint+1),'rotation',rtate(l))
end
%plot([binsize:binsize:3000],[[0 0 0 0 0 0 0 0 0]'; cumsum((fractal_correction_vec.*d(l,11:end))')./sum((fractal_correction_vec.*d(l,11:end))')], ...
%    'Color',colornow, ...
%    'LineWidth',mylinewidth(l-mint+1))
% mystr = strcat(num2str(mypressures(l-mint+1)/1000),' bar')
% if l == 3; mystr = '250 mbar';end
% if l == 4; mystr = '0.5 bar'; end
% if l == 5; mystr = '1 bar'; end
% rtate = [0 0 20 35 40 50 55 55];
% ht(l-mint+1) = text(binsize*max(find(ZI(Ic,:)<0.5)),0.5,mystr,'FontSize',19', ...
%     'FontWeight','bold','Color',colza(a),'BackgroundColor','w')
% set(ht(l-mint+1),'rotation',rtate(l))
hold on
end
box on

colz = parula(size(d,1));
colz = fliplr(colz)

mylocations = [10 20 35 70 90 125]


error_bar_plotting_method = 2;

if error_bar_plotting_method == 1;
    %Notice I very slightly narrow the error bar on the plot for visibility!
    %By 20 mbar to create a small gap between data_opt 1 and data_opt 2. Different colors also help.
    %Plot error bar of bootstrapping on fit to data
    if candidates == 1
    for l = find(atmospoi==ceil(mean(bsl))):(find(atmospoi==ceil(mean(bsu))))
        if (-1+find((atmospoi==ceil(mean(bsl))))) - l < 12; lwx=0.5;else lwx=2;end
        plot([0:binsize:3000],ZI(l,:),'Color',[0.4,0.4,0.4],'LineWidth',lwx)
        hold on
    end
    end

    for l = (find((atmospoi==ceil(mean(do(1).bsl))))):find(atmospoi==ceil(mean(do(1).bsu)))
        if l-(1+find((atmospoi==ceil(mean(do(1).bsl))))) < 12; lwx=0.5;else lwx=2;end
        plot([0:binsize:3000],ZI(l,:),'Color',[0.8,0.8,0.8],'LineWidth',lwx)
        hold on
    end

elseif error_bar_plotting_method == 2; %hachured patches
    if candidates == 1
    mypatch = [ZI(find(atmospoi==bsl),:), ... %non-combined RCMs
        fliplr(ZI(find(atmospoi==bsu),:))];
    myx = [[0:binsize:3000],[3000:-binsize:0]];
    h1=patch(myx,mypatch,[.500 .500 .500]);hPatch1 = findobj(h1, 'Type', 'patch');
    hold on
    hatchfill(hPatch1,'single',90,1)
     end
    
    mypatch = [ZI(find(atmospoi==do(1).bsl),:), ... %combined w/RCMs
        fliplr(ZI(find(atmospoi==do(1).bsu),:))];
    myx = [[0:binsize:3000],[3000:-binsize:0]];
    h1=patch(myx,mypatch,[0.500 .500 .500]);hPatch1 = findobj(h1, 'Type', 'patch');
    hatchfill(hPatch1,'single',0,5) 
    
   
    
    %Overlay boundary lines (different color)
    if candidates ==1
    for l = [find(atmospoi==bsl),(find(atmospoi==bsu))]
        if (-1+find((atmospoi==ceil(mean(bsl))))) - l < 12; lwx=.5;else lwx=.5;end
        plot([0:binsize:3000],ZI(l,:),'Color',[.500 .500 .500],'LineWidth',lwx)
        hold on
    end
    end

    if candidates == 1
    for l = [(find((atmospoi==do(1).bsl))),find(atmospoi==do(1).bsu)]
        if l-(1+find((atmospoi==ceil(mean(do(1).bsl))))) < 12; lwx=.5;else lwx=.5;end
        plot([0:binsize:3000],ZI(l,:),'Color',[.500 .500 .500],'LineWidth',lwx)
        hold on
    end
    end
    
end
   %Plot best fit from bootstrap procedure
   if candidates == 1
    plot([0:binsize:3000],ZI(find(atmospoi==ceil(mean(bsb))),:),'Color',[.25 .25 .25],'LineWidth',2,'LineStyle','--')
   end
    %Plot best fit from bootstrap procedure
    plot([0:binsize:3000],ZI(find(atmospoi==ceil(mean(do(1).bsb))),:),'Color',[.25 .25 .25],'LineWidth',2)
    
    set(gca,'xscale','log','TickDir','out')

    

%SHOWING DATA:
l=1;
plot(cumsum(c(l,2:end)')/sum(c(l,2:end)'),'LineStyle',':','Color','w','LineWidth',1.3)
hold on

    
if strcmp(data_option,'August2012')==1

    plot(cumsum(dih(1,:)')/sum(dih(1,:)'),'k','LineWidth',2)

    if size(d,1)==6
    legend({'6 mbar','125 mbar','250 mbar','0.5 bar','1 bar','2 bar','Zero-atmosphere','Observations'},'FontSize',16)
    elseif size(d,1)==8
    legend({'6 mbar','125 mbar','250 mbar','0.5 bar','1 bar','2 bar','3 bar','5 bar','Zero-atmosphere','Observations'},'FontSize',16)
    end

else %[September or December data options] if strcmp(data_option,'September2012')==1

    
    if show_individual_dtms == 0;    
        plot(cumsum(dicombined(1,:)')/sum(dicombined(1,:)'),'k','LineWidth',4) %Ancient crater edges
        hold on
        if candidates == 1
        plot(cumsum((dicombined(1,:)+dicombined(2,:))')/sum((dicombined(1,:)+dicombined(2,:))'),'k','LineWidth',4,'LineStyle','--') %Ancient crater edges + raised circular mesas
        %plot(cumsum((dicombined(1,:)+dicombined(3,:)+dicombined(2,:))')/sum((dicombined(1,:)+dicombined(3,:)+dicombined(2,:))'),'k','LineWidth',2,'LineStyle','-.') %Ancient crater edges + raised circular mesas + candidate ancient craters
   
        end
    end

    
    if show_individual_dtms == 1
        if show_this_dtm(1) == 1;
        plot(cumsum(dih(1,:)')/sum(dih(1,:)'),'k','LineWidth',4) %Ancient crater edges
        hold on
        plot(cumsum((dih(1,:)+dih(2,:))')/sum((dih(1,:)+dih(3,:))'),'k','LineWidth',4,'LineStyle','--') %Ancient crater edges + raised circular mesas
        plot(cumsum((dih(1,:)+dih(2,:))')/sum((dih(1,:)+dih(2,:))'),'k','LineWidth',4,'LineStyle','-.') %Ancient crater edges + raised circular mesas + candidate ancient craters

        end
        if show_this_dtm(2) == 1

        plot(cumsum(dih2(1,:)')/sum(dih2(1,:)'),'m','LineWidth',4) %Ancient crater edges
        hold on
        plot(cumsum((dih2(1,:)+dih2(3,:))')/sum((dih2(1,:)+dih2(3,:))'),'m','LineWidth',4,'LineStyle','--') %Ancient crater edges + raised circular mesas
        plot(cumsum((dih2(1,:)+dih2(3,:)+dih2(2,:))')/sum((dih2(1,:)+dih2(3,:)+dih2(2,:))'),'m','LineWidth',4,'LineStyle','-.') %Ancient crater edges + raised circular mesas + candidate ancient craters
        end
        
                if size(d,1)==6
        legend({'6 mbar','125 mbar','250 mbar','0.5 bar','1 bar','2 bar','Zero-atmosphere','OBS:CR','OBS:CR+RCM','OBS:CR+RCM+CAND'},'FontSize',16)
        elseif size(d,1)==8
        legend({'205 mbar','410 mbar','820 mbar','1.6 bar','2.5 bar','4.1 bar','Zero-atmosphere','OBS:CR','OBS:CR+RCM','OBS:CR+RCM+CAND'},'FontSize',16)

        %legend({'6 mbar','125 mbar','250 mbar','0.5 bar','1 bar','2 bar','3 bar','5 bar','Zero-atmosphere','OBS:CR','OBS:CR+RCM','OBS:CR+RCM+CAND'},'FontSize',18)
        end
    end
    

    
end
%END OF SHOWI

grid off;xlim([10 200/binsize])
xlim([10 200])
box on;
set(gca,'xminortick','on');set(gca,'yminortick','on')

set(gca,'FontSize',30,'TickDir','out'); xlabel('Crater diameter (m)','FontWeight','bold'); ylabel('Cumulative fraction of craters','FontWeight','bold')

filename = [data_option,'_fractal',fractal_correction]
savefig(filename)

saveas(gcf,[data_option,fractal_correction],'png')


% title({'Comparison of model crater size?frequency distributions to observations. Solid black line corresponds to definite embedded craters. Dashed black', ...
% 'line additionally includes candidates. Stair-stepping in the data curves corresponds to individual craters. Coloured lines show model', ...
% 'predictions for atmospheric filtering of small impactors at different pressures. Grey hatched regions correspond to 2? statistical error', ...
% 'envelopes around the best-fit palaeopressure to the data (best fits shown by thick grey lines). Survey incompleteness leads to overestimates of', ...
% 'median crater size, so best fits are upper limits. Dash-dot lines (if shown) correspond to model predictions with a fractal correction included.'},'fontsize',12)
% %title('Early Mars atmospheric pressure: comparison of model crater size-frequency distributions to observations')






% %CDF WITH THRESHHOLD:
% %colscal = 1/-min(isfinite(log10(sum(m(:,:,12)'>0)./length(m(:,:,12)'))).*log10(sum(m(:,:,12)'>0)./length(m(:,:,12)')));
% 
% for l =1:size(d,1)
% plot([0:1:2999],[[0 0 0 0 0 0 0 0 0]'; cumsum(d(l,11:end)')/sum(d(l,11:end)')], ...
%     'Color',[min(1,1-(1+colscal*log10(sum(m(l,:,12)>0)./length(m(l,:,12))))) max(0,1+colscal*log10(sum(m(l,:,12)>0)./length(m(l,:,12)))) 0], ...
%     'LineWidth',1,'LineStyle','--')
% hold on
% end
% 
% 
% l=1;
% plot(cumsum(c(l,2:end)')/sum(c(l,2:end)'))
% 
% if size(d,1)==6
% legend({'6 mbar','125 mbar','250 mbar','0.5 bar','1 bar','2 bar','Observations','Zero-atmosphere'},'FontSize',18)
% elseif size(d,1)==8
% legend({'6 mbar','125 mbar','250 mbar','0.5 bar','1 bar','2 bar','3 bar','5 bar','Observations','Zero-atmosphere'},'FontSize',18)
% end
% set(gca,'FontSize',20); xlabel('Diameter (m)'); ylabel('Fraction of craters smaller than')


%RESAMPLE *LINEARLY IN CRATER DIAMETER*:
%colscal = 1/-min(isfinite(log10(sum(m(:,:,12)'>0)./length(m(:,:,12)'))).*log10(sum(m(:,:,12)'>0)./length(m(:,:,12)')));


%deprecated
%   plot([binsize:binsize:3000],[[0 0 0 0 0 0 0 0 0]'; cumsum((fractal_correction_vec.*d(l,11:end))')./sum((fractal_correction_vec.*d(l,11:end))')], ...
%SHOWING MODEL OUTPUT *WITHOUT* THE FRACTAL CORRECTIONs

l=1;
%plot(cumsum((([0:binsize:2999.5].*c(l,2:end))'))./sum(([0:binsize:2999.5].*(c(l,2:end)))'))



%for l =1:size(c,1)
%plot(cumsum(c(l,2:end)')/sum(c(l,2:end)'))
%hold on
%end

show_extra_plots = 0;

if show_extra_plots ==1
    figure
    plot(flipud(cumsum(flipud(c'))))
    set(gca,'yscale','log')
    set(gca,'xscale','log')
    grid on

    figure
    plot(flipud(cumsum(flipud(e'))))
    set(gca,'yscale','log')
    set(gca,'xscale','log')
    grid on
end
% 
% 
% figure
% %plot(flipud(cumsum(flipud(d(:,2:end)')))./sum(d(:,2:end)))
% set(gca,'yscale','log')
% set(gca,'xscale','log')
% grid on

clear m


%atmospheric pressure estimate
smallest_crater = 2*min(cr) 
smallest_lowerbound = 2*min(cl(1,:))
smallest_upperbound = 2*min(cl(2,:))

pressure_estimate = exp(-elevation/10700)*(0.5*(atmospoi(atmospoi==ceil(mean(do(1).bsl))) + atmospoi(atmospoi==ceil(mean(do(1).bsu)))))


if strcmp(data_option,'mawrth_2_2018')==1
    save('dihmawrth2.mat','dih')
    save('distributions_mawrth2.mat','bsu','bsl','bsb')
    save('canddformixed_mawrth2.mat','c','d','smallest_lowerbound','smallest_upperbound')

    save('mawrth_2_data.mat','smallest_crater','pressure_estimate')
    xlswrite('mawrth_2_rawcraters.xlsx',cr')
elseif strcmp(data_option,'mawrth_phyllo_2018')==1
    save('dihmawrth_phyllo.mat','dih')
    save('distributions_mawrth_phyllo.mat','bsu','bsl','bsb')
    save('canddformixed_mawrth_phyllo.mat','c','d','smallest_lowerbound','smallest_upperbound')

    save('mawrth_phyllo_data.mat','smallest_crater','pressure_estimate')
    xlswrite('mawrth_phyllo_rawcraters.xlsx',cr')
elseif strcmp(data_option,'meridiani_2_2018')==1
    save('dihmeridiani_2.mat','dih')
    save('distributions_meridiani_2.mat','bsu','bsl','bsb')
    save('canddformixed_meridiani_2.mat','c','d','cl','cr')

    save('mawrth_meridiani_2_data.mat','smallest_crater','pressure_estimate','smallest_lowerbound','smallest_upperbound')
    xlswrite('meridiani_2_rawcraters.xlsx',cr')
end


% save('distributions.mat','bsu','bsl','bsb')
% save('canddformixed.mat','c','d')
% 
% save('mawrth_phyllo_data.mat','smallest_crater','pressure_estimate')