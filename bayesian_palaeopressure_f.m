

%BAYESIAN_PALEOPRESSURE_FOR_SASHA -> MODIFIED FOR MY OWN SINISTER FITTING
%PURPOSES 
%Edwin Kite -> EDIT BY SASHA

%not intended as a standalone script

%purpose: (1) construct a grid of the probability of getting N craters (the observed N craters) 
%in a bin of dR,
%assuming poisson statistics, as a function of Pressure and dummy rate/flux
%parameter F
%
%note: the 'accumulation rate' really is a dummy parameter because the real
%world flux associated with it will change as we change the atmspheric
%pressure
%
%(2) convolve this grid with the observations - 
%(3) output the best fits and the errors.
%this should be equivalent to a bootstrap. the bootstrap is just doing
%poisson statistics.


%take table of model predicted frequencies
%reweight and normalize for the fractal correction (NO NEED. THIS HAS BEEN DONE DURING CONSTRUCTION OF "ZI" - SEE "tll" VARIABLE)

clear ptot plog
load data_option
%Load all the data
    if strcmp(data_option,'mawrth_2_2018')==1
        clear dd mc mc_old
        dd(1,:,:) = xlsread('dch_mawrth2_pass1.xlsx'); %diameter data
        dd(:,3,:) = dd(:,1,:) + dd(:,2,:); %SW: 12/14/18 sum of DEF + CAND 
    else
        %error('Only new data!');
    end
    
    if strcmp(data_option,'mawrth_phyllo_2018')==1
        clear dd mc mc_old
        dd(1,:,:) = xlsread('dch_mawrth_phyllo.xlsx'); %diameter data
        dd(:,3,:) = dd(:,1,:) + dd(:,2,:); %SW: 12/14/18 sum of DEF + CAND 
    else
        %error('Only new data!')
    end
    
    if strcmp(data_option,'meridiani_2_2018')==1
        clear dd mc mc_old
        dd(1,:,:) = xlsread('dch_meridiani_2.xlsx'); %diameter data
        dd(:,3,:) = dd(:,1,:) + dd(:,2,:); %SW: 12/14/18 sum of DEF + CAND 
    else
        %error('Only new data!');
    end

%QUICK AND DIRTY WORKAROUND FOR HOLSAPPLE K1 ISSUE (SEE EMAILS FROM
%KEITH HOLSAPPLE (U. WASHINGTON)
%Commenting this out as we will do a better workaround in quickcheckpressureeffects_for_sasha.m    
%     disp('Warning! Do not double-count the Holsapple K1 correction!')
%     
%     diameter_correction_factor = (0.24/0.132)^(1/3);
%     for dind = 600:-1:2
%         dd(:,:,dind) = dd(:,:,round(dind/diameter_correction_factor));
%     end
%END QUICK AND DIRTY WORKAROUND
    
   
    %LOAD INDIVIDUAL MIXED P FILES HERE:
    

    
    model = ZI;
    warning('To minimize effects of model shot noise, derivative of the model is smoothed. See notebook Mars Climate 5 for rationale for the choice of the smoothing parameter.')
    
    %Get rid of data for very large crater sizes:
    dd(:,:,601:end) = []; model(:,1201:end) = [];
    %disp('inpainting ...')
    %model = 10.^inpaint_nans(log10(model),4);
    %disp('... done.')
    
    %Fill in the model for small values by finding the log-log slope of the
    %cdf
    %at 0.01 and using it to extend the model to very small values.
    if sum(sum(isnan(model(:,1:100))))>0 
        for l = 1:size(model,1)
            if min(find(model(l,:)>1e-10)) >1
                esr = min(find(model(l,:)>1e-2)); %end-slope root
                es = (log10(model(l,esr+5))-log10(model(l,esr))) ...
                    /(log10(esr+5) - log10(esr));%end slope
                fnz = min(find(model(l,:)>1e-10));%first non zero
                dr = log10(fnz) - log10([1:(-1+fnz)]);%meters
                df = dr.*es;
                nv = 10.^(log10(model(l,fnz)) - df);%new values
                model(l,1:(-1+fnz)) = nv;
            end
        end
    end
    %renormalize the model.
    
    %Note: I have confirmed that the histogram bin edges for data and model
    %are aligned (modulo the factor of 2). Edwin Kite 23 Feb 2013.
    
    %experiments with fflist suggest that it makes a truly negligible, tiny difference to
    %the result (though not zero). Therefore I have omitted it in favor of
    %just 1 flux value.
    fflist = 1;%logspace(-0.5,.5,11);%[1/1.5 1/1.1 1 1.1 1.5] %ones(1000);logspace(1e-3,1e3,30);
    nfluxes = max(size(fflist));
    
% psuedo for case dtmall dtm1 dtm2
            for ip = 1:size(ZI,1) %index of Pressure
                %disp(ip)
                mc(ip,:) = [model(ip,1), diff(squeeze(model(ip,:)))]; %Model for this Case; ... diff'd. Prepending the first entry to keep thins even.
                %mc(ip,:) = mc(ip,:)
                for mm = 1:size(mc,2) %Warning: this only works when the zeros are at parts of the cdf that are not very important to the overall fit (e.g., very big craters).
                    if mc(ip,mm)<1e-35;%flood-forward;
                        lastflood = -2+mm+min(find(mc(ip,mm:end)>=1e-35));
                        if isempty(lastflood)==1;lastflood=size(mc,2);end
                        sharingbins = lastflood-mm+2;
                        mc(ip,(mm-1:lastflood)) = mc(ip,mm-1)./sharingbins;
                    end
                end
                
                %I have looked at the mc histogram and verified that 1e-7
                %is seperated from 'real' data by ~3 orders of magnitude
                mc(ip,mc(ip,:)<1e-35) = 1e-35; %e-6; mc(ip,min(find(mc(ip,:)>0)));%diff rounding errors trip problems later on
            
                %Added October 2013 as a workaround for the problems with
                %Bayesian code not itting DEF crater curve correctly
                mc(ip,:) = smooth(mc(ip,:),50);
            end

for dtmcase = 1 %:3

        % pseudo here goeth case specific data extracts
        for data_case = [3,1]%,4] %1 = def only ('conservative'). 3 = def + candidates ('preferred')
             disp(data_case)
             
             dc = squeeze(dd(dtmcase,data_case,:))'; %Data for this Case
                         disp(nansum(dc))

% psuedo for each data type (conservative and preferred)
% pseudo so now we have a data vector
% 
% for each Pressure

                scalefactor = (nansum(dc)./nansum(mc(ip,:))); % ... and scaled;%note summation
            %clear ptot plog ptotp poisson 
            for ip = 1:size(ZI,1) %index of Pressure
                for ffi = 1:nfluxes;%; %flux factor index
                    ff = fflist(ffi);
                    di = linspace(1,size(dc,2),(size(dc,2))); %data index
                    %  what is the probability of Data given Model?
                    k = dc;%(di); %actual number of craters per bin
                    lambda = ff.*(mc(ip,di*2) + mc(ip,di*2-1)).*scalefactor;
                    poissonlog = (log10(lambda.^k.* exp(-lambda)) - log10(factorial(k))); %Sivia 2nd edition equation 5.49.             
                    poissonlog(isinf(poissonlog)) = -NaN;
                    
                    % sum the logs of the probabilities of each Bin (gives p(Data Vector | Model)
                    plog(ip,ffi,:) = poissonlog;
                    %<DEPRECATED> very small crater sizes are NOT included
                    %(unobservable; threshholded in model in some cases)
                    %ptot(ip,ffi)   = (nansum(poissonlog(11:300)));%350)));%end))); %overall LOG10 probability (sum the logs)
                    
                    %Include all crater sizes:
                    ptot(ip,ffi)   = (nansum(poissonlog(1:500)));%was 300 %350)));%end))); %overall LOG10 probability (sum the logs)
                
                end
            end

            %Convolve with log-flat prior in interval (10 mbar --> 5
            %bars].
            pressures = 10:10:5000;
            prior = diff(log10([1:10:5001]));
            %prior = repmat(prior,nfluxes,1)';
            
            arb_shift = -max(max(ptot));
            poisson = 10.^(ptot+arb_shift);
            if size(poisson,2)>1
                ptotp = sum(poisson').*prior; %marginalizing on flux factor (ffi)
            else
                ptotp = poisson'.*prior;
            end
            
                    ptotp(isnan(ptotp)==1)=0;
                    ptotp = cumsum(ptotp)./sum(ptotp);
                    


           %    model(datatype,timearound,aa,tshr,:) = poisson;
            %   lambdam(datatype,timearound,aa,tshr,:) = lambda;
               [BP,IP,JP] = unique(ptotp);
%95.4499736%	4.5500

               lowpress =  interp1(BP,(pressures(IP)),0.0455);%  0.5*0.3173); %log10(1 sigma low accumulation rate)
               highpress = interp1(BP,(pressures(IP)),0.9545);%1-0.5*0.3173); %log10(1 sigma high accumulation rate)
               lowpress1s =  interp1(BP,(pressures(IP)),0.5*0.3173); %log10(1 sigma low accumulation rate)
               highpress1s = interp1(BP,(pressures(IP)),1-0.5*0.3173); %log10(1 sigma high accumulation rate)
               medianp_t = interp1(BP,(pressures(IP)),0.5);
               
             %  accumulation rates:
               errorlow(dtmcase,data_case) = lowpress;
              errorhigh(dtmcase,data_case) = highpress;
             
               errorlow1s(dtmcase,data_case) = lowpress1s;
              errorhigh1s(dtmcase,data_case) = highpress1s;
            bayesmedian(dtmcase,data_case) = 10.^interp1(BP,log10(pressures(IP)),         0.5);
                             [null_var, I] = max(diff(ptotp));
                  ptp(dtmcase,data_case,:) = ptotp;
                bestfit(dtmcase,data_case) = pressures(I);
                medianp(dtmcase,data_case) = medianp_t;

        end
end


    
% end
