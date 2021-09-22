%%% Producing the matrix of parameters (rho and c) that satisfy the
%%% epidemiological criteria for given values of beta_F and beta_E, but 
%%% with the survivor class of individuals having the potential of being 
%%% infectious. This is then used to produce Fig A.21 in the Appendix of 
%%% Chapter 2.

%%% The parameter sweep includes starting at different months of the year,
%%% as an ASF outbreak can occur at any time. We start from a time of 0,
%%% 0.3 and 0.6 years to cover potential starting times throughout the
%%% year.

%%% Uses the function ASFEstoniaINFCHR, which is the model representative 
%%% for Estonia but with variable parameter values, infectious survivors
%%% and matching with two types of survivor discussed in Stahl et al.
%%% (2019) - citation in main paper.

clearvars
options = odeset('Refine',1, 'RelTol',1e-6,'AbsTol',1e-6,'NonNegative',1:11);
        %ode tolerances and restrictions

%% Parameter Sweep starting in 0th month (January)
var = [];
psvar = [];
Tstart = 0;

for rho1 = 15:-1:0
    rho = rho1/100;
    for c1 = 0:12
        c = c1/40;

        %SOLVING THE SYSTEM OF ODES
        y0 = [0.8, 0.1, 0, 0, 0, 1, 0.1, 0, 0, 0, 0];
        [t,y] = ode45(@(t,y)  ASFModelEstoniaINFCHR(t, y, 0, bF, bE, rho, c), [Tstart 10+Tstart], y0, options);
        sus = y(:,1) + y(:,6);  inf = y(:,2) + y(:,4) + y(:,7) + y(:,9);
        sur = y(:,3) + y(:,8);  rec = y(:,5) + y(:,10);
        totalpop = sus + inf + sur + rec;   
        prevalence = inf./totalpop;
        pseudoprev = (inf + sur)./totalpop;

        if min(totalpop) < 0.05
            break
        end

        %CRITERIA

        %Average Total Population Per Year
        avgy = zeros(10,1);
        a = 1;
        for i = 1:10
            b = find(t >= i+Tstart);
            b = b(1);
            avgy(i) = mean(totalpop(a:b));
            a = b;
        end

        if min(avgy) < 0.1  %average cannot go below 95%
            break
        end

        I = find(inf == max(inf)); %finding peak infected time
        Itime = t(I(1));
        D = find(y(:,11) >= 0.02); %finding first noted infected case
        if length(D) >= 1 %virus must be noticed
            Dtime = t(D(1));
            if Itime < (Dtime + 10/12) && Itime > (Dtime + 4/12) && Dtime < 2%Peak 6 months
                Crash = find(t > Itime + 1); %Crash in Pop within year of epidemic
                CrDL = Crash(1);
                if min(totalpop(1:CrDL)) < 0.3 %85%- 95% drop in pop within year of epidemic
                    t24index = find(t >= 2 + Itime);
                    t24index = t24index(1);
                    t38index = find(t <= 38/12 + Itime);
                    t38index = t38index(end);
                    n = 0;
                    for i=t24index:t38index
                        if prevalence(i) > 0.01 && prevalence(i) < 0.03 && min(prevalence(I:i))> 0.005 %within 1-3% after 30-44 months                     
                            var = [var; bF, bE, rho, c];
                            break
                        end
                        if pseudoprev(i) > 0.01 && pseudoprev(i) < 0.03 && min(pseudoprev(I:i))> 0.005 
                            psvar = [psvar; bF, bE, rho, c];
                            break
                        end
                    end
                end
            end
        end
        %END OF CRITERIA
    end
end

save('var0')

%% Parameter Sweep starting in 0.3th year (March/April)
clearvars
var = [];
psvar = [];
Tstart = 0.3;

for rho1 = 15:-1:0
    rho = rho1/100;
    for c1 = 0:12
        c = c1/40;

        %SOLVING THE SYSTEM OF ODES
        y0 = [0.8, 0.1, 0, 0, 0, 1, 0.1, 0, 0, 0, 0];
        [t,y] = ode45(@(t,y)  ASFModelEstoniaINFCHR(t, y, 0, bF, bE, rho, c), [Tstart 10+Tstart], y0, options);
        sus = y(:,1) + y(:,6);  inf = y(:,2) + y(:,4) + y(:,7) + y(:,9);
        sur = y(:,3) + y(:,8);  rec = y(:,5) + y(:,10);
        totalpop = sus + inf + sur + rec;   
        prevalence = inf./totalpop;
        pseudoprev = (inf + sur)./totalpop;

        if min(totalpop) < 0.05
            break
        end

        %CRITERIA

        %Average Total Population Per Year
        avgy = zeros(10,1);
        a = 1;
        for i = 1:10
            b = find(t >= i+Tstart);
            b = b(1);
            avgy(i) = mean(totalpop(a:b));
            a = b;
        end

        if min(avgy) < 0.1  %average cannot go below 95%
            break
        end

        I = find(inf == max(inf)); %finding peak infected time
        Itime = t(I(1));
        D = find(y(:,11) >= 0.02); %finding first noted infected case
        if length(D) >= 1 %virus must be noticed
            Dtime = t(D(1));
            if Itime < (Dtime + 10/12) && Itime > (Dtime + 4/12) && Dtime < 2%Peak 6 months
                Crash = find(t > Itime + 1); %Crash in Pop within year of epidemic
                CrDL = Crash(1);
                if min(totalpop(1:CrDL)) < 0.3 %85%- 95% drop in pop within year of epidemic
                    t24index = find(t >= 2 + Itime);
                    t24index = t24index(1);
                    t38index = find(t <= 38/12 + Itime);
                    t38index = t38index(end);
                    n = 0;
                    for i=t24index:t38index
                        if prevalence(i) > 0.01 && prevalence(i) < 0.03 && min(prevalence(I:i))> 0.005 %within 1-3% after 30-44 months                     
                            var = [var; bF, bE, rho, c];
                            break
                        end
                        if pseudoprev(i) > 0.01 && pseudoprev(i) < 0.03 && min(pseudoprev(I:i))> 0.005 
                            psvar = [psvar; bF, bE, rho, c];
                            break
                        end
                    end
                end
            end
        end
        %END OF CRITERIA
    end
end

save('var03')

%% Parameter Sweep starting in 0.6th year (~July)
clearvars
var = [];
psvar = [];
Tstart = 0.6;

for rho1 = 15:-1:0
    rho = rho1/100;
    for c1 = 0:12
        c = c1/40;

        %SOLVING THE SYSTEM OF ODES
        y0 = [0.8, 0.1, 0, 0, 0, 1, 0.1, 0, 0, 0, 0];
        [t,y] = ode45(@(t,y)  ASFModelEstoniaINFCHR(t, y, 0, bF, bE, rho, c), [Tstart 10+Tstart], y0, options);
        sus = y(:,1) + y(:,6);  inf = y(:,2) + y(:,4) + y(:,7) + y(:,9);
        sur = y(:,3) + y(:,8);  rec = y(:,5) + y(:,10);
        totalpop = sus + inf + sur + rec;   
        prevalence = inf./totalpop;
        pseudoprev = (inf + sur)./totalpop;

        if min(totalpop) < 0.05
            break
        end

        %CRITERIA

        %Average Total Population Per Year
        avgy = zeros(10,1);
        a = 1;
        for i = 1:10
            b = find(t >= i+Tstart);
            b = b(1);
            avgy(i) = mean(totalpop(a:b));
            a = b;
        end

        if min(avgy) < 0.1  %average cannot go below 95%
            break
        end

        I = find(inf == max(inf)); %finding peak infected time
        Itime = t(I(1));
        D = find(y(:,11) >= 0.02); %finding first noted infected case
        if length(D) >= 1 %virus must be noticed
            Dtime = t(D(1));
            if Itime < (Dtime + 10/12) && Itime > (Dtime + 4/12) && Dtime < 2%Peak 6 months
                Crash = find(t > Itime + 1); %Crash in Pop within year of epidemic
                CrDL = Crash(1);
                if min(totalpop(1:CrDL)) < 0.3 %85%- 95% drop in pop within year of epidemic
                    t24index = find(t >= 2 + Itime);
                    t24index = t24index(1);
                    t38index = find(t <= 38/12 + Itime);
                    t38index = t38index(end);
                    n = 0;
                    for i=t24index:t38index
                        if prevalence(i) > 0.01 && prevalence(i) < 0.03 && min(prevalence(I:i))> 0.005 %within 1-3% after 30-44 months                     
                            var = [var; bF, bE, rho, c];
                            break
                        end
                        if pseudoprev(i) > 0.01 && pseudoprev(i) < 0.03 && min(pseudoprev(I:i))> 0.005 
                            psvar = [psvar; bF, bE, rho, c];
                            break
                        end
                    end
                end
            end
        end
        %END OF CRITERIA
    end
end

save('var06')