%%% Producing the matrix of transmission coefficients that satisfy the
%%% epidemiological criteria for given values of rho. This is then used to
%%% produce Fig A.2 in the Appendix of Chapter 2.

%%% The parameter sweep includes starting at different months of the year,
%%% as an ASF outbreak can occur at any time. We start from a time of 0,
%%% 0.3 and 0.6 years to cover potential starting times throughout the
%%% year.

%%% Uses the function ASFEstoniaCRTR, which is the model representative for
%%% Estonia but with variable parameter values, to test the criteria

clearvars
options = odeset('Refine',1, 'RelTol',1e-8,'AbsTol',1e-9,'NonNegative',1:8);

%In Fig A.2, three values for the parameter rho were used, rho = 0.83, 0.85
%and 0.9. This script will need to be run for each value of rho to obtain
%the relevant data matrices.

%%% After running through this script for the desired value of rho, the 
%%% valid parameters for each starting time can be combined using the
%%% script titled 'Trans_Coeff_Combine'

rho = 0.83;


%% Starting simulations from January
var = [];
Tstart = 0;

for b_F = 610:680
    bF = b_F/10;
    minprev = 0;
    for b_E = 10:1:24
        bE = b_E/10;
        %SOLVING THE SYSTEM OF ODES
        y0 = [0.796, 1.2, 0.002, 0.002, 0, 0, 0, 0];
        [t,y] = ode45(@(t,y)  ASFModelEstoniaCRTR(t, y, 0, bF, bE, rho), [Tstart 10+Tstart], y0, options);
        totalpop = y(:, 1) + y(:, 2)+ y(:, 3) + y(:, 4) + y(:, 5) + y(:, 6);
        prevalence = (y(:, 3) + y(:, 4))./totalpop;
        infected = y(:,3) + y(:,4);

        if min(totalpop) < 0.05     %If totalpop has dropped too low, criteria not matched
            break                   %save time and break loop
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

        if min(avgy) < 0.1                      %average cannot go below 95%
            break
        end

        I = find(infected == max(infected));    %finding peak infected time
        Itime = t(I(1));
        D = find(y(:,7) >= 0.02);               %finding first noted infected case
        if length(D) >= 1                       %virus must be noticed
            Dtime = t(D(1));
            if Itime < (Dtime + 10/12) && Itime > (Dtime + 4/12) && Dtime < 2   %Peak 6 months
                Crash = find(t > Itime + 1);                                    %Crash in Pop within year of epidemic
                CrDL = Crash(1);
                if min(totalpop(1:CrDL)) < 0.3 && min(avgy) > 0.1               %85%- 95% drop in pop within year of epidemic
                    t24index = find(t >= 2 + Itime);
                    t24index = t24index(1);
                    t38index = find(t <= 38/12 + Itime);
                    t38index = t38index(end);
                    for i=t24index:t38index
                        if prevalence(i) > 0.01 && prevalence(i) < 0.03 && min(prevalence(I:i))> 0.005      %within 1-3% after 30-44 months                     
                            var = [var; bF, bE];
                            break
                        end
                    end
                end
            end
        end
        %END OF CRITERIA
    end
end
save('var0year')

%% Starting from 3/10ths into the year
var = [];
Tstart = 0.3;
for b_F = 610:680
    bF = b_F/10;
    minprev = 0;
    for b_E = 10:1:24
        bE = b_E/10;
        %SOLVING THE SYSTEM OF ODES
        y0 = [0.796, 1.2, 0.002, 0.002, 0, 0, 0, 0];
        [t,y] = ode45(@(t,y)  ASFModelEstoniaCRTR(t, y, 0, bF, bE, rho), [Tstart 10+Tstart], y0, options);
        totalpop = y(:, 1) + y(:, 2)+ y(:, 3) + y(:, 4) + y(:, 5) + y(:, 6);
        prevalence = (y(:, 3) + y(:, 4))./totalpop;
        infected = y(:,3) + y(:,4);

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

        if min(avgy) < 0.1                      %average cannot go below 95%
            break
        end

        I = find(infected == max(infected));    %finding peak infected time
        Itime = t(I(1));
        D = find(y(:,7) >= 0.02);               %finding first noted infected case
        if length(D) >= 1                       %virus must be noticed
            Dtime = t(D(1));
            if Itime < (Dtime + 10/12) && Itime > (Dtime + 4/12) && Dtime < 2   %Peak 6 months
                Crash = find(t > Itime + 1);                                    %Crash in Pop within year of epidemic
                CrDL = Crash(1);
                if min(totalpop(1:CrDL)) < 0.3 && min(avgy) > 0.1               %85%- 95% drop in pop within year of epidemic
                    t24index = find(t >= 2 + Itime);
                    t24index = t24index(1);
                    t38index = find(t <= 38/12 + Itime);
                    t38index = t38index(end);
                    for i=t24index:t38index
                        if prevalence(i) > 0.01 && prevalence(i) < 0.03 && min(prevalence(I:i))> 0.005      %within 1-3% after 30-44 months                     
                            var = [var; bF, bE];
                            break
                        end
                    end
                end
            end
        end
        %END OF CRITERIA
    end
end

save('var03rdyear')

%% Starting 6/10th into the year

var = [];
Tstart = 0.6;

for b_F = 610:680
    bF = b_F/10;
    minprev = 0;
    for b_E = 10:1:24
        bE = b_E/10;
        %SOLVING THE SYSTEM OF ODES
        y0 = [0.796, 1.2, 0.002, 0.002, 0, 0, 0, 0];
        [t,y] = ode45(@(t,y)  ASFModelEstoniaCRTR(t, y, 0, bF, bE, rho), [Tstart 10+Tstart], y0, options);
        totalpop = y(:, 1) + y(:, 2)+ y(:, 3) + y(:, 4) + y(:, 5) + y(:, 6);
        prevalence = (y(:, 3) + y(:, 4))./totalpop;
        infected = y(:,3) + y(:,4);

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

        if min(avgy) < 0.1                      %average cannot go below 95%
            break
        end

        I = find(infected == max(infected));    %finding peak infected time
        Itime = t(I(1));
        D = find(y(:,7) >= 0.02);               %finding first noted infected case
        if length(D) >= 1                       %virus must be noticed
            Dtime = t(D(1));
            if Itime < (Dtime + 10/12) && Itime > (Dtime + 4/12) && Dtime < 2   %Peak 6 months
                Crash = find(t > Itime + 1);                                    %Crash in Pop within year of epidemic
                CrDL = Crash(1);
                if min(totalpop(1:CrDL)) < 0.3 && min(avgy) > 0.1               %85%- 95% drop in pop within year of epidemic
                    t24index = find(t >= 2 + Itime);
                    t24index = t24index(1);
                    t38index = find(t <= 38/12 + Itime);
                    t38index = t38index(end);
                    for i=t24index:t38index
                        if prevalence(i) > 0.01 && prevalence(i) < 0.03 && min(prevalence(I:i))> 0.005      %within 1-3% after 30-44 months                     
                            var = [var; bF, bE];
                            break
                        end
                    end
                end
            end
        end
        %END OF CRITERIA
    end
end
save('var06thyear')