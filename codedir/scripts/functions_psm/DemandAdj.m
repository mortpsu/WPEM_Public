function out=DemandAdj(IterationCount,demand_chg,region)

    %     demand_chg = demand_chg_zones; % uncomment for testing
    %     region     = 'zones';          % uncomment for testing
    DemandChange = demand_chg(end,:);
    
    %		define current iteration number
    kkkk = IterationCount;
    
    
    if (kkkk == 1)											% iteration = 1
    
        % First time with feedback, usually too large -- explain?
        DemandChangeCGE = DemandChange;
        DemandStep      = DemandChange;
        DemandChange(DemandChange < -0.01) = -0.01;
        DemandChangeADJ = DemandChange;
    
    else											        % iteration > 1
        
        %%%% 		step 1: load previous iteration percentage change in electricity demand (from rem)
        if (strcmp(region,'zones')==1)
            load('DemandChangeData_zones.mat', 'DemandChangeCGE', 'DemandChangeADJ');
            DemandChangeCGE = [DemandChangeCGE; DemandChange];
        else
            load('DemandChangeData_econr.mat', 'DemandChangeCGE', 'DemandChangeADJ');
            DemandChangeCGE = [DemandChangeCGE; DemandChange];
        end
        
        %%%% 		step 2: check if converged is reached across each region
        if ((DemandChangeCGE(kkkk-1,:) == DemandChange(:)))
        
            % DO NOT ADJUST -> CONVERGED           :END
        
        else % ITERATION >= 2 AND NOT YET CONVERGED
        
            %%%		step2a: if not converged, define the difference between iterations
            DemandStep = DemandChangeADJ(kkkk-1,:) - DemandChange;
            
            %%%     step2b: loop through regions and adjust them depending
            %%%             on the size of demand change
            for ii = 1:length(DemandChange)
            
                %%% If econ model wants small adjustment, just pass through
                if (abs(DemandChange(ii)) < 0.001)
                
                    % Do nothing, leave this change as is
                    
                    %%% Is the new demand change and the previous one in the same
                    %%% direction or opposite?
                elseif (sign(DemandChangeADJ(kkkk-1,ii)) == sign(DemandChange(ii)))
                
                    if (abs(DemandStep(ii)) < 0.001)     % do region ii
                    
                        % DO NOT ADJUST -> VERY SMALL CHANGE FROM LAST ITERATION
                    
                    elseif (abs(DemandStep(ii)) <= 0.01) % small difference - split the difference
                        DemandStep(ii) = (DemandStep(ii))/2;
                    elseif ( DemandStep(ii) >  0.3)      % CGE wants a bigger demand reduction, take a step
                        DemandStep(ii) = 0.05;
                    elseif ( DemandStep(ii) >  0.2)      % CGE wants a bigger demand reduction, take a step
                        DemandStep(ii) = 0.03;
                    elseif ( DemandStep(ii) >  0.1)      % CGE wants a bigger demand reduction, take a step
                        DemandStep(ii) = 0.02;
                    elseif ( DemandStep(ii) >  0.01)     % CGE wants a bigger demand reduction, take a step
                        DemandStep(ii) = 0.01;
                    elseif ( DemandStep(ii) <  -0.01)    % CGE wants a bigger demand increase, take a step
                        DemandStep(ii) = -0.01;
                    else                                % moderate change take a small step
                        DemandStep(ii) = DemandStep(ii) / 4;
                    end
                    
                    %		step2b -- define the new pchg in demand based on step
                    DemandChange(ii) = DemandChangeADJ(kkkk-1,ii) - DemandStep(ii);
                    
                    
                    % 		step2c -- check if stuck cycling -- explain?
                    if     ((kkkk > 4) && (DemandChange(ii) == DemandChangeADJ(kkkk-2,ii)) && (DemandChangeADJ(kkkk-1,ii) == DemandChangeADJ(kkkk-3,ii)))
                        DemandChange(ii) = (DemandChangeADJ(kkkk-2,ii) + DemandChangeADJ(kkkk-1,ii))/2;
                    elseif ((kkkk > 10) && (DemandChange(ii) < DemandChangeADJ(kkkk-1,ii)) && (DemandChangeADJ(kkkk-1,ii) > DemandChangeADJ(kkkk-2,ii)) && (DemandChangeADJ(kkkk-2,ii) < DemandChangeADJ(kkkk-3,ii)))
                        DemandChange(ii) = (DemandChangeADJ(kkkk-2,ii) + DemandChangeADJ(kkkk-1,ii))/2;
                    elseif ((kkkk > 10) && (DemandChange(ii) > DemandChangeADJ(kkkk-1,ii)) && (DemandChangeADJ(kkkk-1,ii) < DemandChangeADJ(kkkk-2,ii)) && (DemandChangeADJ(kkkk-2,ii) > DemandChangeADJ(kkkk-3,ii)))
                        DemandChange(ii) = (DemandChangeADJ(kkkk-2,ii) + DemandChangeADJ(kkkk-1,ii))/2;
                    end
                else
                    %%% demand change has opposite signs between iterations
                    
                    if (DemandChange(ii) < -0.01)    
                    
                        DemandChange(ii) = -0.01;   % take a small step in this new direction
                    
                    elseif (DemandChange(ii) > 0.01)  
                    
                        DemandChange(ii) = 0.01;   % take a small step in this new direction
                    else
                        % new demand
                        % DO NOT ADJUST -> VERY SMALL CHANGE
                    end
                end % sign statement
            end % for ii loop
        end % converged
        
        % 		step 3: save requested and adjusted demand changes
        DemandChangeADJ = [DemandChangeADJ; DemandChange];
    
    end %iteration = 1
    
    %		display adjusted percentage change in electricity demand (from rem)
    disp('demand change - Adjusted');
    disp(DemandChange);
    
    %		save percentage change in electricity demand as Matlab data
    if (strcmp(region,'zones')==1)
    save('DemandChangeData_zones.mat', 'DemandChangeCGE', 'DemandChangeADJ');
    else
    save('DemandChangeData_econr.mat', 'DemandChangeCGE', 'DemandChangeADJ');
    end
    
    out = DemandChange;
end