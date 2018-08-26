function predict_postseason

    % Predicts baseball postseason match results based on regular-season
    % Elo ratings, then compares predictions to actual outcomes. Intended
    % use is to tune Elo parameters for maximum predictive accuracy
    
    % Daniel W. Dichter
    % 2018-08-26
    
    %% Constants

    Elo_mean = 1000; % leave Elo_mean = 1000 as constant
        
    %% Inputs
    
    r_p = 2008:2017; % range, predict, for which to predict postseason outcomes
    
% Inputs, Elo:
%   [K, k, h_f_a]
%   [0.08, 11, 0.03]
% Outputs:
%   Correctly predicted: 56 of 82
%   Accuracy: 68.2927%

% Inputs, Elo:
%   [K, k, h_f_a]
%   [0.075, 11.0714, 0.017857]
% Outputs:
%   Correctly predicted: 211 of 301
%   Accuracy: 70.0997%
    
    % For coarse tuning
%     rln = 5; % resolution
%     Elo_K = logspace(-1, 3, rln);
%     Elo_k = logspace(-1, 3, rln);
%     h_f_a = linspace(0.01,0.05,rln);
    
    % For fine tuning
%     var = 0.25; % percentage variation
%     rln = 8; % resolution
%     guess = [0.1, 10, 0.02];
%     Elo_K = linspace(guess(1)*(1-var), guess(1)*(1+var), rln);
%     Elo_k = linspace(guess(2)*(1-var), guess(2)*(1+var), rln);
%     h_f_a = linspace(guess(3)*(1-var), guess(3)*(1+var), rln);
    
    % For single cases
    guess = [0.08, 11, 0.03];
    Elo_K = guess(1);
    Elo_k = guess(2);
    h_f_a = guess(3);
    
    %% Load and initialize
    
    clc
    load('live_ball_game_log_n','G_n','T') % load games
    G_o = G_n; % original data, read-only
    
    case_counter = 1; % case counter
    case_quantity = length(Elo_K)*length(Elo_k)*length(h_f_a); % case quantity
    
    % Display a header
    switch(length(r_p))
        case 1
            disp(['Predicting postseason ' num2str(r_p)])
        otherwise
            disp(['Predicting postseasons ' num2str(min(r_p)) '-' num2str(max(r_p))])
            disp(['Elo_mean: ' num2str(Elo_mean) ' (constant)'])
            disp(['Cases requested: ' num2str(case_quantity)])
            disp(['Estimated runtime: ' num2str(round(0.81*case_quantity/60)) ' min']) % work laptop
    end
    
    %% Iterate through all Elo configurations, assessing accuracy over requested time domain
    accuracy_best = 0;
    for E_K = Elo_K
        for E_k = Elo_k
            for hfa = h_f_a

                tic

                [R, dn] = calculate_elo(E_K,E_k,hfa,0);

                pred_correct    = 0;
                pred_attempted  = 0;

                %% Predict postseason using Elo
                
                for year = r_p

                    G = G_o(G_o(:,1) == year,:);
                    G = G(G(:,9)==2,:); % extract postseason games only

                    % Crop games that do not represent match (series) outcomes
                    T_ = unique(G(:,[5,7])); % teams present in playoffs
                    s_o = zeros(size(G,1),1); % series outcomes
                    for t = 1:length(T_)
                        f_a = max([find(G(:,5)==T_(t)); find(G(:,7)==T_(t))]); % final appearance
                        s_o(f_a) = 1;
                    end
                    G = G(s_o==1,:);

                    if isempty(G)
                        warning('No playoff games found')
                        continue
                    else
            %             disp(['   Year: ' num2str(year)])
            %             disp(['   Games: ' num2str(size(G,1))])
                    end

                    r_i = find(G(1,4)-1 < dn, 1, 'first')-1; % index representing ratings at end of regular season

                    for g = 1:size(G,1) % for each postseason game

                        R_c = [ % ratings, end of regular season
                                R(r_i,G(g,5)),...
                                R(r_i,G(g,7)),...
                              ];

                        % Game outcome format - 1: away win, 2: home win

                        % Generate outcome prediction
                        if R_c(1) >  R_c(2)
                            pred = 1;
                        else
                            pred = 2;
                        end
                        pred_attempted = pred_attempted+1;

                        % Determine actual outcome
                        if G(g,6)>G(g,8)
                            res = 1;
                        else
                            res = 2;
                        end

                        if pred == res % if prediction was correct
                            pred_correct = pred_correct+1;
                        end

                    end
                end
                
                accuracy = pred_correct/pred_attempted;
                
                if accuracy > accuracy_best
                    accuracy_best = accuracy;
                    disp('================================')
                    disp(['Case ' num2str(case_counter) ', ' num2str(floor(case_counter/case_quantity*100)) '% complete'])
                    disp('Inputs, Elo:')
                    disp('  [K, k, h_f_a]')
                    disp(['  [' num2str(E_K) ', ' num2str(E_k) ', ' num2str(hfa) ']'])
                    disp('Outputs:')
                    disp(['  Correctly predicted: ' num2str(pred_correct) ' of ' num2str(pred_attempted) ]) 
                    disp(['  Accuracy: ' num2str(accuracy*100) '%'])
                else
                    disp(['Case ' num2str(case_counter) ', ' num2str(floor(case_counter/case_quantity*100)) '% complete'])
                end

                case_counter = case_counter+1;

%                 toc
            end
        end
    end

end












































