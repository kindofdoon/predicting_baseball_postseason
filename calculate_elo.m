function [R, dn] = calculate_elo(Elo_K,Elo_k,h_f_a,show_plot)

    % Calculates Elo ratings for all MLB teams over time from 1920 (start of live-ball era) to present
    
    % Returns:
        % R: rating matrix; each row represents a day, and each column
            % represents a team; teams are pulled from live_ball_game_log_n
        % dn: datenum array providing a timestamp for each row of R

    % Daniel W. Dichter
    % 2018-08-26

    %% Cosntants
    
    Elo_mean = 1000; % leave Elo_mean = 1000 as constant

    %% Inputs
    
    plot_show_top = 2;
    plot_show_bot = 2;
    
    plot_highlight_quantity = 2;
	include_postseason = 1; % set to 0 to crop postseason games; 1 to include them
    plot_date_range = {'2018-Jul-01','2018-Jul-31'}; % show range, in yyyy-mmm-dd format
    
    % Teams to explicitly highlight in plot
    show_teams = {
                    'New York Yankees'
                    'Boston Red Sox'
                 };

    %% Setup

    load('live_ball_game_log_n','G_n','T')
    G = G_n;

    if ~include_postseason
        G(G(:,9)==2,:) = []; % crop post-season games
    end

    %% Display the input conditions

    if isempty(G)
        error('Elo calculation: no games match specified input conditions')
    else
%         disp('Input games:')
%         disp(['   include_postseason: ' num2str(include_postseason)])
%         disp(['   Years: ' num2str(min(G(:,1))) '-' num2str(max(G(:,1)))])
%         disp(['   Games: ' num2str(size(G,1))])
%         disp(['   Regular: ' num2str(length(find(G(:,9)==1)))])
%         disp(['   Post: ' num2str(length(find(G(:,9)==2)))])
% 
%         disp('Elo calculation input:')
%         disp(['   Mean: ' num2str(Elo_mean)])
%         disp(['   K: ' num2str(Elo_K)])
    end

    %% Calculate game-by-game Elo

    R = zeros(size(G,1)+1,length(T))-1; % Elo matrix, "R" for ratings
    R(1,:) = Elo_mean;

    for g = 1:size(G,1)-1 % for each game played

        g = g+1; % offset to account for pre-season initialization

        i_t = [G(g,5), G(g,7)]; % index, teams

        % Elo calculation: 
            % https://en.wikipedia.org/wiki/Elo_rating_system#Mathematical_details
            % https://en.wikipedia.org/wiki/Logistic_function

        % Expected outcomes
        E(1) = 1 / (1+exp(-Elo_k*(R(g-1,i_t(1))-R(g-1,i_t(2)) )));
        E(2) = 1 - E(1);

        E = E + [-h_f_a, h_f_a]; % apply home field advantage
        
        if G(g,6) > G(g,8) % if away team wins
            S = [1 0];
        else
            S = [0 1];
        end

        % Update scores
        for i = 1:2
            R(g,i_t(i)) = R(g-1,i_t(i)) + Elo_K*(S(i)-E(i));
        end

        p_f = 1:length(T); % indices to propagate forward without change
        p_f(i_t) = [];
        for i = p_f
            R(g,i) = R(g-1,i);
        end

    end

    %% Extract day-by-day Elo

    i_p = zeros(size(G,1),1); % indices to preserve; set to 1 to preserve row

    for d = 1:size(G,1)-1
        if G(d+1,4) > G(d,4) % if this is the last game of the day
            i_p(d) = 1; % preserve this row
        end
    end
    R = R(i_p==1,:);
    dn = G(i_p==1,4);
    
    plot_date_range_numeric = [ % show range, numeric
                datenum(plot_date_range{1},'yyyy-mmm-dd')
                datenum(plot_date_range{2},'yyyy-mmm-dd')
            ];

    %% Plot ratings over time
    
    if show_plot

        % Crop to specified time domain
        R_plot = R;
        dn_plot = dn;
        R_plot (dn<min(plot_date_range_numeric) | dn>max(plot_date_range_numeric),:) = [];
        dn_plot(dn<min(plot_date_range_numeric) | dn>max(plot_date_range_numeric),:) = [];

        % Crop inactive teams
        disp('Inactive teams:')
        for t = length(T):-1:1
            if max(R_plot(:,t)) == min(R_plot(:,t)) % if no changes over specified window
                disp(['  ' T{t}])
                T(t) = [];
                R_plot(:,t) = [];
            end
    %         disp([T{t} ': [' num2str(min(R_plot(:,t))) ' ' num2str(max(R_plot(:,t))) ']' ])
        end

        x = double(dn_plot);
        y = R_plot;

        cols = [ % Source: https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
                    230, 25, 75
                    60, 180, 75
                    255, 225, 25
                    0, 130, 200
                    245, 130, 48
                    145, 30, 180
                    70, 240, 240
                    240, 50, 230
                    210, 245, 60
                    250, 190, 190
                    0, 128, 128
                    230, 190, 255
                    170, 110, 40
    %                 255, 250, 200
                    128, 0, 0
                    170, 255, 195
                    128, 128, 0
    %                 255, 215, 180
                    0, 0, 128
                    128, 128, 128
    %                 255, 255, 255
                    0, 0, 0
               ]/255;

        figure(1)
        clf
        hold on
        set(gcf,'color','white')
        set(gcf,'position',[50 150 800 600])
        grid on

        % Determine which teams to highlight to avoid clutter
        highlight_status = zeros(length(T),1); % highlight status
        R_r = [min(y(end,:)),max(y(end,:))]; % rating range
        d_target = (max(R_r)-min(R_r))/(plot_highlight_quantity-1);

        for i = 1:plot_highlight_quantity
            target = min(R_r) + (i-1)*d_target;
            [~, ind] = min(abs(y(end,:)-target));
            highlight_status(ind) = 1;
        end

        col_ind = 1; % color index

        % Non-highlighted teams first
        for t = 1:length(T)

            % Highlight top and bottom teams
            if length(find(R_plot(end,t)<R_plot(end,:)))<plot_show_top || length(find(R_plot(end,t)>R_plot(end,:)))<plot_show_bot
                highlight_status(t) = 1;
            end

            % Highlight specifically-requested teams
            for s = 1:length(show_teams)
                if strcmp(show_teams{s},T{t})
                    highlight_status(t) = 1;
                end
            end

            if ~highlight_status(t)
                col = 0.8 + zeros(1,3);
                lw = 0.25;
                plot(x,y(:,t),'color',col,'linewidth',lw)
            end

        end

        % Highlighted teams last
        for t = 1:length(T)
            if highlight_status(t)
                col = cols(col_ind,:);
                col_ind = col_ind+1;
                lw = 3;
                text(x(end),y(end,t),[' ' T{t}],'fontsize',9)
                plot(x,y(:,t),'color',col,'linewidth',lw)
            end
        end

        datetick('x', 'yyyy-mm-dd')
        xlim([min(dn_plot),max(dn_plot)])
        ylabel('Elo rating, tailored to postseason strength')
        title(['Postseason strength, ' plot_date_range{1} ' to ' plot_date_range{2}])
        ax_pos = get(gca,'position');
        set(gca,'position',[ax_pos(1) ax_pos(2) ax_pos(3:4)*0.9]) % shift axis left slightly for text labels
        
        %% Print and plot ratings/rankings

        R_sorted = cell(length(T),2);

        R_sorted(:,1) = T;
        for t = 1:length(T)
            R_sorted{t,2} = R_plot(end,t);
        end
        R_sorted = sortrows(R_sorted,-2);

        disp(' ')
        disp('Team ratings at end of requested window:')
        for t = 1:size(R_sorted,1)
            disp([num2str(t) '. ' R_sorted{t,1} ': ' num2str(R_sorted{t,2})])
        end

        R_sorted = sortrows(R_sorted,2);

        figure(2)
        clf
        hold on
        set(gcf,'color','white')
        barh(cell2mat(R_sorted(:,2)),'edgecolor','flat','facecolor',0.7 + zeros(1,3),'barwidth',0.7)
        lims = [min(cell2mat(R_sorted(:,2))) max(cell2mat(R_sorted(:,2)))];
        lims = [lims(1)-diff(lims)*0.1, lims(2)];
        xlim(lims)
        for t = 1:size(R_sorted,1)
            text(min(lims)+diff(lims)*0.02,t,[num2str(size(R_sorted,1)-t+1) ': ' R_sorted{t,1}],'fontsize',9)
        end
        set(gcf,'position',[900 150 800 800])
        grid on
        grid minor
        title(['Postseason strength, ' plot_date_range{2}])
        set(gca,'YTickLabel',[]);
        xlabel('Elo rating, tailored to postseason strength')

    end

end












































