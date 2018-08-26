function home_field_advantage
    
    % Compares team and overall win rates at home versus away over
    % specified time domain, for roughly quantifying home field advantage
    
    % Daniel W. Dichter
    % 2018-08-26
    
    %% Inputs
    
    year_range = 2008:2018;
    include_regular = 1;
    include_postseason = 1;
    
    %% Extract games matching input ranges
    
    clc
    load('live_ball_game_log_n','G_n','T')
    G = G_n;

    if ~include_regular
        G(G(:,9)==1,:) = []; % crop regular games
    end
    if ~include_postseason
        G(G(:,9)==2,:) = []; % crop post-season games
    end
    
    G(G(:,1)<min(year_range),:) = [];
    G(G(:,1)>max(year_range),:) = [];

    %% Display the input conditions

    if isempty(G)
        error('No games match specified input conditions')
    else
        disp('Input games:')
        disp(['   include_regular: ' num2str(include_regular)])
        disp(['   include_postseason: ' num2str(include_postseason)])
        disp(['   Years: ' num2str(min(G(:,1))) '-' num2str(max(G(:,1)))])
        disp(['   Games: ' num2str(size(G,1))])
        disp(['   Regular: ' num2str(length(find(G(:,9)==1)))])
        disp(['   Post: ' num2str(length(find(G(:,9)==2)))])
    end
    
    %% Extract record for each team at home and away, for each year
    
    H.W = zeros(length(year_range),length(T)) -1; % home wins
    H.L = zeros(size(H.W)) -1; % home losses
    
    A.W = zeros(size(H.W)) -1; % away wins
    A.L = zeros(size(H.W)) -1; % away losses
    
    for y = 1:length(year_range)
        
        G_y = G(G(:,1)==year_range(y),:); % match year
        
        for t = 1:length(T)
            
            % Home games
            G_h = G_y(G_y(:,7)==t,:);
            if isempty(G_h)
                H.W(y,t) = 0;
                H.L(y,t) = 0;
            else
                H.W(y,t) = sum(G_h(:,6)<G_h(:,8));
                H.L(y,t) = sum(G_h(:,6)>G_h(:,8));
            end
            
            % Away games
            G_a = G_y(G_y(:,5)==t,:);
            if isempty(G_a)
                A.W(y,t) = 0;
                A.L(y,t) = 0;
            else
                A.W(y,t) = sum(G_a(:,6)>G_a(:,8));
                A.L(y,t) = sum(G_a(:,6)<G_a(:,8));
            end
            
        end
        
    end
    
    %% Crop teams that are partially or completely inactive
    
    G_P = H.W+H.L + A.W+A.L; % games played
    is_active = ~min(G_P)==0;
    disp(['Teams active over range: ' num2str(length(find(is_active==1)))])
    
    %% Plot results
    
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
    
    for f = 1:3
        figure(f)
        clf
        hold on
        set(gcf,'color','white')
        grid on
        ylabel('Win rate, 0-1')
        xlim([min(year_range) max(year_range)])
        a_p = get(gca,'position');
        switch f
            case {1,2}
                set(gca,'position',[a_p(1)-.025 a_p(2:4)]) % shift axis left slightly for text labels
            case 3
                set(gca,'position',[a_p(1)-.025 a_p(2)-0.025 a_p(3:4)]) % shift axis left slightly for text labels
            otherwise
                % Do nothing
        end
    end
    
    c = 1; % color index
    
    for t = 1:length(T)
        
        if ~is_active(t)
            continue
        end
        
        wra = A.W(:,t)./(A.W(:,t)+A.L(:,t)); % win rate, away
        wrh = H.W(:,t)./(H.W(:,t)+H.L(:,t)); % win rate, home
        
        if c <= size(cols,1)
            lw = 2;
            col = cols(c,:);
            c = c+1;
            figure(1)
            text(year_range(end),wra(end),[' ' T{t}],'fontsize',8)
            figure(2)
            text(year_range(end),wrh(end),[' ' T{t}],'fontsize',8)
        else
            lw = 0.5;
            col = zeros(1,3)+ 0.75;
        end
        
        figure(1)
        plot(year_range,wra,'linewidth',lw,'color',col)
        figure(2)
        plot(year_range,wrh,'linewidth',lw,'color',col)
        
    end
    
    figure(1)
    title('Away games')
    figure(2)
    title('Home games')
    
    figure(3)
    
    wra_o = sum(A.W,2)./(sum(A.W,2)+sum(A.L,2)); % win rate, away, overall, all teams
    wrh_o = sum(H.W,2)./(sum(H.W,2)+sum(H.L,2)); % win rate, home, overall, all teams
    
    plot(year_range,wra_o,'linewidth',2,'color',cols(1,:))
    plot(year_range,wrh_o,'linewidth',2,'color',cols(2,:))
    
    text(year_range(end),wra_o(end),' Away')
    text(year_range(end),wrh_o(end),' Home')
    
    title({
            'Home field advantage, all teams'
            ['\rmAverage difference from neutral field: \pm' num2str(round(mean(wrh_o-wra_o)*100000)/1000/2) '%']
         })
    
end












































