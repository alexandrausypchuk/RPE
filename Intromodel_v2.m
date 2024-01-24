function results = Exp2_sim_git(sim)


% This script simulates the experimental paradigms included in the paper by
%Maes et al. entitled:

% Causal evidence supporting the proposal that dopamine transients function 
% as a temporal difference prediction error. https://doi.org/10.1101/520965.


    %
    % USAGE: results = sim_TDRL(sim)
    %
    % INPUTS:
    %   sim - string specifying simulation:
    %           'Blocking' - blocking paradigm with inhibition of DA
    %           'Second_Order' - second order conditioning paradigm with inhibition
    %           of DA

    %
    % OUTPUTS:
    %   results - see linearTDRL.m
    %
    %Matt Gardner November 2018
    
%This provides a transfer function for determining Pavlovian conditioned responses
%from value
CR = @(V,max,Voff,ks) ks*(V - Voff); 
%parameters used:
Voff = .05; %.4 original, tried .8
ks = 5; %was 3 original
%Different maximal responding were used depending on whether the animals
%were connected to the patch fibers or not
MaxFree = 55;
MaxFiber = 35;

%switch sim
    %case 'Blocking'
        
        %Stimuli used in the Pavlovian Blocking Design. Each row of the arrays
        %represents a single time step. Each column represents whether an
        %element was present for the particular stimulus in a binary manner.
        %randTrials adds the ITI between trials
        A =  [1 0 0 0 0; 0 0 0 0 0];
        %B =  [0 0 0 0 1; 0 0 0 0 0];
        AX = [1 1 0 0 0; 0 0 0 0 0];
        %BZ = [0 0 0 1 1; 0 0 0 0 0];                
        X =  [0 1 0 0 0; 0 0 0 0 0];
        %Z =  [0 0 0 1 0; 0 0 0 0 0];

        %Group variables
        
        
        %This sets the alphas for the SR and U learning:
        alpha = 0.05;
        
        R = struct();
            
        %stimuli are randomized within each execution of the TDRL model.
        %This is included for designs in which the order of stimulus
        %presentations affects the value.
        for i = 1:100
            
            %this randomizes trial order for each stage of conditioning
            S{1} = randTrials_2(6,16,A);%Trials: 48 16*6
            S{2} = randTrials_2(12,8,A);%Trials: 6 days with 2 reps of (8 of A and 3 of B)
            S{3} = randTrials_2(18,2,AX);%Trials: 6 days with 3 reps of (2 of AX, 2 of AY, 2 of BZ)
            S{4} = randTrials_2(1,1,X);%Trials: 3 reps of (1 of AX, 1 of AY, 1 of BZ)
            
            %This gets a vector of the stage numbers. This is used incase
            %there are identical stimuli in each stage.
            stg = stagenum(S);
            
            M = cell2mat(S');
            
            %logicals of the first state of each stimulus
            L.a = fstate(A,M) & stg == 1;
            L.a2 = fstate(A,M) & stg == 2;
            %L.b = fstate(B,M) & stg == 2;
            L.ax = fstate(AX,M) & stg == 3;
            %L.ay = fstate(AY,M) & stg == 3; 
            %L.bz = fstate(BZ,M) & stg == 3;
            L.x = fstate(X,M) & stg == 4;
            %L.y = fstate(Y,M) & stg == 4;
            %L.z = fstate(Z,M) & stg == 4;
            
            FN = fieldnames(L);

            %rewarded states
            %r = double(nstate(L.a,2) | nstate(L.a2,2) | nstate(L.ax,2) | nstate(L.ay,2) | nstate(L.bz,2));
            r = double(nstate(L.a,2) | nstate(L.a2,2) | nstate(L.ax,2) ); % | nstate(L.bz,2)
            
            %This implements stimulation during reward after AX
            %trials
            opto = 1*nstate(L.a2,2) + 1*nstate(L.ax,2); %for AXstim
            %opto = 0.5;
            %opto(n) = 1 + opto(n-1);
            
            %Initial values of each of the stimuli. We set cross-modality
            %generalization at 0.2 and generalization of same modality
            %stimuli at 0.7
            ui = [0; 0.2; 0.2; 0.2; 0.2];
            
            results = LinearTDRL2_git(M,r,ui,alpha,opto,'none');
            for k = 1:numel(FN)
                R(1).(FN{k})(i,:) = CR([results(L.(FN{k})).V],MaxFiber,Voff,ks);    
            end
            
            results = LinearTDRL2_git(M,r,ui,alpha,opto,'tonic');
            for k = 1:numel(FN)
                R(2).(FN{k})(i,:) = CR([results(L.(FN{k})).V],MaxFiber,Voff,ks);
            end
            
            results = LinearTDRL2_git(M,r,ui,alpha,opto,'V');
            for k = 1:numel(FN)
                R(3).(FN{k})(i,:) = CR([results(L.(FN{k})).V],MaxFiber,Voff,ks);
            end

        end
           

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        %The following code creates a panel with line plots of each stage of conditioning
        %as a different column, as well as a final column with a bar plot
        %of the behavior during the extinction test
                     
        plotwidths = [0 .16 .16 .12 .35]; %[0 .2 .2 .15 .25]
        plotheights = [.13 .13 .13]; 
        spcx = [0.1 0.01 0.01 0.07];
        spcy = [0.15 0.15 0.15];
        positions = cell(3,4);
        
        %This determines the subplot locations for the panel
        for sx = 1:4
            for sy = 1:3
                
                positions{4-sy,sx} = [sum(spcx(1:sx)) + sum(plotwidths(1:sx)), sum(spcy(1:sy)) + sum(plotheights(1:sy)), plotwidths(sx + 1),plotheights(sy)];
            end
        end
        
        V = zeros(numel(R),2);
        
        for k = 1:numel(R)
            %V(k,:) = mean([R(k).x(:,1) R(k).y(:,1) R(k).z(:,1)]);
            V(k,:) = mean(R(k).x(:,1)); %, R(k).z(:,1)]
        end
        
        
        
        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [.1, .1, .5, .8],'color', 'white'); %[.1, 0.1, .4, .8]
        hold on
        models = {'Controls', 'Error', 'Value'};
        
        for m = 1:3 %for loop to make the graph for each model from model 2 to 3
            

            subplot('Position',positions{m,4})
            

            %b = bar(1:numel(mean([R(m).x(:,1) R(m).z(:,1) R(1).x(:,1) R(1).z(:,1)])),diag(mean([R(m).x(:,1) R(m).z(:,1) R(1).x(:,1) R(1).z(:,1)])),.95, 'stacked','FaceColor', 'c');
            b = bar(1:numel(mean(R(m).x(:,1))),diag(mean(R(m).x(:,1))),.95, 'stacked','FaceColor', 'c'); %R(m).z(:,1) R(m).z(:,1)
            b(1).FaceColor = [25 75 245]./255; %bar 1 
            %b(2).FaceColor = [1 1 1];
            %b(3).FaceColor = [225 55 165]./255;
            %b(4).FaceColor = [1 1 1];


            %legend({'X' 'Z'},'FontSize',15,'Location','North');
            set(gca,'FontSize',11,'XLim', [.3 10.5],'YLim',[0 40],'XTick', [], 'XTickLabel', {'X' 'Y'},'YTick', 5:10:35,'color', 'none'); %XTICK 1 2 3.5 4.5 6 7 8.5 9.5
            %X lim was .4 to 2.6
            %xlabel(models{m},'FontSize',15)
            ylabel('Value','FontSize',11);
            %title('Blocking','FontSize',25,'FontWeight','Bold');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot('Position',positions{m,1})
            title(models{m})
            Resp = [];
            for i = 1:6
                range = (i-1)*6+1:i*6;
                %Resp(i,:) = mean(mean(R(1).a(:,range),2)); %#ok<*AGROW>
                Resp(i,:) = mean(mean(R(m).a(:,range),1));
            end
            plot(1:6,Resp(:,1),'-o', 'color',[25 75 245]./255)
            %hold on
            %plot(1:6,Resp(:,2),'-o','color',[225 55 165]./255)
            
            set(gca,'FontSize',11, 'XTick', 1:6, 'XLim',[0.5 6.5],'YLim',[0 40],'YTick', 5:10:35, 'color', 'none')
            ylabel('Value','FontSize',11);
            %xlabel('Days','FontSize',15);


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot('Position',positions{m,2})
            title(models{m})
            Resp = [];
            for i = 1:6
                range = (i-1)*6+1:i*6;
                Resp(i,:) = mean(mean(R(m).a2(:,range),1)); %mean(R(m).b(:,range),2)  mean(R(1).b(:,range),2)+1
            end
            plot(7:12,Resp(:,1),'-o', 'color',[25 75 245]./255)
            %hold on
            %plot(7:12,Resp(:,2),'-o','color',[25 75 245]./255, 'LineStyle','--')
            %hold on
            %plot(7:12,Resp(:,2),'-o','color',[225 55 165]./255)
            %hold on
            %plot(7:12,Resp(:,4),'-o','color',[225 55 165]./255, 'LineStyle','--')
            set(gca,'FontSize',11, 'XTick', 7:12, 'XLim',[6.5 12.5],'YLim',[0 40],'YTick',[],'color', 'none');%0:10:50)
            %xlabel('Days','FontSize',15);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot('Position',positions{m,3})
            title(models{m})
            Resp = [];
            for i = 1:4
                range = (i-1)*6+1:i*6;
                Resp(i,:) = mean(mean(R(m).ax(:,range),1)); %mean(R(m).ay(:,range),2)+.5 mean(R(m).bz(:,range),2)  mean(R(1).bz(:,range),2)+1
            end
            plot(1:4,Resp(:,1),'-o','MarkerSize',5, 'color',[25 75 245]./255)
            %hold on
            %plot(1:4,Resp(:,2),'-o','MarkerSize',5,'color',[25 75 245]./255, 'LineStyle','--')
            %hold on
            %plot(1:4,Resp(:,2),'-o','MarkerSize',5,'color',[225 55 165]./255)
            %hold on
            %plot(1:4,Resp(:,4),'-o','MarkerSize',5,'color',[225 55 165]./255, 'LineStyle','--')
            
            set(gca,'FontSize',11, 'XTick', 1:4, 'XLim',[0.5 4.5],'YLim',[0 40],'YTick',[],'color', 'none');
            %  ylabel('CR','FontSize',15);
            %xlabel('Days','FontSize',15);

        end



        
end
%end

function S1 = fstate(i,S)
%Provides a logical vector (N x 1) of the first state of 
%stimulus i by template matching i within the state matrix J. This script
%will check up to three rows of stimulus i within the state matrix J.
%Several rows of each stimulus are checked incase the first row or two of
%stimulus i are identical to another stimulus.
    S1 = false(size(S,1),1);
    
    if size(i,1) == 1
        S1 = all(S == i(1,:),2);
    elseif size(i,1) == 2
        S1 = all(S == i(1,:),2) & (all([S(2:end,:); zeros(1,size(i,2))] ==  i(2,:),2));
    elseif size(i,1) > 2  
        S1 = all(S == i(1,:),2) & (all([S(2:end,:); zeros(1,size(i,2))] == i(2,:),2)) & (all([S(3:end,:); zeros(2,size(i,2))] == i(3,:),2));
    end    
end

function Sn = nstate(X,n)
%Provides a logical vector of state n within a given stimulus by shifting
%the logicals of stimulus X by n values. For example, n = 2 provides an N x 1
%logical vector of the second row of the stimulus in X, where X is an N x 1
%logical vector. If n = -1, this returns the logicals for the prior state.

    if n >= 0
        Sn = [false(n-1,1); X(1:end-n+1)];
    else
        Sn = [X(-1*n+1:end); false(-1*n,1)];
    end
end

function Stage = stagenum(S)
%provides an N X 1 vector of the stage number of conditioning for
%the state matrix J. J is a cell array whose size is 1 X the number of
%stages.
    Stage = zeros(size(cell2mat(S'),1),1);
    
    ind_start = 1;
    for i = 1:numel(S)
        ind_end = ind_start + size(S{i},1) - 1;
        
        Stage(ind_start:ind_end) = i;
        ind_start = ind_end + 1;
        
    end    

end