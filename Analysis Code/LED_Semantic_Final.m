folder = uigetdir;
files = dir(folder);

for i = 3:length(files)
    fileNames(i-2).name = dir(fullfile(folder, files(i).name, '*.csv'));
end

subjectNum = ["10"; "11"; "12"; "1"; "2"; "3"; "4"; "5"; "6"; "7"; "8"; "9"];
newFolders = ["Covert B", "Overt B", "Covert B_BB_BBB", "Overt B_BB_BBB", "Covert B_BBB", "Overt B_BBB"];
markers = ["o", "s", "d", "^", "*", "x"];
lineStyles = ["-", "--", "-", "--", "-", "--"];
colors = ['#000000', "#A9A9A9", 'b', "#00CCCC", "r", "#ffc6c4"];
legends = ["Covert B", "Overt B", "Covert B/BB/BBB", "Overt B/BB/BBB", "Covert B/BBB", "Overt B/BBB"];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~INDIVIDUAL PLOTS START HERE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
for j = 1:length(fileNames)
    numFiles = fileNames(j).name;
    mkdir(newFolders(j));
    fileInfo = dir(newFolders(j));
    filePath = fileInfo(1).folder;
    Left_slope = [];
    Left_intercept = [];
    Chi2_left = [];
    Left_slope_sd = [];
    Left_intercept_sd = [];
    Red_Chi2_left = [];
    R_left = [];
    Right_slope = [];
    Right_intercept = [];
    Chi2_right = [];
    Right_slope_sd = [];
    Right_intercept_sd = [];
    Red_Chi2_right = [];
    R_right = [];
    for ji = 1:size(numFiles) %making csv into data input
        newData = [];
        figure(ji);
        fileName = fullfile(numFiles(ji).folder, numFiles(ji).name);
        rawData = table2array(readtable(fileName));
        dist = rawData(:,2);
        uniqueDist = unique(rawData(:,2))';
        RT = rawData(:,3);
        nonZeroRT = find(RT ~= 0);
        newRT = RT(nonZeroRT, :);
        newDist = dist(nonZeroRT, :);
        count = 0;
        for jii = 1:size(newDist) %getting rid of outliers
            if (newRT(jii) > (prctile(newRT, 25) - 1.5*iqr(newRT))) && (newRT(jii) < (prctile(newRT, 75) + ...
                    1.5*iqr(newRT)))
                newData(jii-count).dist = newDist(jii);
                newData(jii-count).RT = newRT(jii);
            else
                count = count + 1;
            end
        end
        angle = [];
        meanRT = [];
        stdDev = [];
        yErr = [];
        for jk = 1:length(uniqueDist) %finding means, stdDev, etc.
            for jki = 1:length(newData)
                if newData(jki).dist == uniqueDist(jk)
                    angle = [angle, newData(jki).RT];
                end
            end
            meanRT = [meanRT, mean(angle)];
            stdDev = [stdDev, std(angle)];
            yErr = [yErr, std(angle)/sqrt(length(angle))];
            angle = [];
        end
        vals(1).uniqueDist = uniqueDist(uniqueDist<=0);
        vals(1).meanRT = meanRT(uniqueDist<=0);
        vals(1).stdDev = stdDev(uniqueDist<=0);
        vals(1).yErr = yErr(uniqueDist<=0);

        vals(2).uniqueDist = uniqueDist(uniqueDist>=0);
        vals(2).meanRT = meanRT(uniqueDist>=0);
        vals(2).stdDev = stdDev(uniqueDist>=0);
        vals(2).yErr = yErr(uniqueDist>=0);
        for jl = 1:2
            weights = (1./vals(jl).yErr).^2;
            f = @(x, xPoints, yPoints, w)sum(w.*((yPoints- ((xPoints.*x(1))+x(2))).^2));
            optFun = @(x)f(x, vals(jl).uniqueDist, vals(jl).meanRT, weights);
            ms = MultiStart;
            OLSFit = polyfit(vals(jl).uniqueDist, vals(jl).meanRT, 1);
            guessParams = [OLSFit(1), OLSFit(2)];
            problem = createOptimProblem('fmincon', 'x0', guessParams, 'objective', optFun, ...
                'lb', [OLSFit(1)-50, OLSFit(2)-500], 'ub', [OLSFit(1)+50, OLSFit(2)+500]);
            params = run(ms, problem, 25);
            slope = params(1);
            intercept = params(2);
            chiSquareFit(jl).slope = slope;
            chiSquareFit(jl).intercept = intercept;
            chi2Val = optFun(params);
            chiSquareFit(jl).chi2Val = chi2Val;
            syms sErr;
            slopeErr = solve(f([sErr, intercept], vals(jl).uniqueDist, vals(jl).meanRT, weights) == chi2Val + 1, sErr);
            chiSquareFit(jl).slopeErr = double(slopeErr(2) - slope);
            syms iErr;
            intErr = solve(f([slope, iErr], vals(jl).uniqueDist, vals(jl).meanRT, weights) == chi2Val+1, iErr);
            chiSquareFit(jl).interceptErr = double(intErr(2)-intercept);
            chiSquareFit(jl).redChiSquare = chi2Val/ (length(vals(jl).uniqueDist) - 2);
    
            R = corrcoef(vals(jl).uniqueDist, vals(jl).meanRT);
            chiSquareFit(jl).R = R(1,2);
            
            grid on;
            scatter(uniqueDist, meanRT, 10, markers(j), 'MarkerEdgeColor', colors(j), 'MarkerFaceColor', colors(j), ...
                'HandleVisibility', 'off');
            hold on;
            errorbar(uniqueDist, meanRT, yErr, '.' , 'color' , colors(j), 'CapSize', 0, 'HandleVisibility', 'off');
            hold on;
            plot(vals(jl).uniqueDist, polyval(params,vals(jl).uniqueDist), 'color', colors(j), 'linestyle', ...
                lineStyles(j), 'HandleVisibility', 'off');
            hold on;
        end
        h = plot(0, 3000, 'Color', colors(j), 'MarkerFaceColor', colors(j), 'MarkerEdgeColor', colors(j), ...
            'LineStyle', lineStyles(j), 'Marker', markers(j), 'DisplayName', legends(j));
        pbaspect([2 3 1]);
        legend('show', 'Location', 'northeastoutside');
        xlabel('Eccentricity (degrees)');
        xticks(uniqueDist(1):15:max(uniqueDist));
        ylabel('Reaction Time (ms)');
        title(['Subject ', subjectNum(ji), ': Reaction Time vs. Eccentricity']);
        xline(0, 'Color', '#C9C9C9', 'HandleVisibility', 'off');
        xlim([(uniqueDist(1) -10) (uniqueDist(end)+10)]);
        ylim([150 700]);
        saveas(gcf, fullfile(filePath, append('Subject', subjectNum(ji), '.png')));
        hold on;
        
        keyParams = table2array(struct2table(chiSquareFit));

        Left_slope = [Left_slope; keyParams(1,1)];
        Left_intercept = [Left_intercept; keyParams(1,2)];
        Chi2_left = [Chi2_left; keyParams(1,3)];
        Left_slope_sd = [Left_slope_sd; keyParams(1,4)];
        Left_intercept_sd = [Left_intercept_sd; keyParams(1,5)];
        Red_Chi2_left = [Red_Chi2_left; keyParams(1,6)];
        R_left = [R_left; keyParams(1,7)];
        Right_slope = [Right_slope; keyParams(2,1)];
        Right_intercept = [Right_intercept; keyParams(2,2)];
        Chi2_right = [Chi2_right; keyParams(2,3)];
        Right_slope_sd = [Right_slope_sd; keyParams(2,4)];
        Right_intercept_sd = [Right_intercept_sd; keyParams(2,5)];
        Red_Chi2_right = [Red_Chi2_right; keyParams(2,6)];
        R_right = [R_right; keyParams(2,7)];
    end
    
    paramTable = table(subjectNum, Right_slope, Right_slope_sd, Left_slope, Left_slope_sd, Right_intercept, ...
        Right_intercept, Right_intercept_sd, Left_intercept, Left_intercept_sd, Chi2_right, ... 
        Red_Chi2_right, Chi2_left, Red_Chi2_left, R_right, R_left);
    
    writetable(paramTable, fullfile(filePath,  append(newFolders(j), '_Param_Info.xlsx')))
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~AGGREGATED PLOTS START HERE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

newFolders = ["Aggregate", "Overt B", "Covert B_BB_BBB", "Overt B_BB_BBB", "Covert B_BBB", "Overt B_BBB"];
markers = ["o", "s", "d", "^", "*", "x"];
lineStyles = ["-", "--", "-", "--", "-", "--"];
colors = ['#000000', "#A9A9A9", 'b', "#00CCCC", "r", "#ffc6c4"];
legends = ["Covert B", "Overt B", "Covert B/BB/BBB", "Overt B/BB/BBB", "Covert B/BBB", "Overt B/BBB"];

%Making our single directory for aggregate and defining its location
mkdir("Aggregate");
fileInfo = dir("Aggregate");
filePath = fileInfo(1).folder;

for jj = 1:length(fileNames)
    numFiles = fileNames(jj).name;
    Left_slope = [];
    Left_intercept = [];
    Chi2_left = [];
    Left_slope_sd = [];
    Left_intercept_sd = [];
    Red_Chi2_left = [];
    R_left = [];
    Right_slope = [];
    Right_intercept = [];
    Chi2_right = [];
    Right_slope_sd = [];
    Right_intercept_sd = [];
    Red_Chi2_right = [];
    R_right = [];
    aggFileNames = [];
    for k = 1:size(numFiles)%aggregating data
        aggFileNames(k).name = fullfile(numFiles(k).folder, numFiles(k).name);
    end
    rawData = [];
    for ki = 1:length(aggFileNames)
        rawData = vertcat(rawData, table2array(readtable(aggFileNames(1,ki).name)));
    end
    newData = [];
    figure(ji+1);
    dist = rawData(:,2);
    uniqueDist = unique(rawData(:,2))';
    RT = rawData(:,3);
    nonZeroRT = find(RT ~= 0);
    newRT = RT(nonZeroRT, :);
    newDist = dist(nonZeroRT, :);
    count = 0;
    for jii = 1:size(newDist) %getting rid of outliers
        if (newRT(jii) > (prctile(newRT, 25) - 1.5*iqr(newRT))) && (newRT(jii) < (prctile(newRT, 75) + 1.5*iqr(newRT)))
            newData(jii-count).dist = newDist(jii);
            newData(jii-count).RT = newRT(jii);
        else
            count = count + 1;
        end
    end
    angle = [];
    meanRT = [];
    stdDev = [];
    yErr = [];
    for jk = 1:length(uniqueDist) %finding means, stdDev, etc.
        for jki = 1:length(newData)
            if newData(jki).dist == uniqueDist(jk)
                angle = [angle, newData(jki).RT];
            end
        end
        meanRT = [meanRT, mean(angle)];
        stdDev = [stdDev, std(angle)];
        yErr = [yErr, std(angle)/sqrt(length(angle))];
        angle = [];
    end
    vals(1).uniqueDist = uniqueDist(uniqueDist<=0);
    vals(1).meanRT = meanRT(uniqueDist<=0);
    vals(1).stdDev = stdDev(uniqueDist<=0);
    vals(1).yErr = yErr(uniqueDist<=0);

    vals(2).uniqueDist = uniqueDist(uniqueDist>=0);
    vals(2).meanRT = meanRT(uniqueDist>=0);
    vals(2).stdDev = stdDev(uniqueDist>=0);
    vals(2).yErr = yErr(uniqueDist>=0);
    for jl = 1:2 %Chi Square and Fit
        weights = (1./vals(jl).yErr).^2;
        f = @(x, xPoints, yPoints, w)sum(w.*((yPoints- ((xPoints.*x(1))+x(2))).^2));
        optFun = @(x)f(x, vals(jl).uniqueDist, vals(jl).meanRT, weights);
        ms = MultiStart;
        OLSFit = polyfit(vals(jl).uniqueDist, vals(jl).meanRT, 1);
        guessParams = [OLSFit(1), OLSFit(2)];
        problem = createOptimProblem('fmincon', 'x0', guessParams, 'objective', optFun, ...
            'lb', [OLSFit(1)-50, OLSFit(2)-500], 'ub', [OLSFit(1)+50, OLSFit(2)+500]);
        params = run(ms, problem, 25);
        slope = params(1);
        intercept = params(2);
        chiSquareFit(jl).slope = slope;
        chiSquareFit(jl).intercept = intercept;
        chi2Val = optFun(params);
        chiSquareFit(jl).chi2Val = chi2Val;
        syms sErr;
        slopeErr = solve(f([sErr, intercept], vals(jl).uniqueDist, vals(jl).meanRT, weights) == chi2Val + 1, sErr);
        chiSquareFit(jl).slopeErr = double(slopeErr(2) - slope);
        syms iErr;
        intErr = solve(f([slope, iErr], vals(jl).uniqueDist, vals(jl).meanRT, weights) == chi2Val+1, iErr);
        chiSquareFit(jl).interceptErr = double(intErr(2)-intercept);
        chiSquareFit(jl).redChiSquare = chi2Val/ (length(vals(jl).uniqueDist) - 2);
    
        R = corrcoef(vals(jl).uniqueDist, vals(jl).meanRT);
        chiSquareFit(jl).R = R(1,2);
        
        grid on;
        scatter(uniqueDist, meanRT, 10, markers(jj), 'MarkerEdgeColor', colors(jj), 'MarkerFaceColor', colors(jj), ...
            'HandleVisibility', 'off');
        hold on;
        errorbar(uniqueDist, meanRT, yErr, '.' , 'color' , colors(jj), 'CapSize', 0, 'HandleVisibility', 'off');
        hold on;
        plot(vals(jl).uniqueDist, polyval(params,vals(jl).uniqueDist), 'color', colors(jj), 'linestyle', ...
            lineStyles(jj), 'HandleVisibility', 'off');
        hold on;
    end
    %Plotting
    h = plot(0, 3000, 'Color', colors(jj), 'MarkerFaceColor', colors(jj), 'MarkerEdgeColor', colors(jj), ...
        'LineStyle', lineStyles(jj), 'Marker', markers(jj), 'DisplayName', legends(jj));
    pbaspect([2 3 1]);
    legend('show', 'Location', 'northeastoutside');
    xlabel('Eccentricity (degrees)');
    xticks(uniqueDist(1):15:max(uniqueDist));
    ylabel('Reaction Time (ms)');
    title(['Aggregate: Reaction Time vs. Eccentricity']);
    xline(0, 'Color', '#C9C9C9', 'HandleVisibility', 'off');
    xlim([(uniqueDist(1) -10) (uniqueDist(end)+10)]);
    ylim([150 625]);
    saveas(gcf, fullfile(filePath, 'Aggregate.png'));
    hold on;
    
    %Gathering parameters
    keyParams = table2array(struct2table(chiSquareFit));

    Left_slope = [Left_slope; keyParams(1,1)];
    Left_intercept = [Left_intercept; keyParams(1,2)];
    Chi2_left = [Chi2_left; keyParams(1,3)];
    Left_slope_sd = [Left_slope_sd; keyParams(1,4)];
    Left_intercept_sd = [Left_intercept_sd; keyParams(1,5)];
    Red_Chi2_left = [Red_Chi2_left; keyParams(1,6)];
    R_left = [R_left; keyParams(1,7)];
    Right_slope = [Right_slope; keyParams(2,1)];
    Right_intercept = [Right_intercept; keyParams(2,2)];
    Chi2_right = [Chi2_right; keyParams(2,3)];
    Right_slope_sd = [Right_slope_sd; keyParams(2,4)];
    Right_intercept_sd = [Right_intercept_sd; keyParams(2,5)];
    Red_Chi2_right = [Red_Chi2_right; keyParams(2,6)];
    R_right = [R_right; keyParams(2,7)];

    %Putting parameters into table
    paramTable = table("Aggregate", Right_slope, Right_slope_sd, Left_slope, Left_slope_sd, Right_intercept, ...
        Right_intercept, Right_intercept_sd, Left_intercept, Left_intercept_sd, Chi2_right, ... 
        Red_Chi2_right, Chi2_left, Red_Chi2_left, R_right, R_left);
    
    writetable(paramTable, fullfile(filePath,  append('Aggregate_Param_Info.xlsx')))
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~STACKED PLOTS START HERE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stackPlot = ["Covert B", "Overt B", "Covert B_BB_BBB", "Overt B_BB_BBB", "Covert B_BBB", "Overt B_BBB"];
colors = ['#000000', "b", "#026440", "#00CCCC", "r", "#ffc6c4", "#90ee90", "m", "#b29700", "#ffa500", "#654321",...
    "#800080", "#003060"];
legends = ["10"; "11"; "12"; "1"; "2"; "3"; "4"; "5"; "6"; "7"; "8"; "9"; "Aggregate"];
titles = ["Covert B", "Overt B", "Covert B/BB/BBB", "Overt B/BB/BBB", "Covert B/BBB", "Overt B/BBB"];

for jjj = 1:length(fileNames)
    numFiles = fileNames(jjj).name;
    %Making our single directory for aggregate and defining its location
    mkdir("Stacked Plots");
    fileInfo = dir("Stacked Plots");
    filePath = fileInfo(1).folder;
    
    figure(ji+1+jjj);
    
    Left_slope = [];
    Left_intercept = [];
    Chi2_left = [];
    Left_slope_sd = [];
    Left_intercept_sd = [];
    Red_Chi2_left = [];
    R_left = [];
    Right_slope = [];
    Right_intercept = [];
    Chi2_right = [];
    Right_slope_sd = [];
    Right_intercept_sd = [];
    Red_Chi2_right = [];
    R_right = [];
    for ji = 1:length(numFiles)+1 %making csv into data input
        newData = [];
        if ji == length(numFiles)+1
            aggFileNames = [];
            for k = 1:size(numFiles)%aggregating data
                aggFileNames(k).name = fullfile(numFiles(k).folder, numFiles(k).name);
            end
            rawData = [];
            for ki = 1:length(aggFileNames)
                rawData = vertcat(rawData, table2array(readtable(aggFileNames(1,ki).name)));
            end
        else
            fileName = fullfile(numFiles(ji).folder, numFiles(ji).name);
            rawData = table2array(readtable(fileName));
        end
        dist = rawData(:,2);
        uniqueDist = unique(rawData(:,2))';
        RT = rawData(:,3);
        nonZeroRT = find(RT ~= 0);
        newRT = RT(nonZeroRT, :);
        newDist = dist(nonZeroRT, :);
        count = 0;
        for jii = 1:size(newDist) %getting rid of outliers
            if (newRT(jii) > (prctile(newRT, 25) - 1.5*iqr(newRT))) && (newRT(jii) < (prctile(newRT, 75) + ...
                    1.5*iqr(newRT)))
                newData(jii-count).dist = newDist(jii);
                newData(jii-count).RT = newRT(jii);
            else
                count = count + 1;
            end
        end
        angle = [];
        meanRT = [];
        stdDev = [];
        yErr = [];
        for jk = 1:length(uniqueDist) %finding means, stdDev, etc.
            for jki = 1:length(newData)
                if newData(jki).dist == uniqueDist(jk)
                    angle = [angle, newData(jki).RT];
                end
            end
            meanRT = [meanRT, mean(angle)];
            stdDev = [stdDev, std(angle)];
            yErr = [yErr, std(angle)/sqrt(length(angle))];
            angle = [];
        end
        vals(1).uniqueDist = uniqueDist(uniqueDist<=0);
        vals(1).meanRT = meanRT(uniqueDist<=0);
        vals(1).stdDev = stdDev(uniqueDist<=0);
        vals(1).yErr = yErr(uniqueDist<=0);

        vals(2).uniqueDist = uniqueDist(uniqueDist>=0);
        vals(2).meanRT = meanRT(uniqueDist>=0);
        vals(2).stdDev = stdDev(uniqueDist>=0);
        vals(2).yErr = yErr(uniqueDist>=0);
        for jl = 1:2 %chi square fit and preliminary plot stuff
            weights = (1./vals(jl).yErr).^2;
            f = @(x, xPoints, yPoints, w)sum(w.*((yPoints- ((xPoints.*x(1))+x(2))).^2));
            optFun = @(x)f(x, vals(jl).uniqueDist, vals(jl).meanRT, weights);
            ms = MultiStart;
            OLSFit = polyfit(vals(jl).uniqueDist, vals(jl).meanRT, 1);
            guessParams = [OLSFit(1), OLSFit(2)];
            problem = createOptimProblem('fmincon', 'x0', guessParams, 'objective', optFun, ...
                'lb', [OLSFit(1)-50, OLSFit(2)-500], 'ub', [OLSFit(1)+50, OLSFit(2)+500]);
            params = run(ms, problem, 25);
            slope = params(1);
            intercept = params(2);
            chiSquareFit(jl).slope = slope;
            chiSquareFit(jl).intercept = intercept;
            chi2Val = optFun(params);
            chiSquareFit(jl).chi2Val = chi2Val;
            syms sErr;
            slopeErr = solve(f([sErr, intercept], vals(jl).uniqueDist, vals(jl).meanRT, weights) == chi2Val + 1, sErr);
            chiSquareFit(jl).slopeErr = double(slopeErr(2) - slope);
            syms iErr;
            intErr = solve(f([slope, iErr], vals(jl).uniqueDist, vals(jl).meanRT, weights) == chi2Val+1, iErr);
            chiSquareFit(jl).interceptErr = double(intErr(2)-intercept);
            chiSquareFit(jl).redChiSquare = chi2Val/ (length(vals(jl).uniqueDist) - 2);
    
            R = corrcoef(vals(jl).uniqueDist, vals(jl).meanRT);
            chiSquareFit(jl).R = R(1,2);
            
            grid on;
            scatter(uniqueDist, meanRT, 10, "o", 'MarkerEdgeColor', colors(ji), 'MarkerFaceColor', colors(ji), ...
                'HandleVisibility', 'off');
            hold on;
            errorbar(uniqueDist, meanRT, yErr, '.' , 'color' , colors(ji), 'CapSize', 0, 'HandleVisibility', 'off');
            hold on;
            plot(vals(jl).uniqueDist, polyval(params,vals(jl).uniqueDist), 'color', colors(ji), 'linestyle', ...
                "-", 'HandleVisibility', 'off');
            hold on;
        end
        %plotting
        h = plot(0, 3000, 'Color', colors(ji), 'MarkerFaceColor', colors(ji), 'MarkerEdgeColor', colors(ji), ...
            'LineStyle', "-", 'Marker', "o", 'DisplayName', legends(ji));
        pbaspect([2 3 1]);
        legend('show', 'Location', 'northeastoutside');
        xlabel('Eccentricity (degrees)');
        xticks(uniqueDist(1):15:max(uniqueDist));
        ylabel('Reaction Time (ms)');
        title(titles(jjj));
        xline(0, 'Color', '#C9C9C9', 'HandleVisibility', 'off');
        xlim([(uniqueDist(1) -10) (uniqueDist(end)+10)]);
        ylim([150 700]);
        saveas(gcf, fullfile(filePath, append('Stacked ', stackPlot(jjj), '.png')));
        hold on;
    end
end
