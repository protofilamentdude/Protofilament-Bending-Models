%%2021_Dec_08_4_PF_Bending_Model_MgDataFit
%version 8 fits with reported forces (not nominal), and weights the data
%with the inverse of the variance for a given condition (Mg and load). 
%version 7 fits with nominalized forces
%version 6 fits the bovine amp vs. force data directly. 
%Version 5 fits the number of segments per PF rather than the PF height 
%version 4 implements a broader fitting function to get more robust error
%estimation for parameters. 

%% Fitting multi-PF model to data: 

% First import bovine load-dep data (+) Mg2+. 
% The data will be fit by using the unloaded amplitude (determined by 
% straight-line-fit to experimental data), and then the spring constant
% will be tuned such that the curves fit the data best. 

%extracting bovine data for input to fitting function:
[amp, CIs, loads, Mg, MgLineFits, allAmp, allLoad, allStd] = BovineCurves(bovine);


% Fitting
for i = 1:length(Mg); %fitting each Mg conc separately

    %estimating zero-force deflection with line fit to Mg data:
    ZeroFest = polyval(MgLineFits(i,:),0);
    %converting zero force amplitude estimates to curl height: 
    ZeroFest = Pulse2curlH(ZeroFest);

    %setting nonlinear fit inputs: 
x0 = [175,3.2]; %guess for torsional stiffness is 175 pN*nm/radian


xdata = allLoad{i}';
ydata = Pulse2curlH(allAmp{i})
weights = 1./(allStd{i}).^2'./min(1./(allStd{i}).^2); %weighting based on inverse of variance, normalized
weights(weights>2)=2; %minimizing weighting to 2-fold:
%squeezing weighting?

[beta, resnorm, J,CovB,~,output] = nlinfit(xdata, ydata,@(x,xdata) LoadCurve(x,xdata), x0,'Weights',weights);
J(isnan(J)) = eps;
paramCIArray{i} = nlparci(beta, resnorm,'jacobian', J); %calculate confidence interval for params
paramArray{i} = beta;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Storing previously determined fit values here: 
%PrevParams = {[129.35,14.02],[182.54,14.43],[161.93,17.21],[174.38,18.74]};
%using previously fitted parameters to generate fit curves (amp vs. load)

params = paramArray;
FtoFit = linspace(0,8,100);
for i = 1:length(Mg);
    CurlDeflections =  LoadCurve(params{i},FtoFit);
    FittedDeflections{i} = CurlH2pulse(CurlDeflections);
end

%plotting data: 
figure; hold on; cmap = lines; %new figure window and est. colormap
for i = 1:length(Mg);
errorbar(loads{i}, amp{i}, CIs{i},'s','MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:),'LineWidth',1.5);
MgAmpLegend{i} = strcat(num2str(Mg(i)), ' mM Mg^2^+');
end

%plotting fitted curves over data:
for i = 1:length(Mg);
      plot(FtoFit,FittedDeflections{i}, 'Color',cmap(i,:));
end

%Plotting parameters and CI's: 
%converting param cell array to array:
paramArrayMat = cell2mat(paramArray');

%plotting torsional stiffness per dimer vs. Mg
    figure; hold on;
    for i = 1:length(Mg);
        temp = paramCIArray{i}; %extracting relevant CI's
        ZAError(i,1:2) = temp(1,:)-paramArrayMat(i,1); %subtracting out parameter value to obtain relative error values
    end
    errorbar(Mg,paramArrayMat(:,1),ZAError(:,1),ZAError(:,2),[],[],'o','MarkerFaceColor','k','MarkerEdgeColor','k');
    xlabel('[Mg^2^+] (mM)');
    ylabel('Stiffness per dimer (pN*nm*)');

%plotting segment number vs. Mg
   figure; hold on;
    for i = 1:length(Mg);
        temp = paramCIArray{i}; %extracting relevant CI's
        kError(i,1:2) = temp(2,:)-paramArrayMat(i,2); %subtracting out parameter value to obtain relative error values
    end
    errorbar(Mg,paramArrayMat(:,2),kError(:,1),kError(:,2),[],[],'o','MarkerFaceColor','k','MarkerEdgeColor','k'); 
    xlabel('[Mg^2^+] (mM)');
    ylabel('Mean dimers per PF');


%% Functions

function [amp, CIs, loads, Mg, MgLineFits, allAmp, allLoad, allStd] = BovineCurves(bovine)
    %takes bovine wave data imported as a cell array, outputs amplitudes at tested 
    %forces (cell array), 95% CI's (student t) for amplitudes (cell array), tested forces, Mg2+
    %concentrations tested (array), the linear fit coefficients for each Mg2+
    %conc (polynomial coefficients) (array), all the amplitudes, all the
    %reported loads, and standard deviations for each requested load
    %force/Mg concentration. 
    %Import yeasty bovine wave data into table without header 'bovine', sorted in date
    %order. Be sure to compensate for column placement between different data
    %spreadsheets (specifically presence of baseline noise column)


    %Creating Index for Quantitative Quality (stdev estimate)
    Thresh = 3; %multiplier for baseline stdev to establish threshold
    %finds wave amplitudes that are greater than Thresh*(baseline stdev)
    BovineQuantQuality = (abs(cell2mat(bovine(:,5))) - Thresh.*cell2mat(bovine(:,6))) > 0;
    BovineQuantQualityIndex = find(BovineQuantQuality); %Finds row indeces of such amplitudes

    %Finds quality bovine waves as evaluated qualitatively
    BovineQualQuality = cell2mat(bovine(:,13));
    BovineQualQualityIndex = find(BovineQualQuality); %finds row indeces for such amplitudes

    %Comparing quantitative and qualitative measures for wave quality
    PercentQualityAgreement = 100.*sum(BovineQuantQuality == BovineQualQuality)/length(BovineQuantQuality);

        %How does varying the Thresh multiplier change percent quality
        %agreement?
        ThreshVar = linspace(1,5,15); %creating vector of threshold multipliers to test
        for j=1:length(ThreshVar);
            BovineQuantQualityVar{j} = (abs(cell2mat(bovine(:,5))) - ThreshVar(j).*cell2mat(bovine(:,6))) > 0;
            PercentQualityAgreementVar(j) = 100.*sum(BovineQuantQualityVar{j} == BovineQualQuality)/length(BovineQuantQuality);
        end
        figure;
        scatter(ThreshVar, PercentQualityAgreementVar, 'ko');
        xlabel('Threshold Multiplier Value');
        ylabel({'Percent Agreement', 'Qualitative vs. Quantitative Wave Quality'});
        set(gcf, 'Position',[12, 205, 370, 414]);


    %% Plotting wave amplitude vs. load force for various [Mg] tested
    %Finding indeces by Mg concentration
    MgTargets = [1, 6, 12, 20]; %set Mg concentrations to target for analysis and plotting. 
    MgConc = unique(cell2mat(bovine(:,12))); %find all Mg concentrations tested
    MgConc = intersect(MgConc, MgTargets); %determines set of Mg that are targeted and included in dataset. 
Mg = MgConc; %FUNCTION OUTPUT

    Preloads = [2,3,4,6]; %Define preloads tested
%loads = Preloads; %FUNCTION OUTPUT
    QualityIndex = BovineQualQualityIndex; %pick your favorite quality index
    cmap = lines; %establishing color map

    figure; hold on;
    xlabel('Load Force (pN)');
    ylabel('Wave Amplitude (nm)');
    MgBump = .15; %establish value for smearing Mg data across load axis to prevent overlap

    for i=1:length(MgConc);

        %finds indeces all waves of [Mg]=i, places in ith index of cell array
        MgIndex{i} = find(cellfun(@(x) MgConc(i) == x, bovine(:,12))); 
        individualAmps = [];
        individualLoads = [];
        individualStd = []; 
        for j=1:length(Preloads)

            PreloadIndex{j} = find(cellfun(@(x) abs((Preloads(j)-x))<0.45, bovine(:,3))); %finds indeces of waves with Preloads(j)+/-1 pN, stores indeces in ith index of cell array

            bovineMgAmp{j} = cell2mat(bovine(  intersect(intersect(QualityIndex, MgIndex{i}),PreloadIndex{j})  ,5)); %finds amplitudes of quality waves at [Mg](i) and Preloads(j), places them in ith index of cell array
            bovinePreloads{j} = cell2mat(bovine(  intersect(intersect(QualityIndex, MgIndex{i}),PreloadIndex{j})  ,3)); %finds preloads of quality waves at [Mg](i) and Preloads(j), places them in ith index of cell array
            bovinePercentQualMg(j) = length(intersect(QualityIndex, MgIndex{i}))./length(MgIndex{i});
            bovineMeanAmp(j) = nanmean(bovineMgAmp{j}); %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            bovineStdev(j) = nanstd(bovineMgAmp{j});
            bovineTscore(j) = tinv(.975,length(bovineMgAmp{j})-1);
            bovineCI(j) = bovineStdev(j).*bovineTscore(j)./sqrt(length(bovineMgAmp{j}));
            bovineMeanPreloads(j) = nanmean(bovinePreloads{j});
            %scatter(bovinePreloads{j} ,bovineMgAmp{j},'o','MarkerEdgeColor', cmap(i,:), 'MarkerFaceColor',cmap(i,:)); %plots 
            individualAmps = [individualAmps; bovineMgAmp{j}];
            individualLoads = [individualLoads; bovinePreloads{j}];
            individualStd = [individualStd; bovineStdev(j).*ones(length(bovineMgAmp{j}),1)];
        end
        loads{i} = bovineMeanPreloads;
        allAmpArray{i} = individualAmps;
        allLoadArray{i} = individualLoads;
        allStdArray{i} = individualStd;
        bovineStats(i,1:2) = {bovineMeanAmp, bovineCI}; %capture means and CI's for all preloads across a particular [Mg2+]

        errorbar(Preloads+MgBump*i, bovineMeanAmp, bovineCI,'s','MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:),'LineWidth',1.5);
        hold on; scatter(allAmpArray{i},allLoadArray{i},'.','MarkerFaceColor',cmap(i,:));
        idx = find(~isnan(bovineMeanAmp)); %determines non-NaN values of Mean Amplitudes for a given Mg concentration
        MgFitLines{i} = polyfit(Preloads(idx), bovineMeanAmp(idx),1); %generates linear fits to amplitude data per Mg concentration
    end
amp = bovineStats(:,1);
CIs = bovineStats(:,2);
MgLineFits = cell2mat(MgFitLines');
allAmp = allAmpArray;
allLoad = allLoadArray;
allStd = allStdArray;

    %Plot baseline noise and 95% CI in amplitude plot, getting baseline from all waves
    bovineBaselineNoise = cell2mat(bovine(:,6));
    bovineMeanBSN = nanmean(bovineBaselineNoise);
    bovineStdevBSN = nanstd(bovineBaselineNoise);
    bovineCI_BSN = 2.*bovineStdevBSN;
    bovineBaselineNoiseSorted = sort(bovineBaselineNoise); %sort baseline noise stdev's
    bovineBottom10BSN = mean(bovineBaselineNoiseSorted(1:round(length(bovineBaselineNoiseSorted)/10))); %take the average of the smallest 10% of baseline noises

    fplot(@(x) 3.*bovineMeanBSN, [min(Preloads) max(Preloads)+2],'k--'); %plots 3*baseline noise mean ?
    fplot(@(x) 3.*bovineBottom10BSN, [min(Preloads) max(Preloads)+2],'k:'); %plots 3*baseline noise mean ? for lowest 10% of noise values



    %Create legend for Amplitude plot
    for i = 1:length(MgConc)
        MgAmpLegend{i} = strcat(num2str(MgConc(i)), ' mM Mg^2^+');
    end
    MgAmpLegend{i+1} = '3*mean BSN \sigma';
    MgAmpLegend{i+2} = '3*mean lowest 10% BSN \sigma';

    legend(MgAmpLegend);
    set(gcf,'Position',[400, 72, 797, 547]);
    xlim([2, max(Preloads+2)]); %set x-axis limits for plot

    %plotting fit lines for select Mg curves:
    MgFit = [1, 6, 12, 20]; %define Mg conc's that we would like to fit lines to
    MgFit = intersect(MgFit, MgConc); %find the fit lines for Mg concentrations desired that are also represented in the dataset

        for k = 1:length(MgConc); %iterates through all Mg concentrations
            if ismember(MgConc(k),MgFit); %triggers following code only if the Mg considered is part of the set we would like to fit lines to
                plot( linspace(min(Preloads), max(Preloads),100)+ MgBump*k , polyval(MgFitLines{k},linspace(min(Preloads), max(Preloads),100))  ,'--','color',cmap(k,:));
            end 
        end
end

function amp = LoadCurve(x, Fdata)
                %inputs: x is a vector; x = [k, segN0, ang];
                %k is the torsional spring constant (pN*nm/rad), segN0 is 
                %the contour length, ang is the relaxed angle per dimer. 
                %ang is optional. If not specified, will be set to 23 deg. 
                
                %outputs: amp = distance of bead from MT wall (nm)
                

                %This model calculates force vs. deflection for a flat surface pushed down
                %onto a radial array of 4 protofilaments. The protofilaments are modeled as
                %series-linked torsion springs and are arranged assuming a 13-pf
                %microtubule. It is assumed there is no radial torsion along the PF axis
                %(that PF's are infinitely stiff out of their plane of bending). 
                
                % A segment is assumed to be 8.2 nm in length
                %to be 8.2 nm. One
                %parameter alpha gives the ratio of the 


                %script and figure setup:
                screensize = get( groot, 'Screensize' );

                %set global parameters
                    kGlobal = x(1); %pN*radian^-1 %declares a global k to call later. this is the torsional spring constant at each node (dimer interface)
                    rGlobal = 8.2; % (nm)
                    if length(x)>=3
                        if x(3)>0 & x(3)<40
                    aGlobal = -x(3)/180*pi; %(radians)
                        else
                            disp('Pick a relaxed angle 0<ang<40 (degrees)')
                        end
                    else
                    aGlobal = -23/180*pi; %(radians)
                    end
                    GlobalnptsF =151 ; %number of forces to try
                    GlobalFlims = [-4 10]; %force range to try

                    %Global parameters for distribution of PF's about MT: 
                    MTn = 13; %number of PF's per MT
                    MTr = 12.5; %(nm)   some expression of the MT radius (assuming each linkage
                    %is positioned at the outer surface of the MT)
                    MTa = 2*pi/MTn; %some expression of the angle between PF's (radians)
                    MTas = MTa.*[2.5 1.5 .5 .5 1.5 2.5]; %array of angles of PF's on a microtubule, assuming 4 PF's in model 
                    MTasAlt = MTa.*[2 1 0 1 2];
                    %MTas = MTasAlt; %temporarily change angles to alternative configuration
                    
                    %creating linearly spaced rotations of MT tip to
                    %average: 
                        for k = 1:11;
                            MTasLinspace{k,:}= abs([-3.5 -2.5 -1.5 -.5 .5 1.5 2.5 3.5]+(k-1).*0.1)
                        end
                    
                    MTasArray = {MTas; MTasAlt}; 
                    
                    FArray = [];

                
                
                %% 1-segment PF's model:

                k = kGlobal; %Nm/radian - arbitrary measure of tortional stiffness
                r = rGlobal % length of single lever arm
                ai = aGlobal; %relaxed angle between segments
                nptsF =GlobalnptsF ; %number of forces to try
                Flims = GlobalFlims; %force range to try
                F = linspace(Flims(1), Flims(2), nptsF);

                %setting non-linear torque balances as two elements of a function vector
                 T = @(a,F,k,r,ai) [r.*F.*cos(a(1)) - k.*(a(1)-ai)
                                      ];

                a0 = [0]; %initial guess for a0

                for i = 1:length(F); %solving the system of nonlinear equations at each force. 
                a1(i,:) = fsolve(@(a) T(a,F(i),k,r,ai), a0); %storing the resulting angles 
                end

                   subsetF1 = GlobalFlims;
                    groupedVar1 = [F' a1];
                    subsetVals1 = groupedVar1(F<subsetF1(2) & F>subsetF1(1),:);
                    %picking density of solutions to plot: 
                    rowStep =  1;
                    subsetVals1 = subsetVals1((1:floor(size(subsetVals1,1)./rowStep)).*rowStep,:); %sampling SubsetVals by interval = rowstep
                    %reassigning selected values to variables
                    F = subsetVals1(:,1);
                    a1 = subsetVals1(:,2);

                    for i = 1:length(F); %iterate through each force, plotting the linkage positions
                        linkage1(i,:) = [0 r.*cos(a1(i,1))   0 -r.*sin(a1(i,1)) ];
                        %plot(linkage2(i,1:3), linkage2(i,4:6),'Color',Fcmap(i,:));
                    end

                %% Simple model of two-lever tortional spring with initial angles of 23 deg: 

                %From the torque balance equations, we will find the vector of angles that
                %minimizes the sum of the squares of the torque equations, given a specific
                %force and initial angle estimate. 

                k = kGlobal; %Nm/radian - arbitrary measure of tortional stiffness
                r = rGlobal % length of single lever arm
                ai = aGlobal; %relaxed angle between segments
                nptsF =GlobalnptsF ; %number of forces to try
                Flims = GlobalFlims %force range to try
                F = linspace(Flims(1), Flims(2), nptsF);

                %setting non-linear torque balances as two elements of a function vector
                 T = @(a,F,k,r,ai) [r.*F.*cos(a(1)) + r.*F.*cos(a(1)+a(2)) - k.*(a(1)-ai)
                     r.*F.*cos(a(1)+a(2)) - k.*(a(2)-ai)];

                a0 = [0, 0]; %initial guess for a0

                for i = 1:length(F); %solving the system of nonlinear equations at each force. 
                a2(i,:) = fsolve(@(a) T(a,F(i),k,r,ai), a0); %storing the resulting angles 
                end

                %plot the numerical solution for angle one over the analytical solution to
                %test agreement; Don't need to inspect solutions for andgles for this
                %effort. 
                %scatter(F, a2(:,1)','*k');
                %legend('Analytical solution, 1st angle','Numerical solution, 1st angle');

                %     %figure; hold on;
                %     plot(F, a2(:,1)','k');
                %     plot(F, a2(:,2)','g');
                %     xlabel('Force (au)');
                %     ylabel('angle (radians)')
                %     legend('angle 1', 'angle 2');
                %     title('angles of series-linked tortion springs, numerical, n=2');
                %     set(gcf, 'Position', [screensize(3)*(1/5+1/4)+2.*18 screensize(4)/2 screensize(3)/4 screensize(4)/2.42]);

                %plotting the range of motion of the linkages: %Not necessary for model
                %fitting. 
                    %figure; hold on;
                    %cmap = colormap(jet);

                    %picking subset of forces to plot linkage positions for:
                    subsetF2 = GlobalFlims;
                    groupedVar2 = [F' a2];
                    subsetVals2 = groupedVar2(F<subsetF2(2) & F>subsetF2(1),:);
                    %picking density of solutions to plot: 
                    rowStep =  1;
                    subsetVals2 = subsetVals2((1:floor(size(subsetVals2,1)./rowStep)).*rowStep,:); %sampling SubsetVals by interval = rowstep
                    %reassigning selected values to variables
                    F = subsetVals2(:,1);
                    a2 = subsetVals2(:,2:3);
                    %assign colormap based on absolute value of force
                    maxF2 = max(abs(F)); minF2 = min(abs(F)); %finding max and min of force range
                    for i = 1:length(F);
                    %Fcmap(i,:) = cmap(floor(abs(F(i))./(maxF2-minF2).*(length(cmap)-1))+1,:)
                    end

                    for i = 1:length(F); %iterate through each force, plotting the linkage positions
                        linkage2(i,:) = [0 r.*cos(a2(i,1))  r.*cos(a2(i,1))+r.*cos(a2(i,1)+a2(i,2)) 0 -r.*sin(a2(i,1))  -r.*sin(a2(i,1))-r.*sin(a2(i,1)+a2(i,2))];
                        %plot(linkage2(i,1:3), linkage2(i,4:6),'Color',Fcmap(i,:));
                    end


                    %xlim([0 6]); ylim([-2 6]); set(gca,'DataAspectRatio',[1 1 1]); 
                    %xlabel('x-position of linkage (au)');
                    %ylabel({'y-position of linkage'; '(segment lengths)'});
                    %title({'Positions of 2-bar linkage'; 'various forces'});

                    %set(gcf, 'Position', [screensize(3)*(1/5+2/4)+3.*18 screensize(4)/2 screensize(3)/7 screensize(4)/2.42]);
                %% Simple model of three-lever tortional spring with initial angles of 23 deg (PART 2, numerical method):

                %Same strategy as above:

                k = kGlobal; %Nm/radian - arbitrary measure of tortional stiffness
                r = rGlobal % length of single lever arm
                ai = aGlobal; %relaxed angle between segments ~23 deg
                nptsF =GlobalnptsF ; %number of forces to try
                Flims = GlobalFlims %force range to try
                F = linspace(Flims(1), Flims(2), nptsF);

                %setting non-linear torque balances as three elements of a function vector
                 T = @(a,F,k,r,ai) [r.*F.*(cos(a(1)) + cos(a(1)+a(2)) + cos(a(1)+a(2)+a(3))) - k.*(a(1)-ai)
                     r.*F.*(cos(a(1)+a(2)) + cos(a(1)+a(2)+a(3))) - k.*(a(2)-ai)
                     r.*F.*cos(a(1)+a(2)+a(3)) - k.*(a(3)-ai)];

                a0 = [0, 0, 0]; %initial guess for a0

                for i = 1:length(F); %solving the system of nonlinear equations at each force. 
                a3(i,:) = fsolve(@(a) T(a,F(i),k,r,ai), a0); %storing the resulting angles 
                end

                %figure; hold on;
                for j = 1:length(a0);
                %plot(F, a3(:,j)');
                end
                %xlabel('Force (au)');
                %ylabel('angle (radians)')
                %legend('angle 1', 'angle 2', 'angle 3');
                %title('angles of series-linked tortion springs, numerical, n=3');
                %set(gcf, 'Position', [1 57 screensize(3)/4 screensize(4)/2.85]);

                %plotting the range of motion of the linkages: 
                   % figure; hold on;
                    %cmap = colormap(jet);

                    %picking subset of forces to plot linkage positions for:
                    subsetF3 = GlobalFlims;
                    groupedVar3 = [F' a3];
                    subsetVals3 = groupedVar3(F<subsetF3(2) & F>subsetF3(1),:);
                    %picking density of solutions to plot: 
                    rowStep =  1;
                    subsetVals3 = subsetVals3((1:floor(size(subsetVals3,1)./rowStep)).*rowStep,:); %sampling SubsetVals by interval = rowstep
                    %reassigning selected values to variables
                    F = subsetVals3(:,1);
                    a3 = subsetVals3(:,2:end);
                    %assign colormap based on absolute value of force
                    maxF3 = max(abs(F)); minF3 = min(abs(F)); %finding max and min of force range
                    for i = 1:length(F);
                     %Fcmap3(i,:) = cmap(floor(abs(F(i))./(maxF3-minF3).*(length(cmap)-1))+1,:)
                    end


                    for i = 1:length(F); %iterate through each force, plotting the linkage positions
                        xpos = [0 r.*cos(a3(i,1))  r.*cos(a3(i,1))+r.*cos(a3(i,1)+a3(i,2)) r.*cos(a3(i,1))+r.*cos(a3(i,1)+a3(i,2))+r.*cos(a3(i,1)+a3(i,2)+a3(i,3))];
                        ypos = [0 r.*sin(a3(i,1))  r.*sin(a3(i,1))+r.*sin(a3(i,1)+a3(i,2)) r.*sin(a3(i,1))+r.*sin(a3(i,1)+a3(i,2))+r.*sin(a3(i,1)+a3(i,2)+a3(i,3))];
                        linkage3(i,:) = [xpos -ypos];
                        %plot(linkage3(i,1:4), linkage3(i,5:8),'Color',Fcmap3(i,:));
                    end
                    %xlabel('x-position of linkage (au)');
                   % ylabel({'y-position of linkage'; '(segment lengths)'});
                    %title({'Positions of 3-bar linkage'; 'various forces.'})
                   % xlim([0 6]); ylim([-2 6]); set(gca,'DataAspectRatio',[1 1 1]);

                    %set(gcf, 'Position', [screensize(3)/4+18 57 screensize(3)/7 screensize(4)/2.85]);

                %% Simple model of four-lever tortional spring with initial angles of 23 deg 
                %Same strategy as above:

                k = kGlobal; %Nm/radian - arbitrary measure of tortional stiffness
                r = rGlobal % length of single lever arm
                ai = aGlobal; %relaxed angle between segments ~23 deg
                nptsF =GlobalnptsF ; %number of forces to try
                Flims = GlobalFlims; %force range to try
                F = linspace(Flims(1), Flims(2), nptsF);

                %setting non-linear torque balances as three elements of a function vector
                 T = @(a,F,k,r,ai) [r.*F.*(cos(a(1)) + cos(a(1)+a(2)) + cos(a(1)+a(2)+a(3)) + cos(a(1)+a(2)+a(3)+a(4))) - k.*(a(1)-ai)
                     r.*F.*(cos(a(1)+a(2)) + cos(a(1)+a(2)+a(3))+cos(a(1)+a(2)+a(3)+a(4))) - k.*(a(2)-ai)
                     r.*F.*(cos(a(1)+a(2)+a(3))+cos(a(1)+a(2)+a(3)+a(4))) - k.*(a(3)-ai)
                     r.*F.*(cos(a(1)+a(2)+a(3)+a(4))) - k.*(a(4)-ai)];

                a0 = [0, 0, 0, 0]; %initial guess for a0

                for i = 1:length(F); %solving the system of nonlinear equations at each force. 
                a4(i,:) = fsolve(@(a) T(a,F(i),k,r,ai), a0); %storing the resulting angles 
                end

                %plotting angles at each node: (don't need this part of the model plotted
                    % %here
                    % figure; hold on;
                    % for j = 1:length(a0);
                    % plot(F, a4(:,j)');
                    % end
                    % xlabel('Force (au)');
                    % ylabel('angle (radians)')
                    % legend('angle 1', 'angle 2', 'angle 3', 'angle 4');
                    % title('angles of series-linked tortion springs, numerical, n=4');
                    % set(gcf, 'Position', [screensize(3)*(1/4+1/7) 57 screensize(3)/4 screensize(4)/2.85]);

                %plotting the range of motion of the linkages: 
                    %figure; hold on;
                    %cmap = colormap(jet);

                    %picking subset of forces to plot linkage positions for:
                    subsetF4 = GlobalFlims;
                    groupedVar4 = [F' a4];
                    subsetVals4 = groupedVar4(F<subsetF4(2) & F>subsetF4(1),:);
                    %picking density of solutions to plot: 
                    rowStep =  1;
                    subsetVals4 = subsetVals4((1:floor(size(subsetVals4,1)./rowStep)).*rowStep,:); %sampling SubsetVals by interval = rowstep
                    %reassigning selected values to variables
                    F = subsetVals4(:,1);
                    a4 = subsetVals4(:,2:end);
                    %assign colormap based on absolute value of force
                    maxF4 = max(abs(F)); minF4 = min(abs(F)); %finding max and min of force range
                    for i = 1:length(F);
                %    Fcmap4(i,:) = cmap(floor(abs(F(i))./abs(maxF4-minF4).*(length(cmap)-1)+1),:)
                    end


                    for i = 1:length(F); %iterate through each force, plotting the linkage positions
                        xpos = [0 
                                r.*cos(a4(i,1))  
                                r.*cos(a4(i,1))+r.*cos(a4(i,1)+a4(i,2))
                                r.*cos(a4(i,1))+r.*cos(a4(i,1)+a4(i,2))+r.*cos(a4(i,1)+a4(i,2)+a4(i,3))
                                r.*cos(a4(i,1))+r.*cos(a4(i,1)+a4(i,2))+r.*cos(a4(i,1)+a4(i,2)+a4(i,3))+r.*cos(a4(i,1)+a4(i,2)+a4(i,3)+a4(i,4))
                                ];
                        ypos = [0 
                                r.*sin(a4(i,1))  
                                r.*sin(a4(i,1))+r.*sin(a4(i,1)+a4(i,2))
                                r.*sin(a4(i,1))+r.*sin(a4(i,1)+a4(i,2))+r.*sin(a4(i,1)+a4(i,2)+a4(i,3))
                                r.*sin(a4(i,1))+r.*sin(a4(i,1)+a4(i,2))+r.*sin(a4(i,1)+a4(i,2)+a4(i,3))+r.*sin(a4(i,1)+a4(i,2)+a4(i,3)+a4(i,4))
                                ];
                        linkage4(i,:) = [xpos' -ypos'];
                        %plot(linkage4(i,1:5), linkage4(i,6:10),'Color',Fcmap4(i,:));
                    end
                %     xlabel('x-position of linkage (au)');
                %     ylabel({'y-position of linkage'; '(segment lengths)'});
                %     title({'Positions of 4-bar linkage'; 'various forces.'})
                %     xlim([0 6]); ylim([-2 6]); set(gca,'DataAspectRatio',[1 1 1]);
                %     
                %     set(gcf, 'Position', [screensize(3)/4+18 57 screensize(3)/7 screensize(4)/2.85]);
                %% Simple model of five-lever tortional spring with initial angles of zero 
                %Same strategy as above:

                k = kGlobal; %Nm/radian - arbitrary measure of tortional stiffness
                r = rGlobal % length of single lever arm
                ai = aGlobal; %relaxed angle between segments ~23 deg
                nptsF =GlobalnptsF ; %number of forces to try
                Flims = GlobalFlims; %force range to try
                F = linspace(Flims(1), Flims(2), nptsF);

                %setting non-linear torque balances as three elements of a function vector
                 T = @(a,F,k,r,ai) [r.*F.*(cos(a(1)) + cos(a(1)+a(2)) + cos(a(1)+a(2)+a(3)) + cos(a(1)+a(2)+a(3)+a(4))+cos(a(1)+a(2)+a(3)+a(4)+a(5))) - k.*(a(1)-ai)
                     r.*F.*(cos(a(1)+a(2)) + cos(a(1)+a(2)+a(3))+cos(a(1)+a(2)+a(3)+a(4))+cos(a(1)+a(2)+a(3)+a(4)+a(5))) - k.*(a(2)-ai)
                     r.*F.*(cos(a(1)+a(2)+a(3))+cos(a(1)+a(2)+a(3)+a(4))+cos(a(1)+a(2)+a(3)+a(4)+a(5))) - k.*(a(3)-ai)
                     r.*F.*(cos(a(1)+a(2)+a(3)+a(4))+cos(a(1)+a(2)+a(3)+a(4)+a(5))) - k.*(a(4)-ai)
                     r.*F.*(cos(a(1)+a(2)+a(3)+a(4)+a(5))) - k.*(a(5)-ai)]

                a0 = [0, 0, 0, 0, 0]; %initial guess for a0

                for i = 1:length(F); %solving the system of nonlinear equations at each force. 
                a5(i,:) = fsolve(@(a) T(a,F(i),k,r,ai), a0); %storing the resulting angles 
                end

                % %plotting angles at each node: (don't need this part of the model plotted
                    % figure; hold on;
                    % for j = 1:length(a0);
                    % plot(F, -a5(:,j)');
                    % end
                    % xlabel('Force (au)');
                    % ylabel('angle (radians)')
                    % legend('angle 1', 'angle 2', 'angle 3', 'angle 4', 'angle 5');
                    % title('angles of series-linked tortion springs, numerical, n=5');
                    % set(gcf, 'Position', [screensize(3)*(1/4+1/7) 57 screensize(3)/4 screensize(4)/2.85]);

                %plotting the range of motion of the linkages: 
%                     figure; hold on;
                     cmap = colormap(jet);

                    %picking subset of forces to plot linkage positions for:
                    subsetF5 = GlobalFlims;
                    groupedVar5 = [F' a5];
                    subsetVals5 = groupedVar5(F<subsetF5(2) & F>subsetF5(1),:);
                    %picking density of solutions to plot: 
                    rowStep =  1;
                    subsetVals5 = subsetVals5((1:floor(size(subsetVals5,1)./rowStep)).*rowStep,:); %sampling SubsetVals by interval = rowstep
                    %reassigning selected values to variables
                    F = subsetVals5(:,1);
                    a5 = subsetVals5(:,2:end);
                    %assign colormap based on absolute value of force
                    maxF5 = max(abs(F)); minF5 = min(abs(F)); %finding max and min of force range
                    for i = 1:length(F);
                    Fcmap5(i,:) = cmap(floor(abs(F(i))./(maxF5-minF5).*(length(cmap)-1)+1),:);
                    end


                    for i = 1:length(F); %iterate through each force, plotting the linkage positions
                        xpos = [0 
                                r.*cos(a5(i,1))  
                                r.*cos(a5(i,1))+r.*cos(a5(i,1)+a5(i,2))
                                r.*cos(a5(i,1))+r.*cos(a5(i,1)+a5(i,2))+r.*cos(a5(i,1)+a5(i,2)+a5(i,3))
                                r.*cos(a5(i,1))+r.*cos(a5(i,1)+a5(i,2))+r.*cos(a5(i,1)+a5(i,2)+a5(i,3))+r.*cos(a5(i,1)+a5(i,2)+a5(i,3)+a5(i,4))
                                r.*cos(a5(i,1))+r.*cos(a5(i,1)+a5(i,2))+r.*cos(a5(i,1)+a5(i,2)+a5(i,3))+r.*cos(a5(i,1)+a5(i,2)+a5(i,3)+a5(i,4))+r.*cos(a5(i,1)+a5(i,2)+a5(i,3)+a5(i,4)+a5(i,5))
                                ];
                        ypos = [0 
                                r.*sin(a5(i,1))  
                                r.*sin(a5(i,1))+r.*sin(a5(i,1)+a5(i,2))
                                r.*sin(a5(i,1))+r.*sin(a5(i,1)+a5(i,2))+r.*sin(a5(i,1)+a5(i,2)+a5(i,3))
                                r.*sin(a5(i,1))+r.*sin(a5(i,1)+a5(i,2))+r.*sin(a5(i,1)+a5(i,2)+a5(i,3))+r.*sin(a5(i,1)+a5(i,2)+a5(i,3)+a5(i,4))
                                r.*sin(a5(i,1))+r.*sin(a5(i,1)+a5(i,2))+r.*sin(a5(i,1)+a5(i,2)+a5(i,3))+r.*sin(a5(i,1)+a5(i,2)+a5(i,3)+a5(i,4))+r.*sin(a5(i,1)+a5(i,2)+a5(i,3)+a5(i,4)+a5(i,5))
                                ];
                        linkage5(i,:) = [xpos' -ypos'];
%                         plot(linkage5(i,1:6), linkage5(i,7:12),'Color',Fcmap5(i,:));
                    end
                %     xlabel('x-position of linkage (au)');
                %     ylabel({'y-position of linkage'; '(segment lengths)'});
                %     title({'Positions of 5-bar linkage'; 'various forces.'})
                %     xlim([0 6]); ylim([-2 6]); set(gca,'DataAspectRatio',[1 1 1]);
                %     
                %     set(gcf, 'Position', [screensize(3)/4+18 57 screensize(3)/7 screensize(4)/2.85]);

                    %Calculating and plotting the energy in each torsion spring: 

                    subsetEF5 = [0 subsetF5(2)]; %selecting the subset of values for positive forces. 
                    subsetEVals5 = groupedVar5(groupedVar5(:,1)<subsetEF5(2) & groupedVar5(:,1)>subsetEF5(1),:);
                    a5 = -subsetEVals5(:,2:end); %selects angles only. change angle sign to reflect model sketch. 
                    ai = -ai; %change sign of ai as well to reflect sketch
                    for j = 1:size(a5,2);
                        Fval(:,j) = -k.*(a5(:,j)-ai);
                        E5(:,j) = cumtrapz(-(a5(:,j)-ai), Fval(:,j) ); %calculating the energy for each angle. 
                    end

                %     figure; hold on; 
                %     for i = 1:size(a5,2);
                %     plot(subsetEVals5(:,1), E5(:,i));
                %     end
                %     
                %     xlabel('Force (au)');
                %     ylabel('Energy (au)')
                %     legend('angle 1', 'angle 2', 'angle 3', 'angle 4', 'angle 5');
                %     title('Energy of of series-linked tortion springs, numerical, n=5');
                %     set(gcf, 'Position', [screensize(3)*(1/4+1/7) 57 screensize(3)/4 screensize(4)/2.85]);

                    %% Plotting y-position (amplitude) vs. force for a variety of segment lengths
                    %Ensure the unloaded angles are set to the same value across the
                    %segment-length models. 23 deg seems to be a consensus value
                    %(Gudimchuck et al 2020)

                    %For now we're just going to manually plot from the different models. 
%                     figure; hold on; %plotting y ampliude vs. force: 
%                     %finding segment-terminal y-positions for only positive force values
%                     %(no pulling up on segments). 
% 
%                     plot(subsetVals1(subsetVals1(:,1)>=0,1), linkage1(subsetVals1(:,1)>=0, end));
%                     plot(subsetVals2(subsetVals2(:,1)>=0,1), linkage2(subsetVals2(:,1)>=0, end));
%                     plot(subsetVals3(subsetVals3(:,1)>=0,1), linkage3(subsetVals3(:,1)>=0, end));
%                     plot(subsetVals4(subsetVals4(:,1)>=0,1), linkage4(subsetVals4(:,1)>=0, end));
%                     plot(subsetVals5(subsetVals5(:,1)>=0,1), linkage5(subsetVals5(:,1)>=0, end));
%                     xlabel('Force (pN)');
%                     ylabel({'Protofilament deflection'; 'in direction of force (nm)'});
%                     title({'Protofilament deflection vs. Force','for single PF''s of various lengths'});
%                     legend('1-segment','2-segment', '3-segment', '4-segment', '5-segment');
%                     xlim([0 12]); ylim([0 4.*rGlobal]);
%                     set(gcf, 'Position', [924   530   298   282]);
                %% 4-PF model of deflection: 
                    %We will be using the force-deflection relationship calculated earlier
                    %for single protofilaments in the geometrical context of a 13-PF MT to
                    %generate the force-deflection profile for the MT tip at large. 
            for ij = 1:size(MTasArray,1); %calculating across different MT angles
                MTas = MTasArray{ij};
                    %1-segment PF's
                        %first find the unloaded height of the PF relative to the MT axis: 
                        ul1 = interp1(subsetVals1(:,1), linkage1(:, end),0);
                        ul1H = cos(MTas).*(MTr + ul1); %defines an array containing the y-axis height of each PF relative to the MT axis 
                        %create vector of y-axis heights to probe: 
                        Bh1 = linspace(max(ul1H), MTr, 100); %range from furthest PF tip to MT lattice

                        %calculating force at each bead height:
                        for i = 1:length(Bh1)
                            for j = 1:length(ul1H) %iterating over number of PF's being modeled
                                F1subunit(i,j) = [
                                    heaviside(-(Bh1(i)-ul1H(j))).*interp1(linkage1(:,end), subsetVals1(:,1), ( Bh1(i)-MTr.*cos(MTas(j)) )./cos(MTas(j))   ,'linear',0) .*cos(MTas(j));
                                    ];
                            end
                        end
                        F1Tot = sum(F1subunit,2);



                    %2-segment PF's
                            %first find the unloaded height of the PF relative to the MT axis: 
                        ul2 = interp1(subsetVals2(:,1), linkage2(:, end),0);
                        ul2H = cos(MTas).*(MTr + ul2); %defines an array containing the y-axis height of each PF relative to the MT axis 
                        %create vector of y-axis heights to probe: 
                        Bh2 = linspace(max(ul2H), MTr, 100); %range from furthest PF tip to MT lattice

                        %calculating force at each bead height:
                        for i = 1:length(Bh2)
                            for j = 1:length(ul2H) %iterating over number of PF's being modeled
                                F2subunit(i,j) = [
                                    heaviside(-(Bh2(i)-ul2H(j))).*interp1(linkage2(:,end), subsetVals2(:,1), ( Bh2(i)-MTr.*cos(MTas(j)) )./cos(MTas(j))  ,'linear',0 ) .*cos(MTas(j));
                                    ];
                            end
                        end
                        F2Tot = sum(F2subunit,2);
                    %3-segment PF's
                                %first find the unloaded height of the PF relative to the MT axis: 
                        ul3 = interp1(subsetVals3(:,1), linkage3(:, end),0);
                        ul3H = cos(MTas).*(MTr + ul3); %defines an array containing the y-axis height of each PF relative to the MT axis 
                        %create vector of y-axis heights to probe: 
                        Bh3 = linspace(max(ul3H), MTr, 100); %range from furthest PF tip to MT lattice

                        %calculating force at each bead height:
                        for i = 1:length(Bh3)
                            for j = 1:length(ul3H) %iterating over number of PF's being modeled
                                F3subunit(i,j) = [
                                    heaviside(-(Bh3(i)-ul3H(j))).*interp1(linkage3(:,end), subsetVals3(:,1), ( Bh3(i)-MTr.*cos(MTas(j)) )./cos(MTas(j))   ,'linear',0) .*cos(MTas(j));
                                    ];
                            end
                        end
                        F3Tot = sum(F3subunit,2);
                    %4-segment PF's
                                %first find the unloaded height of the PF relative to the MT axis: 
                        ul4 = interp1(subsetVals4(:,1), linkage4(:, end),0);
                        ul4H = cos(MTas).*(MTr + ul4); %defines an array containing the y-axis height of each PF relative to the MT axis 
                        %create vector of y-axis heights to probe: 
                        Bh4 = linspace(max(ul4H), MTr, 100); %range from furthest PF tip to MT lattice

                        %calculating force at each bead height:
                        for i = 1:length(Bh4)
                            for j = 1:length(ul4H) %iterating over number of PF's being modeled
                                F4subunit(i,j) = [
                                    heaviside(-(Bh4(i)-ul4H(j))).*interp1(linkage4(:,end), subsetVals4(:,1), ( Bh4(i)-MTr.*cos(MTas(j)) )./cos(MTas(j))   ,'linear',0) .*cos(MTas(j));
                                    ];
                            end
                        end
                        F4Tot = sum(F4subunit,2);
                    %5-segment PF's
                                %first find the unloaded height of the PF relative to the MT axis: 
                        ul5 = interp1(subsetVals5(:,1), linkage5(:, end),0);
                        ul5H = cos(MTas).*(MTr + ul5); %defines an array containing the y-axis height of each PF relative to the MT axis 
                        %create vector of y-axis heights to probe: 
                        Bh5 = linspace(max(ul5H), MTr, 100); %range from furthest PF tip to MT lattice

                        %calculating force at each bead height:
                        for i = 1:length(Bh5)
                            for j = 1:length(ul5H) %iterating over number of PF's being modeled
                                F5subunit(i,j) = [
                                    heaviside(-(Bh5(i)-ul5H(j))).*interp1(linkage5(:,end), subsetVals5(:,1), ( Bh5(i)-MTr.*cos(MTas(j)) )./cos(MTas(j))   ,'linear',0) .*cos(MTas(j));
                                    ];
                            end
                        end
                        F5Tot = sum(F5subunit,2);
                %Plotting force vs. deflection for all 4-PF models: 
                
                %%%%%%%%%%% Not plotting for purpose of fitting. 
%                 figure; hold on; 
%                     plot(F1Tot, Bh1-MTr);
%                     plot(F2Tot, Bh2-MTr);
%                     plot(F3Tot, Bh3-MTr);
%                     plot(F4Tot, Bh4-MTr);
%                     plot(F5Tot, Bh5-MTr);
% 
%                     xlabel('Force (pN)');
%                     ylabel({'Bead Height from MT wall (nm)'});
%                     title({'Deflection vs. Force for PF curls at MT tips','with PF''s of 1-5 segments'});
%                     legend('1-segment PFs','2-segment PFs', '3-segment PFs', '4-segment PFs', '5-segment PFs');
%                     xlim([0 12]); ylim([0 4.*rGlobal]);
%                     set(gcf,'Position', [1240 530 297 283]);
                  FArray = cat(3,FArray,[F1Tot'; F2Tot'; F3Tot'; F4Tot'; F5Tot']);
            end
            BhArray = [Bh1; Bh2; Bh3; Bh4; Bh5];
            FMean = mean(FArray,3);
            
            %Plotting force vs. deflection for mean of all 4-PF models considered: 
                           %%%%%%%%%%% Not plotting for purpose of fitting. 
                figure; hold on; 
                    plot(FMean(1,:), Bh1-MTr);
                    plot(FMean(2,:), Bh2-MTr);
                    plot(FMean(3,:), Bh3-MTr);
                    plot(FMean(4,:), Bh4-MTr);
                    plot(FMean(5,:), Bh5-MTr);

                    xlabel('Force (pN)');
                    ylabel({'Bead Height from MT wall (nm)'});
                    title({'Deflection vs. Force for PF curls at MT tips','with PF''s of 1-5 segments'});
                    legend('1-segment PFs','2-segment PFs', '3-segment PFs', '4-segment PFs', '5-segment PFs');
                    xlim([0 12]); ylim([0 4.*rGlobal]);
                    set(gcf,'Position', [1240 530 297 283]);
    

%get deflections for two curves with zeroF amplitude nearest input zeroF
    %data I'm looking for is in FMean and Bh<n>. pass these through dummy
    %variables to allow for selection of curves other than mean btw MT
    %rotations.
    Ffit = FMean; %getting Force points from model
    
    %getting zero-force deflection for each curve using interp: 
    for i = 1:size(Ffit,1)
        ZeroAmp(i) = interp1(Ffit(i,:),BhArray(i,:),0,'pchip')-MTr;
    end
    disp('Zero Amps are:'); disp(ZeroAmp);%%%%%%%%%%%%%%%%%%%
    %set alpha as ratio of higher curve to lower curve
        %find indeces of theoretical unloaded amps above and below query
        %unloaded amp:
        segNumber = x(2); %sets segment number to input parameter 2
        %how to get from segment number to zero amp query? 
        ZeroAmpQuery = interp1(1:5,ZeroAmp,segNumber);
        FitBoundIndeces = [max(find(ZeroAmp<ZeroAmpQuery)) min(find(ZeroAmp>ZeroAmpQuery))];
        disp('The fit bound index is:'); disp(FitBoundIndeces); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        alpha = (ZeroAmpQuery-ZeroAmp(FitBoundIndeces(1)))./(ZeroAmp(FitBoundIndeces(2))-ZeroAmp(FitBoundIndeces(1)));
        %construct a weighted sum with alpha

        %inverting independent and dependent variables to allow us to
        %average across the same force vector: 
        %constructing a force vector: 
        TheoreticalForces = linspace(0,min([max(Ffit(FitBoundIndeces(1),:)), max(Ffit(FitBoundIndeces(2),:))]),100);
        ResampledDispHigh = interp1(Ffit(FitBoundIndeces(2),:), BhArray(FitBoundIndeces(2),:)-MTr,TheoreticalForces);
        ResampledDispLow = interp1(Ffit(FitBoundIndeces(1),:), BhArray(FitBoundIndeces(1),:)-MTr,TheoreticalForces);
        
        TheoreticalDeflections = alpha.*ResampledDispHigh + (1-alpha).*ResampledDispLow; 
amp = interp1(TheoreticalForces, TheoreticalDeflections,Fdata);
plot(TheoreticalForces, TheoreticalDeflections,'r');
end

function height = Pulse2curlH(p);
%This script generates an analytical model of the wave assay,
%showing the relationship between curl height (h), bead size (r), tether
%length (tau), and measured wave amplitude (a) 

%The geometry of loaded bead, tether, and microtubule define a leverage
%system with the following rlationship between variables: 

%sqrt(tau*(tau+2*r)) = (h^2 + 2*r*h + a^2)/(2*a)
         g = @(tau, r, h, a) sqrt(tau.*(tau+2.*r))-(h.^2 + 2.*r.*h + a^2)./(2.*a);

        tau = 36; %tether length
        r = 220; %bead radius
        height = zeros(1,length(p)); %initializing height vector



    for i = 1:length(p);
        if isnan(p(i));
            dummy(i) = NaN;
        
        else
        h0 = p(i); %initial guess for numeric solve for h
        a = p(i); %set input
        dummy(i) = fsolve(@(h) g(tau, r, h, a),[h0]);
        end
    end
    height = dummy;
end

function pulse = CurlH2pulse(height);
%This script generates an analytical model of the wave assay,
%showing the relationship between curl height (h), bead size (r), tether
%length (tau), and measured wave amplitude (a) 

%The geometry of loaded bead, tether, and microtubule define a leverage
%system with the following rlationship between variables: 

%sqrt(tau*(tau+2*r)) = (h^2 + 2*r*h + a^2)/(2*a)
         g = @(tau, r, h, a) sqrt(tau.*(tau+2.*r))-(h.^2 + 2.*r.*h + a^2)./(2.*a);

        tau = 36; %tether length
        r = 220; %bead radius
        pulse = zeros(1,length(height)); %initializing height vector



    for i = 1:length(height);
        if isnan(height(i));
            dummy(i) = NaN;
        
        else
        a0 = height(i); %initial guess for numeric solve for a
        h = height(i); %set input, reasonable guess is h
        dummy(i) = fsolve(@(a) g(tau, r, h, a),[a0]);
        end
    end
    pulse = dummy;
end
