%% HyPASS OBJECTIVE 1 ONLY: exposed FFPE tissue thickness versus SWIR signal
% Direct block-face illumination: exposed paraffin-INFILTRATED tissue, NO top
% paraffin layer; horizontal tissue-to-bulk-paraffin interface below tissue.
% 2016-reference versus higher-scattering FFPE SENSITIVITY scenarios.
% This file runs independently with HyPASS_Wavelength_Optical_Properties.csv.
% MATLAB R2020b+; no toolboxes needed; NO illumination-angle simulation.
%
% Source: Zhang et al., J Biomed Opt 2016, 21(12), 126006,
% https://doi.org/10.1117/1.JBO.21.12.126006
% Reference scenario uses Fig. 9 optical inputs PREVIOUSLY adopted from
% rat-skin literature; NOT measurements of FFPE tissue or chicken tissue.
% The reference mu_s' values below are discrete approximate reported inputs;
% verify against the paper before publication.
%
% Model limitations: plane-parallel infinite tissue layer, assumed optical
% properties (including n), no lateral tissue/paraffin heterogeneity, detector
% modeled as normal cone NA, no microscope PSF / spatial resolution modeling.

clear; close all; clc;
rng(61,'twister');
HERE=fileparts(mfilename('fullpath')); if isempty(HERE), HERE=pwd; end

% Output folder outside Dropbox
OUT = fullfile(tempdir, 'HyPASS_MC_Results', datestr(now,'yyyymmdd_HHMMSS_FFF'));
% Unique run folder avoids overwriting PDFs that may still be open.

if ~exist(OUT, 'dir')
    mkdir(OUT);
end

fprintf('Results will be saved to:\n%s\n', OUT);
csvfile=fullfile(HERE,'HyPASS_Wavelength_Optical_Properties.csv');
assert(isfile(csvfile),'Missing HyPASS_Wavelength_Optical_Properties.csv');
S=readtable(csvfile);

%% PARAMETERS — edit these
p.lambda_nm = [900 1300 1500 1600]; % exact reference knots; 900–1600 nm only
p.thickness_um = [5 10 25 50 75 100 150 200 300 500 750 1000 1500 2000 3500 5000];
p.N = 500000;                    % incident photons PER wavelength/thickness/scenario
p.maxEvents = 15000;            % per-photon interaction safety cap
p.NA = 0.20;                    % assumed camera collection NA
p.paraffinDepth_mm = 10;        % layer beneath exposed tissue
p.nAir = 1.00;
p.nTissue = 1.44;              % illustrative paraffin-infiltrated tissue
p.nParaffin = 1.46;            % illustrative bulk embedding paraffin
p.gTissue = 0.90;              % common value in BOTH scenarios to isolate mu_s'
p.gParaffin = 0.85;
p.ffpeWaterFraction = 0.10;     % assumed, NOT measured
p.ffpeBackgroundAbs_cm = 0.20;
p.paraffinAbs_cm = 0.30;        % assumed, NOT measured
p.paraffinScatterPrime_cm = 0.50; % assumed, NOT measured
p.highScatterA500_cm = 20;     % HIGH-SCATTERING sensitivity comparator
p.highScatterExponent = 1.20;

% Coefficients used in Zhang et al. (2016), Fig. 9, units cm^-1.
% Refer to paper for provenance; these are not direct FFPE measurements.
refLambda = [900 1300 1500 1600];
refMuSprime = [1.30 1.07 0.95 0.85];

assert(all(p.lambda_nm>=900 & p.lambda_nm<=1600), ...
 'Reference interpolation is restricted to 900–1600 nm; do not extrapolate to 600–800 nm.');
assert(all(p.thickness_um>0),'Thickness must be positive for log X-axis.');
assert(all(p.NA>0 & p.NA<1),'NA must be between zero and one.');
waterMuA = interp1(S.Wavelength_nm,S.Measured_water_muA_cm,p.lambda_nm,'pchip',NaN);
assert(all(isfinite(waterMuA)),'Missing water absorption values at one or more wavelengths.');

% Same provisional FFPE absorption in both scenarios to isolate scattering.
muAT = p.ffpeBackgroundAbs_cm + p.ffpeWaterFraction * waterMuA;
muAP = p.paraffinAbs_cm * ones(size(p.lambda_nm));
muSP = (p.paraffinScatterPrime_cm/(1-p.gParaffin))*ones(size(p.lambda_nm));
refSprime = interp1(refLambda,refMuSprime,p.lambda_nm,'pchip',NaN);
highSprime = p.highScatterA500_cm*(p.lambda_nm/500).^(-p.highScatterExponent);
scenarioSprime = [refSprime(:), highSprime(:)];
scenarioName = {'Zhang2016_reference','Higher_scattering_assumption'};

%% OBJECTIVE 1: thickness sweep only — TWO scattering scenarios
% Model is at NORMAL incidence: no angle sweep in this file.
rows=zeros(numel(p.lambda_nm)*numel(p.thickness_um)*2,14);
r=0;
for iscenario=1:2
    for iw=1:numel(p.lambda_nm)
        muST=scenarioSprime(iw,iscenario)/(1-p.gTissue);
        for j=1:numel(p.thickness_um)
            dT=p.thickness_um(j)*1e-4; % um -> cm
            x=mcLayered(p.N,0,dT,p.paraffinDepth_mm/10,...
                muAT(iw),muST,p.gTissue,p.nTissue,...
                muAP(iw),muSP(iw),p.gParaffin,p.nParaffin,p);
            r=r+1;
            rows(r,:)=[iscenario,p.lambda_nm(iw),p.thickness_um(j),...
                scenarioSprime(iw,iscenario),x];
        end
        fprintf('Scenario %d/2, %d nm completed\n',iscenario,p.lambda_nm(iw));
    end
end
T=array2table(rows,'VariableNames',{'ScenarioID','Wavelength_nm','TissueThickness_um',...
    'Tissue_muSprime_cm','TotalDiffuse','CameraDiffuse',...
    'CameraDiffuse_EnteredParaffin','EscapedBottom','AbsorbedInTissue',...
    'AbsorbedInParaffin','TopSurfaceSpecular',...
    'MedianMaximumDepthDetected_um','P90MaximumDepthDetected_um','EnergyClosure'});
T.Scenario=categorical(reshape(scenarioName(T.ScenarioID),[],1));
T.CameraDiffuseSE=sqrt(max(0,T.CameraDiffuse.*(1-T.CameraDiffuse))./p.N);
writetable(T,fullfile(OUT,'HyPASS_Objective1_Thickness_Comparison.csv'));

Optics=table(p.lambda_nm(:),muAT(:),refSprime(:),highSprime(:),...
 'VariableNames',{'Wavelength_nm','Tissue_muA_assumed_cm','Reference_muSprime_cm','Higher_muSprime_cm'});
writetable(Optics,fullfile(OUT,'HyPASS_Objective1_Optical_Inputs.csv'));

%% FIGURE 1: X=log thickness; Y=LINEAR camera-collected diffuse fraction
fig=figure('Color','w','Position',[50 50 1100 500]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
c=lines(numel(p.lambda_nm));
for iscenario=1:2
    nexttile; hold on;
    for iw=1:numel(p.lambda_nm)
        idx=T.ScenarioID==iscenario & T.Wavelength_nm==p.lambda_nm(iw);
        semilogx(T.TissueThickness_um(idx),T.CameraDiffuse(idx),'-o',...
            'Color',c(iw,:),'LineWidth',1.6,'MarkerSize',4,...
            'DisplayName',sprintf('%d nm',p.lambda_nm(iw)));
    end
    set(gca,'XScale','log','YScale','linear','Layer','top');
    xlim([min(p.thickness_um),max(p.thickness_um)]);
    ymax=max(T.CameraDiffuse(T.ScenarioID==iscenario))*1.15;
    if ymax<=0,ymax=1/p.N;end
    ylim([0 ymax]);
    % Faint shaded band, conventional histologic section thickness
    patch([5 10 10 5],[0 0 ymax ymax],[0.6 0.6 0.6],...
        'FaceAlpha',0.14,'EdgeColor','none','HandleVisibility','off');
    xlabel('Exposed FFPE tissue thickness (\mum)');
    ylabel('Camera-collected diffuse photon fraction');
    if iscenario==1,title('2016 reference scattering inputs');
    else,title('Higher-scattering assumption');end
    grid on;legend('Location','best');
end
sgtitle('HyPASS Objective 1: tissue thickness over embedding paraffin');
hideAxesToolbars(fig);
exportgraphics(fig,fullfile(OUT,'HyPASS_Objective1_Thickness.png'),'Resolution',300);
exportgraphics(fig,fullfile(OUT,'HyPASS_Objective1_Thickness.pdf'),'ContentType','vector');

%% FIGURE 2: detected photons that entered the underlying paraffin
fig2=figure('Color','w','Position',[80 80 900 520]);hold on;
for iscenario=1:2
 for iw=1:numel(p.lambda_nm)
    idx=T.ScenarioID==iscenario & T.Wavelength_nm==p.lambda_nm(iw);
    vv=T.CameraDiffuse_EnteredParaffin(idx);
    if iscenario==1,ls='-';else,ls='--';end
    semilogx(T.TissueThickness_um(idx),vv,'LineStyle',ls,...
        'Color',c(iw,:),'LineWidth',1.6,...
        'DisplayName',sprintf('%d nm, %s',p.lambda_nm(iw),scenarioName{iscenario}));
 end
end
set(gca,'XScale','log','YScale','linear');
xlabel('Exposed FFPE tissue thickness (\mum)');
ylabel('Collected fraction that entered underlying paraffin');
title('Contribution from photons visiting underlying paraffin');grid on;
legend('Location','eastoutside');
hideAxesToolbars(fig2);
exportgraphics(fig2,fullfile(OUT,'HyPASS_Objective1_Paraffin_Contribution.png'),'Resolution',300);
exportgraphics(fig2,fullfile(OUT,'HyPASS_Objective1_Paraffin_Contribution.pdf'),'ContentType','vector');

save(fullfile(OUT,'HyPASS_Objective1_Settings.mat'),'p','Optics');
fprintf('\nEnergy closure: min=%.6f, median=%.6f\n',min(T.EnergyClosure),median(T.EnergyClosure));
if any(T.EnergyClosure<0.999)
 warning('Some photon histories did not terminate before maxEvents. Check EnergyClosure.');
end
fprintf('Finished Objective 1. Outputs saved in %s\n',OUT);
% In MATLAB, use: winopen(OUT)  to open the saved results folder.

%% Local functions: analog layered MC transport
function result=mcLayered(N,theta,dT,dP,muaT,musT,gT,nT,muaP,musP,gP,nP,p)
% Two physical layers: exposed FFPE tissue [0,dT], paraffin [dT,dT+dP].
% Infinite lateral extent. Normal camera collects top-exiting photons
% within NA in air. Fresnel reflections and Snell refractions at both interfaces.
% Analog MC: photon absorbed OR scattered (not partial weight).
% Result = [Diffuse,Collected,CollectedAfterParaffin,Bottom,AbsT,AbsP,...
%           SpecularTop,MedianMaxDepth,90thMaxDepth,Closure]
tallies=zeros(1,7);
sampleDepth=zeros(N,1); collected=0;
naAngle=asin(min(p.NA,0.999));
muT=[muaT+musT, muaP+musP]; g=[gT,gP]; mua=[muaT,muaP];
n=[nT,nP];
Rs=fresnel(p.nAir,nT,cosd(abs(theta)));
for i=1:N
    if rand<Rs
        tallies(7)=tallies(7)+1; continue;
    end
    thetaT=asind(sind(theta)*p.nAir/nT);
    u=[sind(thetaT),0,cosd(thetaT)];
    z=0; layer=1; visitedP=false; maxDepth=0;
    done=false;
    for ev=1:p.maxEvents
        % Distance to next interaction, versus planar boundary.
        if u(3)>0
            if layer==1, zb=dT; else, zb=dT+dP; end
        elseif u(3)<0
            if layer==1, zb=0; else, zb=dT; end
        else
            zb=inf;
        end
        if isfinite(zb), dist=max(0,(zb-z)/u(3)); else, dist=inf; end
        s=-log(max(realmin,rand))/max(muT(layer),eps);
        if s >= dist
            % Encounter planar interface before volume interaction
            z=zb; maxDepth=max(maxDepth,z);
            fromLayer=layer;
            if z==0 && u(3)<0
                n1=nT; n2=p.nAir;
            elseif z==dT && u(3)>0
                n1=nT; n2=nP;
            elseif z==dT && u(3)<0
                n1=nP; n2=nT;
            else
                n1=nP; n2=p.nAir;
            end
            cos1=abs(u(3));
            R=fresnel(n1,n2,cos1);
            if rand<R
                u(3)=-u(3);
            else
                cos2=sqrt(max(0,1-(n1/n2)^2*(1-cos1*cos1)));
                u(1:2)=(n1/n2)*u(1:2);
                u(3)=sign(u(3))*cos2;
                u=u/norm(u);
                if z==0 && u(3)<0
                    tallies(1)=tallies(1)+1;
                    if acos(min(1,max(-1,-u(3))))<=naAngle
                        tallies(2)=tallies(2)+1;
                        if visitedP, tallies(3)=tallies(3)+1; end
                        collected=collected+1;sampleDepth(collected)=maxDepth;
                    end
                    done=true;
                elseif z==dT+dP && u(3)>0
                    tallies(4)=tallies(4)+1;done=true;
                elseif z==dT
                    if fromLayer==1,layer=2;visitedP=true;else,layer=1;end
                end
            end
            if done,break;end
            % prevent floating point sticking on interface
            z=max(0,z+1e-11*sign(u(3)));
        else
            z=z+s*u(3);maxDepth=max(maxDepth,z);
            if rand < mua(layer)/muT(layer)
                if layer==1,tallies(5)=tallies(5)+1;
                else,tallies(6)=tallies(6)+1;end
                done=true;break;
            end
            u=rotateHG(u,g(layer));
        end
    end
    if ~done
        % Hard capped photon histories counted as unclosed; inspect closure.
        % High closure deficit indicates maxEvents or material depth too small.
    end
end
if collected>0
    v=sort(sampleDepth(1:collected));
    medianDepth=v(max(1,ceil(.5*collected)));
    p90Depth=v(max(1,ceil(.9*collected)));
else,medianDepth=NaN;p90Depth=NaN;end
result=[tallies/N,medianDepth*1e4,p90Depth*1e4,sum(tallies([1 4 5 6 7]))/N];
end

function un=rotateHG(u,g)
r=rand;
if abs(g)<1e-12,ct=2*r-1;
else,ct=(1+g^2-((1-g^2)/(1-g+2*g*r))^2)/(2*g);end
ct=max(-1,min(1,ct));st=sqrt(max(0,1-ct^2));
phi=2*pi*rand;
if abs(u(3))>.99999
    un=[st*cos(phi),st*sin(phi),sign(u(3))*ct];
else
    h=sqrt(1-u(3)^2);
    un=[st*(u(1)*u(3)*cos(phi)-u(2)*sin(phi))/h+u(1)*ct,...
        st*(u(2)*u(3)*cos(phi)+u(1)*sin(phi))/h+u(2)*ct,...
        -st*h*cos(phi)+u(3)*ct];
end
un=un/norm(un);
end

function R=fresnel(n1,n2,c)
c=max(0,min(1,abs(c)));
st=(n1/n2)*sqrt(max(0,1-c^2));
if st>=1,R=1;return;end
ct=sqrt(max(0,1-st^2));
rs=(n1*c-n2*ct)/(n1*c+n2*ct);
rp=(n2*c-n1*ct)/(n2*c+n1*ct);
R=.5*(rs^2+rp^2);
end


function hideAxesToolbars(f)
% Remove MATLAB's interactive toolbar from exported figures.
ax=findall(f,'Type','axes');
for k=1:numel(ax)
    try
        ax(k).Toolbar.Visible='off';
    catch
        % Some older MATLAB releases may not expose axes.Toolbar.
    end
end
end
