% plots and figures for the covariance estimation paper


classdef plots
    
    properties(Constant)
        figPath = '~/cov/figures/src/'
    end
    
    methods(Static)
        
        function fig3
            
            useQuad = false;
            
            if useQuad
                fname = 'comparison-quad';
                pairs = {
                    };
            else
                fname = 'comparison';
                pairs = [
                    10  0  .09  .05
                    30  0  .09  .05
                    40  0  .09  .05
                    70  0  .09  .05
                    
                    30 10  .015  .01
                    40 10  .015  .01
                    70 10  .015  .01
                    
                    40 30  .015  .01
                    70 30  .015  .01
                    
                    70 40  .015  .01
                    ];
            end
            nbins = 24;
            
            c = covest.CovMatrix & 'nfolds>1';
            c1 = pro(c, 'method->m1','cv_loss->l1');
            c2 = pro(c, 'method->m2','cv_loss->l2');
            
            for i=1:size(pairs,1)
                r = struct('m1',pairs(i,1),'m2',pairs(i,2));
                %[l1,l2] = fetchn(c1*c2 & r, 'l1', 'l2');
                [l1,l2] = fetchn(covest.ActiveCells, c1*c2 & r, 'avg(l1)->ll1','avg(l2)->ll2');
                
                fig = Figure(1, 'size', [40 30]);
                rr = pairs(i,3);
                bins = linspace(-rr,rr,nbins);
                x = l2-l1;
                %                xarrow = fetch1(r & arrowKey, 'l1-l2->d');
                a = hist(min(1,x),bins);
                bar(bins,a)
                p = signrank(x);
                fprintf('%2d-%2d Significance %1.0e, median difference = %g\n', ...
                                     pairs(i,2),pairs(i,1),p, median(x))
                xlim([-1 1]*rr*(1+1.1/nbins))
                ylim([-.5 14])
                hold on
                plot([0 0], ylim,'k')
                plot([0 0]+median(x), ylim, ':','Color',[.6 0 0])
                hold off
                box off
                set(gca, 'YColor',[1 1 1]*0.99,'YTick',[])
                ticks = (-1:1)*pairs(i,4);
                colormap gray
                
                set(gca,'XTick',ticks,'XTickLabel',nozero(ticks))
                set(gca,'YTick',10)
                set(gca,'Position', [.05    .2    .8    .8])
                hold on
                PlotAxisAtOrigin([1 1])
                hold off
                
                fig.cleanup
                fig.save(fullfile(covest.plots.figPath,sprintf('%s-%02d-vs-%02d.eps',fname,pairs(i,1),pairs(i,2))))
            end
            
            function s=nozero(f)
                % remove leading zeros in decimal fractions
                if isscalar(f)
                    s = sprintf('%g',f);
                    if strncmp(s,'-0.',3)
                        s(2)='';
                    elseif strncmp(s,'0.',2)
                        s(1)='';
                    end
                else
                    s = arrayfun(@(f) nozero(f), f, 'uni', false);
                end
            end
        end
        
        
        
        function compareCorrMethods
            % compare correlation coefficients with variance estimated in
            % each bin or globally
            
            for key = fetch(covest.ActiveCells & 'preprocess_method_num=5' & 'high_repeats')'
                clf
                [X,evokedBins,sel,xyz] = fetch1(covest.ActiveCells*covest.Traces & key, ...
                    'trace_segments','evoked_bins', 'selection', 'cell_xyz');
                evokedBins = min(evokedBins,size(X,1));
                
                % compare with original binned traces
                %X = X(1:min(end,evokedBins),:,:,sel);
                X = X(:,:,:,sel);
                
                [nBins,nConds,nTrials,nCells] = size(X);
                
                % subtract mean response
                if size(X,1)<=evokedBins
                    M = nanmean(X,3);
                else
                    M1 = nanmean(X(1:evokedBins,:,:,:),3);   % binwise mean of evoked response
                    M2 = reshape(nanmean(reshape(X(evokedBins+1:end,:,:,:),[],nCells)),1,1,1,nCells);  % common mean in intertrial periods
                    M2 = repmat(M2,size(X,1)-evokedBins,nConds);
                    M = cat(1,M1,M2);
                end
                
                X = bsxfun(@minus, X, M);
                
                % method 1: common variance
                V = reshape(nanvar(reshape(X,[],nCells)), 1,1,1,nCells);
                c1 = getCorr(X,V);
                
                %                 % method 2: condition-specific variance
                %                 V = reshape(nanvar(reshape(permute(X, [1 3 2 4]),[],nConds,nCells)), 1,nConds,1,nCells);
                %                 c2 = getCorr(X,V);
                
                % method 3: bin-specific variance
                if size(X,1)<=evokedBins
                    V = nanvar(X,[],3);
                else
                    V1 = nanvar(X(1:evokedBins,:,:,:),[],3);   % binwise mean of evoked response
                    V2 = reshape(nanvar(reshape(X(evokedBins+1:end,:,:,:),[],nCells)),1,1,1,nCells);  % common mean in intertrial periods
                    V2 = repmat(V2,size(X,1)-evokedBins,nConds);
                    V = cat(1,V1,V2);
                end
                c3 = getCorr(X,V);
                
                
                cc1 = c1;
                cc2 = c3;
                clear c0 c1 c2 c3
                
                
                subplot 431, imagesc(cc1,[-1 1]*.2), axis image off, title 'common variance'
                colorbar
                subplot 432, imagesc(cc2,[-1 1]*.2), axis image off, title 'binned variance'
                colorbar
                colormap(covest.lib.doppler)
                
                n = size(cc1,1);
                [i,j] = meshgrid(1:n,1:n);
                cc1 = cc1(i<j);
                cc2 = cc2(i<j);
                mm1 = mean(cc1);
                mm2 = mean(cc2);
                
                r = [-.1 .3];
                subplot 434, hist(cc1,linspace(r(1),r(2),100)), xlim(r), ylim_ =ylim; xlabel correlations, hold on, plot(mm1*[1 1],ylim_,'r'), hold off, grid on
                subplot 435, hist(cc2,linspace(r(1),r(2),100)), xlim(r), ylim(ylim_), xlabel correlations, hold on, plot(mm2*[1 1],ylim_,'r'), hold off, grid on
                
                ix = linspace(r(1)/2,r(2)/2,50);
                
                subplot (4,3,[3 6]), densityPlot(ix, ix, cc1, cc2), xlabel 'common variance', ylabel 'binned variance'
                hold on, plot(mm1,mm2,'r+','MarkerSize',40), set(refline(1),'Color','r'), hold off
                
                if count(covest.OriTuning & setfield(key,'high_repeats',0)) %#ok<SFLD>
                    [ori,p] = fetch1(covest.OriTuning & setfield(key,'high_repeats',0), 'von_pref','von_p_value'); %#ok<SFLD>
                    ori(p<0.05) = nan;
                    ori = ori(sel)*180/pi;
                    [ori1,ori2] = meshgrid(ori,ori);
                    oDiff = oriDiff(ori1(i<j),ori2(i<j));
                    lowerLimits = [0 15 45];
                    
                    bins = sum(bsxfun(@ge, oDiff(~isnan(oDiff)), lowerLimits),2);
                    meanCorrs1 = accumarray(bins, cc1(~isnan(oDiff)), size(lowerLimits'), @mean);
                    meanCorrs2 = accumarray(bins, cc2(~isnan(oDiff)), size(lowerLimits'), @mean);
                    
                    subplot 437, bar(0.5:length(lowerLimits), meanCorrs1), ylim_=ylim*1.2;  ylim([0 ylim_(2)])
                    set(gca,'XTick',0:length(lowerLimits),'XTickLabel',[lowerLimits 90]), grid on, xlim([0 length(lowerLimits)])
                    subplot 438, bar(0.5:length(lowerLimits), meanCorrs2), ylim([0 ylim_(2)])
                    set(gca,'XTick',0:length(lowerLimits),'XTickLabel',[lowerLimits 90]), grid on, xlim([0 length(lowerLimits)])
                    xlabel '\DeltaOri (degrees)'
                end
                
                dist = arrayfun(@(i,j) norm(xyz(i,:)-xyz(j,:)), i,j);
                dist = dist(i<j);
                lowerLimits = [0 25 75 150];
                bins = sum(bsxfun(@ge, dist, lowerLimits),2);
                meanCorrs1 = accumarray(bins, cc1, size(lowerLimits'), @mean);
                meanCorrs2 = accumarray(bins, cc2, size(lowerLimits'), @mean);
                
                subplot(4,3,10), bar(0.5:length(lowerLimits), meanCorrs1), ylim_=ylim*1.2; ylim([0 ylim_(2)])
                set(gca,'XTick',0:length(lowerLimits)-1,'XTickLabel',lowerLimits), grid on, xlim([0 length(lowerLimits)])
                subplot(4,3,11), bar(0.5:length(lowerLimits), meanCorrs2), ylim([0 ylim_(2)])
                set(gca,'XTick',0:length(lowerLimits)-1,'XTickLabel',lowerLimits), grid on, xlim([0 length(lowerLimits)])
                xlabel 'distance (\mum)'
                
                
                set(gcf,'PaperSize',[12 12],'PaperPosition',[0 0 12 12])
                print('-dpng','-r400', sprintf('~/dev/temp/corrComp5-%04d',mod(key.aod_scan_start_time,1e4)))
            end
            
            
            function c = getCorr(X,V)
                z = reshape(bsxfun(@rdivide,X,sqrt(V+eps)),[],size(X,4));
                c = corrcov(nancov(z));
            end
            
        end
        
    end
end


function densityPlot(x,y,c1,c2)
img = hist3([c2 c1],{x,y});
img = 255*img/max(img(:));
cmap = [1 1 1; jet(255)];
%(1-bsxfun(@power, gray(256), [.2 .3 .4])
subimage(x, y, img, cmap)
set(gca,'YDir','normal')
grid on
end


function d=oriDiff(ori1,ori2)
% compute the absolute difference between orientations ori1 and ori2 (in degrees)
ori1 = mod(ori1, 180);
ori2 = mod(ori2, 180);
b1 = min(ori1,ori2);
b2 = max(ori1,ori2);
d = min(b2-b1,b1+180-b2);
d(isnan(ori1))=nan;
d(isnan(ori2))=nan;
end