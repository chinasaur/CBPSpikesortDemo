function DisplaySortedSpikes(data_pp, spike_times, spike_amps, init_waveforms, ...
                             snippets, recon_snippets, ...
                             params, name)

nchan=size(data_pp.data,1);

%** Add residuals

figure(params.plotting.first_fig_num); clf
set(gcf, 'MenuBar', 'none', 'Name', name);
subplot(2,1,1);  
inds = params.plotting.dataPlotInds;
plotChannelOffset = 6*ones(length(inds),1)*([1:nchan]-1); %**magic number
plot((inds-1)*data_pp.dt, data_pp.data(:,inds)' + plotChannelOffset, 'k');
axis tight
yrg = get(gca,'Ylim');   xrg = get(gca,'Xlim');
title(sprintf('Data, filtered & whitened, nChannels=%d, %.1fkHz', nchan, 1/(1000*data_pp.dt)));

subplot(2,1,2);
bandHt = 0.12;
yinc = bandHt*(yrg(2)-yrg(1))/length(init_waveforms);
clrs = hsv(length(init_waveforms));
patch([xrg'; xrg(2); xrg(1)], [yrg(2)*[1;1]; (yrg(1)+(1+bandHt)*(yrg(2)-yrg(1)))*[1;1]], ...
      0.9*[1 1 1], 'EdgeColor', 0.9*[1 1 1]);
set(gca,'Ylim', [yrg(1), yrg(2)+bandHt*(yrg(2)-yrg(1))]);
hold on
%** Should do proper interpolation
midChan = ceil(nchan/2);
for n=1:length(init_waveforms)
    spkInds = (spike_times{n} > inds(1)) & (spike_times{n} < inds(end));
    tInds = spike_times{n}(spkInds);
    plot((tInds-1)*data_pp.dt, (yrg(2)+(n-0.5)*yinc)*ones(1,length(tInds)), '.', 'Color', clrs(n,:));
    trace = zeros(length(inds),nchan);   
    trace(round(tInds)-inds(1)+1,midChan) = spike_amps{n}(spkInds)';
    trace = conv2(trace, reshape(init_waveforms{n},[],nchan), 'same');
    plot((inds-1)*data_pp.dt, trace + plotChannelOffset, 'Color', clrs(n,:));
    plot((inds-1)*data_pp.dt, plotChannelOffset, 'k');
end
hold off
xlabel('time (sec)');
title('Recovered spikes');

% Residual Histograms
existingFig = ishghandle(params.plotting.first_fig_num+2);
figure(params.plotting.first_fig_num+1); clf
set(gcf, 'Name', 'Residual histograms');
if (~existingFig) % only do this if we've created a new figure (avoid clobbering user changes)
    set(gcf, 'ToolBar', 'none');
end
resid = cell2mat(cellfun(@(c,cr) c-cr, snippets, recon_snippets, 'UniformOutput', false));
subplot(2,1,1); 
%mx = max(cellfun(@(c) max(abs(c(:))), snippets));
mx = max(abs(data_pp.data(:)));
[N, Xax] = hist(resid, mx*[-50:50]/101);
plot(Xax,N); set(gca,'Yscale','log'); rg=get(gca,'Ylim');
hold on
gh=plot(Xax, max(N(:))*exp(-(Xax.^2)/2), 'r', 'LineWidth', 2); 
plot(Xax,N); set(gca,'Ylim',rg); set(gca, 'Xlim', [-mx mx]); 
hold off; 
if (nchan < 1.5)
    title('Histogram, filtered/whitened data with spikes removed');
else
    title(sprintf('Histograms, filtered/whitened data with spikes removed (%d channels)', nchan));
end
legend(gh, 'univariate Gaussian');
subplot(2,1,2); 
mx = max(sqrt(sum(data_pp.data.^2,1)));
[N,Xax] = hist(sqrt(sum(resid.^2, 2)), mx*[0:100]/100);
chi = 2*Xax.*chi2pdf(Xax.^2, nchan);
bar(Xax,N); set(gca,'Yscale','log'); yrg= get(gca, 'Ylim'); 
hold on;
ch= plot(Xax, (max(N)/max(chi))*chi, 'r', 'LineWidth', 2);
hold off; set(gca, 'Ylim', yrg); set(gca, 'Xlim', [0 mx]);
title('Histogram, magnitude over filtered/whitened channel(s), with spikes removed');
legend(ch, 'chi-distribution, univariate Gaussian');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0)
  % PC scatter plot.  ***not clear this is useful...
  figure(params.plotting.first_fig_num+5); clf
  wfs = cell2mat(cellfun(@(c) c(:), init_waveforms, 'UniformOutput', false));
  thresh= params.clustering.spike_threshold;
    
%  if (size(wfs,2)<3), wfs = [wfs, PCs(:, 3-size(wfs,2))]; end
%  [ax,s] = svd(wfs);  ax = ax(:,[1,2]);
  ax = PCs(:,[1,2]);

  proj_wfs = ax'*wfs;

  % YUK
  cluster_pars = params.clustering;
  if isempty(cluster_pars.window_len), cluster_pars.window_len = params.general.waveform_len; end
  cluster_pars.align_mode = data_pp.polarity;

  proj_snippets = []; snippet_ids = []; snippet_dist2wf = [];
  for n=1:size(wfs,2)
    %***Shouldn't have to round the times!
    snippets = ax'*ConstructSnippetMatrix(data_pp.data, round(spike_times{n}), cluster_pars);
    proj_snippets = [proj_snippets, snippets];
    distances = sqrt(sum((snippets-repmat(proj_wfs(:,n),1,size(snippets,2))).^2))';
    snippet_dist2wf = [snippet_dist2wf; distances];
    snippet_ids = [snippet_ids; n*ones(length(distances), 1)];
  end
  
  %** need to fix display of snippets with two or more spikes...
  hold on
  for n=1:size(wfs,2) %plot central cluster first
    sn = proj_snippets(:, ((snippet_ids==n)&(snippet_dist2wf<thresh)));
    plot(sn(1,:), sn(2,:), '.', 'Color', 0.5*clrs(n,:)+0.5*[1 1 1]);
  end
  for n=1:size(wfs,2) %then plot outliers
    sn = proj_snippets(:, ((snippet_ids==n)&(snippet_dist2wf>=thresh)));
    plot(sn(1,:), sn(2,:), '.', 'Color', 0.5*clrs(n,:)+0.5*[1 1 1]);
    plot(proj_wfs(1,n), proj_wfs(2,n), 'o', 'MarkerSize', 9, 'LineWidth', 2,...
	 'MarkerEdgeColor', 'black', 'MarkerFaceColor', clrs(n,:));
  end
  xl = get(gca, 'XLim'); yl = get(gca, 'YLim');
  plot([0 0], yl, '-', 'Color', 0.8 .* [1 1 1]);
  plot(xl, [0 0], '-', 'Color', 0.8 .* [1 1 1]);
  th=linspace(0, 2*pi, 64);
  nh= plot(thresh*sin(th),thresh*cos(th), 'k', 'LineWidth', 2);
  legend(nh,sprintf('spike threshold = %.1f',thresh));
  axis equal
  hold off
  xlabel('PC 1'); ylabel('PC 2');  title('CBP results');
end
%   Fig6: projection into PC space of segments, with spike assignments (as in paper)
