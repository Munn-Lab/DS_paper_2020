function [freq] = partialisation(cfg, freq, field)

crsspctrm = getfield(freq, field);
Nrpt      = size(crsspctrm,1);
Nsgn      = size(crsspctrm,2);
Nfrq      = size(crsspctrm,3);
Ntim      = size(crsspctrm,4);
%create signalcombination-matrix, each non-nan entry corresponding
%with the index of the signalcombination in the cross-spectrum
%negative entries should be conjugated before further processing.
%the underlying idea is to temporarily convert the cross-spectra into 
%a matrix-shape, so that the partialisation can take place efficiently
%this mainly involves a lot of careful book-keeping
fprintf('indexing the channels and channel combinations\n');
sgnindx    = [1:length(freq.label)]';
cmbindx = zeros(size(freq.labelcmb));
for i=1:size(freq.labelcmb,1)
  cmbindx(i,1) = strmatch(freq.labelcmb{i,1}, freq.label, 'exact');
  cmbindx(i,2) = strmatch(freq.labelcmb{i,2}, freq.label, 'exact');
end
sgncmbmat  = nan*zeros(length(sgnindx),length(sgnindx));
for i=1:size(cmbindx,1)
  sgncmbmat(cmbindx(i,2),cmbindx(i,1)) = -i;
  sgncmbmat(cmbindx(i,1),cmbindx(i,2)) = i;
end

fprintf('checking whether partialisation of the requested channel(s) is possible\n');
prtindx = [];
for j = 1:length(cfg.partchan)
  prtindx = [prtindx; strmatch(cfg.partchan{j}, freq.label, 'exact')];
end
if isempty(prtindx), warning('no partialisation will be performed, requested channel(s) not present in data'); end;

prtindx     = sort(prtindx);
sgnindx     = setdiff(sgnindx, prtindx);
prtsgnmat   = sgncmbmat(prtindx,prtindx);
rstsgnmat   = sgncmbmat(sgnindx,sgnindx);
Nprtsgn     = length(prtindx);
Nrstsgn     = length(sgnindx);  
prtsgnlabel = freq.label(prtindx); 
label       = freq.label(sgnindx);

%the reformatting of the cross-spectral densities leads to a re-ordering of the signal-combinations
[ind1, ind2] = find(~isnan(rstsgnmat));
labelcmb     = [label(ind1) label(ind2)];
cmbindx      = [ind1 ind2];
for j = 1:size(cmbindx,1)
  duplicate(j,1) = find([cmbindx(:,1)==cmbindx(j,2)] .* [cmbindx(:,2)==cmbindx(j,1)]);
end %but we end up with a redundant amount of cross-spectra
powindx = find(duplicate==duplicate(duplicate)); %indices of auto-spectra
sel     = setdiff(1:size(cmbindx,1),powindx);
duplicate(powindx) = nan;
for j = 1:length(duplicate)
  if ~isnan(duplicate(j)), duplicate(find(duplicate==j)) = nan; end
end

%these cross-spectral density matrices should be complete
if any(isnan(prtsgnmat)),    error('partialisation of the requested channel(s) is not possible'); end;
prtrstsgnmat = sgncmbmat(prtindx,sgnindx);
if any(isnan(prtrstsgnmat)), error('partialisation of the requested channel(s) is not possible'); end;
rstprtsgnmat = sgncmbmat(sgnindx,prtindx);
if any(isnan(rstprtsgnmat)), error('partialisation of the requested channel(s) is not possible'); end;

%allocate memory
csdrr = nan+complex(zeros(Nrstsgn.^2,Nfrq,Ntim)      , zeros(Nrstsgn.^2,Nfrq,Ntim));
csdpp = nan+complex(zeros(Nprtsgn.^2,Nfrq,Ntim)      , zeros(Nprtsgn.^2,Nfrq,Ntim));
csdrp = nan+complex(zeros(Nrstsgn.*Nprtsgn,Nfrq,Ntim), zeros(Nrstsgn.*Nprtsgn,Nfrq,Ntim));
csdpr = nan+complex(zeros(Nprtsgn.*Nrstsgn,Nfrq,Ntim), zeros(Nprtsgn.*Nrstsgn,Nfrq,Ntim));

progress('init', cfg.feedback, 'partialising out the requested channel(s)')
for m = 1:size(crsspctrm,1)
%  progress(m/size(crsspctrm,1));
  %a bit of hocus-pocus with the cross-spectral densities
  csdrr(find(~isnan(rstsgnmat)),:,:)  = squeeze(crsspctrm(m,abs(rstsgnmat(find(~isnan(rstsgnmat)))),:,:));
  csdpp(:) = squeeze(crsspctrm(m,abs(prtsgnmat),:,:));  
  csdrp(:) = squeeze(crsspctrm(m,abs(rstprtsgnmat),:,:));
  csdpr(:) = squeeze(crsspctrm(m,abs(prtrstsgnmat),:,:));
  
  %take care of the conjugates
  csdrr(find(rstsgnmat<0),:,:)     = conj(csdrr(find(rstsgnmat<0),:,:));
  csdpp(find(prtsgnmat<0),:,:)     = conj(csdpp(find(prtsgnmat<0),:,:));
  csdrp(find(rstprtsgnmat<0),:,:)  = conj(csdrp(find(rstprtsgnmat<0),:,:));
  csdpr(find(prtrstsgnmat<0),:,:)  = conj(csdpr(find(prtrstsgnmat<0),:,:));
  
  %perform the partialisation
  for j = 1:Nfrq
    for k = 1:Ntim
      rr  = reshape(squeeze(csdrr(:,j,k)),[Nrstsgn Nrstsgn]);
      rp  = reshape(squeeze(csdrp(:,j,k)),[Nrstsgn Nprtsgn]);
      pp  = reshape(squeeze(csdpp(:,j,k)),[Nprtsgn Nprtsgn]);
      pr  = reshape(squeeze(csdpr(:,j,k)),[Nprtsgn Nrstsgn]);
      rrp = rr - rp*pinv(pp)*pr;    
      csdp(m,:,j,k) = rrp(find(~isnan(rrp)));
    end
  end
end
progress('close');

freq.label     = label;
freq.labelcmb  = labelcmb([powindx; find(~isnan(duplicate))], :);
freq           = setfield(freq, field, csdp(:, [powindx; find(~isnan(duplicate))], :, :));



