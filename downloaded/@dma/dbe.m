%From paper: Wang et al., Automatically determining the number of clusters in unlabelled data sets, TKDE, 2009
%
%I      : a nxn VAT matrix (scaled or unscaled, but with smaller values denoting cluster formations)
%alpha  : a value within [0.01, 0.03] that corresponds to avg intracluster distances 
%c      : the number of detected clusters
%details: interim results

%JY Goulermas, 2010

function [c, details] = dbe(I, alpha, display)

  n = size(I,1);
  if ~exist('alpha', 'var') || ~isfinite(alpha)
    alpha = 0.03;
  end
  
  if ~exist('display', 'var')
    display = false;
  end

  mx = max( I(:) );
  if mx ~= 1.0, I = I ./ mx; end
 
  %contrast enhancement
  s  = graythresh(I);
  I1 = 1 - exp(-I./s);
  
  %binarisation
  s  = graythresh(I1);
  I2 = im2bw(I1, s); %s is often higher than it needs for small n
  
  %morphological cleaning
  l1 = alpha * n;
  se = strel('line', l1, 0);
  I3 = imopen(I2, se);
  
  %distance transformation & projection
  I4 = bwdist(I3, 'euclidean');
  H1 = ProjectMainDiag(I4); %also approximated as: sum(imrotate(I4, 45),1);
  l2 = max(3, 2 * alpha * n); %better increase this more for low a & n
  H2 = smooth(H1, l2, 'moving');
  H3 = diff(H2,1);
  
  %peak-valley detection
  H       = sign(H3);
  H(H>0)  = +1;
  H       = diff(H);
  peaks   = find(H<0) +1; %also could use: [values, peaks] = findpeaks(H2)
  valleys = find(H>0) +1;
  
  %weak peak removal
  valleys    = [1; valleys; numel(H3)]; %assert( numel(peaks) == numel(valleys)-1 )
  l3         = 2 * alpha * n;
  idx        = find( diff(valleys) < l3 );
  peaks(idx) = [];

  %cluster calculation & return
  c                  = numel(peaks);
  details.peaks      = peaks;
  details.binary     = I3;
  details.projection = H1;
   
  if display
    PlotResults();
  end
  
  %----------------------------------------------
  %PlotResults: to show all interim calculations in one figure
  function [] = PlotResults()
    figure
    colormap gray
    subplot(2,4,1), imagesc(I), title('original D'), axis square
    subplot(2,4,2), imagesc(I1), title('contrast enhanced'), axis square
    subplot(2,4,3), imhist(I1), title('histogram of enhanced')
    subplot(2,4,4), imagesc(I2), title('binarised'), axis square
    subplot(2,4,5), imagesc(I3), title('opened binary'), axis square
    subplot(2,4,6), imagesc(I4), title('distance transform'), axis square
    subplot(2,4,7), plot(H2, 'k'), hold on, plot(H1, 'r--'), title('raw & smoothed projection')
    subplot(2,4,8), plot(H3), hold on,
                    plot(peaks, H3(peaks), 'sr'),
                    plot(valleys, H3(valleys), 'ok'), title('peaks & valleys')
  end
  %----------------------------------------------
  
end


%----------------------------------------------
%ProjectMainDiag:
%A is alpha square matrix
%v is the sum of all antidiagonals
function [v] = ProjectMainDiag(A)
  n = size(A,1);
  A = fliplr(A);
  v = [];
  for k = (n-1) : -1 : -(n-1)
    v = [v; sum( diag(A, k) ) ];
  end
end
%----------------------------------------------