clc
rng(0)

%a small example
D = [ 0, .73, .19, .71, .16
  .73 0 .59 .12 .78
  .19 .59 0 .55 .19
  .71 .12 .55 0 .74
  .16 .78 .19 .74 0];

D_ready=false; %dont chnage this

switch(1) 

case 1
c = [0.2,0.35;
     0.1, 0.7;
     0.7, 0.4];
n = [30, 50, 80]*   1;
s = [0.04, 0.04, 0.02]*1.9;
X = [];
for i = 1 : numel(n)
Y = repmat( c(i,:), n(i), 1);
N = randn( n(i), size(c,2) )*s(i) ;
%N = randn( n(i), size(c,2) )*s(i) - s(i)/2;
X = [X; Y+N];
end
  
case -1
c = [0.2,0.35;
     0.1, 0.7;];
n = [30, 50];
s = [0.025, 0.02]*2.9;
X = [];
for i = 1 : numel(n)
Y = repmat( c(i,:), n(i), 1);
N = randn( n(i), size(c,2) )*s(i) ;
%N = randn( n(i), size(c,2) )*s(i) - s(i)/2;
X = [X; Y+N];
end

case 2
c = [0,0;  3,4;  6,0];
s = ones(1,3) * 0.83162;
n = round( [0.15,0.35,0.5] * 50 );
X = [];
for i = 1 : numel(n)
Y = repmat( c(i,:), n(i), 1);
N = randn( n(i), size(c,2) )*s(i) ;
%N = randn( n(i), size(c,2) )*s(i) - s(i)/2;
X = [X; Y+N];
end

case 44
c = [0,0;  3,4;  6,0; 11,11];
s = ones(1,4) * 1.23162;
n = round( [0.15,0.35,0.5,.5]* 50 );
X = [];
for i = 1 : numel(n)
Y = repmat( c(i,:), n(i), 1);
N = randn( n(i), size(c,2) )*s(i) ;
%N = randn( n(i), size(c,2) )*s(i) - s(i)/2;
X = [X; Y+N];
end

case 3
k=10;
[x,y] = meshgrid( linspace(0,1,k), linspace(0,1.8,k) );
X = [x(:), y(:)];

case -3
k=10;
[x,y] = meshgrid(1:1:k, 1:1:k);
X = [x(:), y(:)];

case 4
a1=linspace(0,2*pi, 50)';
a2=linspace(0,2*pi, 70)';
X1 = [ cos(a1), sin(a1) ] * 1;
X2 = [ cos(a2), sin(a2) ] * 2;
X = [X1;X2];

case -4
a=linspace(0,2*pi,20)';
X = [ cos(a), sin(a) ]*1e88;

case 5
a = linspace(0,2*pi * 5,150)';
r = 1:1:numel(a); 
r=r.^1.5;
X = [ cos(a), sin(a) ].*repmat(r',1,2) ;

case 22
C = rand(7,2) * 200;
%C = [1,1; 2,2; 1,2; 2,1; 1.5, 1.5] * 99;
n = 30;
X = [];
for i = 1 : size(C,1)
  X = [ X; mvnrnd( C(i,:), [1,0;0,1] *5, n ) ];
end
   
    
case 111
  D_ready = true;
%  X = load('iris.mat');
 % X = X.Iris;
%  X = rand(255,5)*1000-500;
  %X = load('/Volumes/yannisdisk/D.mat'); D=X.D;X=D;  
  %X = load('/Volumes/yannisdisk/workspace/work/papers/seriation datasets/alex/mds2.mat'); X = X.MDS2; D=X;
  
  %X = load('/Users/yannis/Desktop/FispinCorr.mat');
  %X = X.C;
  %D = 1-X;
  
  %D = (D+D')/2;

  load('/Volumes/yannisdisk/workspace/work/papers/_submitted ones/ppc/_resources/data textual from Tingting/Reuters_5topic.mat')
  D=R;  X=D;

case 999
  X = load('/Volumes/yannisdisk/workspace/work/papers/datasets/toy/spiral.txt');
  X = X(1:1:end,1:2);
  
  %X = rm(X, 'mode', 'minkowski', 'mink_p', 1).run().R; %assume that X is not known; just some distance info of it...
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%messy demo dma class code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
colormap jet

subplot(231)
plot( X(:,1), X(:,2), 'k.', 'markersize', 8 )
axis equal
title('X data points')

subplot(232)
if ~D_ready
[n,m] = size(X);
rel   = rm();
%D     = rel.set( X, 'mode', 'minkowski', 'mink_p', 2.3).run().R;
D     = rel.set( X, 'mode', 'euclidean').run().R;
   %D = rel.set(X, 'mode', 'knn', 'knn_type', 'undirected', 'knn_k', 7, 'knn_weighting', 'metric').run().R;
   %D(isinf(D)) = 100*max(max(D(~isinf(D))));
%D    = 1-rel.set(X, 'mode', 'correlation').run().R; %used,for instance, in dma.corder()
%D = D./sum(sum(D)); %make some objectives easier to see
D=real(D);
end
imagesc(D)
axis square
title('S original')

subplot(233)
if ~D_ready
rp = randperm(n);
else
  n = size(D,1);
  rp = 1:n;
end
D = D(rp,rp);
imagesc(D)
title('S shuffled')
axis square

subplot(234)
plot(X(:,1), X(:,2), 'k-.')
title('linked scatterplot')
axis equal

drawnow
time = cputime;


o2 = dma(D,'r2e'); ext_perm = o2.run().perm;
%o2= dma(D,'hc'); o2.hc_args= {'method', 'average', 'ordering', 'olo'}; ext_perm = o2.run().perm;
%ext_perm = 1:n; %[1:2:n, 2:2:n]; %silly order for testing
%
o2=dma(D, 'vat'); hybs{1} = o2.run().perm'; %o2=dma(D, 'mds'); hybs{2} = o2.run().perm';
%
o               = dma(D);
o.svd_args      = {'subset', [1,2], 'gap', false};
o.fpca_args     = {'dims', 2};
o.svat_args     = {'clusters', 15, 'sample' 0.6};
o.revat_sample  = 0.6;
o.spin_args     = {'method', 'sts', 'exact', true};
o.mds_args      = {'method', 'classical', 'dims', 1};
o.hc_args       = {'method', 'average', 'ordering', 'olo', 'ext_perm', ext_perm };
o.spectral_args = {'method', 'ding', 'ext_mix', 1.0, 'ext_perm', ext_perm };
o.tsp_args      = {'method', 'heuristic', 'starts', inf, 'precision', 2};
%o.ga_args       = {'objective', 'are', 'W', dma.template(n,'linear') };
o.ga_args       = {'objective', 'qap', 'W', dma.template(n,'banded', 11), 'pop_size', 300, 'tol_fun', 1e-7, 'hybrids', 'off' };
%o.ga_args       = {'objective', 'arc', 'pop_size', 100};
o.qap_args      = {'method', 'fw', 'max_iters', 555, 'W', dma.template(n,'linear') };
o.mode          = 1;
o.algo          = 'gncr'; %<<<<*****<<<<<<<<<<<<<********
o               = o.run();

fprintf('\ntime to run seriation: %2.2fsecs \n', cputime-time )

subplot(2,3,5)
imagesc(o.Q)
title('S ordered')
axis square

if ~D_ready
  subplot(2,3,6)
  X = X(o.perm,:);
  plot(X(:,1), X(:,2), 'k-*')
  %for q = 1 : size(X,1), hold on, plot( X(q,1), X(q,2), 'r*') , pause,  end
  hold on, plot(X(1,1), X(1,2), 'r.', 'markersize', 8), plot(X(end,1), X(end,2), 'g.', 'markersize', 8), 
  title('linked scatterplot ordered')
  axis equal
end



%keyboard

if strcmp( o.algo, 'revat') 
  o.plot_revat_profiles();
  pause, close(2:2)
end
%[c, details] = dma.dbe(o.Q, nan, true), pause, close


if 000 %stepwise paused graph
  figure(2)
  plot(X(:,1), X(:,2), 'm-')
  axis equal
  hold on
  for r = 1 : size(X,1)
    plot(X(r,1),X(r,2), 'ko')
    pause
  end
  close(2)
end

%o.error_window=25;
o = o.analyse();
fprintf('Original D:\n')
disp( o.scores.D )
fprintf('\nReordered Q:\n')
disp( o.scores.Q )
fprintf('\n\n\n')

if (false)
  struct2mat = @(x) cell2mat(struct2cell(x));
  diff = struct2mat(o.scores.Q) - struct2mat(dma.measures_naive(o.Q));
  max(abs(diff))
end










