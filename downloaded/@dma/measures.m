%Calculates various measures to evaluate the quality of seriated matrices.
%D can be a (square) symmetric distance matrix (one-mode), or a general perhaps non-square (two-mode) matrix.
%Anti-Robinson behaviour is measured, so that smaller values are near the main diagonal for one-mode D.
%Note1 : some measures only use the upper triangle of the (square) D.
%Note2 : some measures not meaningful to both one- and two-mode cases.
%m     : a structure with different measures (see initilisation section below).

%A Kostopoulos & JY Goulermas, 2012

function [m] = measures(D)

  [n,d] = size(D);

  %Square matrices:
  m.hp    = 0;            %hamiltonian path                             (loss)  %Caraux and Pinloche 2005 & Hahsler 2008 
  m.in    = 0;            %inertia criterion                            (merit) %Caraux and Pinloche 2005 & Hahsler 2008 
  m.ls    = 0;            %least squares criterion                      (loss)  %Caraux and Pinloche 2005 & Hahsler 2008 
  m.arc   = 0;            %anti-robinson compatibility                  (merit) %Hubert et al 2001        & Hahsler 2008
  m.warc  = 0;            %weighted anti-robinson compatibility         (merit) %Hubert et al 2001        & Hahsler 2008
  m.are   = 0;            %anti-robinson events violation               (loss)  %Chen 2002                & Hahsler 2008 & Streng 1991
  m.ware  = 0;            %weighted anti-robinson events violation      (loss)  %Chen 2002                & Hahsler 2008 & Streng 1991
  m.dware = 0;            %doubly weighted anti-robinson violation      (loss)  %Chen 2002                & Hahsler 2008 & Streng 1991
  m.g     = zeros(n-2,1); %generalised anti-robinson violation          (loss)  %Wu, Tien & Chen 2010
  m.rg    = zeros(n-2,1); %relative generalised anti-robinson violation (loss)  %Wu, Tien & Chen 2010


  %Rectangular matrices:
  m.moe  = 0; %measure of effectiveness          (merit) %McCormick et al 1972 & Hahsler 2008
  m.mstr = 0; %stress with moore neighbourhood   (loss)  %Niermann 2005        & Hahsler 2008
  m.nstr = 0; %stress with neumann neighbourhood (loss)  %Niermann 2005        & Hahsler 2008
  m.bcc  = 0; %bertin classification criterion   (loss)  %Pilhoefer 2012


  %m.moe:
  shifts = [zeros(n,1), D(:,1:end-1)] + [zeros(1,d); D(1:end-1,:)]; %equivalent to D*R + U'*D
  m.moe  = sum(sum(D .* shifts));
  %also:
  %R     = diag(ones(d-1,1), 1); %D*R shifts D right
  %U     = diag(ones(n-1,1), 1); %U*D shifts D up
  %H     = R' + R;
  %V     = U + U';
  %m.moe = sum(sum( D.* (D*H + V*D) )) /2; %or
  %m.moe = trace(D'*D*H) /2 + trace(D*D'*V) /2;


  %m.m/nstr:
  horz   = sum(sum( ( D(1:n  ,1:d-1) - D(1:n  ,2:d) ).^2 ));
  vert   = sum(sum( ( D(1:n-1,1:d)   - D(2:n  ,1:d) ).^2 ));
  diag1  = sum(sum( ( D(1:n-1,1:d-1) - D(2:n  ,2:d) ).^2 ));
  diag2  = sum(sum( ( D(2:n  ,1:d-1) - D(1:n-1,2:d) ).^2 ));
  m.mstr = 2 * (horz + vert + diag1 + diag2);
  m.nstr = 2 * (horz + vert);


  %m.bcc:
  cs    = cumsum( cumsum(D(1:n-1,d:-1:2), 1), 2 );
  m.bcc = sum(sum( D(2:n,1:d-1) .* cs(1:end,end:-1:1) ));



  if n ~= d %for two-mode matrices, ignore all measures below (no 2-mode square then:-)
    return;
  end


  %m.hp:
  m.hp = sum( diag(D,1) );


  %m.in, m.ls:
  [c,r] = meshgrid(1:n);
  W     = abs((r-c));
  m.in  = sum(sum( D .* W.^2 )); %same as: trace(D * W.^2);
  m.ls  = sum(sum( (D-W).^2 ));  %same as: k=sum(sum(D.^2 + W.^2));  %m.ls = k-2*trace(D*W);


  %m.robinson measures:
  D = triu(D,1);
  for i = 2 : n-1
    diff_hor = triu( D(1:n-i,2:n+1-i) - D(1:n-i,i+1:n) );
    diff_ver = triu( D(i:n-1,i+1:n)   - D(1:n-i,i+1:n) );

    m.arc  = m.arc  - sum(sum( sign(diff_hor) + sign(diff_ver) ));
    m.warc = m.warc - sum(sum( diff_hor       + diff_ver ));

    sign_hor = diff_hor > 0;
    sign_ver = diff_ver > 0;

    m.are   = m.are   + sum(sum( sign_hor                 + sign_ver ));
    m.ware  = m.ware  + sum(sum( sign_hor.*diff_hor       + sign_ver.*diff_ver ));
    m.dware = m.dware + sum(sum( sign_hor.*diff_hor*(i-1) + sign_ver.*diff_ver*(i-1) ));

    %tr  = tril(ones(n-i),window-i); %old code
    %m.g = m.g + sum(sum( sign_hor.*tr + sign_ver.*tr ));
    for j = 1 : size(sign_hor,1)
      w        = n-1-size(sign_hor,1) + j;
      m.g(w-1) = m.g(w-1) + sum(sum( diag(sign_hor,j-1) + diag(sign_ver,j-1) ) );
    end

  end%i-loop

  %m.rg = m.g / ( 2 * ( n * window*(window-1)/2 - (window^3-window)/3 ) ); %old code
  m.g     = cumsum(m.g);
  w_range = (2:n-1)';
  w_norm  = n * w_range .* (w_range-1)  -  2/3*(w_range.^3 - w_range);
  m.rg    = m.g ./ w_norm;

end%func





