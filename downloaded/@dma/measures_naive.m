%Same as measures.m but implemented in a naive way (e.g. with for-loops),
%for the sole purpose of cross-checking the results of measures.m
%e.g., try for some D:
%struct2mat= @(x) cell2mat(struct2cell(x))
%diff= struct2mat(dma.measures(D)) - struct2mat(dma.measures_naive(D))

function [m] = measures_naive(D)

  [n,d] = size(D);

  m.hp    = 0;
  m.in    = 0;
  m.ls    = 0;
  m.arc   = 0;
  m.warc  = 0;
  m.are   = 0;
  m.ware  = 0;
  m.dware = 0;
  m.g     = zeros(n-2,1);
  m.rg    = zeros(n-2,1);
  
  m.moe   = 0;
  m.mstr  = 0;
  m.nstr  = 0;
  m.bcc   = 0;
  
  
  
  for i=1:n
    for j=1:d
      if j==1, l=0; else l=D(i,j-1); end
      if j==d, r=0; else r=D(i,j+1); end
      if i==1, t=0; else t=D(i-1,j); end
      if i==n, b=0; else b=D(i+1,j); end
      m.moe= m.moe + D(i,j) * (l+r+t+b);
    end
  end
  m.moe = m.moe / 2;
  
  
  
  for i=1:n
    for j=1:d
      for k=max(1,i-1):min(n,i+1)
        for l=max(1,j-1):min(d,j+1)
          
          m.mstr = m.mstr + (D(i,j) - D(k,l))^2;
          
        end
      end
    end
  end
  
  
  
  for i=1:n
    for j=1:d
        
      for k=max(1,i-1):min(n,i+1)
        m.nstr = m.nstr + (D(i,j) - D(k,j))^2;
      end
      
      for l=max(1,j-1):min(d,j+1)
        m.nstr = m.nstr + (D(i,j) - D(i,l))^2;    
      end
      
    end
  end
  
  
  
  for i=2:n
    for j=1:d-1
      for k=1:i-1
        for l=j+1:d
          m.bcc = m.bcc + D(i,j)*D(k,l);
        end
      end
    end
  end
  %also:
%   m.bcc2 = 0;
%   for r = 2 : n
%     for c = 1 : d-1
%       m.bcc2 = m.bcc2 + D(r,c) * sum(sum( D(1:r-1,c+1:d) ));
%     end
%   end
  



  if n ~= d, return;   end
  
  
  div = zeros(n-2,1);
  for i=1:n-2
    for j=i+2:n
      for k=i+1:j-1
        Dh = D(i,j)-D(i,k);
        Dv = D(i,j)-D(k,j);
        
        m.arc  = m.arc  + sign(Dh) + sign(Dv);
        m.warc = m.warc +      Dh  +      Dv;
        
        if Dh < 0
          m.are   = m.are   - sign(Dh);
          m.ware  = m.ware  - Dh;
          m.dware = m.dware - Dh*abs(j-k);
        end
        if Dv < 0
          m.are   = m.are   - sign(Dv);
          m.ware  = m.ware  - Dv;
          m.dware = m.dware - Dv*abs(i-k);
        end

      end
    end
    
    for w=2:n-1
      for j=i+2:min(i+w,n)
        for k=i+1:j-1
          Dh  = D(i,k) - D(i,j);
          Dv  = D(k,j) - D(i,j);
          if Dh > 0
            m.g(w-1) = m.g(w-1) + sign(Dh);
          end
          if Dv > 0
            m.g(w-1) = m.g(w-1) + sign(Dv);
          end
          div(w-1) = div(w-1) + 2;
        end
      end
    end
    
  end
  
  for w=2:n-1
    m.rg(w-1) = m.g(w-1) / div(w-1);
  end
  
  
  for i=1:n-1
    m.hp = m.hp + D(i,i+1);
  end
  
  
  
  for i = 1:n
    for j = 1:n
      m.ls = m.ls + (D(i,j) - abs(i-j))^2;
    end
  end
  
  
  
  for i = 1:n
    for j = 1:n
      m.in = m.in + D(i,j)*(i-j)^2;
    end
  end
  
  
  
end%func




