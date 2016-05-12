function [l2, sh1, h1, peakdiff] = errorNLP (basename, range, vmin, vmax)
  
  l2  = zeros(size(range));
  sh1 = zeros(size(range));
  h1  = zeros(size(range));
  
  peakdiff = zeros(size(range));
  
  for i = 1:length(range)    
    filename = sprintf (basename, range(i));
    out = csvread (filename, 1, 0, 'emptyvalue', NaN);
    
    valid = ! isnan (out(:, 1));
    [Vexp, Cexp, dCexp]  = deal (out(valid, 1),
                                 out(valid, 2),
                                 out(valid, 3));
    
    [Vnum, Cnum, dCnum]  = deal (out(:, 4),
                                 out(:, 5),
                                 out(:, 6));
    
    if (! exist ('vmin'))
      vmin = max ([min(Vexp) min(Vnum)]);
    endif                  
    
    if (! exist ('vmax'))
      vmax = min ([max(Vexp) max(Vnum)]);
    endif                  
    
    V = linspace (vmin, vmax, 1e4);
    
    Cexp = interp1 (Vexp, Cexp, V, 'extrap');
    dCexp = interp1 (Vexp, dCexp, V, 'extrap');
    
    Cnum = interp1 (Vnum, Cnum, V, 'extrap');
    dCnum = interp1 (Vnum, dCnum, V, 'extrap');
    
#     [~, Mexp] = max(dCexp);
#     [~, Mnum] = max(dCnum);
#     
#     Vshift = V(Mnum) - V(Mexp);
    Vshift = 0;
    Vcomp = V - Vshift;
    
#     plot(V, dCexp, Vcomp, dCnum); pause(0.1);
    
    keep = (Vcomp > vmin) & (Vcomp < vmax);
    Vcomp = Vcomp(keep);
    
    l2(i)  = sqrt (trapz (Vcomp, (Cexp - Cnum)(keep).^2));
    sh1(i) = sqrt (trapz (Vcomp, (dCexp - dCnum)(keep).^2));
    h1(i)  = sqrt (l2(i).^2 + sh1(i).^2);
    
    peakdiff(i) = abs(max(dCnum) - max(dCexp));
  endfor

endfunction
