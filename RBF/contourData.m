function C = contourData(nContNodes,flag_geom,nCorners)
% Generate the integral boundary and quadrature nodes;
% Only implement "star";
%
%

if strcmp(flag_geom,'star')
  r        = 0.3; 
  k        = nCorners;   
  tt       = linspace(0,2*pi*(1 - 1/nContNodes),nContNodes);
  C        = zeros(2,nContNodes);
  C(1,:)   =   1.5*cos(tt) + (r/2)*            cos((k+1)*tt) + (r/2)*            cos((k-1)*tt);   
  C(2,:)   =       sin(tt) + (r/2)*            sin((k+1)*tt) - (r/2)*            sin((k-1)*tt);   
else
  fprintf(1,'This option for the geometry is not implemented.\n')
end
C = C';
end