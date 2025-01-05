function C = contourData(nContNodes,flag_geom,nCorners)
% Generate the integral boundary and quadrature nodes;
% Only implement "star";
%
% nContNodes:number of contour nodes/quadrature nodes
% flag_geom: "star"
% 
%
% Second-order derivative only for double layer potential.
%

if strcmp(flag_geom,'star')
  r        = 0.3; 
  k        = nCorners;   
  tt       = linspace(0,2*pi*(1 - 1/nContNodes),nContNodes);
  C        = zeros(6, nContNodes);
  % in order: G_1(t),G_1'(t),G_1''(t), G_2(t),G_2'(t),G_2''(t)
  C(1,:)   =   1.5*cos(tt) + (r/2)*            cos((k+1)*tt) + (r/2)*            cos((k-1)*tt);   
  C(2,:)   = - 1.5*sin(tt) - (r/2)*(k+1)*      sin((k+1)*tt) - (r/2)*(k-1)*      sin((k-1)*tt);  
  C(3,:)   = - 1.5*cos(tt) - (r/2)*(k+1)*(k+1)*cos((k+1)*tt) - (r/2)*(k-1)*(k-1)*cos((k-1)*tt);    
  C(4,:)   =       sin(tt) + (r/2)*            sin((k+1)*tt) - (r/2)*            sin((k-1)*tt);   
  C(5,:)   =       cos(tt) + (r/2)*(k+1)*      cos((k+1)*tt) - (r/2)*(k-1)*      cos((k-1)*tt);  
  C(6,:)   = -     sin(tt) - (r/2)*(k+1)*(k+1)*sin((k+1)*tt) + (r/2)*(k-1)*(k-1)*sin((k-1)*tt);   
else
  fprintf(1,'This option for the geometry is not implemented.\n')
end

end