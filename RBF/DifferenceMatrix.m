% Gregory E. Fasshauer, page 342, program 37.2

function DM = DifferenceMatrix( dataXcoord, centerXcoord )

% dataXcoord:  [ A B C D ], M-by-1
% centerXcoord:[ a b c ],   N-by-1
% return a matrix of M-by-N.

% XX=       YY=
% A A A     a b c
% B B B     a b c
% C C C     a b c
% D D D     a b c

[ XX, YY ] = ndgrid( dataXcoord(:), centerXcoord(:) );
DM = XX - YY;