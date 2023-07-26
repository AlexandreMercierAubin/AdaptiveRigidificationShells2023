clear;
dimension = 3;
F = sym( 'F', [dimension,dimension], 'real' );
Eprev = sym( 'Eprev', [dimension,dimension], 'real' );
h = sym( 'h', 'real'  );

E = 0.5.*(F'*F - eye(dimension));

Edot = (E-Eprev)./h;
Edot = simplify(Edot);
ccode( Edot, 'File', 'mexEdiff.c' );