function StrategyI( w0, Q, name )
     
     %Return values in kOhm and uF
     scale_C = 10^(-7);
     R1 = 1;
     R2 = 4*Q^2;
     C = 1/(2*Q);
     
     kf = w0;
     km = 1/kf*C*(1/scale_C);
     
     R1 = km*R1;
     R2 = km*R2;
     C1 = 1/(kf*km)*C;
     C2 = C1;
     
     k = 2*Q^2;
     
     results = struct ('R1', R1, 'R2', R2,'C1', C1, 'C2', C2, 'k', k);
     assignin('base', name, results);

end

