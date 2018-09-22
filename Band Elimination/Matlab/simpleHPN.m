function simpleHPN( w0, wz, Q, name )
     
     %Return values in Ohm and uF
     scale_C = 10^(-6);
     k1 = (w0^2)/(wz^2)-1;
     k2 = ((2+k1)*(Q^2))/((2+k1)*(Q^2)+1);
     k = k2*(w0^2)/(wz^2);
     R1 = 1;
     R2 = (Q^2)*((k1+2)^2);
     R3 = 1;
     R4 = (Q^2)*(k1+2);
     C = 1/(Q*(2+k1));
     
     kf = w0;
     km = 1/kf*C*(1/scale_C);
     
     R1 = km*R1;
     R2 = km*R2;
     R3 = km*R3;
     R4 = km*R4;
     
     C1 = 1/(kf*km)*C*(1/scale_C);
     C2 = C1;
     C3 = 1/(kf*km)*k1*C*(1/scale_C);
     
     results = struct ('R1', R1, 'R2', R2, 'R3', R3, 'R4', R4, 'C1', C1, 'C2', C2, 'C3', C3, 'k', k);
     assignin('base', name, results);

end

