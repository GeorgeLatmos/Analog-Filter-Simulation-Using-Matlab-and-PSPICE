function BoctorLPN( w0, wz, Q, name )
        
        %Return values in kOhm and uF
        C = 10^(-6);

        k1 = 0.98;
        Wz = wz/w0;

        R1 = 2/(k1*Wz^2-1);
        R2 = 1/(1-k1);
        R3 = 1/2*(k1/(Q^2)+k1*Wz^2-1);
        R4 = 1/k1;
        R5 = 1;
        R6 = R5;
        C1 = k1/(2*Q);
        C2 = 2*Q;
        k = 1/(k1/(Q^2)+k1*Wz^2+1);

        kf = w0;
        km = 1/kf*C2*(1/C);

        R1 = km*R1/1000;
        R2 = km*R2/1000;
        R3 = km*R3/1000;
        R4 = km*R4/1000;
        R5 = km*R5/1000;
        R6 = km*R6/1000;
        C1 = 1/(km*kf)*C1*10^6;
        C2 = 1/(km*kf)*C2*10^6;
        
        results = struct ('R1', R1, 'R2', R2, 'R3', R3, 'R4', R4, 'R5', R5, 'R6', R6, 'C1', C1, 'C2', C2, 'k', k);
        assignin('base', name, results);

end

