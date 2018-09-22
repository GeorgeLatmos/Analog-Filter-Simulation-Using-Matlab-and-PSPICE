function Strategy_Enhancement( w0, Q, name )
     
     scale_C = 10^-7;
     beta = 1;
     
     C1 = 1;
     C2 = 1;
     R1 = 1/sqrt(beta);
     R2 = sqrt(beta);
     k = (Q*(beta+2)-sqrt(beta))/(2*Q-sqrt(beta));
     
     kf = w0;
     km = 1/kf*C1*(1/scale_C);
     
     R1 = km*R1;
     R2 = km*R2;
     C1 = 1/(kf*km)*C1;
     C2 = 1/(kf*km)*C2;
     
     RA = 1000;
     RB = (k-1)*RA;
     
     w0_new = 1/sqrt(R1*R2*C1*C2);
     Q_new = sqrt(R2/R1)/(sqrt(C2/C1)+sqrt(C1/C2)-1/(k-1*R2/R1*sqrt(C1/C2)));
     
     TF = tf([-2*Q_new*w0_new 0],[1 w0_new/Q_new w0_new^2]);
     
     w0 = w0_new;
     Q = Q_new;
     
     s = j*w0;
     gain_at_w0 = norm(evalfr(TF, s));
     
     results = struct ('R1', R1, 'R2', R2, 'C1', C1, 'C2', C2, 'RA', RA, 'RB', RB, 'w0', w0, 'Q', Q, 'TF', TF, 'Gain_w0', gain_at_w0);
     assignin('base', name, results);

end
