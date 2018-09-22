function StrategyII( w0, Q, name )
     
     scale_C = 10^-8;
     C1 = 1;
     C2 = 1;
     R1 = 1/(2*Q);
     R2 = 2*Q;
     
     kf = w0;
     km = 1/kf*C1*(1/scale_C);
     
     R1 = km*R1;
     R2 = km*R2;
     C1 = 1/(kf*km)*C1;
     C2 = 1/(kf*km)*C2;
     
     w0_new = 1/sqrt(R1*R2*C1*C2);
     Q_new = 1/(sqrt((R1*C2)/(R2*C1))+sqrt((R1*C1)/(R2*C2)));
     
     TF = tf([1 0 0],[1 w0_new/Q_new w0_new^2]);
     
     results = struct ('R1', R1, 'R2', R2, 'C1', C1, 'C2', C2, 'w0', w0_new,'Q', Q_new,'TF', TF);
     assignin('base', name, results);

end

