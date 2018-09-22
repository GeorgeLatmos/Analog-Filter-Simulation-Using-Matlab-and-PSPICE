function PassiveFilter( w0,p,name )

    scale_C = 10^-6;
    R1 = 1;
    C1 = 1/p;
    
    kf = w0;
    km = 1/kf*C1*(1/scale_C);
    
    R1 = km*R1;
    C1 = 1/(kf*km)*C1;
    
    TF = tf(1/(C1*R1), [1 1/(C1*R1)]);
    
    results = struct ('R1', R1, 'C1', C1, 'TF', TF);
    assignin('base', name, results);


end

