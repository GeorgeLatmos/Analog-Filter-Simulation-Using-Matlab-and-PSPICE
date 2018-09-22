function PassiveFilter( w0, name )

    scale_C = 10^-8;
    R1 = 1;
    C1 = 1;
    
    kf = w0;
    km = 1/kf*C1*(1/scale_C);
    
    R1 = km*R1;
    C1 = 1/(kf*km)*C1;
    
    TF = tf([1 0], [1 1/(C1*R1)]);
    
    results = struct ('R1', R1, 'C1', C1, 'TF', TF);
    assignin('base', name, results);


end

