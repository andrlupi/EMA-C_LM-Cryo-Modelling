function Ts = Tfinal(R1,R2,Tc,Te,qb)
    a = (R1*R2*qb)+(Tc*R2)+(Te*R1);

    b = (R1+R2);

    Ts = a/b;
end