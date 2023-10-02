function SS304 = SS304(T)
%Calculates the Specific Heat SH, Thermal Condutivity TC, Young Moduli YM
%and Length Expansion LE coefficients given a Temperature T of the
%StainlessSteal 304, output is a vector
SS304_nist = readtable("Material Properties/SS304-TC_SH_YM_YM_LE.tsv", "FileType", "text", 'Delimiter', '\t');


%Polinomial Expansion coefficients
C = table2array(SS304_nist(1:9,2:6));

%Polinomial Expansion terms
n=0:1:8;
Lx = log10(T).^n;
Tx = T.^n;
%initializing results vector in order TC_SH_YM1_YM2_LE
SS304 = zeros(5,1); 

for i = 1:5
    if i<=2
        SS304(i) = 10^dot(C(:,i),Lx);
    else
        SS304(i) = dot(C(:,i),Tx);
    end
end
