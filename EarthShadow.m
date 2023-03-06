function bShadow=EarthShadow(Si,Ri)
bShadow=0; Rlim=6379;
cosfai= dot(Si,Ri)/norm(Si)/norm(Ri);

if cosfai<0
    sinfai=sqrt(1-cosfai^2);
    if  sinfai<Rlim/norm(Ri); 
        bShadow=1;
    end
end