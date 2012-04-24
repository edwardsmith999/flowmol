nsurfacebins=	2*( gnbins(1)   * gnbins(2) ...
    +  (gnbins(3)-2)* gnbins(2) ...
    +  (gnbins(3)-2)*(gnbins(1)-2));
n = 1;
for kcell=1:gnbins(3)
    for jcell=1:gnbins(2)
        for icell=1:gnbins(1)
            
            %Remove inner part of domain
            if((icell > (1) && icell < (gnbins(1))) && ...
               (jcell > (1) && jcell < (gnbins(2))) && ...
               (kcell > (1) && kcell < (gnbins(3))))
                continue
            end
            %Remove outer cells leaving only 1 layer of surface cells
            if((icell < (1) || icell > (gnbins(1))) || ...
               (jcell < (1) || jcell > (gnbins(2))) || ...
               (kcell < (1) || kcell > (gnbins(3))))
                continue
            end
            
            surfacebins(n,1)=icell;
            surfacebins(n,2)=jcell;
            surfacebins(n,3)=kcell;
            n = n + 1
            
        end
    end
end