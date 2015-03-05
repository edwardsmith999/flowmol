nsurfacebins=	2*( globalnbins(1)   * globalnbins(2) ...
    +  (globalnbins(3)-2)* globalnbins(2) ...
    +  (globalnbins(3)-2)*(globalnbins(1)-2));
n = 1;
for kcell=1:globalnbins(3)
    for jcell=1:globalnbins(2)
        for icell=1:globalnbins(1)
            
            %Remove inner part of domain
            if((icell > (1) && icell < (globalnbins(1))) && ...
               (jcell > (1) && jcell < (globalnbins(2))) && ...
               (kcell > (1) && kcell < (globalnbins(3))))
                continue
            end
            %Remove outer cells leaving only 1 layer of surface cells
            if((icell < (1) || icell > (globalnbins(1))) || ...
               (jcell < (1) || jcell > (globalnbins(2))) || ...
               (kcell < (1) || kcell > (globalnbins(3))))
                continue
            end
            
            surfacebins(n,1)=icell;
            surfacebins(n,2)=jcell;
            surfacebins(n,3)=kcell;
            n = n + 1
            
        end
    end
end