
module grid_export

        include "param.inc"

        real         suxix  (0:ngx+1,0:ngy  ),suxiy  (0:ngx+1,0:ngy  )   &
                ,svetax (0:ngx  ,0:ngy+1),svetay (0:ngx  ,0:ngy+1)   &
                ,spz    (0:ngx  ,0:ngy  ),swz    (0:ngx  ,0:ngy  )   &
                ,surxix (0:ngx+1,0:ngy  ),surxiy (0:ngx+1,0:ngy  )   &
                ,svretax(0:ngx  ,0:ngy+1),svretay(0:ngx  ,0:ngy+1)   &
                ,swrz   (0:ngx  ,0:ngy  ),vp     (0:ngx  ,0:ngy  )
        real        zpp(ngz-1),xpp(ngx-1,ngy-1),ypp(ngx-1,ngy-1),   &
                 zpu(ngz-1),xpu(ngx  ,ngy-1),ypu(ngx  ,ngy-1),   &
                 zpv(ngz-1),xpv(ngx-1,ngy  ),ypv(ngx-1,ngy  ),   &
                 zpw(ngz  ),xpw(ngx-1,ngy-1),ypw(ngx-1,ngy-1),   &
                 zpg(ngz  ),xpg(ngx  ,ngy  ),ypg(ngx  ,ngy  )


end module

