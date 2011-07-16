subroutine arrayzero
!cccccccccccccccccccccccccc
!   initialize arrays
!cccccccccccccccccccccccccc
        use data_export
        use mesh_export
   
        !-- Pencil Layout ---
        p  = 0.0
        dp = 0.0

        uc = 0.0
        vc = 0.0
        wc = 0.0

        u = 0.0
        v = 0.0
        w = 0.0

        conx = 0.0
        cony = 0.0
        conz = 0.0
  
        rhs = 0.0
   
        !-- Transpose Layout ---
        phatr = 0.0
        !COMPAQ qT    = 0.0                !Initialized before allocation==> crash on compaq
   
        !-- Other variables ---
        cpu_wavet = 0.0

        vkz    = 0.0
        omegak = 0.0

        return
end

