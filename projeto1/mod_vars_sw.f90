module vars_sw
    use vars
    implicit none
    real(8),dimension(:),allocatable             :: eig
    real(8),dimension(:,:),allocatable           :: diag_pos, diag_neg
    real(8),dimension(:,:),allocatable           :: inv_t_xi, t_xi
    real(8),dimension(:,:),allocatable           :: A_pos, A_neg
    real(8),dimension(:,:),allocatable           :: inv_teta, teta
    real(8),dimension(:,:),allocatable           :: B_pos, B_neg
    real(8),dimension(:,:),allocatable           :: aux_mult
    real(8),dimension(:,:),allocatable           :: B_sys, aux_deltaQ, aux_deltaQ_til
    real(8),dimension(:,:,:),allocatable         :: flux_residue, Identy
    real(8),dimension(:,:,:),allocatable         :: A_pos_sys, A_neg_sys
    real(8),dimension(:,:,:),allocatable         :: B_pos_sys, B_neg_sys
    real(8),dimension(:,:,:),allocatable         :: deltaQ, deltaQ_til
    real(8),dimension(:,:,:),allocatable         :: main_sw,lower_sw,upper_sw
    integer(4)                                   :: index, index_i, index_j
    real(8)                                      :: eig_pos, eig_neg
    real(8)                                      :: rho_t
contains
    
    subroutine allocate_vars_sw
        implicit none
        allocate(eig(dim))
        allocate(diag_pos(dim,dim),diag_neg(dim,dim))
        allocate(flux_residue(imax,jmax,dim))
        allocate(deltaQ(imax,jmax,dim))
        allocate(deltaQ_til(imax,jmax,dim))
        allocate(aux_mult(dim,dim))
    end subroutine allocate_vars_sw

    subroutine allocate_vars_sys_ksi
        implicit none
        allocate(main_sw(dim,dim,imax-2),lower_sw(dim,dim,imax-2),upper_sw(dim,dim,imax-2))
        allocate(inv_t_xi(dim,dim),t_xi(dim,dim))
        allocate(A_pos(dim,dim),A_neg(dim,dim))
        allocate(A_pos_sys(dim,dim,imax),A_neg_sys(dim,dim,imax),Identy(dim,dim,imax))
        allocate(aux_deltaQ_til(dim,imax-2))
        allocate(B_sys(dim,imax-2))
        diag_neg       = 0.0d0
        diag_pos       = 0.0d0
        lower_sw       = 0.0d0
        main_sw        = 0.0d0 
        upper_sw       = 0.0d0
        B_sys          = 0.0d0
        aux_deltaQ_til = 0.0d0
        deltaQ_til     = 0.0d0
        deltaQ         = 0.0d0  
        A_pos_sys      = 0.0d0
        A_pos          = 0.0d0
        A_neg_sys      = 0.0d0
        A_neg          = 0.0d0
        Identy         = 0.0d0
        inv_t_xi       = 0.0d0
        t_xi           = 0.0d0
        flux_residue   = 0.0d0
        aux_mult       = 0.0d0
    end subroutine allocate_vars_sys_ksi

    subroutine deallocate_vars_ksi
        implicit none
        deallocate(main_sw,lower_sw,upper_sw)
        deallocate(inv_t_xi,t_xi)
        deallocate(A_pos,A_neg)
        deallocate(A_pos_sys,A_neg_sys)
        deallocate(aux_deltaQ_til)
        deallocate(B_sys)
        deallocate(Identy)
    end subroutine deallocate_vars_ksi

    subroutine allocate_vars_sys_eta
        implicit none
        allocate(aux_deltaQ(dim,jmax-2))
        allocate(main_sw(dim,dim,jmax-2),lower_sw(dim,dim,jmax-2),upper_sw(dim,dim,jmax-2))
        allocate(inv_teta(dim,dim),teta(dim,dim))
        allocate(B_pos(dim,dim),B_neg(dim,dim),Identy(dim,dim,jmax))
        allocate(B_pos_sys(dim,dim,jmax),B_neg_sys(dim,dim,jmax))
        allocate(B_sys(dim,jmax-2))

        Identy     = 0.0d0
        aux_deltaQ = 0.0d0
        lower_sw   = 0.0d0
        main_sw    = 0.0d0
        upper_sw   = 0.0d0
        B_sys      = 0.0d0
        B_pos_sys  = 0.0d0
        B_neg_sys  = 0.0d0  
        teta       = 0.0d0
        inv_teta   = 0.0d0 
    end subroutine allocate_vars_sys_eta

    subroutine deallocate_vars_left
        implicit none
        deallocate(aux_deltaQ)
        deallocate(main_sw,lower_sw,upper_sw)
        deallocate(B_sys,B_pos,B_neg)
        deallocate(B_pos_sys,B_neg_sys)
        deallocate(Identy)
        
        deallocate(inv_teta,teta)
        deallocate(eig)
        deallocate(diag_pos,diag_neg)
        deallocate(flux_residue)
        deallocate(aux_mult)
        deallocate(deltaQ)
        deallocate(deltaQ_til)        
    end subroutine deallocate_vars_left
end module vars_sw