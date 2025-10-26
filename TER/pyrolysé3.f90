program pyro
    implicit none

    integer,parameter ::pr=8
    !integer :: Nx, Ncoeff, Ny 
    integer,parameter :: Nx=50
    integer,parameter :: Ny=50
    integer,parameter :: Ncoeff=(Nx*Ny)**2
    real(pr),parameter :: lambdav=1.
    real(pr),parameter :: lambdap=1.
    real(pr),parameter :: Lm=3*10**(6)
    integer,parameter :: cpp=1000
    integer,parameter :: cpv=1000
    integer :: i,j,n,kmax,Nmax
    real(pr)::var,test,eps, rhop, rhov, Ta, Acoeff,Temp,deltat, deltatx, deltaty,tps, conduct
    ! real(pr), dimension(:,:), allocatable :: D ,A
    ! real (pr), dimension(:), allocatable :: T,b,rho,rcp,rho_star,rcp_star, Aval, Mval
    ! integer, dimension (:), allocatable :: Mi, Mj , Ai, Aj
    !real(pr), dimension(Nx**2,Nx**2):: A
    real(pr),dimension(:),allocatable :: maille1,maille2,maille3,temp1,temp2,temp3,A
    real (pr), dimension(:),allocatable:: T,b,rho,rcp,rho_star,rcp_star
    real (pr), dimension(:),allocatable :: Aval, Mval
    integer,dimension(:),allocatable  :: Mi, Mj , Ai, Aj
    real(pr) :: res, Neumann,T0 
    integer :: k
    
    ! character(len=100) :: ligne1, ligne2 
    ! integer :: position_equal1, position_equal2, Nx1, Nx2, ios

    ! ! Ouvrir le fichier en mode lecture
    ! open(unit=10, file='data.txt', status='old', action='read', iostat=ios)

    ! ! Vérifier si l'ouverture du fichier s'est déroulée avec succès
    ! if (ios /= 0) then
    !     write(*,*) 'Erreur lors de l''ouverture du fichier.'
    !     stop
    ! end if

    ! ! Lire la première ligne du fichier
    ! read(10, '(A)') ligne1

    ! ! Trouver la position du signe égal (=) dans la première ligne
    ! position_equal1 = INDEX(ligne1, '=')

    ! ! Lire la valeur après le signe égal dans la première ligne
    ! read(ligne1(position_equal1+1:), *) Nx1

    ! ! Afficher la valeur de Nx 
    ! write(*,*) 'La valeur de Nx est :', Nx1

    ! ! Lire la deuxième ligne du fichier
    ! read(10, '(A)') ligne2

    ! ! Trouver la position du signe égal (=) dans la deuxième ligne
    ! position_equal2 = INDEX(ligne2, '=')

    ! ! Lire la valeur après le signe égal dans la deuxième ligne
    ! read(ligne2(position_equal2+1:), *) Nx2

    ! ! Fermer le fichier
    ! close(unit=10)

    ! ! Afficher la valeur de Ny 
    ! write(*,*) 'La valeur de Ny est :', Nx2
    
    ! Ncoeff=Nx*(5*Nx-4)


    ! allocate(D(Nx, Nx),A(Nx**2,Nx**2))
    ! allocate(T(Nx*Nx))
    ! allocate(b(Nx*Nx))
    ! allocate(rho(Nx*Nx))
    ! allocate(rcp(Nx*Nx))
    ! allocate(rho_star(Nx*Nx))
    ! allocate(rc p_star(Nx*Nx))
    ! allocate(Aval(Ncoeff), Ai(Ncoeff),Aj(Ncoeff), Mval (Ncoeff))
    ! allocate(Mi(Ncoeff), Mj(Ncoeff))
    Nmax = 100
    allocate(maille1(Nmax),maille2(Nmax),maille3(Nmax), temp1(Nmax),temp2(Nmax),temp3(Nmax))
    allocate(T(Nx*Ny),b(Nx*Ny),rho(Nx*Ny),rcp(Nx*Ny),rho_star(Nx*Ny),rcp_star(Nx*Ny),A(Nx*Ny))
    allocate (Aval(Ncoeff) , Mval(Ncoeff) )
    allocate(Mi(Ncoeff), Mj(Ncoeff) , Ai(Ncoeff), Aj(Ncoeff))

    deltat=1
    deltatx=1./Nx
    deltaty=1./Ny
        !
    rhov=1500.
    rhop=1000.

    !A=1.
    Ai=1
    Aj=1
    Aval=1.
    


    b=1.
    T0=293.
    eps=0.0001
    kmax = 200000


    

   

   
!----------------------------------------------------LES BAILS---------------------------------------------------------------------------------
    !initialisation
    T=T0
    rho=rhov
    rcp=rhov*cpv
    rho_star=rhov
    rcp_star=rhov*cpv
    



    call sauv(0,T,rho)
    !open(unit=15,file='matrice')

    do n=1,Nmax
        print*,real(n)/real(Nmax)*100.,"%"
        tps = (n+1)*deltat
        Neumann = calcule_flux(tps)
        !calcul des valeurs intermédiares
        print*,"etape val inter"
        do i =1,Nx
            do j=1,Ny
                rho_star((i-1)*Ny+j) = calcule_rho(i,j,tps, rhop, rhov, Acoeff, T, Nx,Ny)
                rcp_star((i-1)*Ny+j) = rho_cp(i,j,deltatx,deltaty,tps,T)
            enddo
        enddo
        
        !Systeme linéaire
        print*,"remplissage creux"
        call remplissage_creu(deltatx,deltaty,deltat,tps,T,Ai,Aj,Aval)
       !call  dense(Ai,Aj,Aval)
    
        call remplissage_B(T0,rcp,rcp_star,rho,rho_star,Lm,Neumann,B,deltat,deltaty)
        !print*,B
        !calcul de Tn+1 et rhon+1
        T = gradC(Ai,Aj,Aval,T,B,eps,kmax)
        do i=1,Nx
            do j =1,Ny
                rho((i-1)*Nx+j) = calcule_rho(i,j,tps, rhop, rhov, Acoeff, T,Nx,Ny)
                rcp((i-1)*Nx+j) = rho_cp(i,j,deltatx,deltaty,tps,T)
            enddo
        enddo
        maille1(n+1) = rho(4)
        maille2(n+1) = rho(5)
        maille3(n+1) = rho(76)

        temp1(n+1) = T(1)
        temp2(n+1) = T(361)
        temp3(n+1) = T(721)
        
        if ( mod(n,1) == 0) then
            call sauv(n,T,rho)
            ! A=0.
            ! do k=1,Ncoeff
            !     A(Ai(k),Aj(k))=Aval(k)
            ! end do
            ! do k=1,Nx**2
            !     write(15,*) A(k,:)
            ! end do
            ! write(15,*) '---------------------------------------------------------------------'


        endif
    end do 
close(15)
print*,"50%",int(Ny*0.50)
print*,"5%",int(Ny*0.05)
open(unit=5,file='rho_mailles')
do n=1,Nmax
    write(5,*)maille1(n),maille2(n),maille3(n)
enddo
close(5)

open(unit=6,file='temp_mailles')
do n=1,Nmax
    write(6,*)temp1(n),temp2(n),temp3(n)
enddo
close(6)
    !-------------test
    
    !Mi=(/1.,1.,1.,1.,2.,2.,2.,2.,3.,3.,3.,3.,4.,4.,4.,4./)
    !Mj=(/1.,2.,3.,4.,1.,2.,3.,4.,1.,2.,3.,4.,1.,2.,3.,4./)
    !Mval=(/1.,3.,3.,3.,3.,1.,3.,3.,3.,3.,1.,3.,3.,3.,3.,1./)
    !B = (/1.,2.,3.,4./)
    !T = 0.
    !T = gradC(Mi,Mj,Mval,T,B,eps,kmax)
    !print*,T
    !print*,prodMV(Mi,Mj,Mval,B)
    !-----------------fin de test


    deallocate(maille1,maille2,maille3, temp1, temp2, temp3,T,b,rho,rcp,rho_star,rcp_star,Aval, Mval, Mi, Mj , Ai, Aj)
    ! deallocate(D,A)
    ! deallocate(T,b, rho)
    ! deallocate(rcp, rho_star, rcp_star)
    ! deallocate(Aval, Ai, Aj, Mval)
    ! deallocate(Mi, Mj)
   

!--------------------------------------------------FIN DES BAILS---------------------------------------------------------------------------------------------
    contains


   
    subroutine remplissage_creu(deltatx,deltaty,deltat,tps,T,Ai,Aj,Aval)

        real(pr),intent(in)                           :: deltatx,deltaty,deltat,tps
        real(pr),dimension(:),intent(in)          :: T
        real(pr),dimension(Ncoeff),intent(out)        :: Aval
        integer,dimension(Ncoeff),intent(out)         :: Ai, Aj

        integer                                       :: k,m1,m2,indice

        

        indice = 1

        do m1=1,Nx
            do m2=1,Nx
                if (m1==m2) then !D general
                    do k=1,Ny-1
                        Ai(indice)=(m1-1)*Ny+k
                        Aj(indice)=(m1-1)*Ny+k
                        Aval(indice) = alpha(m2,k,deltatx,deltaty,deltat,tps,T)
                        indice = indice+1
                        !print*,alpha(m2,k,deltatx,deltaty,deltat,tps,T)
                        ! if (m2==1 .and. k==1) then
                        !     print*,"Coin NO"
                        ! endif

                        Ai(indice) = (m1-1)*Ny+k
                        Aj(indice) = (m1-1)*Ny+k+1
                        Aval(indice) = -phi(m2,k+1,deltatx,deltaty,deltat,tps,T)
                        indice = indice+1

                        Ai(indice)=(m1-1)*Ny+k+1
                        Aj(indice)=(m1-1)*Ny+k
                        Aval(indice) = beta(m2,k,deltatx,deltaty,deltat,tps,T)
                        indice = indice+1
                    enddo
                                
                    Ai(indice)=m1*Ny
                    Aj(indice)=m1*Ny
                    Aval(indice) = alpha(m2,Nx,deltatx,deltaty,deltat,tps,T)
                    indice = indice+1
                    
                    

                else if (m2==m1+1) then !matrice B
                    do k=1,Ny
                        Ai(indice)=(m1-1)*Ny+k
                        Aj(indice)=(m2-1)*Ny+k
                        Aval(indice) = -gamma_coeff(m2,k,deltatx,deltaty,deltat,tps,T)
                        indice = indice+1
                    enddo

                else if (m2==m1-1) then !matrice C
                    do k=1,Ny
                        Ai(indice)=(m1-1)*Ny+k
                        Aj(indice)=(m2-1)*Ny+k
                        Aval(indice) = lambda(m2,k,deltatx,deltaty,deltat,tps,T)
                        indice = indice+1
                    enddo

                end if
            enddo
        enddo

    end subroutine remplissage_creu
    
    
    subroutine remplissage_B(T0,rcp,rcp_star,rho,rho_star,Lm,Neumann,B,deltat,deltaty)
    integer::i,j
    real(pr),intent(in)::T0,Lm,Neumann, deltat, deltaty
    real(pr),dimension(:),intent(in)::rcp,rcp_star,rho,rho_star
    real(pr),dimension(Nx*Ny),intent(out)::B

    do i =1,Nx
        do j=1,Ny

            B((i-1)*Ny+j) = T0 + (T((i-1)*Ny+j)-T0)*(rcp((i-1)*Ny+j)/rcp_star((i-1)*Ny+j)) &
            & + Lm*(rho_star((i-1)*Ny+j)-rho((i-1)*Ny+j))/(rcp_star((i-1)*Ny+j))

            if (i==1) then
                B((i-1)*Ny+j) = B((i-1)*Ny+j) + deltat*Neumann/deltaty
            endif
        enddo
    enddo


    end subroutine remplissage_B

    !produit scalaire entre 2 vecteurs
    function prodVV(v1,v2) result(res)
        real(pr),dimension(:),intent(in)::v1,v2
        real(pr)::res
        integer::i
        res=0._pr

        do i=1,size(v1)
            res = res + v1(i)*v2(i)
        enddo
    end function prodVV

    !produit matriciel entre une matrice et un vecteur, retourne un vecteur non creu
    function prodMV(Mi,Mj,Mval,v1) result(v2)
        real(pr),dimension(:),intent(in)::Mval
        integer,dimension(:),intent(in)::Mi,Mj
        real(pr),dimension(:),intent(in)::v1
        real(pr),dimension(Nx*Ny)::v2
        integer::i,k

        v2 = 0.

        do i=1,Nx*Ny
            do k=1,Ncoeff
                if (Mi(k)==i) then
                    v2(i) = v2(i) + Mval(k)*v1(Mj(k))
                endif
            enddo
        enddo

    end function prodMV

    function diag(Ai,Aj,Aval)result(M)
        real(pr),dimension(Ncoeff)::Aval
        integer,dimension(Ncoeff)::Ai,Aj
        real(pr),dimension(Nx*Ny)::M
        integer ::k

        do k=1,Ncoeff
            if (Ai(k)==Aj(k)) then
                M(Ai(k)) = Aval(k)
            endif
        enddo
    end function diag

    !Gradient Conjugué, retourne le vecteur x du systeme Ax=b
   function gradC(Ai, Aj, Aval, x0, b, eps, kmax) result(xk)
        real(pr), dimension(:), intent(in) :: Aval
        integer, dimension(:), intent(in) :: Ai, Aj
        real(pr), dimension(Nx*Ny) :: M
        real(pr), dimension(:), intent(in) :: x0, b
        real(pr), dimension(Nx*Ny) :: xk, xk1, rk, rk1, p, Ap, z
        integer :: k, kmax
        real(pr) :: eps, rh, gam, delt, alph, rh0

        ! Preconditioneur (Jacobi)
        M = diag(Ai, Aj, Aval)  
        do k=1,Nx*Ny
            if (M(k)==0) then
                print*,"DIV0"
            endif
        enddo

        ! Initialization
        xk = x0
        xk1 = x0
        rk = b - prodMV(Ai, Aj, Aval, x0)
        rk1 = rk
        k = 0

        do while (sqrt(prodVV(rk, rk)/prodVV(b,b)) > eps .and. k < kmax)
            z = rk
            rh0=1.

            ! Apply Jacobi preconditioner
            z = z !/ M

            rh = prodVV(rk, z)
            if (k == 0) then
                p = z
            else
                gam = rh / rh0
                p = gam * p + z
            endif

            Ap = prodMV(Ai, Aj, Aval, p)
            delt = prodVV(p, Ap)
            alph = rh / delt
            xk1 = xk + alph * p
            rk1 = rk - alph * Ap
            rh0 = rh
            k = k + 1

            xk = xk1
            rk = rk1
        end do

        !print *, "Number of iterations:", k
    end function gradC


    function calcule_rho(i,j,tps, rhop, rhov, Acoeff, T,Nx,Ny) result(rho)
        real(pr), intent(in)                          :: tps
        real(pr)                                      ::rhop, rhov, Acoeff
        integer, intent(in)                           ::Nx,Ny,i,j
        real(pr),dimension(:), intent(in)         ::T
        real(pr)                                      ::rho
        
  
        Acoeff=1000._pr
        Ta=6000._pr
        rhov=1500._pr
        rhop=1000._pr
        rho = rhop + (rhov - rhop) * exp((rhov*Acoeff*exp(-(Ta/T((i-1)*Ny+j))) * tps)/(rhop - rhov))
        !print*,rho,(rhov*Acoeff*exp(-(Ta/T((i-1)*Nx+j))) * tps)/(rhop - rhov)
    end function calcule_rho

    function avancement(i,j,dx,dy,tps,T) result(res)

        integer                          :: i,j
        real(pr)                         :: dx,dy,tps,r
        real(pr),dimension(:),intent(in)       :: T
        real(Pr)::res

        r = calcule_rho(i,j,tps, rhop, rhov, Acoeff, T,Nx,Ny)
        res = (rhov-r / (rhov-rhop))
        !print*,res
    end function avancement

    function rho_cp(i,j,dx,dy,tps,T) result(res)
        integer                          :: i,j
        real(pr)                         :: dx,dy,tps,r
        real(pr),dimension(:),intent(in)       :: T
        real(Pr)::res

        res = (1._pr-avancement(i,j,dx,dy,tps,T))*rhov*cpv + avancement(i,j,dx,dy,tps,T)*rhop*cpp
        !print*,res
    end function rho_cp

    function conductivite(i,j,dx,dy,tps,T) result(res)

        integer                          :: i,j
        real(pr)                         :: dx,dy,tps,r
        real(pr),dimension(:),intent(in)       :: T
        real(Pr)::res
        res = (1._pr-avancement(i,j,dx,dy,tps,T))*lambdav + avancement(i,j,dx,dy,tps,T)*lambdap

    end function conductivite


    function calcule_flux(tps) result(Neumann)
        real(pr),intent(in)::tps
        real(pr) ::Neumann

        if (tps<= 50.) then
            Neumann = 10000.*tps
        else
            Neumann = 500000. -9000.*(tps-50)
        endif
    end function calcule_flux

    !sauvegarde le vecteur T dans un .txt nommé d'après le temps
    subroutine sauv(n, T, rho)
    implicit none
    integer, intent(in) :: n
    real(pr), dimension(:), intent(in) :: T, rho
    character(50) :: filename
    character(4) :: extension

    ! Save T to file
    extension = ".txt"
    write(filename, '(I0,A2, A4)') n,"_T", extension
    open(unit=10, file=trim(filename), status='replace', action='write')
    do i = 1, size(T, 1)
        write(10, *) T(i)
    end do
    close(10)

    ! Save rho to file
    write(filename, '(I0, A4,A4)') n,"_RHO", extension  ! Use A7 to include space for the extension
    open(unit=11, file=filename, status='replace', action='write')
    do i = 1, size(rho, 1)
        write(11, *) rho(i)
    end do
    close(11)
    end subroutine sauv





    function Diff_ip1demi(i,j,deltatx,deltaty,tps,T) result(res)

        integer                          :: i,j
        real(pr)  ,intent(in)                        :: deltatx,deltaty,tps
        real(pr),dimension(:),intent(in)       :: T
        real(Pr)::res
        
        res = (conductivite(i+1,j,deltatx,deltaty,tps,T)+conductivite(i,j,deltatx,deltaty,tps,T)) &
        & / (2*rho_cp(i,j,deltatx,deltaty,tps,T))

    end function Diff_ip1demi

    function Diff_im1demi(i,j,deltatx,deltaty,tps,T) result(res)

        integer                          :: i,j
        real(pr)  ,intent(in)                        :: deltatx,deltaty,tps
        real(pr),dimension(:),intent(in)       :: T
        real(Pr)::res
        res = (conductivite(i,j,deltatx,deltaty,tps,T)+conductivite(i-1,j,deltatx,deltaty,tps,T)) &
        &/ (2*rho_cp(i,j,deltatx,deltaty,tps,T))

    end function Diff_im1demi

    function Diff_jp1demi(i,j,deltatx,deltaty,tps,T) result(res)

        integer                          :: i,j
        real(pr)  ,intent(in)                        :: deltatx,deltaty,tps
        real(pr),dimension(:),intent(in)       :: T
        real(Pr)::res

        res = (conductivite(i,j+1,deltatx,deltaty,tps,T)+conductivite(i,j,deltatx,deltaty,tps,T)) &
        &/ (2*rho_cp(i,j,deltatx,deltaty,tps,T))

    end function Diff_jp1demi

    function Diff_jm1demi(i,j,deltatx,deltaty,tps,T) result(res)

        integer                          :: i,j
        real(pr)  ,intent(in)                        :: deltatx,deltaty,tps
        real(pr),dimension(:),intent(in)       :: T
        real(Pr)::res

        res = (conductivite(i,j,deltatx,deltaty,tps,T)+conductivite(i,j-1,deltatx,deltaty,tps,T)) &
        & / (2*rho_cp(i,j,deltatx,deltaty,tps,T))
        ! print*,res,i,j
        
    ! CL flux non nul

    end function Diff_jm1demi

    function cl_flux(tps) result(res)
        real(pr) ,intent(in)  :: tps
        real(pr)   :: res

        if ( tps <= 50._PR ) then
            res = 10000*tps
        else
            res = 500000-9000*(tps-50.)
        end if

    end function cl_flux

    ! coeffs pour la matrice A de A*T = b

    function alpha(i,j,deltatx,deltaty,deltat,tps,T) result(res)

        integer                          :: i,j
        real(pr)  ,intent(in)                        :: deltat,deltatx,deltaty,tps
        real(pr),dimension(:),intent(in)       :: T
        real(Pr)::res
        if ( i == 1 .and. j/=1 .and. j/=Nx ) then !bord Nord
            res = 1. + (deltat/(deltatx*deltatx))*(Diff_ip1demi(i,j,deltatx,deltaty,tps,T))&
                    & + (deltat/(deltaty*deltaty))*(Diff_jp1demi(i,j,deltatx,deltaty,tps,T) &
                    &- Diff_jm1demi(i,j,deltatx,deltaty,tps,T))
        else if ( i == Nx .and. j/=1 .and. j/=Nx) then !bord Sud
            res = 1. - (deltat/(deltatx*deltatx))*(Diff_im1demi(i,j,deltatx,deltaty,tps,T))&
                    & + (deltat/(deltaty*deltaty))*(Diff_jp1demi(i,j,deltatx,deltaty,tps,T)&
                    &-Diff_jm1demi(i,j,deltatx,deltaty,tps,T))
        else if ( j == 1 .and. i/=1 .and. i/=Nx) then !bord Ouest
            res = 1. + (deltat/(deltatx*deltatx))*(Diff_ip1demi(i,j,deltatx,deltaty,tps,T)&
                    &- Diff_im1demi(i,j,deltatx,deltaty,tps,T))&
                    & + (deltat/(deltaty*deltaty))*(Diff_jp1demi(i,j,deltatx,deltaty,tps,T))
        else if ( j == Nx .and. i/=1 .and. i/=Nx) then !bord Est
            res = 1. + (deltat/(deltatx*deltatx))*(Diff_ip1demi(i,j,deltatx,deltaty,tps,T)&
                &-Diff_im1demi(i,j,deltatx,deltaty,tps,T)) -(deltat/(deltaty*deltaty))*Diff_jm1demi(i,j,deltatx,deltaty,tps,T)
        else if ( i == 1  .and. j == 1) then !coin NO
            res = 1. + (deltat/(deltatx*deltatx))*Diff_ip1demi(i,j,deltatx,deltaty,tps,T) &
                & + (deltat/(deltaty*deltaty))*Diff_jp1demi(i,j,deltatx,deltaty,tps,T)
        else if ( i == 1 .and. j == Nx ) then !coin NE
            res = 1. + (deltat/(deltatx*deltatx))*Diff_ip1demi(i,j,deltatx,deltaty,tps,T) &
                & - (deltat/(deltaty*deltaty))*Diff_jm1demi(i,j,deltatx,deltaty,tps,T)
            !print*,res
        else if ( i == Nx .and. j == 1 ) then !coin SO
            res = 1. - (deltat/(deltatx*deltatx))*Diff_im1demi(i,j,deltatx,deltaty,tps,T) &
                & + (deltat/(deltaty*deltaty))*Diff_jp1demi(i,j,deltatx,deltaty,tps,T)
        else if ( i == Nx .and. j == Nx ) then !coin SE
            res = 1. - (deltat/(deltatx*deltatx))*Diff_im1demi(i,j,deltatx,deltaty,tps,T) &
                & - (deltat/(deltaty*deltaty))*Diff_jm1demi(i,j,deltatx,deltaty,tps,T)
        
        else !intérieur du maillage
            res = 1. + (deltat/(deltaty*deltaty))*(Diff_jp1demi(i,j,deltatx,deltaty,tps,T)&
                &-Diff_jm1demi(i,j,deltatx,deltaty,tps,T)) &
                & + (deltat/(deltatx*deltatx))*(Diff_ip1demi(i,j,deltatx,deltaty,tps,T)&
                &-Diff_im1demi(i,j,deltatx,deltaty,tps,T))
        end if

    end function alpha

    function beta(i,j,deltatx,deltaty,deltat,tps,T) result(res)

        integer                          :: i,j
        real(pr)  ,intent(in)                        :: deltat,deltatx,deltaty,tps
        real(pr),dimension(:),intent(in)       :: T
        real(Pr)::res

    
        res = deltat/(deltaty*deltaty)*Diff_jm1demi(i,j,deltatx,deltaty,tps,T)
        

    end function beta

    function lambda(i,j,deltatx,deltaty,deltat,tps,T) result(res)

        integer                          :: i,j
        real(pr)  ,intent(in)                        :: deltat,deltatx,deltaty,tps
        real(pr),dimension(:),intent(in)       :: T
        real(Pr)::res

        
        res = deltat/(deltatx*deltatx)*Diff_im1demi(i,j,deltatx,deltaty,tps,T)

    end function lambda


    function gamma_coeff(i,j,deltatx,deltaty,deltat,tps,T) result(res)
        integer                          :: i,j
        real(pr)  ,intent(in)                        :: deltat,deltatx,deltaty,tps
        real(pr),dimension(:),intent(in)       :: T
        real(Pr)::res

       
        res = deltat/(deltatx*deltatx)*Diff_ip1demi(i,j,deltatx,deltaty,tps,T)
        

    end function gamma_coeff

    function phi(i,j,deltatx,deltaty,deltat,tps,T) result(res)

        integer                          :: i,j
        real(pr)  ,intent(in)                        :: deltat,deltatx,deltaty,tps
        real(pr),dimension(:),intent(in)       :: T
        real(Pr)::res

        res = deltat/(deltaty*deltaty)*Diff_jp1demi(i,j,deltatx,deltaty,tps,T)

    end function phi

    subroutine dense(Ai,Aj,Aval)
        real(pr),dimension(:),intent(in)::Aval
        integer,dimension(:),intent(in)::Ai,Aj
        real(pr),dimension(Nx*Ny,Nx*Ny)::A
        integer :: k
        A = 0.

        do k=1,Ncoeff
            A(Ai(k),Aj(k)) = Aval(k)

        enddo
        do k=1,Nx*Ny
            print*,A(k,:)
        end do
    end subroutine dense 


end program pyro
