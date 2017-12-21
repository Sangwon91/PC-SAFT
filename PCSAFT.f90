module mod_PCSAFT

    implicit none

    public ! for global variables

    integer, parameter :: isoT = 1, isoP = 2, fitPara = 3
    integer            :: calType      ! isoT, isoP, fitPara

    integer, parameter :: maxComp = 11 ! max component
    integer, parameter :: maxSite = 5  ! max component
    integer, parameter :: nSite   = 2  ! for cross association model
    integer            :: n            ! number of component

    real(kind=8), parameter :: NA = 6.022e23_8, kb = 13.806044503487215 ! avorgadro number, boltzmann constant<MPa*A^3/K>
    real(kind=8), parameter :: PI = 3.141592653589793

    real(kind=8), dimension(0:6) :: a0, a1, a2, b0, b1, b2 ! dispersion parameters

    real(kind=8) :: T, P ! temper, pressure
    real(kind=8) :: rho  ! density in unit A^-3

    character(len=31), dimension(maxComp) :: mName                 ! molecule name
    real(kind=8), dimension(maxComp) :: x, y, z                    ! liquid fraction, vapor fraction, input fraction
    real(kind=8), dimension(maxComp) :: mass, eps, sig, seg        ! molecule weight, epsilon, sigma, segment number
    real(kind=8), dimension(maxComp) :: deAB, dkAB                 ! epsilon, sigma, segment number
    real(kind=8), dimension(maxComp,maxComp) :: eAB, kAB           ! association well edpth, association volume
    real(kind=8), dimension(maxComp,maxComp) :: kij                ! binary interaction parameters
    real(kind=8), dimension(maxComp,maxSite) :: Xs, XsOld, difX    ! Xs(i,j) is the mole fraction of component i that is not bonded at site j

    real(kind=8), dimension(maxComp) :: d
    real(kind=8), dimension(maxComp) :: giihs
    real(kind=8), dimension(maxComp,maxComp) :: gijhs

    real(kind=8) :: m
    real(kind=8) :: ze0, ze1, ze2, ze3, ze3_1, pi6 = PI / 6.0_8    ! zeta, zeta^n and 1 - zeta
    real(kind=8) :: eta, e1, e2, e3, e4, e_1, e_2
    real(kind=8), dimension(0:6) :: ea ! eta array ea(n) = eta^n
    real(kind=8) :: c1, c2
    real(kind=8) :: m2es3, m2e2s3
    real(kind=8) :: I1, I2
    real(kind=8) :: ahs, ahc, adisp, anon, aasc, atot, ares
    real(kind=8), dimension(maxComp,maxComp) :: dAB        ! delAB

    real(kind=8), dimension(:,:), allocatable :: Lpq, iLpq ! Lmabdapq : Radosz method
    real(kind=8), dimension(50) :: psi, dXs
    logical :: singularQ, criQ
    integer :: dp, dq, dim_

    ! for x derivation (notation : df/dx = fx)

    real(kind=8), dimension(maxComp,0:3) :: zex ! d zeta / dxi
    real(kind=8), dimension(maxComp,maxComp) :: giihsx
    real(kind=8), dimension(maxComp,maxComp,maxComp) :: gijhsx
    real(kind=8), dimension(maxComp) :: c1x, c2x
    real(kind=8), dimension(maxComp) :: m2es3x, m2e2s3x
    real(kind=8), dimension(maxComp) :: I1x, I2x
    real(kind=8), dimension(maxComp) :: ahsx, ahcx, adispx, anonx, aascx, atotx, aresx
    real(kind=8), dimension(maxComp,maxComp,maxComp) :: dABx
    real(kind=8), dimension(maxComp,maxComp,maxSite) :: Xsx
    real(kind=8), dimension(maxComp) :: mures, phi

    ! for rho derivation (notation : df/d(rho) = fd)

    real(kind=8), dimension(maxComp,maxComp) :: gijhsd ! d indicate rho(density)
    real(kind=8), dimension(maxComp,maxComp) :: dABd
    real(kind=8), dimension(maxComp,maxSite) :: Xsd
    real(kind=8) :: zhs, zhc, zdisp, znon, zasc, zres, ztot, ptot
    real(kind=8) :: deI1, deI2
    real(kind=8) :: anond

    real(kind=8) :: dmax, pmax, dmin, pmin

    real(kind=8) :: x0, x1, dx, T0, T1, dT

    ! mixture calculation

    real(kind=8), dimension(maxComp) :: fl, fv, Ke   ! Ke = equilibium constant
    real(kind=8), dimension(maxComp) :: srhol, srhov ! if Cal Psat save density
    real(kind=8), dimension(maxComp) :: gam          ! activity coefficient

end module mod_PCSAFT
!----------------------------------------------------------------------------------------------------
program main

    use mod_PCSAFT

    implicit none

    integer :: c
    character(len=61) :: filename

    call system("title PC-SAFT equation of state calculator")
    call system("mode con: cols=120 lines=2500")
    call system("color 2f")

    call system("date/t")
    call system("time/t")

    c = 1
    write(*,"('--- calculation start ---')")
    open(99,file='node.txt')
     do while (.true.)
      read(99,"(A)",end=1) filename
      call GetInf(filename)
!      go to 2
      select case (caltype)
       case (0); call PsatAuto ('output_'//filename,c)
       case (1); call Pxy ('output_'//filename)
       case (2); call Txy ('output_'//filename)
       case (3); call PrintPeos ('output_'//filename, x, T)
       case (4)
         call bublP
         call dewP
       case (5); call FindBinary
       case default;
      end select
!2     call calGam
     end do
    close(99)
    deallocate(Lpq,iLpq)
1   write(*,"('--- calculation end ---')")
    call system('pause')

end program main
!----------------------------------------------------------------------------------------------------
subroutine GetInf(filename)

    use mod_PCSAFT

    implicit none

    character(len=61), intent(in) :: filename

    integer :: i, j
    logical :: error

    error = .true.

    open(1, file=filename,err=2)

        read(1,*) n

        if (n==0) stop 'error : zero component'
        do i = 1, n
          read(1, *,err=1) mName(i), mass(i), seg(i), sig(i), eps(i), dkAB(i), deAB(i)
        end do
        if (n==2) then
         do i = 1, n
          read(1,*,err=1) kij(i,1:n)  ! read binary interaction parameter matrix
         end do
         read(1,*,err=1) calType
         if (calType == 5   ) write(*,"('fitting')")
         if (calType == isoT) read(1,*,err=1) T
         if (calType == isoP) read(1,*,err=1) P
         read(1,*,err=1) x0
         read(1,*,err=1) x1
         read(1,*,err=1) dx
        else if (n==1) then
         calType = 0
         read(1,*,err=1) T0
         read(1,*,err=1) T1
         read(1,*,err=1) dT
        else
         calType = 4 ! for test
         do i = 1, n
          read(1,*,err=1) kij(i,1:n)  ! read binary interaction parameter matrix
         end do
         read(1,*,err=1) T
         read(1,*,err=1) P
         read(1,*,err=1) z(1:n)
        end if

        error = .false.

    close(1)

    do i=1, n; do j=1, n;
      if (i==j) kij(i,j) = 0.0
      if (i<j)  kij(j,i) = kij(i,j)
    end do; end do

    ! check information
    write(*,"('component information =>')")
    write(*,"()")
    do i = 1, n
      write(*,"('com',I2,' = ',A31)") i, mName(i)
    end do
    write(*,"()")

    write(*,"('noc = ', I10)") n
    write(*,"('com = ', 20I10)") (i,i=1,n)
    write(*,"('eps = ', 20F10.4)") eps(1:n)
    write(*,"('sig = ', 20F10.4)") sig(1:n)
    write(*,"('seg = ', 20F10.4)") seg(1:n)
    write(*,"('eAB = ', 20F10.4)") deAB(1:n)
    write(*,"('kAB = ', 20F10.4)") dkAB(1:n)

    write(*,"('kij = ')")
    do i = 1, n
        write(*,"('      ',11F10.4)") kij(i,1:n)
    end do
    write(*,"('calType = ',I3)") calType
    if (calType==1.or.calType==2) then
     write(*,"('x0  = ',F10.4)") x0
     write(*,"('x1  = ',F10.4)") x1
     write(*,"('dx  = ',F10.4)") dx
    else if (calType==0) then
     write(*,"('T0  = ',F10.4)") T0
     write(*,"('T1  = ',F10.4)") T1
     write(*,"('dT  = ',F10.4)") dT
    else if (calType==4) then
     write(*,"('T   = ',F10.4)") T
     write(*,"('P   = ',F10.4)") P
     write(*,"('zi  = ',20F10.4)") z(1:n)
     write(*,"('sum of z = ', F10.4)") sum(z(1:n))
    end if

    ! use binary mixing rule
    eAB = 0.0
    kAB = 0.0
    do i=1,n; do j=1,n
        eAB(i,j) = 0.5*(deAB(i)+deAB(j))
        kAB(i,j) = (dkAB(i)*dkAB(j))**0.5 * (2.0*(sig(i)*sig(j))**0.5/(sig(i)+sig(j)))**3
    end do; end do
    if (allocated(Lpq)) deallocate(Lpq,iLpq)
    allocate(Lpq(n*2,n*2),iLpq(n*2,n*2))

1   if (error) stop 'input read error'
2   if (error) stop 'file open faild'

end subroutine GetInf
!----------------------------------------------------------------------------------------------------
subroutine SetNormalTerms (f,rho_,T_) ! mole fraction, density, temper

    use mod_PCSAFT

    implicit none

    real(kind=8), intent(in) :: rho_, T_
    real(kind=8), dimension(maxComp) :: f ! mol fraction
    integer :: i, ii, j, jj
    real(kind=8), dimension(maxComp) :: dgiihsd
    real(kind=8), dimension(0:6)     :: ind  ! indexer
    real(kind=8) :: hm                   ! harmonic mean
    real(kind=8) :: del, sum1, sum2      ! for convergence test

    ! constant for dispersion

    a0(0) = 0.9105631445
    a0(1) = 0.6361281449
    a0(2) = 2.6861347891
    a0(3) =-26.547362491
    a0(4) = 97.759208784
    a0(5) =-159.59154087
    a0(6) = 91.297774084

    a1(0) =-0.3084016918
    a1(1) = 0.1860531159
    a1(2) =-2.5030047259
    a1(3) = 21.419793629
    a1(4) =-65.255885330
    a1(5) = 83.318680481
    a1(6) = -33.746922930

    a2(0) =-0.0906148351
    a2(1) = 0.4527842806
    a2(2) = 0.5962700728
    a2(3) =-1.7241829131
    a2(4) =-4.1302112531
    a2(5) = 13.776631870
    a2(6) =-8.6728470368

    b0(0) = 0.7240946941
    b0(1) = 2.2382791861
    b0(2) =-4.0025849485
    b0(3) =-21.003576815
    b0(4) = 26.855641363
    b0(5) = 206.55133841
    b0(6) =-355.60235612

    b1(0) =-0.5755498075
    b1(1) = 0.6995095521
    b1(2) = 3.8925673390
    b1(3) =-17.215471648
    b1(4) = 192.67226447
    b1(5) =-161.82646165
    b1(6) =-165.20769346

    b2(0) = 0.0976883116
    b2(1) =-0.2557574982
    b2(2) =-9.1558561530
    b2(3) = 20.642075974
    b2(4) =-38.804430052
    b2(5) = 93.626774077
    b2(6) =-29.666905585

    ! calculation part

    d(1:n) = sig(1:n)*(1.0_8 - 0.12_8*dexp(-3.0_8*eps(1:n)/T_))
    m      = sum(f(1:n)*seg(1:n))
    ze0    = pi6 * rho_ * sum(f(1:n)*seg(1:n))
    ze1    = pi6 * rho_ * sum(f(1:n)*seg(1:n)*(d(1:n)**1))
    ze2    = pi6 * rho_ * sum(f(1:n)*seg(1:n)*(d(1:n)**2))
    ze3    = pi6 * rho_ * sum(f(1:n)*seg(1:n)*(d(1:n)**3))
    ze3_1  = 1.0 - ze3 ! many usage
    if (ze3_1<=0.0) stop 'error : 1 - ze3 < 0'

    giihs(1:n) = 1.0/ze3_1 + 0.5*d(1:n)*3.0*ze2/(ze3_1**2) + 0.25*(d(1:n)**2)*2.0*ze2*ze2/(ze3_1**3)
    if (any(giihs(1:n)<=0.0)) stop 'error : giihs < 0.0'

    eta = ze3

    e1  = eta
    e2  = e1 * eta
    e3  = e2 * eta
    e4  = e3 * eta
    e_1 = 1.0 - e1
    e_2 = 2.0 - e1

    c1  = 1.0 + m*2.0*e1*(4.0-e1)/e_1**4 + (1.0-m)*(20.0*e1-27.0*e2+12.0*e3-2.0*e4)/(e_1*e_2)**2
    c1  = 1.0 / c1

    c2  = -c1*c1*(m*(-4.0*e2+20.0*e1+8.0)/e_1**5 + (1-m)*(2.0*e3+12.0*e2-48.0*e1+40.0)/(e_1*e_2)**3)

    m2es3 = 0.0
    do i = 1, n; do j = 1, n ! double roop
     m2es3 = m2es3 + f(i)*f(j)*seg(i)*seg(j)*(eps(i)*eps(j))**0.5*(1.0-kij(i,j))/T_ *(0.5*(sig(i)+sig(j)))**3
    end do; end do           ! i and j

    m2e2s3 = 0.0
    do i = 1, n; do j = 1, n ! double roop
     m2e2s3 = m2e2s3 + f(i)*f(j)*seg(i)*seg(j)*(eps(i)*eps(j))*(1.0-kij(i,j))**2/T_/T_*(0.5*(sig(i)+sig(j)))**3
    end do; end do           ! i and j

    ea(0) = 1.0
    ea(1) = eta
    ea(2) = ea(1)*eta
    ea(3) = ea(2)*eta
    ea(4) = ea(3)*eta
    ea(5) = ea(4)*eta
    ea(6) = ea(5)*eta

    I1 = sum((a0 + (m-1.0)/m * a1 + (m-1.0)*(m-2.0)/m/m * a2)*ea)
    I2 = sum((b0 + (m-1.0)/m * b1 + (m-1.0)*(m-2.0)/m/m * b2)*ea)

    ahs   = 1.0/ze0*( 3.0*ze1*ze2/ze3_1 + ze2**3 / ze3 / ze3_1**2 + (ze2**3/ze3**2 - ze0)*log(ze3_1) )
    ahc   = m * ahs - sum(f(1:n)*(seg(1:n)-1.0)*log(giihs(1:n)))
    adisp = -2.0*PI*rho_*I1*m2es3 - PI*rho_*m*c1*I2*m2e2s3
    anon  = ahc + adisp

    do i = 1, n; do j = 1, n
     hm = d(i)*d(j)/(d(i)+d(j))
     gijhs(i,j) = 1/ze3_1 + hm*3.0*ze2/ze3_1**2 + hm*hm*2.0*ze2*ze2/ze3_1**3
    end do; end do
    if (any(gijhs(1:n,1:n)<0.0)) stop 'error : gijsh < 0'
    do i = 1, n; do j = 1, n
     dAB(i,j) = gijhs(i,j)*(dexp(eAB(i,j)/T_)-1.0)*(0.5*(sig(i)+sig(j)))**3*kAB(i,j)
    end do; end do

    ! matrix iteration method

    Xs  = 0.5
    del = 1.0

    do while (del > 1.0e-9)
     XsOld = Xs
     do i = 1, n
      do ii = 1, nSite     ! ii = site : only A and B type
       sum1 = 0.0
       do j = 1, n
        sum2 = 0.0
        do jj = 1, nSite
         if (ii==jj) cycle
         sum2 = sum2 + XsOld(j,jj)*dAB(i,j)
        end do
        sum1 = sum1 + f(j)*sum2
       end do
       Xs(i,ii) = 1.0 / (1.0 + rho_*sum1)
       difX(i,ii) = dabs((Xs(i,ii)-XsOld(i,ii))/Xs(i,ii))
      end do
     end do
     del = maxval(difX)
    end do
    if (any(Xs(1:n,:) < 0.0)) stop 'Xs < 0.0'

    aasc = 0.0
    do i = 1, n
     sum1 = 0.0
     do ii = 1, 2
      sum1 = sum1 + log(Xs(i,ii)) - 0.5*Xs(i,ii) + 0.5
     end do
     aasc = aasc + f(i)*sum1
    end do

    ares = anon + aasc

    ! compressibility factor

    do i = 1, n; do j = 1, n
     hm = d(i)*d(j)/(d(i)+d(j))
     gijhsd(i,j) = 1.0/rho_*(ze3/ze3_1/ze3_1 + hm*(3.0*ze2/ze3_1**2 + 6.0*ze2*ze3/ze3_1**3) + hm*hm*(4.0*ze2**2/ze3_1**3 + 6.0*ze2**2*ze3/ze3_1**4))
    end do; end do

    do i=1, n
     dgiihsd(i) = gijhsd(i,i)
    end do

    do i = 1, n; do j = 1, n
     dABd(i,j) = gijhsd(i,j)*(dexp(eAB(i,j)/T_)-1.0)*(0.5*(sig(i)+sig(j)))**3*kAB(i,j)
    end do; end do

    zhs = ze3/ze3_1 + 3*ze1*ze2/ze0/(ze3_1**2) + (ze2**3)*(3.0-ze3)/ze0/(ze3_1**3)
    zhc = m*zhs - rho_*sum(f(1:n)*(seg(1:n)-1.0)/giihs(1:n)*dgiihsd(1:n))

    forall(i=0:6)
     ind(i) = i + 1
    end forall

    deI1 = sum((a0 + (m-1.0)/m * a1 + (m-1.0)*(m-2.0)/m/m * a2)*ind*ea)
    deI2 = sum((b0 + (m-1.0)/m * b1 + (m-1.0)*(m-2.0)/m/m * b2)*ind*ea)

    zdisp = -2.0*PI*rho_*deI1*m2es3 - PI*rho_*m*(c1*deI2 + c2*eta*I2)*m2e2s3

    znon  = zhc + zdisp

    ! calculate the elements of matrix lambdapq

    Lpq = 0.0

    dp = 0; dq = 0;
    do ii = 1, nSite
     do i = 1, n
      dp = dp + 1
      do jj = 1, nSite
       do j = 1, n
        dq = dq + 1
        if (ii==jj) then
         if (i==j) then
          Lpq(dp,dq) = 1.0
         else
          Lpq(dp,dq) = 0.0
         end if
         cycle
        end if
        Lpq(dp,dq) = f(j)*rho_*dAB(i,j)*Xs(i,ii)**2 ! dAB(i,ii,j,jj) in normal
       end do
      end do
      dq = 0
     end do
    end do

    dim_ = n * 2 ! component * site number
    call inverse (Lpq, iLpq, singularQ, dim_) ! make inverse matrix

    if (singularQ) stop 'singular case orrur'

    ! calculate the elements of matrix psi of rho

    psi = 0.0
    dp = 0
    do ii = 1, nSite
     do i = 1, n
      sum1 = 0.0
      do j = 1, n
       sum2 = 0.0
       do jj = 1, nSite
        if(ii==jj) cycle
        sum2 = sum2 + Xs(j,jj)*dABd(i,j)
       end do
       sum1 = sum1 + f(j)*sum2
      end do
      dp = dp + 1
      psi(dp) = -(Xs(i,ii))**2*(rho_*sum1 + (1.0/Xs(i,ii)-1.0)/rho_)
     end do
    end do

    ! calculate (inverse lambda) * psi

    dXs = 0.0
    do dp = 1, dim_
     do dq = 1, dim_
      dXs(dp) = dXs(dp) + iLpq(dp,dq)*psi(dq)
     end do
    end do

    ! reshape Xsd(p) to Xsd(i, ii)

    dp = 0
    do ii = 1, nSite
     do i = 1, n
      dp = dp + 1
      Xsd(i,ii) = dXs(dp)
     end do
    end do

    zasc = 0.0
    do i = 1, n
     do ii = 1, nSite
      zasc = zasc + f(i)*(1.0/Xs(i,ii)-0.5)*Xsd(i,ii)
     end do
    end do

    zasc = rho_*zasc

    zres = znon + zasc
    ztot = zres + 1.0
    ptot = rho_*kb*T_*ztot

end subroutine SetNormalTerms
!----------------------------------------------------------------------------------------------------
subroutine SetdxTerms (f,rho_,T_) ! mole fraction, density, temper

    use mod_PCSAFT

    implicit none

    real(kind=8), intent(in)         :: rho_, T_
    real(kind=8), dimension(maxComp) :: f ! mol fraction
    integer :: i, j, ii, jj, c
    real(kind=8), dimension(maxComp) :: dsum     ! duumy array for 2 dimension summation
    real(kind=8), dimension(0:6)     :: dd       ! duumy array for 2 dimension summation
    real(kind=8), dimension(maxComp,0:6) :: ddx  ! duumy array for 2 dimension summation
    real(kind=8) :: hm                           ! harmonic mean
    real(kind=8) :: sum1, sum2, sum3

    ! call SetNonTerms ! call total term to use calculation for derivation of a for x
    ! assumption : used after call SetNonTerms

    do i = 0, 3
     zex(1:n,i) = pi6*rho_*seg(1:n)*d(1:n)**i
    end do

    ! component index must be first variables
    if (ze3_1<0.0) stop 'ze3 > 1.0'

    ahsx(1:n) = - 1.0/ze0/ze0 * (log(ze3_1)*(-ze0 + ze2**3/ze3**2) + 3.0*ze1*ze2/ze3_1 + ze2**3/ze3_1/ze3_1/ze3)*zex(1:n,0) + &
    &             1.0/ze0*(3.0*ze2*zex(1:n,1)/ze3_1 + 3.0*ze1*zex(1:n,2)/ze3_1 + 3.0*ze2**2*zex(1:n,2)/ze3_1**2/ze3 + &
    &             3.0*ze1*ze2*zex(1:n,3)/ze3_1**2 - (-ze0 + ze2**3/ze3**2)*zex(1:n,3)/ze3_1 - ze2**3*zex(1:n,3)/(ze3_1*ze3)**2 + &
    &             2.0*ze2**3*zex(1:n,3)/ze3_1**3/ze3 + log(ze3_1)*((-zex(1:n,0)+3.0*ze2**2*zex(1:n,2)/ze3**2)-2.0*ze2**3*zex(1:n,3)/ze3**3))

    do i = 1, n
     giihsx(1:n,i) = d(i)**2*ze2*zex(1:n,2)/ze3_1**3 + 1.5*d(i)*zex(1:n,2)/ze3_1/ze3_1 + 1.5*d(i)**2*ze2**2*zex(1:n,3)/ze3_1**4 + &
    &                3.0*d(i)*ze2*zex(1:n,3)/ze3_1**3 + zex(1:n,3)/ze3_1**2
    end do

    dsum = 0.0
    do i = 1, n
     dsum(1:n) = dsum(1:n) + f(i)*(seg(i)-1.0)/giihs(i)*giihsx(1:n,i)
    end do

    ahcx(1:n) = seg(1:n)*ahs + m*ahsx(1:n) - dsum -(seg(1:n)-1.0)*log(giihs(1:n))

    c1x(1:n) = c2*zex(1:n,3) - c1*c1*(seg(1:n)*2.0*e1*(4.0-e1)/e_1**4-seg(1:n)*(20.0*e1-27.0*e2+12.0*e3-2.0*e4)/(e_1*e_2)**2)

    dsum = 0.0
    do i = 1, n
     dsum(1:n) = dsum(1:n) + f(i)*seg(i)*(eps(1:n)*eps(i))**0.5*(1.0-kij(1:n,i))/T_*(0.5*(sig(1:n)+sig(i)))**3
    end do
    m2es3x(1:n) = 2.0*seg(1:n)*dsum(1:n)

    dsum = 0.0
    do i = 1, n
     dsum(1:n) = dsum(1:n) + f(i)*seg(i)*(eps(1:n)*eps(i))*(1.0-kij(1:n,i))**2/T_/T_*(0.5*(sig(1:n)+sig(i)))**3
    end do
    m2e2s3x(1:n) = 2.0*seg(1:n)*dsum(1:n)

    dd = a0 + (m-1.0)/m*a1 + (m-1.0)*(m-2.0)/m/m*a2
    do i = 1, n
     ddx(i,:) = seg(i)/m/m*a1 + seg(i)/m/m*(3.0-4.0/m)*a2
    end do
    I1x(1:n) = 0.0
    do i = 0, 6
     I1x(1:n) = I1x(1:n) + dd(i)*real(i)*zex(1:n,3)*ea(i)/eta + ddx(1:n,i)*ea(i)
    end do

    dd = b0 + (m-1.0)/m*b1 + (m-1.0)*(m-2.0)/m/m*b2
    do i = 1, n
     ddx(i,:) = seg(i)/m/m*b1 + seg(i)/m/m*(3.0-4.0/m)*b2
    end do
    I2x(1:n) = 0.0
    do i = 0, 6
     I2x(1:n) = I2x(1:n) + dd(i)*real(i)*zex(1:n,3)*ea(i)/eta + ddx(1:n,i)*ea(i)
    end do

    do i = 1, n; do j = 1, n
     hm = d(i)*d(j)/(d(i)+d(j))
     gijhsx(1:n,i,j) = 4.0*hm*hm*ze2*zex(1:n,2)/ze3_1**3 + 3.0*hm*zex(1:n,2)/ze3_1**2 + 6.0*hm*hm*ze2**2*zex(1:n,3)/ze3_1**4 + &
     &                 6.0*hm*ze2*zex(1:n,3)/ze3_1**3 + zex(1:n,3)/ze3_1**2
    end do; end do

    do i = 1, n; do j = 1, n
     dABx(1:n,i,j) = gijhsx(1:n,i,j)*(dexp(eAB(i,j)/T_)-1.0)*(0.5*(sig(i)+sig(j)))**3*kAB(i,j)
    end do; end do

    adispx(1:n) = - 2.0*PI*rho_*(I1x(1:n)*m2es3 + I1*m2es3x(1:n)) - PI*rho_*((seg(1:n)*c1*I2 + m*c1x(1:n)*I2 + m*c1*I2x(1:n))*m2e2s3 + m*c1*I2*m2e2s3x(1:n))

    anonx(1:n) = ahcx(1:n) + adispx(1:n)

   ! calculate the elements of matrix psi of xm for all component

   psi = 0.0
   do c = 1, n ! c : component index
    dp = 0
    do ii = 1, nSite
     do i = 1, n
      sum1 = 0.0
      do j = 1, n
       sum2 = 0.0
       sum3 = 0.0
       do jj = 1, nSite
        if (ii==jj) cycle
        sum2 = sum2 + Xs(j,jj)*dABx(c,i,j)
        if (c==j) then
         sum3 = sum3 + Xs(c,jj)*dAB(i,c)
        end if
       end do
       sum1 = sum1 + f(j)*sum2 + sum3
      end do
      dp = dp + 1
      psi(dp) = -(Xs(i,ii))**2*rho_*sum1
     end do
    end do

    ! calculate (inverse lambda) * psi for all component

    dXs = 0.0
    do dp = 1, dim_
     do dq = 1, dim_
      dXs(dp) = dXs(dp) + iLpq(dp,dq)*psi(dq)
     end do
    end do

    ! reshape Xsx(p) to Xsx(c,i, ii) for all component

    dp = 0
    do ii = 1, nSite
     do i = 1, n
      dp = dp + 1
      Xsx(c,i,ii) = dXs(dp)
     end do
    end do
   end do ! c

   aasc = 0.0
   sum1 = 0.0
   sum2 = 0.0
   sum3 = 0.0

   do c = 1, n
   sum1 = 0.0
    do i = 1, n
     sum2 = 0.0
     sum3 = 0.0
     do ii = 1, nSite
      sum2 = sum2 + (1.0/Xs(i,ii)-0.5)*Xsx(c,i,ii)
      if (i==c) then
       sum3 = sum3 + log(Xs(c,ii))-0.5*Xs(c,ii)+0.5
      end if
     end do
     sum1 = sum1 + sum2*f(i) + sum3
    end do
    aascx(c) = sum1
   end do

   aresx = anonx + aascx

   dsum = sum(f(1:n)*aresx(1:n))
   mures(1:n) = ares + (ztot-1.0) + aresx(1:n) - dsum
   phi(1:n)   = dexp(mures(1:n))/ztot

   !write(*,"('phis = ',2F30.6)") phi(1:n)

end subroutine SetdxTerms
!----------------------------------------------------------------------------------------------------
subroutine Inverse (mat, iMat, singularQ, dim)

    implicit none

    integer, intent(in):: DIM ! dimension of matrix
    integer            :: pivot, i, j
    logical :: singularQ
    real(kind=8), dimension(DIM,DIM), intent(in)  :: mat
    real(kind=8), dimension(DIM,DIM), intent(out) :: iMat
    real(kind=8), dimension(DIM,DIM) :: dMat ! dummy matrix
    real(kind=8), dimension(DIM)   :: temp, itemp
    real(kind=8) :: assumpZero, maxComp, coeff ! maxComp : max component, coeff : coefficient of elimination

    ! 초기화
    singularQ = .false.
    dMat = mat
    iMat = 0.0 ! identity matrix
    do i = 1, DIM
        iMat(i,i) = 1.0
    end do
    assumpZero = maxval(dMat) * 10e-8

    do i = 1, DIM
      ! pivoting

     pivot   = maxloc(abs(dMat(i:DIM,i)),1) + i - 1 ! pivot과 maximum value를 찾는데 내장 함수 사용
     maxComp = dMat(pivot,i)

     if (pivot /= i) then
       temp = dMat(i,:)           ! 행 교환
      itemp = iMat(i,:)
      dMat(i,:) = dMat(pivot,:)
      iMat(i,:) = iMat(pivot,:)
      dMat(pivot,:) =  temp
      iMat(pivot,:) = itemp
     end if

     if (abs(maxComp) < assumpZero) then;
      singularQ = .true.; iMat = 0.0_8; return
     end if

     do j = i+1, DIM    ! 가우스 소거법 적용
      coeff     = dMat(j,i)
      dMat(j,:) = dMat(j,:) - dMat(i,:) * coeff/maxComp
      iMat(j,:) = iMat(j,:) - iMat(i,:) * coeff/maxComp
     end do
    end do

    do i = DIM, 1, -1      ! 후진 대입법으로 계산
     coeff     = dMat(i,i)

     dMat(i,:) = dMat(i,:) / coeff
     iMat(i,:) = iMat(i,:) / coeff

     do j = i-1, 1, -1
      coeff     = dMat(j,i)

      dMat(j,:) = dMat(j,:) - coeff*dMat(i,:)
      iMat(j,:) = iMat(j,:) - coeff*iMat(i,:)
     end do
    end do

end subroutine Inverse
!----------------------------------------------------------------------------------------------------
subroutine FindMinMax (dmax_, pmax_, dmin_, pmin_, x_, T_, critiQ_)

    use mod_PCSAFT

    implicit none

    integer :: i
    integer, parameter :: maxi = 200
    logical :: miniQ
    logical, intent(out) :: critiQ_ ! if super critical, critiQ = true

    real(kind=8), intent(out) :: dmax_, pmax_, dmin_, pmin_
    real(kind=8), dimension(maxComp), intent(in) :: x_
    real(kind=8), intent(in) :: T_
    real(kind=8) :: rhomax, drho, pnew, pold, rhon, rho_
    real(kind=8) :: dold, dnew, ddp, ddrho

    real(kind=8) :: Pcal, dPcal

    dmax = 0.0
    pmax = 0.0
    dmin = 0.0
    pmin = 0.0
    critiQ_ = .false.

    d(1:n) = sig(1:n)*(1.0_8 - 0.12_8*dexp(-3.0_8*eps(1:n)/T_))
    rhomax = 0.8*0.7405/(pi6*sum(x_(1:n)*seg(1:n)*d(1:n)**3))
    call SetNormalTerms (x_, rhomax, T_)
    drho   = rhomax/real(maxi)
    rho_   = 1.0e-8
    pold   = Pcal(x_,rho_,T_)
    miniQ  = .true.
    !call system('pause')

rough : &
    do i = 1, maxi
     rhon = rho_ + drho
     pnew = Pcal(x_,rho_,T_)

     if (miniQ) then
      if (pnew < pold) then
       pmax_ = pnew
       dmax_ = rhon
       miniQ = .false.
       !write(*,*) 'me1'
      end if
     else
      if (pnew > pold) then
       pmin_ = pnew
       dmin_ = rhon
       !write(*,*) 'me2'
       exit rough
      end if
     end if

     pold = pnew
     rho_ = rhon
    end do rough

    if (i-1==maxi) then
     critiQ_ = .true.
     return
    end if

    ddrho = 1.0e-7
    dold  = dmax

    ddp   = (dPcal (x_,dold+ddrho,T_) - dPcal (x_,dold,T_))/ddrho

    do i = 1, maxi
     dnew = dold - dPcal(x_,dold,T_)/ddp
     if (dnew<0.0) exit
     if (dabs((dnew-dold)/dnew)<10e-4) exit
     dold = dnew
     ddp  = (dPcal (x_,dold+ddrho,T_) - dPcal (x_,dold,T_))/ddrho
    end do

    dmax = dnew
    pmax = Pcal (x_,dmax,T_)

    dold  = dmin
    ddp   = (dPcal (x_,dold+ddrho,T_) - dPcal (x_,dold,T_))/ddrho
    do i = 1, maxi
     dnew = dold - dPcal(x_,dold,T_)/ddp
     if (dnew<0.0) exit
     if (dabs((dnew-dold)/dnew)<10e-4) exit
     dold = dnew
     ddp   = (dPcal (x_,dold+ddrho,T_) - dPcal (x_,dold,T_))/ddrho
    end do
    dmin = dnew
    pmin = Pcal (x_,dmin,T_)

end subroutine FindMinMax
!----------------------------------------------------------------------------------------------------
subroutine PrintPeos (filename_, x_, T_)

    use mod_PCSAFT

    implicit none

    integer :: i
    integer, parameter :: ut = 10
    character(len=61), intent(in) :: filename_
    real(kind=8), intent(in) :: T_
    real(kind=8), dimension(maxComp), intent(in) :: x_
    real(kind=8) :: rhomax, drho, pnew, pold, rhon, rho_

    real(kind=8) :: Pcal

    d(1:n) = sig(1:n)*(1.0_8 - 0.12_8*dexp(-3.0_8*eps(1:n)/T_))
    rhomax = 0.8*0.7405/(pi6*sum(x_(1:n)*seg(1:n)*d(1:n)**3))
    drho   = rhomax/200.0
    rho_   = 1.0e-8
    open(ut, file=filename_)
    write(*,"(60('-'))")
    write(ut,"(60('-'))")
    write(*,"(2A30)")  'dens<N/A^3>', 'pressure<MPa>'
    write(ut,"(2A30)") 'dens<N/A^3>', 'pressure<MPa>'
    write(*,"(60('-'))")
    write(ut,"(60('-'))")
    do i = 1, 200
     write(*,"(2F30.10)")  rho_, Pcal (x_,rho_,T_)
     write(ut,"(2F30.10)") rho_, Pcal (x_,rho_,T_)
     rho_ = rho_ + drho
    end do
    close(ut)

end subroutine PrintPeos
!----------------------------------------------------------------------------------------------------
function rhoCal (x_, T_, P_, rhoi)

    use mod_PCSAFT

    implicit none

    integer :: i
    real(kind=8), intent(in) :: T_, P_, rhoi ! rhoi = initial rho
    real(kind=8), dimension(maxComp) :: x_
    real(kind=8) :: drho, dpdr, rhoo, rhon, ptmp ! dpdr = d P / d rho

    real(kind=8) :: rhoCal
    real(kind=8) :: Pcal

    drho = 1.0e-10
    rhoo = rhoi

    dpdr = (Pcal(x_, rhoo+drho, T_) - Pcal(x_, rhoo, T_))/drho
    ptmp = Pcal (x_, rhoo, T_)

    do i = 1, 100
     rhon = rhoo - (ptmp - P_)/dpdr
     if (dabs((rhon-rhoo)/rhon)<10e-10) exit
     rhoo = rhon
     dpdr = (Pcal(x_, rhoo+drho, T_) - Pcal(x_, rhoo, T_))/drho
     ptmp = Pcal (x_, rhoo, T_)

    end do

    rhoCal = rhon

    rhoo = rhon
    dpdr = (Pcal(x_, rhoo+drho, T_) - Pcal(x_, rhoo, T_))/drho
    if (dpdr<0.0) stop 'error in cal rho : center root'
    if (i==100)   stop 'error in cal rho : diverge'

end function rhoCal
!----------------------------------------------------------------------------------------------------
function Pcal (x_,rho_,T_)

    use mod_PCSAFT

    implicit none

    real(kind=8), intent(in) :: T_, rho_
    real(kind=8), dimension(maxComp) :: x_

    real(kind=8) :: Pcal

    call SetNormalTerms (x_,rho_,T_)

    Pcal = ptot

end function Pcal
!----------------------------------------------------------------------------------------------------
function dPcal (x_,rho_,T_)

    use mod_PCSAFT

    implicit none

    real(kind=8), intent(in) :: T_, rho_
    real(kind=8), dimension(maxComp) :: x_
    real(Kind=8) :: drho

    real(kind=8) :: dPcal

    drho = 10e-10
    call SetNormalTerms (x_,rho_,T_)
    dPcal = ptot
    call SetNormalTerms (x_,rho_+drho,T_)
    dPcal = ptot - dPcal
    dPcal = dPcal / drho

end function dPcal
!----------------------------------------------------------------------------------------------------
subroutine Getphi (x_,rho_,T_)

    use mod_PCSAFT

    implicit none

    real(kind=8), intent(in) :: T_, rho_
    real(kind=8), dimension(maxComp), intent(in) :: x_

    call SetNormalTerms (x_,rho_,T_)
    call SetdxTerms (x_,rho_,T_)

end subroutine Getphi
!----------------------------------------------------------------------------------------------------
function Psat(c,T_)

    use mod_PCSAFT

    implicit none

    integer :: i
    integer, parameter :: imax = 100
    integer,      intent(in) :: c   ! component index
    real(kind=8), intent(in) :: T_  ! temperature
    real(kind=8), dimension(maxComp) :: tf ! if (i==c) tf = 1  else tf = 0 True or False
    real(kind=8) :: rhomax, drho, pnew, pold, rhon, rhoi, rhov, rhol
    real(kind=8) :: phiv, phil

    real(kind=8) :: Psat
    real(kind=8) :: rhoCal, Pcal

    tf = 0.0
    do i = 1, n
     if (i==c) tf(i) = 1.0
    end do; i = 0

    d(1:n) = sig(1:n)*(1.0_8 - 0.12_8*dexp(-3.0_8*eps(1:n)/T_))
    rhomax = 0.8*0.7405/(pi6*sum(tf(1:n)*seg(1:n)*d(1:n)**3))

    call FindMinMax (dmax, pmax, dmin, pmin, tf, T_, criQ)

    if (criQ) then ! super critical
     Psat = Pcal (tf, 0.2*rhomax, T_)
     return
    end if

    rhol = 1.1*dmin

    if (pmin<0) then
     pold = 0.5*pmax
    else
     pold = 0.5*(pmax+pmin)
    end if

    do i = 1, imax
     rhoi = pold / T_ / kb ! assump ideal gas
     rhov = rhoCal (tf, T_, pold, rhoi)
     rhoi = rhol
     rhol = rhoCal (tf, T_, pold, rhoi)
     if (dabs((rhov-rhol)/rhol)<1.0e-9) stop 'rhov = rhol'
     if (rhov>rhol) stop 'density fine filed'

     call Getphi(tf, rhov, T_)
     phiv = phi(c)
     call Getphi(tf, rhol, T_)
     phil = phi(c)

     pnew = pold * phil / phiv

     if (dabs((pnew-pold)/pnew)<1.0e-5) exit

     pold = pnew
    end do

    Psat = pnew
    srhol(c) = rhol
    srhov(c) = rhov

    if (i==imax) stop 'error in Past'

end function Psat
!----------------------------------------------------------------------------------------------------
subroutine InvPsat (c,Pin,Tout)

    use mod_PCSAFT

    implicit none

    integer, parameter :: imax = 100
    integer :: i
    integer,      intent(in)  :: c   ! component index
    real(kind=8), intent(in)  :: Pin
    real(kind=8), intent(out) :: Tout
    real(kind=8), dimension(maxComp) :: tf ! if (i==c) tf = 1  else tf = 0 True or False
    real(kind=8) :: Ttmp, Tnew, dPdT, delT, Ptmp

    real(kind=8) :: Psat

    tf = 0.0
    do i = 1, n
     if (i==c) tf(i) = 1.0
    end do; i = 0

    delT = 0.0001

    Tout = 0.0
    Ttmp = 300
    call FindMinMax (dmax, pmax, dmin, pmin, tf, Ttmp, criQ)

    do while (criQ)
      call FindMinMax (dmax, pmax, dmin, pmin, tf, Ttmp, criQ)
      !write(*,*) Ttmp
      Ttmp = Ttmp - 20.0
    end do

    do i = 1, imax
      Ptmp = Psat (c,Ttmp)
      dPdT = (Psat (c,Ttmp+delT) - Ptmp)/delT
      Tnew = Ttmp - (Ptmp - Pin)/dPdT
      if (dabs((Tnew-Ttmp)/Tnew)<1.0e-5) exit
      Ttmp = Tnew
    end do

    Tout = Tnew

    if (i==imax) Tout = 0.0

end subroutine InvPsat
!----------------------------------------------------------------------------------------------------
subroutine PsatAuto (filename, c)

    use mod_PCSAFT

    implicit none

    integer :: i, iT, niT
    integer, parameter :: imax = 100
    integer,      intent(in)      :: c   ! component index
    character(len=61), intent(in) :: filename

    real(kind=8), dimension(maxComp) :: tf ! if (i==c) tf = 1  else tf = 0 True or False
    real(kind=8) :: rhomax, drho, pnew, pold, rhon, rhoi, rhov, rhol
    real(kind=8) :: phiv, phil, tmp

    real(kind=8) :: rhoCal, Pcal

    niT = int((T1-T0)/dT) + 1

    tf = 0.0
    do i = 1, n
     if (i==c) tf(i) = 1.0
    end do; i = 0

    rhomax = 0.8*0.7405/(pi6*sum(tf(1:n)*seg(1:n)*d(1:n)**3))
    call FindMinMax (dmax, pmax, dmin, pmin, tf, T0, criQ)

    rhol = 1.1*dmin

    if (pmin<0) then
     pold = 0.5*pmax
    else
     pold = 0.5*(pmax+pmin)
    end if

    open(1,file=filename)
     write(*,"(100('-'))")
     write(1,"(100('-'))")
     write(*,"(5A20)") 'temper<K>', 'rhov<N/A^3>', 'rhol<N/A^3>', 'pv<MPa>', 'pl<MPa>'
     write(1,"(5A20)") 'temper<K>', 'rhov<N/A^3>', 'rhol<N/A^3>', 'pv<MPa>', 'pl<MPa>'
     write(*,"(100('-'))")
     write(1,"(100('-'))")
    tmp = T0
    do while (.true.)
     if (tmp>T1) exit
     in : do i = 1, imax
       rhoi = pold / tmp / kb ! assump ideal gas
       rhov = rhoCal (tf, tmp, pold, rhoi)
       rhoi = rhol
       rhol = rhoCal (tf, tmp, pold, rhoi)
       if (dabs((rhov-rhol)/rhol)<1.0e-9) stop 'rhov = rhol'
       if (rhov>rhol) stop 'density fine filed'

       call Getphi(tf, rhov, tmp)
       phiv = phi(c)
       call Getphi(tf, rhol, tmp)
       phil = phi(c)

       pnew = pold * phil / phiv

       if (dabs((pnew-pold)/pnew)<1.0e-5) exit in
       pold = pnew
     end do in
     write(*,"(F20.5, 4F20.10)") tmp, rhov, rhol, Pcal(tf,rhov,tmp), Pcal(tf,rhol,tmp)
     write(1,"(F20.5, 4F20.10)") tmp, rhov, rhol, Pcal(tf,rhov,tmp), Pcal(tf,rhol,tmp)
     tmp = tmp + dT
     pold= pnew
    end do

    close(1)

end subroutine PsatAuto
!----------------------------------------------------------------------------------------------------
subroutine Pxy (filename) ! only for binary

    use mod_PCSAFT

    implicit none

    integer :: i
    logical :: sumQ, fracQ, iniQ
    character(len=61), intent(in) :: filename
    real(kind=8), dimension(maxComp) :: dpsat, ynew, xnew
    real(kind=8) :: pold, pnew, rhoi, rhov, rhol, maxv, sumy

    real(kind=8) :: rhoCal, Psat

    iniQ = .true.

    x(1) = x0
    x(2) = 1.0 - x(1) ! binary assumption

    dpsat = 0.0
    do i = 1, n
     dpsat(i) = Psat(i,T)  ! T = system temperature
    end do

    write(*,"('Pi =', 11F10.5)") dpsat(1:n)

    pold = sum(x(1:n)*dpsat(1:n))
    Ke(1:n) = dpsat(1:n)/pold

    call FindMinMax (dmax, pmax, dmin, pmin, x, T, criQ)

    if (pmin<0) then
     pold = 0.5*pmax
    else
     pold = 0.5*(pmax+pmin)
    end if
    rhol = 1.1*dmin
    write(*,"('temper = ', F10.3)") T
    open(1,file=filename)
    write(*,"(60('-'))")
    write(1,"(60('-'))")
    write(*,"(5A20)") 'x1', 'y1', 'p<MPa>'
    write(1,"(5A20)") 'x1', 'y1', 'p<MPa>'
    write(*,"(60('-'))")
    write(1,"(60('-'))")
    do while (.true.)
      if (dx>0.0) then
        if (x(1) > x1+10e-10) exit
      else
        if (x(1) < x1-10e-10) exit
      end if
      sumQ   = .false.
      y(1:n) = Ke(1:n)*x(1:n)
      do while (.not.sumQ)
        rhoi = rhol
        rhol = rhoCal (x, T, pold, rhoi)
        call Getphi (x, rhol, T)
        fl(1:n) = x(1:n)*phi(1:n)
        fracQ = .false.
        do while (.not.fracQ)
          rhoi = pold/kb/T
          rhov = rhoCal (y, T, pold, rhoi)
          call Getphi (y, rhov, T)
          fv(1:n) = y(1:n)*phi(1:n)
          ynew(1:n) = y(1:n)*fl(1:n)/fv(1:n)
          if (iniQ) then ! this logic is necessary only in the beginning
            where (ynew(1:n)>1.0) ynew = ynew/sum(ynew(1:n))
          end if
          maxv = maxval(dabs((ynew(1:n)-y(1:n))/ynew(1:n)))
          if (maxv < 1.0e-5) fracQ = .true.
          y(1:n) = ynew(1:n)
          !if (any(y(1:n)<0.0)) stop 'failed'
        end do
        sumy = sum(y(1:n))
        if (dabs(1.0-sumy)<1.0e-5) sumQ = .true.
        pold = pold * sumy
      end do
      iniQ = .false. ! check the beginning
      !write(*,"(3F20.6)") x(1), y(1), pold
      write(1,"(3F20.6)") x(1), y(1), pold
      Ke(1:n) = y(1:n)/x(1:n)
      x(1) = x(1) + dx
      x(2) = 1.0 - x(1)
    end do

    close(1)

end subroutine Pxy
!----------------------------------------------------------------------------------------------------
subroutine Txy (filename) ! only for binary

    use mod_PCSAFT

    implicit none

    integer :: i, j, c
    logical :: sumQ, fracQ, overQ, iniQ
    character(len=61), intent(in) :: filename
    real(kind=8), dimension(maxComp) :: dpsat, ynew, xnew, dtold
    real(kind=8) :: Told, Tnew, rhoi, rhov, rhol, maxv, sumy, pguess, Prto, pold
    real(kind=8), dimension(maxComp) :: tf ! if (i==c) tf = 1  else tf = 0 True or False
    real(kind=8) :: delT

    real(kind=8) :: rhoCal, Psat, Pcal

    iniQ = .true.

    x(1) = x0
    x(2) = 1.0 - x(1) ! binary assumption

    do c = 1, n
      call InvPsat (c, P, Told)
      dtold(c) = Told
    end do

    Told = 0.0
    do i = 1, n
      do j = 1, n
        Told = Told + x(i)*x(j)*(dtold(i)*dtold(j))**0.5
      end do
    end do
    delT = Told / 10.0

    dpsat = 0.0
    do i = 1, n
     dpsat(i) = Psat(i,Told)  ! T = system temperature
    end do

    pguess  = sum(x(1:n)*dpsat(1:n))
    Ke(1:n) = dpsat(1:n)/pguess
    y(1:n)  = Ke(1:n)*x(1:n)
    overQ   = any(y(1:n) > 1.0001)
    if (overQ) stop 'some yi > 1.0'
    call FindMinMax (dmax, pmax, dmin, pmin, x, Told, criQ)
    if (criQ) stop 'super critical'

    rhol = 1.1*dmin
    write(*,"('p = ', F10.3)") P
    open(1,file=filename)
    write(*,"(60('-'))")
    write(1,"(60('-'))")
    write(*,"(5A20)") 'x1', 'y1', 'T<K>'
    write(1,"(5A20)") 'x1', 'y1', 'T<K>'
    write(*,"(60('-'))")
    write(1,"(60('-'))")
    do while (.true.)
      if (dx>0.0) then
        if (x(1) > x1+10e-10) exit
      else
        if (x(1) < x1-10e-10) exit
      end if
      sumQ  = .false.
      y(1:n) = Ke(1:n)*x(1:n)
      do while (.not.sumQ)
        rhoi = rhol
        rhol = rhoCal (x, Told, P, rhoi)
        call Getphi (x, rhol, Told)
        fl(1:n) = x(1:n)*phi(1:n)
        fracQ = .false.
        do while (.not.fracQ)
          rhoi = P/kb/Told
          rhov = rhoCal (y, Told, P, rhoi)
          call Getphi (y, rhov, Told)
          fv(1:n) = y(1:n)*phi(1:n)
          ynew(1:n) = y(1:n)*fl(1:n)/fv(1:n)
          if (iniQ) then ! this logic is necessary only in the beginning
            where (ynew(1:n)>1.0) ynew = ynew/sum(ynew(1:n))
          end if
          maxv = maxval(dabs((ynew(1:n)-y(1:n))/ynew(1:n)))
          if (maxv < 1.0e-5) fracQ = .true.
          y(1:n) = ynew(1:n)
        end do
        sumy = sum(y(1:n))
        if (dabs(1.0-sumy)<1.0e-5) sumQ = .true.
        if (sumy < 0.0) stop 'error sumy < 0.0'
        Told = Told -delT*log(sumy)
      end do
      iniQ = .false. ! check the first time
      write(*,"(3F20.6)") x(1), y(1), Told
      write(1,"(3F20.6)") x(1), y(1), Told
      Ke(1:n) = y(1:n)/x(1:n)
      x(1) = x(1) + dx
      x(2) = 1.0 - x(1)
    end do

    close(1)

end subroutine Txy
!----------------------------------------------------------------------------------------------------
subroutine dewP   ! 초기값 생성 알고리즘 불안정.

    use mod_PCSAFT

    implicit none

    integer :: i, j
    logical :: sumQ, fracQ, iniQ
    !character(len=61), intent(in) :: filename
    real(kind=8), dimension(maxComp) :: dpsat, ynew, xnew, phiv, phil, tf, phip ! phi pure
    real(kind=8) :: pold, pnew, rhoi, rhov, rhol, maxv, sumx

    real(kind=8) :: rhoCal, Psat, Pcal

    iniQ = .true.

    y = z ! assumption for calculate dew P

    dpsat = 0.0
    do i = 1, n
     dpsat(i) = Psat(i,T)  ! and get density
    end do

    pold    = 1.0 / sum(y(1:n)/dpsat(1:n))
    Ke(1:n) = dpsat(1:n)/pold
    x(1:n)  = y(1:n)/Ke(1:n)
    !  write(*,"(11F10.6)") x(1:n)

    call FindMinMax (dmax, pmax, dmin, pmin, y, T, criQ)

    if (pmin<0) then
      if (pold < pmax) pold = 0.5*pmax
    else
      if (pmin < pold .and. pold < pmax) pold = 0.5*(pmax+pmin)
    end if
    call FindMinMax (dmax, pmax, dmin, pmin, x, T, criQ)
    rhol = 1.1*dmin

    sumQ   = .false.
    do while (.not.sumQ)
      rhoi = pold/kb/T
      rhov = rhoCal (y, T, pold, rhoi)
      call Getphi (y, rhov, T)
      fv(1:n) = y(1:n)*phi(1:n)
      fracQ = .false.
      do while (.not.fracQ)
        !write(*,"('rho l = ', 2F20.6)") rhol, Pcal (x,rhol,T)
        rhoi = rhol
        rhol = rhoCal (x, T, pold, rhoi)
        call Getphi (x, rhol, T)
        fl(1:n) = x(1:n)*phi(1:n)
        xnew(1:n) = x(1:n)*fv(1:n)/fl(1:n)
        !write(*,"('fv/fl = ',11F10.6)") fv(1:n)/fl(1:n), sum(xnew(1:n))
        !write(*,"(12F10.6)") x(1:n), sum(x(1:n))
        if (iniQ) then ! this logic is necessary only in the beginning
          where (xnew(1:n)>1.0) xnew = xnew/sum(xnew(1:n))
        end if
        maxv = maxval(dabs((xnew(1:n)-x(1:n))/xnew(1:n)))
        if (maxv < 1.0e-5) fracQ = .true.
        x(1:n) = xnew(1:n)
        !write(*,"(12F10.6)") x(1:n), sum(x(1:n))
        !if (any(x(1:n)<0.0)) stop 'failed'
      end do
      sumx = sum(x(1:n))
      if (dabs(1.0-sumx)<1.0e-5) sumQ = .true.
      pold = pold / sumx
      !write(*,"(F10.6)") pold
    end do

    write(*,"('xi     =',11F10.6)") x(1:n)
    write(*,"('dew P = ', F10.6)") pold

end subroutine dewP
!----------------------------------------------------------------------------------------------------
subroutine bublP

    use mod_PCSAFT

    implicit none

    integer :: i
    logical :: sumQ, fracQ, iniQ
    !character(len=61), intent(in) :: filename
    real(kind=8), dimension(maxComp) :: dpsat, ynew, xnew
    real(kind=8) :: pold, pnew, rhoi, rhov, rhol, maxv, sumy

    real(kind=8) :: rhoCal, Psat

    iniQ = .true.

    x = z ! assumption for calculate buble P

    dpsat = 0.0
    do i = 1, n
     dpsat(i) = Psat(i,T)  ! T = system temperature
    end do

    !write(*,"('Pi =', 11F10.5)") dpsat(1:n)

    pold = sum(x(1:n)*dpsat(1:n))
    Ke(1:n) = dpsat(1:n)/pold

    call FindMinMax (dmax, pmax, dmin, pmin, x, T, criQ)

    if (pmin<0) then
     pold = 0.5*pmax
    else
     pold = 0.5*(pmax+pmin)
    end if
    rhol = 1.1*dmin

    sumQ   = .false.
    y(1:n) = Ke(1:n)*x(1:n)
    do while (.not.sumQ)
      rhoi = rhol
      rhol = rhoCal (x, T, pold, rhoi)
      call Getphi (x, rhol, T)
      fl(1:n) = x(1:n)*phi(1:n)
      fracQ = .false.
      do while (.not.fracQ)
        rhoi = pold/kb/T
        rhov = rhoCal (y, T, pold, rhoi)
        call Getphi (y, rhov, T)
        fv(1:n) = y(1:n)*phi(1:n)
        ynew(1:n) = y(1:n)*fl(1:n)/fv(1:n)
        if (iniQ) then ! this logic is necessary only in the beginning
          where (ynew(1:n)>1.0) ynew = ynew/sum(ynew(1:n))
        end if
        maxv = maxval(dabs((ynew(1:n)-y(1:n))/ynew(1:n)))
        if (maxv < 1.0e-5) fracQ = .true.
        y(1:n) = ynew(1:n)
        !if (any(y(1:n)<0.0)) stop 'failed'
      end do
      sumy = sum(y(1:n))
      if (dabs(1.0-sumy)<1.0e-5) sumQ = .true.
      pold = pold * sumy
    end do

    write(*,"('yi     =',11F10.6)") y(1:n)
    write(*,"('bubl P = ', F10.6)") pold

end subroutine bublP
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
subroutine OBF (datafile, dkij, error) ! only for binary

    use mod_PCSAFT

    implicit none

    logical :: phaseI, phaseII
    integer :: i, ndata, temp, mid
    logical :: sumQ, fracQ, iniQ, middleQ
    character(len=61), intent(in)  :: datafile
    real(kind=8)     , intent(in)  :: dkij
    real(kind=8)     , intent(out) :: error

    real(kind=8), dimension(100)   :: pexp, texp, xexp, yexp
    real(kind=8), dimension(maxComp) :: dpsat, ynew, xnew
    real(kind=8) :: pold, pnew, rhoi, rhov, rhol, maxv, sumy, Told, midx, midy, sumx

    real(kind=8) :: rhoCal, Psat

    phaseI  = .true.
    phaseII = .false.

    midx  = 0.001
    midy  = 0.001

    error = 0.0

    ndata = 1

    kij = 0.0
    kij(1,2) = dkij
    kij(2,1) = dkij ! off diagonal element

    open(10, file=datafile)

      do while(.true.)
        read(10, *, err=1,end=2) pexp(ndata), texp(ndata), xexp(ndata), yexp(ndata)
        ndata = ndata + 1
      end do
1     continue
2     continue

      ndata = ndata - 1

    close(10)

    if (phaseI) then

    !============== phase 1 calculation (x, P) start ========================

    write(*,"('phase I')")
    write(*,"(7A15)") 'Pexp' ,'Pcal', 'T', 'x1', 'dy1', 'delt error', 'error'
    do temp = 1, ndata

      middleQ = .false.

      if ( dabs(Told - texp(temp)) > 10e-4 ) then

        iniQ    = .true.


        Told = texp(temp)

        x(1) = xexp(temp)
        x(2) = 1.0 - x(1)

        dpsat = 0.0
        do i = 1, n
          dpsat(i) = Psat(i,texp(temp))  ! T = system temperature
        end do

        pold = sum(x(1:n)*dpsat(1:n))
        Ke(1:n) = dpsat(1:n)/pold

        call FindMinMax (dmax, pmax, dmin, pmin, x, texp(temp), criQ)

        if (pmin<0) then
          pold = 0.5*pmax
        else
          pold = 0.5*(pmax+pmin)
        end if
        rhol = 1.1*dmin

      end if
      !write(*,"('temper = ', F10.3)") T

      x(1) = xexp(temp)
      x(2) = 1.0 - x(1)

3     continue

      sumQ   = .false.
      y(1:n) = Ke(1:n)*x(1:n)
      do while (.not.sumQ)
        rhoi = rhol
        rhol = rhoCal (x, Texp(temp), pold, rhoi)
        call Getphi (x, rhol, Texp(temp))
        fl(1:n) = x(1:n)*phi(1:n)
        fracQ = .false.
        do while (.not.fracQ)
          rhoi = pold/kb/Texp(temp)
          rhov = rhoCal (y, Texp(temp), pold, rhoi)
          call Getphi (y, rhov, Texp(temp))
          fv(1:n) = y(1:n)*phi(1:n)
          ynew(1:n) = y(1:n)*fl(1:n)/fv(1:n)
          if (iniQ) then ! this logic is necessary only in the beginning
            where (ynew(1:n)>1.0) ynew = ynew/sum(ynew(1:n))
          end if
          maxv = maxval(dabs((ynew(1:n)-y(1:n))/ynew(1:n)))
          if (maxv < 1.0e-5) fracQ = .true.
          y(1:n) = ynew(1:n)
          !if (any(y(1:n)<0.0)) stop 'failed'
        end do
        sumy = sum(y(1:n))
        if (dabs(1.0-sumy)<1.0e-5) sumQ = .true.
        pold = pold * sumy
      end do
      iniQ = .false. ! check the beginning
      Ke(1:n) = y(1:n)/x(1:n)
      ! cal P & y1

      if ( .not.middleQ ) then

        !error = error + abs((pold-pexp(temp))/pexp(temp))  + 3.0*abs((y(1)-yexp(temp))/yexp(temp))
        error = error + abs((pold-pexp(temp))/pexp(temp))
        !error = error + abs((y(1)-yexp(temp))/yexp(temp))

        write(*, "(7F15.6)") pexp(temp), pold, texp(temp), xexp(temp), &
        & dabs((y(1)-yexp(temp))/yexp(temp)), abs((pold-pexp(temp))/pexp(temp)), error

      end if
      middleQ = .true.
      if ( abs(texp(temp)-texp(temp+1)) < 10e-4 ) then ! same temper data
        if ( xexp(temp+1) > xexp(temp) ) then
          if ( xexp(temp+1) > x(1) + midx) then !
            x(1) = x(1) + midx
            x(2) = 1.0  - x(1)
            !write(*,"('middle')")
            goto 3
          end  if
        else if ( xexp(temp+1) < xexp(temp) ) then
          if ( xexp(temp+1) < x(1) - midx) then
            x(1) = x(1) - midx
            x(2) = 1.0  - x(1)
            !write(*,"('(x,y) = ', 2F10.3, 'x(i+1), x(i) = ', 2F10.3)") x(1), y(1), xexp(temp+1), xexp(temp)
            goto 3
          end  if
        end if
      end if

    end do

    !============== phase 1 calculation (x, P) end ========================
    end if

    if (phaseII) then

    !============== phase 2 calculation (x, P) start ======================

    write(*,"('phase II')")
    write(*,"(7A15)") 'Pexp' ,'Pcal', 'T', 'dx1', 'y1', 'delt error', 'error'
    do temp = 1, ndata

      middleQ = .false.

      if ( dabs(Told - texp(temp)) > 10e-4 ) then

        iniQ    = .true.

        Told = texp(temp)

        y(1) = yexp(temp)
        y(2) = 1.0 - y(1)

        dpsat = 0.0
        do i = 1, n
          dpsat(i) = Psat(i,Told)  ! and get density
        end do

        pold    = 1.0 / sum(y(1:n)/dpsat(1:n))
        Ke(1:n) = dpsat(1:n)/pold
        x(1:n)  = y(1:n)/Ke(1:n)

        call FindMinMax (dmax, pmax, dmin, pmin, y, texp(temp), criQ)

        if (pmin<0) then
          if (pold < pmax) pold = 0.5*pmax
        else
          if (pmin < pold .and. pold < pmax) pold = 0.5*(pmax+pmin)
        end if
        call FindMinMax (dmax, pmax, dmin, pmin, x, texp(temp), criQ)
        rhol = 1.1*dmin

      end if
      !write(*,"('temper = ', F10.3)") T

      y(1) = yexp(temp)
      y(2) = 1.0 - y(1)

4     continue

      sumQ   = .false.
      x(1:n) = y(1:n)/Ke(1:n)
        do while (.not.sumQ)
          rhoi = pold/kb/texp(temp)
          rhov = rhoCal (y, texp(temp), pold, rhoi)
          call Getphi (y, rhov, texp(temp))
          fv(1:n) = y(1:n)*phi(1:n)
          fracQ = .false.
          do while (.not.fracQ)
            rhoi = rhol
            rhol = rhoCal (x, texp(temp), pold, rhoi)
            call Getphi (x, rhol, texp(temp))
            fl(1:n) = x(1:n)*phi(1:n)
            xnew(1:n) = x(1:n)*fv(1:n)/fl(1:n)
            if (iniQ) then ! this logic is necessary only in the beginning
              where (xnew(1:n)>1.0) xnew = xnew/sum(xnew(1:n))
            end if
            maxv = maxval(dabs((xnew(1:n)-x(1:n))/xnew(1:n)))
            if (maxv < 1.0e-5) fracQ = .true.
            x(1:n) = xnew(1:n)
          end do
          sumx = sum(x(1:n))
          if (dabs(1.0-sumx)<1.0e-5) sumQ = .true.
          pold = pold / sumx
        end do
      iniQ = .false. ! check the beginning
      Ke(1:n) = y(1:n)/x(1:n)

      if ( .not.middleQ ) then

        !error = error + abs((pold-pexp(temp))/pexp(temp))  + 3.0*abs((x(1)-xexp(temp))/xexp(temp))
        error = error + abs((pold-pexp(temp))/pexp(temp))
        !error = error + abs((x(1)-xexp(temp))/xexp(temp))

        write(*, "(7F15.6)") pexp(temp), pold, texp(temp), abs(x(1)-xexp(temp))/xexp(temp), &
        & yexp(temp), abs((pold-pexp(temp))/pexp(temp)), error

      end if
      middleQ = .true.
      if ( abs(texp(temp)-texp(temp+1)) < 10e-4 ) then ! same temper data
        if ( yexp(temp+1) > yexp(temp) ) then
          if ( yexp(temp+1) > y(1) + midy) then !
            y(1) = y(1) + midy
            y(2) = 1.0  - y(1)
            !write(*,"('(x,y) = ', 2F10.3, 'y(i+1), y(i) = ', 2F10.3)") x(1), y(1), yexp(temp+1), yexp(temp)
            goto 4
          end  if
        else if ( yexp(temp+1) < yexp(temp) ) then
          if ( yexp(temp+1) < y(1) - midy) then
            y(1) = y(1) - midy
            y(2) = 1.0  - y(1)
            !write(*,"('(x,y) = ', 2F10.3, 'y(i+1), y(i) = ', 2F10.3)") x(1), y(1), yexp(temp+1), yexp(temp)
            goto 4
          end  if
        end if
      end if

    end do

    !============== phase 2 calculation (y, P) end ======================
    end if
    error = error / ndata * 100.0

end subroutine OBF
!----------------------------------------------------------------------------------------------------
subroutine FindBinary

    implicit none

    real(kind=8) :: dkij, dkijn, error, errorn, delk, direc

    !read(*,*) dkij

    delk  = 0.0001

    dkij = 0.00
    call OBF('expdata.txt', dkij, error)
    write(*,"('error = ', F10.3)") error

    dkijn = dkij + delk

    call OBF('expdata.txt', dkijn, errorn)
    write(*,"('kij =  0.001   error = ', F10.3)") errorn

    if ( errorn < error ) then
      direc = 1.0
      write(*,"('positive direction')")
    end if

    dkijn = dkij - delk

    !call OBF('expdata.txt', dkijn, errorn)
    write(*,"('kij = -0.001   error = ', F10.3)") errorn

    if ( errorn < error ) then
      direc = -1.0
      write(*,"('negative direction')")
    end if

    !direc = 1.0
    !dkijn = 0.004
    do while(.true.)
     open(2,file='binary.txt')
       write(2,"(2F10.6)") dkij, error
     close(2)
     error = errorn
     dkij  = dkijn
     dkijn = dkij + direc*delk

     call OBF('expdata.txt', dkijn, errorn)
     write(*,"('kij & error = ', 2F10.6)") dkijn, errorn
     if ( errorn > error ) call system("pause")

    end do

end subroutine FindBinary
