module test_kummer_subroutines

use utils
use chebyshev
use odesolve
use kummer

implicit double precision (a-h,o-z)

double precision :: dnu

contains




subroutine legen0(n,val,der)
implicit double precision (a-h,o-z)
double precision    :: val, der
real*16             :: pi, cosval,sinval, gammaratio,dnu
!
!  Evaluate the Legendre function of the first kind of integer degree n
!  and its derivative at the point 0 to extended precision accuracy.
!
!  Input parameters:
!     dnu - the degree of the polynomial to evaluate
!
!  Output parameters:
!     val - the value of the polynomial at 0
!     der - the value of its derivative at 0 
!
data pi / 3.14159265358979323846264338327950288d0 /

val = 0
der = 0
dnu = n

! Compute cosval = cos(pi/2 * dnu) and sinval = sin(pi/2*dnu)


if ( mod(n,4) == 0) cosval=1
if ( mod(n,4) == 1) cosval=0
if ( mod(n,4) == 2) cosval=-1
if ( mod(n,4) == 3) cosval=0

if ( mod(n,4) == 0) sinval=0
if ( mod(n,4) == 1) sinval=1
if ( mod(n,4) == 2) sinval=0
if ( mod(n,4) == 3) sinval=-1

! Compute the ratio of gamma functions by brute force if n is small, and use an asymptotic
! expansion if not.
gammaratio =  gamma(dnu/2+0.5d0)/gamma(dnu/2+1.0d0)


if (dnu .ge. 100) then
gammaratio = -0.104957025613225462481750387218029248d16/dnu**30-  &
0.3435274821674695538341969877678682d13/dnu**29+  &
0.137012508230061862467393693767518715d14/dnu**28+  &
0.52184449177744945865946227168534838d11/dnu**27-  &
0.20802210228419220675266715332806065d12/dnu**26-  &
0.93358732709156442705307398873761416d9/dnu**25+  &
0.37189644043206042515191168099049809d10/dnu**24+  &
0.19959724959938103205809132306364276d8/dnu**23-  &
0.79435712783055360854592176681049474d8/dnu**22-  &
0.51897923979290808045377229245471765d6/dnu**21+  &
0.20627383965581267895358219780481014d7/dnu**20+  &
0.16766118747003365058835187265431109d5/dnu**19-  &
0.66511566853934577951098546692891475d5/dnu**18-  &
0.6912167061162564516424966534975484d3/dnu**17+  &
0.27339931645089473537297049513979147d4/dnu**16+  &
0.37638891605797570836668516338219077d2/dnu**15-  &
0.14815719515134617540873450267244743d3/dnu**14-  &
0.28341526372214554549138737809124132d1/dnu**13+  &
0.11064038264157344432164622228873395d2/dnu**12+  &
0.31450601605351481364931293755092798d0/dnu**11-  &
0.12103480338672938011211980036999511d1/dnu**10-  &
0.5638860579751321228004007809251394d-1/dnu**9+  &
0.21215037666443619840288699752634574d0/dnu**8+  &
0.18752313014255059774912529012118952d-1/dnu**7-  &
0.6888076310874815972557053234370966d-1/dnu**6-  &
0.14501213286052244152751691019728349d-1/dnu**5+  &
0.5524271728019902534381596578944133d-1/dnu**4+  &
0.4419417382415922027505277263155306d-1/dnu**3-  &
0.3535533905932737622004221810524245d0/dnu**2+  &
0.14142135623730950488016887242096981d1/dnu
gammaratio = gammaratio *sqrt(dnu)
endif

if (n .le. -100) then
dnu = -dnu
gammaratio = 0.104957025613225462481750387218029248d16/dnu**30-  &
0.343527482167469553834196987767868205d13/dnu**29-  &
0.137012508230061862467393693767518715d14/dnu**28+  &
0.521844491777449458659462271685348384d11/dnu**27+  &
0.208022102284192206752667153328060649d12/dnu**26-  &
0.933587327091564427053073988737614164d9/dnu**25-  &
0.37189644043206042515191168099049809d10/dnu**24+  &
0.199597249599381032058091323063642763d8/dnu**23+  &
0.794357127830553608545921766810494744d8/dnu**22-  &
0.518979239792908080453772292454717655d6/dnu**21-  &
0.206273839655812678953582197804810138d7/dnu**20+  &
0.167661187470033650588351872654311091d5/dnu**19+  &
0.66511566853934577951098546692891475d5/dnu**18-  &
0.691216706116256451642496653497548362d3/dnu**17-  &
0.273399316450894735372970495139791469d4/dnu**16+  &
0.376388916057975708366685163382190772d2/dnu**15+  &
0.148157195151346175408734502672447426d3/dnu**14-  &
0.283415263722145545491387378091241324d1/dnu**13-  &
0.110640382641573444321646222288733953d2/dnu**12+  &
0.314506016053514813649312937550927984d0/dnu**11+  &
0.121034803386729380112119800369995106d1/dnu**10-  &
0.563886057975132122800400780925139408d-1/dnu**9-  &
0.212150376664436198402886997526345737d0/dnu**8+  &
0.187523130142550597749125290121189519d-1/dnu**7+  &
0.688807631087481597255705323437096598d-1/dnu**6-  &
0.145012132860522441527516910197283494d-1/dnu**5-  &
0.552427172801990253438159657894413312d-1/dnu**4+  &
0.44194173824159220275052772631553065d-1/dnu**3+  &
0.35355339059327376220042218105242452d0/dnu**2+  &
0.141421356237309504880168872420969808d1/dnu
gammaratio = -gammaratio*sqrt(dnu) * sinval/cosval

endif


val        = 1.0d0/sqrt(pi) * cosval * gammaratio
der        = 2.0d0/sqrt(pi) * sinval * 1.0d0/gammaratio

end subroutine


subroutine lege0(n,x,pol)
implicit double precision (a-h,o-z)
integer n
double precision x,pol
!
!  Evaluate the Legendre polynomial of degree n at the point x using the
!  3-term recurrence relation.  As is well-known, this is somewhat inaccurate
!  when x is near the points +- 1.
!
!  Input parameters:
!
!    n - the degree of the polynomial to evaluate
!    x - the point at which the polynomial is to be evaluated
!
!  Output parameters:
!
!    pol - the value of P_n(x)
!

if (n == 0) then
pol = 1.0d0
else if (n == 1) then
pol = x
else
p1 = 1
p2 = x

do j=2,n
   p  = ((2*j-1)*x*p2-(j-1)*p1)/j
   p1 = p2
   p2 = p
end do

pol = p
endif

end subroutine

subroutine qlegendre(t,val)
implicit double precision (a-h,o-z)
double precision, intent(in)   :: t
double precision, intent(out)  :: val
val = 1.0d0/(1-t**2)**2 + dnu*(dnu+1)/(1-t**2)
end subroutine

subroutine qtest(t,val)
implicit double precision (a-h,o-z)
double precision, intent(in)   :: t
double precision, intent(out)  :: val
val = dnu**2 * (1-t**2*cos(3*t)) 
end subroutine

subroutine odetest(t,y,yp,f,dfdy,dfdyp)
implicit double precision (a-h,o-z)
double precision, intent(in)   :: t,y,yp
double precision, intent(out)  :: f,dfdy, dfdyp
f     =  -dnu**2 * (1-t**2*cos(3*t)) * y
dfdy  =  -dnu**2 * (1-t**2*cos(3*t))  
dfdyp = 0
end subroutine

end module


program test_kummer

use utils
use odesolve
use kummer
use test_kummer_subroutines

implicit double precision (a-h,o-z)

double precision, allocatable :: xscheb(:),whtscheb(:),chebintl(:,:),chebintr(:,:), &
   ucheb(:,:),vcheb(:,:)

double precision, allocatable :: ab(:,:),alpha(:,:),alphap(:,:),alphapp(:,:),alphappp(:,:)
double precision, allocatable :: alphainv(:,:),alphainvp(:,:),abinv(:,:)

double precision, allocatable :: alpha_coefs(:,:),alphainv_coefs(:,:),alphap_coefs(:,:)
double precision, allocatable :: ts(:),vals(:),vals0(:),errs(:)

double precision              :: abin(2,1)
double precision, allocatable :: ys(:,:), yders(:,:), yder2s(:,:), ab2(:,:), vals000(:)

print *,"Compute Legendre functions over [0,0.9]"
print *,""

eps0   = epsilon(0.0d0)
pi     = acos(-1.0d0)
eps    = 1.0d-12
k      = 16
ifleft = 1
nn     = 1000




!
!  Fetch the Chebyshev quadrature and related matrices.
!

call chebexps(k,xscheb,whtscheb,ucheb,vcheb,chebintl,chebintr)


! if (eps0 .lt. 1.0d-16) then
! iw = 101
! open(iw,FILE="sol.dat")
! do ii=0,20
! a         = -1
! b         =  1
! dnu       = 2**ii
! nintsin   = 1
! abin(1,1) = a
! abin(2,1) = b
! val0      = 1
! der0      = dnu
! eps       = 1.0d-15

! call ode_solve_ivp_adap(ier,eps,nintsin,abin,k,xscheb,chebintl,ucheb,nints2,ab2,  &
!    ys,yders,yder2s,odetest,val0,der0)
! call chebpw_eval(nints2,ab2,k,xscheb,ys,b,val)
! write (iw,"(F44.36)") val
! end do
! stop
! else
! iw = 101
! open(iw,FILE="sol.dat")
! allocate(vals000(21))
! do ii=0,20
! read (iw,"(F44.36)") vals000(ii)
! end do

! endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LEGENDRe

write (*,"(A,A,A)")  "nu          ",  "time ","          error"

do iii=0,20
norder   = 2**iii
dnu      = norder
a        =  0.0d0
b        =  0.9d0


call legen0(norder,val0,der0)

!
!  Average the time required to construt the phase functions and solve
!  the IVP.
!

call elapsed(t1)
do ii=1,nn

!
!  Construct the phase functions
!
call kummer_adap(eps,a,b,qlegendre,k,xscheb,chebintl,chebintr,ucheb, &
  nints,ab,alphap,alphapp)

call  kummer_phase(ifleft,k,xscheb,chebintl,chebintr,ucheb, &
   nints,ab,alpha,alphap,alphapp)

!
!  Solve the IVP
!

call kummer_coefs(ifleft,nints,ab,k,xscheb,alphap,alphapp,val0,der0,c1,c2)
end do


call elapsed(t2)
t_phase = (t2-t1)/nn

!
!  Check the error of the Legendre function at the right-hand side
!

call kummer_eval(nints,ab,k,xscheb,alpha,alphap,c1,c2,b,val)
val=val/sqrt(1-b**2)
call lege0(norder,b,val0)
derr = abs(val-val0)

write (*,"(I10.10,D15.7,D15.7)")      norder,t_phase,derr

end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print *,""
write (*,"(A,A,A)")  "  nu          ",  "           time ","          error"


do iii=0,20
dnu  = 2**iii
a    = -1.0d0
b    =  1.0d0
val0 = 1
der0 = dnu

!
!  Average the time required to construt the phase functions and solve
!  the IVP.
!


call elapsed(t1)
do ii=1,nn

!
!  Construct the phase functions
!
call kummer_adap(eps,a,b,qtest,k,xscheb,chebintl,chebintr,ucheb, &
  nints,ab,alphap,alphapp)

call  kummer_phase(ifleft,k,xscheb,chebintl,chebintr,ucheb, &
   nints,ab,alpha,alphap,alphapp)

!
!  Solve the IVP
!

call kummer_coefs(ifleft,nints,ab,k,xscheb,alphap,alphapp,val0,der0,c1,c2)
end do

call elapsed(t2)
t_phase = (t2-t1)/nn

!
!  Check the error at the right-hand side
!

nintsin   = 1
abin(1,1) = a
abin(2,1) = b

call ode_solve_ivp_adap(ier,eps,nintsin,abin,k,xscheb,chebintl,ucheb,nints2,ab2,  &
   ys,yders,yder2s,odetest,val0,der0)

call chebpw_eval(nints2,ab2,k,xscheb,ys,b,val0)
call kummer_eval(nints,ab,k,xscheb,alpha,alphap,c1,c2,b,val)

! val0 = vals000(iii)
derr = abs(val-val0)

write (*,"(D15.7,D15.7,D15.7)")      dnu,t_phase,derr

end do

end program
