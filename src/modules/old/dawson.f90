function daw ( xx )

!*****************************************************************************80
!
!! DAW evaluates Dawson's integral function.
!	Subroutine taken from the package "specfun" found on the internet.
!
!  Discussion:
!
!    This routine evaluates Dawson's integral,
!
!      F(x) = exp ( - x * x ) * Integral ( 0 <= t <= x ) exp ( t * t ) dt
!
!    for a real argument x.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    William Cody
!
!  Reference:
!
!    William Cody, Kathleen Paciorek, Henry Thacher,
!    Chebyshev Approximations for Dawson's Integral,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 171-178.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XX, the argument of the function.
!
!    Output, real ( kind = 8 ) DAW, the value of the function.
!
  implicit none

  real    ( kind = 8 ) daw
  real    ( kind = 8 ) frac
  real    ( kind = 8 ) half
  integer ( kind = 4 ) i
  real    ( kind = 8 ) one
  real    ( kind = 8 ) one225
  real    ( kind = 8 ) p1(10)
  real    ( kind = 8 ) p2(10)
  real    ( kind = 8 ) p3(10)
  real    ( kind = 8 ) p4(10)
  real    ( kind = 8 ) q1(10)
  real    ( kind = 8 ) q2(9)
  real    ( kind = 8 ) q3(9)
  real    ( kind = 8 ) q4(9)
  real    ( kind = 8 ) six25
  real    ( kind = 8 ) sump
  real    ( kind = 8 ) sumq
  real    ( kind = 8 ) two5
  real    ( kind = 8 ) w2
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xlarge
  real    ( kind = 8 ) xmax
  real    ( kind = 8 ) xsmall
  real    ( kind = 8 ) xx
  real    ( kind = 8 ) y
  real    ( kind = 8 ) zero
!
!  Mathematical constants.
!
  data zero / 0.0D+00 /
  data half / 0.5D+00 /
  data one / 1.0D+00 /
  data six25 /6.25D+00 /
  data one225 /12.25d0 /
  data two5 /25.0d0/
!
!  Machine-dependent constants
!
  data xsmall /1.05d-08/
  data xlarge /9.49d+07/
  data xmax /2.24d+307/
!
!  Coefficients for R(9,9) approximation for  |x| < 2.5
!
  data p1/-2.69020398788704782410d-12, 4.18572065374337710778d-10, &
          -1.34848304455939419963d-08, 9.28264872583444852976d-07, &
          -1.23877783329049120592d-05, 4.07205792429155826266d-04, &
          -2.84388121441008500446d-03, 4.70139022887204722217d-02, &
          -1.38868086253931995101d-01, 1.00000000000000000004d+00/
  data q1/ 1.71257170854690554214d-10, 1.19266846372297253797d-08, &
           4.32287827678631772231d-07, 1.03867633767414421898d-05, &
           1.78910965284246249340d-04, 2.26061077235076703171d-03, &
           2.07422774641447644725d-02, 1.32212955897210128811d-01, &
           5.27798580412734677256d-01, 1.00000000000000000000d+00/
!
!  Coefficients for R(9,9) approximation in J-fraction form
!  for  x in [2.5, 3.5)
!
  data p2/-1.70953804700855494930d+00,-3.79258977271042880786d+01, &
           2.61935631268825992835d+01, 1.25808703738951251885d+01, &
          -2.27571829525075891337d+01, 4.56604250725163310122d+00, &
          -7.33080089896402870750d+00, 4.65842087940015295573d+01, &
          -1.73717177843672791149d+01, 5.00260183622027967838d-01/
  data q2/ 1.82180093313514478378d+00, 1.10067081034515532891d+03, &
          -7.08465686676573000364d+00, 4.53642111102577727153d+02, &
           4.06209742218935689922d+01, 3.02890110610122663923d+02, &
           1.70641269745236227356d+02, 9.51190923960381458747d+02, &
           2.06522691539642105009d-01/
!
!  Coefficients for R(9,9) approximation in J-fraction form
!  for  x in [3.5, 5.0]
!
  data p3/-4.55169503255094815112d+00,-1.86647123338493852582d+01, &
          -7.36315669126830526754d+00,-6.68407240337696756838d+01, &
           4.84507265081491452130d+01, 2.69790586735467649969d+01, &
          -3.35044149820592449072d+01, 7.50964459838919612289d+00, &
          -1.48432341823343965307d+00, 4.99999810924858824981d-01/
  data q3/ 4.47820908025971749852d+01, 9.98607198039452081913d+01, &
           1.40238373126149385228d+01, 3.48817758822286353588d+03, &
          -9.18871385293215873406d+00, 1.24018500009917163023d+03, &
          -6.88024952504512254535d+01,-2.31251575385145143070d+00, &
           2.50041492369922381761d-01/
!
!  Coefficients for R(9,9) approximation in J-fraction form
!  for 5.0 < |x|.
!
  data p4/-8.11753647558432685797d+00,-3.84043882477454453430d+01, &
          -2.23787669028751886675d+01,-2.88301992467056105854d+01, &
          -5.99085540418222002197d+00,-1.13867365736066102577d+01, &
          -6.52828727526980741590d+00,-4.50002293000355585708d+00, &
          -2.50000000088955834952d+00, 5.00000000000000488400d-01/
  data q4/ 2.69382300417238816428d+02, 5.04198958742465752861d+01, &
           6.11539671480115846173d+01, 2.08210246935564547889d+02, &
           1.97325365692316183531d+01,-1.22097010558934838708d+01, &
          -6.99732735041547247161d+00,-2.49999970104184464568d+00, &
           7.49999999999027092188d-01/

  x = xx

  if ( xlarge < abs ( x ) ) then

    if ( abs ( x ) <= xmax ) then
      daw = half / x
    else
      daw = zero
    end if

  else if ( abs ( x ) < xsmall ) then

    daw = x

  else

    y = x * x
!
!  ABS(X) < 2.5.
!
    if ( y < six25 ) then

      sump = p1(1)
      sumq = q1(1)
      do i = 2, 10
        sump = sump * y + p1(i)
        sumq = sumq * y + q1(i)
      end do

      daw = x * sump / sumq
!
!  2.5 <= ABS(X) < 3.5.
!
    else if ( y < one225 ) then

      frac = zero
      do i = 1, 9
        frac = q2(i) / ( p2(i) + y + frac )
      end do

      daw = ( p2(10) + frac ) / x
!
!  3.5 <= ABS(X) < 5.0.
!
    else if ( y < two5 ) then

      frac = zero
      do i = 1, 9
        frac = q3(i) / ( p3(i) + y + frac )
      end do

      daw = ( p3(10) + frac ) / x

    else
!
!  5.0 <= ABS(X) <= XLARGE.
!
      w2 = one / x / x

      frac = zero
      do i = 1, 9
        frac = q4(i) / ( p4(i) + y + frac )
      end do
      frac = p4(10) + frac

      daw = ( half + half * w2 * frac ) / x

    end if

  end if

  return
end
