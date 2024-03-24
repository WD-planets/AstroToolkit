Module aovsub
      Use aovconst
      Implicit None
!
!      Interface sortx
!         Module Procedure sortx_m_r4, sortx_m_r8
!      End Interface sortx
!
      Private
      Public amulthw, aovmhw, aovw, cracow, ferraz_glsw, givens, givensacc, &
         lombw, norm_obs, ortint, ortprj, powspw, sortx
!
!
Contains
!
      Pure Subroutine aovmhw (t, v, w, npar, ncov, fr0, frs, nfr, th)
!! Routine for LSQ fit of Fourier harmonics via projection on orthogonal
!! polynomialsfor. Weights are applied to observations.
!
!!Purpose:<br>
!!AOVMH & AOVMHW  fit unevenly sampled time series with a multiharmonic
!!Fourier series, respectively for no weights and weighted observations.
!!They return spectrum for a range of Fourier frequencies. For NH harmonics
!!input NPAR = 2NH+1. For NPAR=3 and SPTYPE=SP_MODEL, AOVMH & AOVMHW results
!!are identical to Ferraz-Mello(1981) generalized Lomb-Scargle (GLS) spectrum.
!!NCOV is provided for conformity but it is ignored.
!!
!!Method: <br>
!!By Szego recurrence AOVMH & AOVMHW generate orthogonal trigonometric
!!polynomials and next by orthogonal projection expand uneven sampled time
!!series in terms of these polynomials. Such an algorithm constitutes the
!!fastest implementation of Fourier series least squares fit. For NPAR>>3
!!these algorithms are particularly suitable to detect smooth non-sinusoidal
!!oscillations involving sharp features (e.g. eclipses and narrow pulses).
!!For detection of sine signals in statistical sense (sensitivity) AOVMH &
!!AOVMHW never perform worse than either Lomb-Scargle or Power Spectrum.
!!For non-sinusoidal signals AOVMH & AOVMHW advantage becomes significant.
!!For identical NPAR, AOVMH & AOVMHW are about 10 times slower than more crude
!!AOV/AOVW routines.
!
!!Copyrights and Distribution:<br>
!!This routines are subject to copyrights by its author,
!!(C) A. Schwarzenberg-Czerny 1998-2005,  alex@camk.edu.pl
!!Its distribution is free, except that distribution of modifications
!!is prohibited.
!!Please quote A.Schwarzenberg-Czerny, 1996, Astrophys. J.,460, L107
!!while using AOVMH & AOVMHW.
         Real (TIME), Intent (In) :: t (:)!! times of observations
         Real (SP), Intent (In) :: v (:)!! values of observations
         Real (SP), Intent (In) :: w (:)!! weights of observations
         Integer, Intent (In) :: npar !! number of parameters (2*nh+1)
         Integer, Intent (In) :: ncov !! dummy
         Real (TIME), Intent (In) :: fr0 !! start frequency
         Real (TIME), Intent (In) :: frs !! frequency step
         Integer, Intent (In) :: nfr !! number of frequencies
         Real (SP), Intent (Out) :: th (nfr)!! result periodogram (th(1)<0 for errors)
!
         Logical :: qform
         Integer :: l, n, nn2, nh, no
         Real (SP) :: sn, eps2, eps3
         Real (TIME) :: ph (size(t))
         Complex (CP), Dimension (size(t)) :: cf, p, q, z, zn
         Complex (CP) :: al
!
         l = ncov ! dummy use
         no = size (t)
         th = 0._SP
         If (npar < 3 .Or. Mod(npar, 2) /= 1) th (1) = - 1.
         If (no <= npar) th (1) = - 2.
         If (th(1) < 0._SP) Return
         nh = npar / 2
         nn2 = npar - 1
         eps2 = Sqrt (epsilon(eps2))
         eps3 = epsilon (eps3) * Sqrt (Real(no, SP))
!
         Do l = 1, nfr
            qform = .True.
            ph = (frs*(l-1)+fr0) * t
            ph = PI2 * (ph-floor(ph))
            z = Exp (cmplx(0._SP, Real(ph, SP), kind=CP))
            ph = ph * nh
            zn = 1._SP
            cf = v * Exp (cmplx(0._SP, Real(ph, SP), kind=CP))
            p = w
!
            Do n = 0, nn2 ! new compact loop
               q = zn * conjg (p)
               sn = 0._SP
               If (qform) sn = sum (Real(q, SP))! (1)
! In some parallel environment one may permanently replace (1) with (2)
               If (sn < eps2) Then
                  qform = .False.
                  sn = sum (Abs(p)**2/w)! (2)
               End If
               sn = 1._SP / Max (sn, eps3)
! NOTE: dot_product(a,b)=sum(conjg(a)*b)
               th (l) = th (l) + Abs (dot_product(p, cf)) ** 2 * sn
               If (n == nn2) Exit
               p = p * z
               al = sum (p) * sn
               p = p - al * q
               zn = zn * z
            End Do
         End Do
      End Subroutine aovmhw
!  End of AOVMH & AOVMHW routines for multiharmonic Fourier spectrum
      Pure Subroutine amulthw (t, v, w, npar, ncov, fr0, frs, nfr, th)
!! Routine for LSQ fit of Fourier harmonics via projection on orthogonal
!! polynomials. Weights are applied to observations.
!
!!Purpose:<br>
!!AOVMH & AOVMHW  fit unevenly sampled time series with a multiharmonic
!!Fourier series, respectively for no weights and weighted observations.
!!They return spectrum for a range of Fourier frequencies. For NH harmonics
!!input NPAR = 2NH+1. For NPAR=3 and SPTYPE=SP_MODEL, AOVMH & AOVMHW results
!!are identical to Ferraz-Mello(1981) generalized Lomb-Scargle (GLS) spectrum.
!!NCOV is provided for conformity but it is ignored.
!!
!!Method: <br>
!!By Szego recurrence AOVMH & AOVMHW generate orthogonal trigonometric
!!polynomials and next by orthogonal projection expand uneven sampled time
!!series in terms of these polynomials. Such an algorithm constitutes the
!!fastest implementation of Fourier series least squares fit. For NPAR>>3
!!these algorithms are particularly suitable to detect smooth non-sinusoidal
!!oscillations involving sharp features (e.g. eclipses and narrow pulses).
!!For detection of sine signals in statistical sense (sensitivity) AOVMH &
!!AOVMHW never perform worse than either Lomb-Scargle or Power Spectrum.
!!For non-sinusoidal signals AOVMH & AOVMHW advantage becomes significant.
!!For identical NPAR, AOVMH & AOVMHW are about 10 times slower than more crude
!!AOV/AOVW routines.
!
!!Copyrights and Distribution:<br>
!!This routines are subject to copyrights by its author,
!!(C) A. Schwarzenberg-Czerny 1998-2005,  alex@camk.edu.pl
!!Its distribution is free, except that distribution of modifications
!!is prohibited.
!!Please quote A.Schwarzenberg-Czerny, 1996, Astrophys. J.,460, L107
!!while using AOVMH & AOVMHW.
         Real (TIME), Intent (In) :: t (:)!! times of observations
         Real (SP), Intent (In) :: v (:)!! values of observations
         Real (SP), Intent (In) :: w (:)!! weights of observations
         Integer, Intent (In) :: npar (:)!! number of parameters (2*nh+1)
         Integer, Intent (In) :: ncov !! dummy
         Real (TIME), Intent (In) :: fr0 !! start frequency
         Real (TIME), Intent (In) :: frs !! frequency step
         Integer, Intent (In) :: nfr !! number of frequencies
         Real (SP), Intent (Out) :: th (nfr)
            !! result periodogram (th(1,1)<0 for errors)
!
         Logical :: qform
         Integer :: l, n, k, j, no, npm, npx, nh2 (size(npar)), nc
         Real (SP) :: sn, eps2, eps3, vr, thl (0:size(npar)), f
         Real (SP), Dimension (size(npar)) :: ev, sig, fc
         Real (TIME) :: ph (size(t))
         Complex (CP), Dimension (size(t)) :: p, q, z, zn
         Complex (CP) :: al, cf (size(npar), size(t))
!
         l = ncov ! dummy use
         no = size (t)
         npm = size (npar)
         npx = maxval (npar) / 2 * 2
         nh2 = npar / 2 * 2
         If (minval(nh2) < 2) th (1) = - 1.
         If (no+1 <= npx) th (1) = - 2.
         If (npm > 1 .And. minval(nh2(2:)-nh2(:npm-1)) <= 0) th (1) = - 3.
         If (th(1) < 0._SP) Return
!
         vr = sum (w*v**2)
         eps2 = Sqrt (epsilon(eps2))
         eps3 = epsilon (eps3) * Sqrt (Real(no, SP))
! C.f. Eadie et al., 1971, p. 70 and 239
         Forall (k=1:npm)
            thl (k) = (nh2(k)+1)
            sig (k) = (no-nh2(k)-2)
            ev (k) = (-(1./thl(k)**2-1./sig(k)**2)/3.+1./thl(k)-1./sig(k)) / 2.
            sig (k) = Sqrt ((1./thl(k)**3+1./sig(k)**3)/3.+ &
               (1./thl(k)**2+1./sig(k)**2+1./thl(k)+1./sig(k))/2.)
         End Forall
!
         Do l = 1, nfr
            thl = 0.
            thl (0) = vr
            qform = .True.
            ph = (frs*(l-1)+fr0) * t
            ph = PI2 * (ph-floor(ph))
            z = Exp (cmplx(0._SP, Real(ph, SP), kind=CP))
            zn = 1._SP
            Do k = 1, npm
               j = nh2 (k) / 2
               cf (k, :) = v * Exp (cmplx(0._SP, Real(j*ph, SP), kind=CP))
            End Do
            p = w
!
            nc = 1
            Do n = 0, npx ! new compact loop
               q = zn * conjg (p)
               sn = 0._SP
               If (qform) sn = sum (Real(q, SP))! (1)
! In some parallel environment one may permanently replace (1) with (2)
               If (sn < eps2) Then
                  qform = .False.
                  sn = sum (Abs(p)**2/w)! (2)
               End If
               sn = 1._SP / Max (sn, eps3)
! NOTE: dot_product(a,b)=sum(conjg(a)*b)
               Forall (k=nc:npm, n <= nh2(k)) thl (k) = thl (k) + &
                  Abs (dot_product(p, cf(k, :))) ** 2 * sn
! Abandon if no improvement
               If (n == nc) Then
                  f = (thl(nc)-thl(nc-1)) * (no-nh2(nc)-2) / &
                     (Max(vr-thl(nc), epsilon(thl))) * &
                     (nh2(nc)-nh2(nc-1))
                  If (f < fc(nc) .Or. n == npx) Exit
                  nc = nc + 1
               End If
               p = p * z
               al = sum (p) * sn
               p = p - al * q
               zn = zn * z
            End Do
! dodac poprawke na korelacje
            Forall (k=1:nc) thl (k) = &
               (Log(thl(k)*(no-nh2(k)-2)/(Max(vr-thl(k), &
                epsilon(thl))*(nh2(k)+1)))-ev(k)) / sig (k)
            th (l) = maxval (thl)
         End Do
      End Subroutine amulthw
!  End of AOVMH & AOVMHW routines for multiharmonic Fourier spectrum
!
      Subroutine aovw (t, vin, w, np, ncov, fr0, frs, nfr, th, iflex)
!! Routine for phase folding and binning of observations using weights
!
!! Purpose:<br>
!! The state-of-art routines for either period search by phase folding and
!! binning of observations (for NP>0) or for search for periodic transits/pulses
!! of short duty cycle, with constat signal elsewhere (for NP<0). In the latter
!! case data are fitted with just two unequal phase bins: in- and out-of transit.
!
!! Method:<br>
!! In any case observations are phase folded and binned into NP equal phase
!! bins starting at NCOV sub-bin phases. Subsequently the algorithm
!! either sums squares of bin averages (for NP>0). For transit search
!! (NP<0) algorithm searches for the maximum value phase bin (considered
!! in-transit for magnitude units). Its squared value and identities stemming
!! from Parceval theorem are employed to evaluate sum of squares for two
!! phase-bins in- and out-of transit. Note: here transits are defined as MAXIMA
!! (i.e. in magnitude units). For flux units (minima) reverse sign of observations.
!! If algorithm detects poor phase coverage (bins with <CTMIN observations)
!! it rebins observations into bins of near equal occupancy, by sorting
!! observations in phase. This is a rare situation hence amounts to no
!! significant slowing down. Thanks to this feature, for null frequency
!! these routines return useful variability indicator as long as the sort
!! routine is stable (e.g. merge sort).
!!
!! Advantage:<br>
!! By checking and repairing phase coverage these routines perform robust and
!! near-optimal period search, returnind spectrum of known classical
!! probability distribution, equally well for small and large data. Slight
!! statistical inefficiency stems from step-like function yielding slightly
!! worse fit that equal size Fourier fit. This may be remedied by calling
!! AOVMH/AOVMHV instead, at the speed cost. In terms of speed the present
!! AOV/AOVW implementation is state-of-art as it employs pre-binning into
!! NP*NCOV sub-bins and only afterwords combines subbins into suitably
!! phased bins. See B. Enoch et al., 2012, A&A 548, A48
!! for independent evaluation of the present sub-bin algorithm.
!
!! Copyrights and Distribution:<br>
!! This package is subject to copyrights by its author,
!! (C) A. Schwarzenberg-Czerny 2003-2018,  alex@camk.edu.pl
!! Its distribution is free, except that distribution of modifications
!! is prohibited.
!! Please quote A. Schwarzenberg-Czerny & J.-Ph. Beaulieu, 2005, MNRAS 365, 165,
!! while using AOV or AOVW.
         Real (TIME), Intent (In) :: t (:)!! times of observations
         Real (SP), Intent (In) :: vin (:)!! values of observations
         Real (SP), Intent (In) :: w (:)!! weights of observations
         Integer, Intent (In) :: np !! number of phase bins
         Integer, Intent (In) :: ncov !! number of phase coverages
         Real (TIME), Intent (In) :: fr0 !! start frequency
         Real (TIME), Intent (In) :: frs !! frequency step
         Integer, Intent (In) :: nfr !! number of frequencies
         Real (SP), Intent (Out) :: th (nfr)!! result periodogram (th(1)<0 for errors)
         Integer, Intent (Out) :: iflex !! number of poor phase coverages
!
         Integer :: ind (size(t)), i, ibin, ip, nbc, ifr, npar, no
         Real (SP) :: ph (size(t)), v (size(t)), sav, ctw, sw, &
        & nct ((Abs(np)+1)*ncov), ave ((Abs(np)+1)*ncov)
         Real (TIME) :: fr, dph
         Logical :: transit
!
         no = size (t)
         npar = Abs (np)
         transit = np < 0
!
         th (1) = 0._SP
         If (npar < 2) th (1) = - 1.
         If (no <= npar+npar .Or. no < CTMIN*npar) th (1) = - 2.
         If (th(1) .Lt. 0._SP) Return
!
         sw = sum (w)
         v = vin * w
!
         nbc = npar * ncov
         iflex = 0
         ctw = CTMIN * sw / no
         Do ifr = 1, nfr ! Loop over frequencies
            fr = Real (ifr-1, kind=TIME) * frs + fr0
            Do ip = 0, 1
               ave = 0._SP
               nct = 0._SP
               If (ip == 0) Then ! Try default fixed bins ...
                  Do i = 1, no ! MOST WORK HERE
                     dph = t (i) * fr ! real(TIME) :: dph, t, fr
                     dph = dph - floor (dph)
                     ph (i) = Real (dph, SP)
! note: mod accounts for rounding up due to double-->real conversion
                     ibin = Int (dph*nbc) + 1
                     ave (ibin) = ave (ibin) + v (i)
                     nct (ibin) = nct (ibin) + w (i)
                  End Do
               Else !... elastic bins, if necesseary
!
                  iflex = iflex + 1 ! sort index ind using key ph
                  Call sortx (ph, ind)! NRf77 indexx would do
                  Do i = 1, no
                     ibin = i * nbc / no + 1
                     ave (ibin) = ave (ibin) + v (ind(i))
                     nct (ibin) = nct (ibin) + w (ind(i))
                  End Do
               End If
  ! counts: sub-bins=>bins
               nct (1+nbc:ncov+nbc) = nct (1:ncov)
               sav = 0._SP
               Do ibin = ncov + nbc, 1, - 1
                  sav = sav + nct (ibin)
                  nct (ibin) = sav
               End Do
               Do ibin = 1, nbc
                  nct (ibin) = nct (ibin) - nct (ibin+ncov)
               End Do
               If (minval(nct(1:nbc)) >= ctw) Exit
            End Do
  ! data: sub-bins=>bins */
            ave (1+nbc:ncov+nbc) = ave (1:ncov)
            sav = 0._SP
            Do ibin = ncov + nbc, 1, - 1
               sav = sav + ave (ibin)
               ave (ibin) = sav
            End Do
            Do ibin = 1, nbc
               ave (ibin) = ave (ibin) - ave (ibin+ncov)
            End Do
!
            If (transit) Then ! AoV statistics for transits ...
               ave (1:nbc) = ave (1:nbc) / nct (1:nbc)
               ibin = maxloc (ave(1:nbc), dim=1)
               sav = ave (ibin)
               sav = sav * sav * nct (ibin) * sw / (sw-nct(ibin))
            Else ! the same for the classical AoV statistics:
               sav = sum (ave(1:nbc)**2/nct(1:nbc)) / ncov
            End If
            th (ifr) = sav
         End Do
         i = npar
         If (transit) i = 2
!
      End Subroutine aovw
!
!  C R A C O W
!
!  Purpose:
!  To find least squares solution of normal equations.
!
!  Input:
!  A(n1,n1) - a symmetric matrix. If matrix A represents normal
!    equations, put the sum of squares of RHS into element A(n1,n1)
!    for calculation of the chi2 and covariance matrix. Only A(i,j)
!    for i>=j are requires for input.
!  N1       - size of matrix A (number of unknowns) +1
!  M        - M>=N1 - number of condition equations
!           - M< N1 - only triangle root and its triangle inverse are
!             returned
!
!  Output:
!  M>=N1:
!    A(i,n1) - unknowns (i=1,...,n)
!    A(i,i)  - their variances
!    A(n1,n1)- variance of an observation of unit weight (Chi2/D.F.)
!    A(i,j)  - for j>=i half of the covariance matrix of unknowns
!
!  M<N1:
!    A(i,j)  - for j>=i inverse triangle root matrix (its square is A**(-1))
!    A(i,j)  - for j<i contains square root matrix except of diagonal.
!              The diagonal contains arithmetic inverses of its elements.
!
!  Method:
!  Banachiewicz-Cholesky decomposition of the normal equations.
!  This subroutine is particularly suitable for well posed problems
!  with near-orthogonal base functions.
!    This method should not be used for near singular matrices,
!  since it sufferes from large numerical errors then. Use the singular
!  value decomposition method in such cases.
!
!  This package is subject to copyrights by its author,
!  (C) A. Schwarzenberg-Czerny 1987-2011,  alex@camk.edu.pl
!
      Function cracow (a, m)
         Real (SP) :: cracow
         Integer, Intent (In) :: m
         Real (SP), Intent (Inout) :: a (:, :)
!
         Integer :: n1, i, j
         Real (SP) :: aii
!
!  Find Banachiewicz-Cholesky triangle root and its inverse Q, so that
!  Q*Q=A**(-1). The unknowns -a(i,n1) are a byproduct. So...
         n1 = size (a, dim=1)
         cracow = 0._SP
         Do i = 1, n1 - 1
            If (a(i, i) <= 0._SP) Go To 99
            aii = Sqrt (a(i, i))
!  ... complete previous i step, ...
            a (1:n1, i) = a (1:n1, i) / aii
!  ... complete next row of the triangle matrix, ...
            Do j = i + 1, n1
               a (i, i+1) = a (i, i+1) - dot_product (a(j, 1:i), a(i+1, 1:i))
            End Do
!  ... compute next column of the inverse matrix.
            a (i, i) = 1._SP / aii
            Do j = 1, i
               a (j, i+1) = - dot_product (a(j, j:i), a(i+1, j:i))
            End Do
         End Do
!
!  Compute covariance matrix and reverse sign of unknowns
         If (m >= n1) Then
            a (n1, n1) = a (n1, n1) / (m-n1+1)
            Do i = 1, n1 - 1
               a (i, n1) = - a (i, n1)
               Do j = 1, i
                  a (j, i) = dot_product (a(i, i:n1-1), a(j, i:n1-1)) * a (n1, n1)
               End Do
            End Do
         End If
         cracow = 1._SP
99       Return
      End Function cracow
!
      Pure Subroutine ferraz_glsw (t, vin, w, npar, ncov, fr0, frs, nfr, th)
!
!! Purpose:
!! ferraz_glsw - returns Ferraz-Mello generalized least squares spectrum for
!! weighted unevenly sampled observations. Asuming mean subtracted data the
!! algorithm fits them with a sinusoid f(t)=Asin(om*t)+Bcos(om*t) and returns
!! weighted squared norm of residuals ||x-f||^2.
!!
!! Reference: S. Ferraz-Mello, 1981, Astronomical Journal 86, 619.
!!
!!  Copyrights and Distribution:
!!  This routine is subject to copyrights by its author,
!!  (C) A. Schwarzenberg-Czerny 2018,  alex@camk.edu.pl
!!  Its distribution is free, except that distribution of modifications is
!!  prohibited.
         Real (TIME), Intent (In) :: t (:)!! times of observations
         Real (SP), Intent (In) :: vin (:)!! values of observations
         Real (SP), Intent (In) :: w (:)!! weights of observations
         Integer, Intent (In) :: npar !! dummy
         Integer, Intent (In) :: ncov !! dummy
         Real (TIME), Intent (In) :: fr0 !! start frequency
         Real (TIME), Intent (In) :: frs !! frequency step
         Integer, Intent (In) :: nfr !! number of frequencies
         Real (SP), Intent (Out) :: th (nfr)!! result periodogram (th(1)<0 for errors)
!
         Integer :: l
         Real (SP), Dimension (size(t)) :: c, s, v
         Real (SP) :: sw, vr, c1, s1, sx, cx, ss, cc, cs, a1, a2, cc1, cc2
         Real (TIME) :: dph (size(t))
!
         l = ncov + npar + size (w)! dummy use
         th (1) = 0._SP
         If (size(t) .Lt. 4) th (1) = - 1.
         If (th(1) .Lt. 0._SP) Return
!
         v = vin * w
         sw = sum (w)
         vr = sum (v*vin)
         Do l = 1, nfr
            dph = (frs*(l-1)+fr0) * t
            s = Real (PI2*(dph-floor(dph)), SP)
            c = Cos (s)
            s = Sin (s)
            c1 = sum (w*c)
            s1 = sum (w*s)
            cx = sum (c*v)
            sx = sum (s*v)
            cc = sum (w*c**2)
            ss = sum (w*s**2)
            cs = sum (w*c*s)
            a1 = Sqrt (1._SP/Max(cc-c1**2/sw, epsilon(a1)*sw))
            cc2 = a1 * (cs-c1*s1/sw)
            a2 = 1._SP / Max (ss-s1**2/sw-cc2**2, epsilon(a2)*sw)
            cc1 = a1 * cx
            th (l) = (cc1**2+a2*(sx-cc1*cc2)**2) / vr
         End Do
      End Subroutine ferraz_glsw
!
      Function givens (r, idf, detnth)
! (C) Alex Schwarzenberg-Czerny,2000      alex@camk.edu.pl
         Implicit None
         Real (SP) :: givens
         Integer, Intent (In) :: idf
         Real (SP), Intent (Inout) :: r (:, :)
         Real (SP), Intent (Out) :: detnth
!
! Purpose: GIVENS calculates the least squares solution of an
! overdetermined linear system by accurate Givens rotations.
! Prior to GIVENS call m-times GIVENSA in order to build the
! triangle system R
! Input:
! r(n1,n1) - triangle system R and its R.H.S
! idf - >0 solve normal equations for idf degrees of freedom
!       <0 and r(n1,n1)=1: inverse triangle matrix R
! Output:
! givens - condition ratio (min eigenv/max eigenv)
! detnth=Det(en)^(1/n) - n-th root of determinant
! For idf>0:
! r(n1,:)       - unknowns
! (r(k,k),k=1,n) variances
! (r(:k,k),k=1,n) half of covariance matrix
! (r(k+1:,k),k=2,n) half of inverse triangle root Q
! r(n1,n1) - chi^2/idf (=RHS variance=std.dev^2)
! For idf<0:
! ((r(j,k),j=1,k),k=1,n) inverse triangle root  Q
! (C) Alex Schwarzenberg-Czerny, Dec.2000,     alex@camk.edu.pl
!
! tested OK
!
         Integer :: n1, n, i, j
         Real (SP), Parameter :: tol = 30._SP * epsilon (1._SP)
         Real (SP) :: rii, riimx, riimn
!
         n1 = size (r, dim=1)
         n = n1 - 1
         r (n1, n1) = (r(n1, n1)-dot_product(r(1:n, n1), r(1:n, n1))) / Max (1, idf)
         Do i = 1, n
            rii = 0._SP
            If (Abs(r(i, i)) > tol) rii = 1._SP / r (i, i)
            r (i, i) = - 1._SP
            Do j = 1, i
               r (i, j) = r (i, j) * rii
               r (i+1, j) = - dot_product (r(j:i, j), r(j:i, i+1))
! en(j,i+1)=-dot_product(en(j,k:i),en(i+1,k:i))
            End Do
         End Do
! Already solved, square Q to get  inverse/covariance
         detnth = 0.
         riimx = Abs (r(1, 1))
         riimn = riimx
         Do i = 1, n
            rii = Abs (r(i, i))
            If (idf .Gt. 0) Then
               Do j = 1, i
                  r (j, i) = dot_product (r(i:n, i), r(i:n, j)) * r (n1, n1)
               End Do
            End If
            If (rii > 0._SP) detnth = detnth - Log (rii)
            riimx = Max (rii, riimx)
            riimn = Min (rii, riimn)
         End Do
         givens = riimn / riimx
         detnth = Exp (2._SP*detnth/n)
!
      End Function givens
!
      Subroutine givensacc (a, r)
! (C) Alex Schwarzenberg-Czerny,2000      alex@camk.edu.pl
         Implicit None
         Real (SP), Intent (In) :: a (:)
         Real (SP), Intent (Inout) :: r (:, :)
! Purpose: GIVENSACC prepares the least squares solution (LSQ) of an
! overdetermined linear system by Givens rotations, a state-of-art
! algorithm. Routine takes current observation equation from a and
! accumulates it into r:
! Input:
! a(n1)    - an equation with its RHS as a(n1) element,
!            multiplied by sqrt of weight.
!            Note: sum of all weights must equall no.observations,
!            lest GIVENS returns wrong covariance estimates.
! r(n1,n1) - triangle system R and its R.H.S
!            set r = 0 before te first call of GIVENSACC.
! Output:
! r(n1,n1) - updated system accounting of a
!
! tested OK
         Integer :: n1, n, i, j
         Real (SP) :: e (size(a)), p, q, s, rii, ei
         n1 = size (a)
         n = n1 - 1
         e = a
         r (n1, n1) = r (n1, n1) + e (n1) * e (n1)
         Do i = 1, n
            rii = r (i, i)
            ei = e (i)
            s = rii * rii + ei * ei
            If (s <= 0._SP) Then
               p = 1._SP
               q = 0._SP
            Else
               s = sign (Sqrt(s), rii)
               p = rii / s
               q = ei / s
            End If
            r (i, i) = s
            e (i) = 0._SP
            Do j = i + 1, n1
               s = q * e (j) + p * r (i, j)
               e (j) = p * e (j) - q * r (i, j)
               r (i, j) = s
            End Do
         End Do
!
      End Subroutine givensacc
!
      Pure Subroutine lombw (t, vin, w, npar, ncov, fr0, frs, nfr, th)
!
!!Purpose:
!!Lombw - returns Lomb least squares spectrum for weighted unevenly sampled
!!observations. Asuming mean subtracted data the algorithm fits them with
!!a sinusoid f(t)=Asin(om*t)+Bcos(om*t) and returns its weighted squared norm
!! ||f||^2.
!!Reference: Lomb, N.R. 1976, Ap&SS, 39, 447.
!!
!!  Copyrights and Distribution:
!!  This routine is subject to copyrights by its author,
!!  (C) A. Schwarzenberg-Czerny 2018,  alex@camk.edu.pl
!!  Its distribution is free, except that distribution of modifications is
!!  prohibited.
         Real (TIME), Intent (In) :: t (:)!! times of observations
         Real (SP), Intent (In) :: vin (:)!! values of observations
         Real (SP), Intent (In) :: w (:)!! weights of observations
         Integer, Intent (In) :: npar !! dummy
         Integer, Intent (In) :: ncov !! dummy
         Real (TIME), Intent (In) :: fr0 !! start frequency
         Real (TIME), Intent (In) :: frs !! frequency step
         Integer, Intent (In) :: nfr !! number of frequencies
         Real (SP), Intent (Out) :: th (nfr)!! result periodogram (th(1)<0 for errors)
!
         Integer :: l
         Real (SP), Dimension (size(t)) :: c, s, v
         Real (SP) :: ys, yc, ss, cc, cs
         Real (TIME) :: dph (size(t))
!
         l = ncov + npar ! dummy use
         th (1) = 0._SP
         If (size(t) .Lt. 4) th (1) = - 1.
         If (th(1) .Lt. 0._SP) Return
!
         v = vin * w
         Do l = 1, nfr
            dph = (frs*(l-1)+fr0) * t
            s = Real (PI2*(dph-floor(dph)), SP)
            c = Cos (s)
            s = Sin (s)
            yc = sum (c*v)
            ys = sum (s*v)
            cc = sum (w*c**2)
            ss = sum (w*s**2)
            cs = sum (w*c*s)
            th (l) = (ys**2*cc+yc**2*ss-2._SP*yc*ys*cs) / &
               Max (ss*cc-cs**2, epsilon(ss)*(ss**2+cc**2)*0.5_SP)
         End Do
      End Subroutine lombw
!
      Subroutine norm_obs (tin, vin, win, t, v, w, epoch, av, vr)
         Real (TIME), Intent (In) :: tin (:)
         Real (SP), Intent (In) :: vin (:), win (:)
         Real (TIME), Intent (Out) :: t (:)
         Real (SP), Intent (Out) :: v (:), w (:)
         Real (TIME), Intent (Out) :: epoch
         Real (SP), Intent (Out) :: av, vr
!
         Logical :: mask (size(tin))
!
         epoch = (minval(tin)+maxval(tin)) * 0.5
         mask = win > 0._SP
         t = pack (tin, mask)
         v = pack (vin, mask)
         w = 1. / pack (win, mask) ** 2
         av = sum (v*w) / sum (w)
         v = v - av
         vr = sum (v**2*w)
      End Subroutine norm_obs
!
      Subroutine ortint (t, x, no, fr, nh2, a, c)
! (C) Alex Schwarzenberg-Czerny, 1999-2007     alex@camk.edu.pl
!
!   ORTINT - Trigonometric polynomial interpolation from known coefficients
!   NOT intended for general use: asumes normalized data, (i.e. mean subtracted
!   from times and values)
!
!
         Integer, Intent (In) :: no, nh2
         Real (TIME), Intent (In) :: t (no)! times
! NOTE: as t() in ortprj and ortint may differ, is up to user to subtract
! the same EPOCH as used in ortprj, (t_min+t_max)/2
         Real (SP), Intent (Out) :: x (no)! values
         Real (TIME), Intent (In) :: fr ! frequency
!
! recurence & expansion coefficients
         Complex (CP), Intent (In) :: a (0:nh2-1), c (0:nh2-1)
         Integer :: n, nn, nn2
         Real (TIME) :: ph (no)
         Complex (CP), Dimension (no) :: cf, p, z, zn
!
         nn = nh2 / 2
         nn2 = nn + nn
!
         ph = fr * t
         ph = PI2 * (ph-floor(ph))
         z = cmplx (Cos(ph), Sin(ph), kind=CP)
         zn = 1._SP
         p = zn
         cf = 0._SP
         Do n = 0, nn2
            cf = cf + c (n) * p
            p = z * p - a (n) * zn * conjg (p)
            zn = zn * z
         End Do
         ph = - nn * ph
         x = Real (cf*cmplx(Cos(ph), Sin(ph), kind=CP), kind=SP)
      End Subroutine ortint
!
      Subroutine ortprj (t, x, w, nfr, frs, fr0, nh2, a, c, per, no)
!  ORTPRJ - Orthogonal projection onto trigonometric polynomials
!     It may be used to produce a whole periodogram or expansion
!     coefficients for just one (last) frequency. The routine asumes data
!     are normalized so that <x>=0 and <t> is small.
!
!   NOT for general use: this is back-end routine performing no data
!     consistency tests, designed for minimum overhead in calls.
!     Normalization and checks should be performed in the calling
!     driver routine.
!   Input:
!     t[no]-time, should be centered;
!     f[no]-data*rw data weighted by sqrt weight
!     rw[no]-sqrt(w) weight square root (all must be positive-no invalid ones)
!   Output:
!     per[nf]-squared norm of model=sum(cr[i]**2+ci[i]**2) i=0,nh2
!     a(r,i)[nh2+1]-complex polynomial recurrence coefficients
!     c(r,i)[nh2+1]-complex expansion coefficients (for last frequency)
!   Output Power normalization:
!     For FFT data: Per=(N/2)A^2, where A-semi amplitude., N=Size(t)
!
!  (C) Alex Schwarzenberg-Czerny, 1999-2007     alex@camk.edu.pl
!
!   Modiffied: Oct 2010: new formula for alpha
!
!  Input:
         Integer, Intent (In) :: no, nfr, nh2
         Real (TIME), Intent (In) :: t (no)! times
         Real (SP), Intent (In) :: x (no), w (no)
!            x-values, w-weights (all >0)
         Real (TIME), Intent (In) :: fr0, frs ! frequency start & step
!  Output:
         Complex (CP), Intent (Out) :: a (0:nh2-1), c (0:nh2-1)
!            recurence (a) & expansion (c) coefficients for the LAST frequency
         Real (SP), Intent (Out) :: per (:)! periodogram,
!            the frequency range is: <fr0,(Size(per)-1)*frs+fr0>
!
         Logical :: qform
         Integer :: n, nn, nn2, ifr
         Real (SP) :: sn, eps2, eps3
         Real (TIME) :: ph (no)! work array
         Complex (CP), Dimension (no) :: cf, p, z, zn, q ! work arrays
         Complex (CP) :: al, sc
!
         nn = nh2 / 2
         nn2 = nn + nn
         per = - 1._SP
         If (nfr <= 0 .Or. no < nn2+2) Then
            Print *, 'ortprj:error: wrong dimensions'
            Return
         End If
!
         eps2 = Sqrt (epsilon(eps2))
         eps3 = epsilon (eps3) * Sqrt (Real(no, SP))
         per = 0._SP
         Do ifr = 1, nfr ! loop over frequencies
            qform = .True.
            ph = (frs*(ifr-1)+fr0) * t
!   (f,g)=sum(w*conjg(f)*g) -definition
!   z=cmplx(cos(ph),sin(ph))
!   f=z^n*x
!   p_0=1
            ph = PI2 * (ph-floor(ph))
            z = cmplx (Cos(ph), Sin(ph), kind=CP)
            ph = ph * nn
            cf = x * cmplx (Cos(ph), Sin(ph), kind=CP)
            zn = 1._SP
            p = w
!
!   c_k=(p_k,f)/(p_k,p_k)
!   alpha_k=(1,zp_k)/(p_k,p_k)
            Do n = 0, nn2
               q = zn * conjg (p)
               sn = 0._SP
               sn = sum (Real(q, SP))! (1)
! In some parallel environment one may permanently replace (1) with (2)
               If (sn < eps2) Then
                  qform = .False.
                  sn = sum (Abs(p)**2/w)! (2)
               End If
               sn = 1._SP / Max (sn, eps3)
! NOTE: dot_product(a,b)=sum(conjg(a)*b)
               sc = dot_product (p, cf)
               per (ifr) = per (ifr) + Abs (sc) ** 2 * sn
               p = p * z
               al = sum (p) * sn
               p = p - al * q
               zn = zn * z
               a (n) = al
               c (n) = sc * sn
!   norm2=(f,f)=sum(abs(p_k,f)^2/(f,f))
!   p_(k+1)=z*p_k-alpha_k*z^k*conjg(p_k)
            End Do
         End Do
      End Subroutine ortprj
!
      Pure Subroutine powspw (t, vin, w, npar, ncov, fr0, frs, nfr, th)
!
!!Purpose:<br>
!!POWSP - returns discrete power spectrum (DPS) for unevenly sampled
!!observations. For obs.v==1 calculates window function (WF).
!!NPAR and NCOV are ignored. The normalization factor is 2/NOBS,
!!corresponding to Lomb spectrum for FFT sampling.<br>
!!Reference: T.Deeming, 1975, Ap&SS 36, 137.
!!
!!  Copyrights and Distribution:<br>
!!  This routine is subject to copyrights by its author,
!!  (C) A. Schwarzenberg-Czerny 1998-2018,  alex@camk.edu.pl
!!  Its distribution is free, except that distribution of modifications is
!!  prohibited.
         Real (TIME), Intent (In) :: t (:)!! times of observations
         Real (SP), Intent (In) :: vin (:)!! values of observations
         Real (SP), Intent (In) :: w (:)!! dummy
         Integer, Intent (In) :: npar !! dummy
         Integer, Intent (In) :: ncov !! dummy
         Real (TIME), Intent (In) :: fr0 !! start frequency
         Real (TIME), Intent (In) :: frs !! frequency step
         Integer, Intent (In) :: nfr !! number of frequencies
         Real (SP), Intent (Out) :: th (nfr)!! result periodogram (th(1)<0 for errors)
!
         Integer :: l, no
         Real (SP) :: fac, v (size(t))
         Real (TIME) :: ph (size(t))
!
         l = ncov + npar + size (w)! dummy use
         no = size (t)
         th (1) = 0._SP
         If (no .Lt. 4) th (1) = - 1.
         If (th(1) .Lt. 0._SP) Return
!
         v = vin * w
         fac = 2 * no / sum (w) ** 2
!
         Do l = 1, nfr
            ph = (frs*(l-1)+fr0) * t
            ph = PI2 * (ph-floor(ph))
            th (l) = Real (sum(v*Cos(ph))**2+sum(v*Sin(ph))**2, SP) * fac
         End Do
      End Subroutine powspw
!
      Subroutine sortx (dat, ina)
! (C) Alex Schwarzenberg-Czerny, 2005      alex@camk.edu.pl
! Stable index sorting by merge sort. Requires extra memory (inb)
         Implicit None
         Real (SP), Intent (In) :: dat (:)
         Integer, Intent (Out) :: ina (:)
!
         Integer :: inb ((size(dat)+1)/2), i, n
!
         n = size (dat)
         If (size(ina) < n) Then
            Write (*,*) 'sortx: wrong dimension'
            Stop
         End If
!
         Forall (i=1:n) ina (i) = i
!
         Call sort_part (1, n)
!
      Contains
!
         Recursive Subroutine sort_part (low, up)
            Integer, Intent (In) :: low, up
!
            Integer :: med
!
            If (low < up) Then
               med = (low+up) / 2
               Call sort_part (low, med)
               Call sort_part (med+1, up)
               Call merge_part (low, med, up)
            End If
         End Subroutine sort_part
!
         Subroutine merge_part (low, med, up)
            Integer, Intent (In) :: low, med, up
!
            Integer :: i, j, r
!
            inb (1:med-low+1) = ina (low:med)
!
            i = 1
            r = low
            j = med + 1
            Do
               If (r >= j .Or. j > up) Exit
               If (dat(inb(i)) <= dat(ina(j))) Then
                  ina (r) = inb (i)
                  i = i + 1
               Else
                  ina (r) = ina (j)
                  j = j + 1
               End If
               r = r + 1
            End Do
!
            ina (r:j-1) = inb (i:i+j-1-r)
         End Subroutine merge_part
!
      End Subroutine sortx
!
End Module aovsub
!
!
!
!
!
!
