Module aov
!
! Fortran 95 routines for period anlysis to be called from python
! (C) Alex Schwarzenberg-Czerny 1998-2011, alex@camk.edu.pl
!
! for debugging one could compile/run them with:
! f2py -m aov -h aov.pyf aovconst.f90 aovsub.f90 aov.f90
! f2py -m aov -c aovconst.f90 aovsub.f90 aov.f90
!
! Data policy decision: all routines callable from aov.f90
!   assume 3-rd column are errors and invalid ones are negative.
! All aovsub.f90 routines assume 3rd column are weights and none is negative.
!   They shoud be called through aov.f90 routine.
#define OMP
!#undef OMP
#ifdef OMP
      Use omp_lib
#endif
!
      Use aovconst
      Use aovsub
      Implicit None
!
      Integer :: NTHREADS = 99
      Logical :: VERBOSE = .True.
      Private
  ! needed for f2py cross-calls
      Public aovdrv, covar, fgrid, fitcor, fouw, peak, prew, simul, &
     & spectrum, totals
!
Contains
!
      Subroutine aovdrv (method, tin, vin, win, npar, ncov, fr0, frs, &
     & nfr, sptype, th, np0, stat, vr, nobs)
!! AOVDRV driver routine for `AovPer` frequency spectrum routines
!! Purpose: example calling of spectrum routines.
         Integer, Intent (In) :: nobs !! number of observations
         Character (Len=*), Intent (In) :: method !! periodogram method
         Real (TIME), Intent (In) :: tin (nobs)!! times
         Real (SP), Intent (In) :: vin (nobs)!! values
         Real (SP), Intent (In) :: win (nobs)!! weights
         Integer, Intent (In) :: npar !! number of parameters
         Integer, Intent (In) :: ncov !! number of phase coverages
         Real (TIME), Intent (In) :: fr0 !! low frequency
         Real (TIME), Intent (In) :: frs !! frequency step
         Integer, Intent (In) :: nfr !! number of frequences
         Character (Len=*), Intent (In) :: sptype !! periodogram method
         Real (SP), Intent (Out) :: th (nfr)!! frequency spectrum/periodogram
         Integer, Intent (Out) :: np0 !! number of parameters
         Real (SP), Intent (Out) :: vr !! variance of observations
         Integer, Intent (Out) :: stat !! status
!
         Integer :: np, ntfr, ith, i0, nth, lnfr, lflex, iflex, n
         Logical, Save :: first = .True.
         Real (TIME), Allocatable :: t (:)
         Real (SP), Allocatable :: lth (:), v (:), w (:)
         Real (SP) :: av
         Real (TIME) :: lfr0, epoch
!
         ith = 0
         n = count (win > 0._SP)
         Select Case (trim(method))
         Case ('AMHW')
            If (n <= npar/2*2) ith = 1
         Case ('AOVW', 'ATRW')
            If (n < npar*CTMIN) ith = 1
         Case ('F_MW')
            If (n <= 3) ith = 1
         Case ('LOMW')
            If (n <= 2) ith = 1
         Case ('PSPW')
            If (n <= 1) ith = 1
         Case Default
            stat = 4
            If (VERBOSE) Print *, 'aovdrv: Wrong METHOD ' // trim &
           & (method)
            Return
         End Select
         If (ith /= 0) Then
            stat = 4
            If (VERBOSE) Print *, 'aovdrv: too few observations or' // &
           & ' npar*CTMIN large in ' // trim (method)
            Return
         End If
         Allocate (t(n), v(n), w(n))
         Call norm_obs (tin, vin, win, t, v, w, epoch, av, vr)
!
         stat = 0
         nth = 0
         ith = 0
         iflex = 0
         ntfr = nfr
         np = Max (2, npar)
! If no OpenMP remove from here ...
#ifdef OMP
         If (NTHREADS > 0) nth = Min (omp_get_num_procs(), NTHREADS)
         Call omp_set_num_threads (nth)
         ntfr = ntfr / Max (nth, 1) + 1
!$omp parallel  private(np,ith,i0,lfr0,lnfr,lth,lflex) reduction(+:iflex)
         ith = omp_get_thread_num ()
! ... to there
#endif
         If (first .And. ith == 0) Then
            Print *, "Aovdrv: number of OMP threads = ", nth
            first = .False.
         End If
         i0 = ith * ntfr
         lnfr = Min (ntfr, nfr-i0)
         lfr0 = fr0 + i0 * frs
         lflex = 0
         Allocate (lth(lnfr))
! get spectrum/periodogram
! all routines assume epoch subtracted from times, all weights positive and 
! weighted mean subtracted from values 
         Select Case (trim(method))
         Case ('AOVW')
            Call aovw (t, v, w, npar, ncov, lfr0, frs, lnfr, lth, &
           & lflex)
            If (ith == 0) np = npar
         Case ('ATRW')
            Call aovw (t, v, w,-npar, ncov, lfr0, frs, lnfr, lth, &
           & lflex)
            If (ith == 0) np = 2
         Case ('AMHW')
            np = npar / 2 * 2 + 1
            Call aovmhw (t, v, w, np, ncov, lfr0, frs, lnfr, lth)
            If (lth(1) < 0.) Print *, "error:", ith, lth (1)
         Case ('F_MW')
            Call ferraz_glsw (t, v, w, np, ncov, lfr0, frs, lnfr, lth)
            If (ith == 0) np = 3
         Case ('LOMW')
            Call lombw (t, v, w, np, ncov, lfr0, frs, lnfr, lth)
            If (ith == 0) np = 2
         Case ('PSPW')
            Call powspw (t, v, w, np, ncov, lfr0, frs, lnfr, lth)
            If (ith == 0) np = 2
         Case Default
            If (ith == 0) stat = 4
            If (VERBOSE .And. ith == 0) Print *, 'aovdrv: warning: Wron&
           &g METHOD, ignored'
         End Select
         If (ith == 0) np0 = np
         iflex = iflex + lflex
         If (lth(1) < 0._SP) stat = Max (stat, Int(Abs(lth(1))+0.5_SP))
         th (i0+1:i0+lnfr) = lth
         Deallocate (lth)
#ifdef OMP
!$omp end parallel
#endif
         If (VERBOSE .And. iflex > 0) Print *, 'AOV/ATR: warning: poor &
            &phase coverage at ', iflex, ' frequencies'
         If (VERBOSE .And. stat /= 0) Print *, 'aovdrv: warning: wrong &
            &periodogram parameters, ignored'
         Call spectrum (size(t), np0, vr, sptype, th, nfr)
      Contains
!
      End Subroutine aovdrv
!
      Subroutine covar (t1, v1, e1, t2, v2, e2, eps, iscale, ifunct, &
         lav, lmi, lmx, cc, cmi, cmx, ilag, n1, n2, nlag)
! TO BE DONE:
! -more flexibility in choosing the lag grid
! Purpose:
!   Calculates cross-correlation function(CCF) of two time series sampled
!   unevenly in time. Note that for good result each hump in the time series
!   should be covered by several observations (over-sampling).
!
         Integer, Intent (In) :: n1, n2, nlag
         Integer, Intent (In) :: iscale, ifunct
         Real (SP), Intent (In) :: eps, v1 (n1), e1 (n1), v2 (n2), e2 (n2)
         Real (TIME), Intent (In) :: t1 (n1), t2 (n2)
         Real (SP), Intent (Out) :: cc (nlag), cmi (nlag), cmx (nlag), &
            lav (nlag), lmi (nlag), lmx (nlag)
         Integer, Intent (Out) :: ilag
! f2py Integer,Intent(IN),optional   :: iscale=0,ifunct=0
! f2py Real(SP),Intent(IN),optional  :: eps=0.
!
!   Input:
!   n1      - number of observations in 1st set
!   n2      - number of observations in 2nd set
!   t1,v1,e1(n1)  - 1st set of observations: time, value, error
!   t2,v2,e2(n2)  - 2nd set of observations: time, value, error,
!          negative errors indicate invalid observations
!   iscale  - (ignored)output scale: iscale/=1 -linear; iscale=1 -logarythmic
!   ifunct  - (ignored)output function: ifunct/=1 -correlation; ifunct=1 -structure
!   nlag    - Size of lag arrays
!
!   Output:
!   lav,lmi,lmx(nlag)- average & range of lags and ...
!   cc,cmi,cmx(nlag)-       ... average & range of correlations
!   ilag    - actual number of computed lags
!
! Method/Reference:
! Alexander, T., 1997, Astronomical Time Series,
! Eds. D. Maoz, A. Sternberg, and E.M. Leibowitz,
! (Dordrecht: Kluwer),218, 163
!
!  Copyrights and Distribution:
!  This package is subject to copyrights by its author,
!  (C) A. Schwarzenberg-Czerny 2011,  alex@camk.edu.pl
!
! n1*n2 could be large dimension
         Integer :: ind (n1*n2), mdt, mct, idt, ict, i, i1, i2, nct, &
            ib1 (2*n1*n2/nlag), ib2 (2*n1*n2/nlag)
         Real (SP) :: dt (n1*n2), e, z, s, r, r2, ad1, ad2, sd, sw
         Real (SP), Allocatable :: db1 (:), db2 (:), wb1 (:), wb2 (:), &
            wb (:), cb (:), lags (:)
!
         mdt = iscale + ifunct ! use dummy variables to prevent warnings
         mdt = n1 * n2
         nct = mdt / nlag
         mct = 2 * nct
         lav = 0._SP
         lmi = 0._SP
         lmx = 0._SP
         cc = 0._SP
         cmi = 0._SP
         cmx = 0._SP
!
         Forall (i1=1:n1, i2=1:n2) dt ((i1-1)*n2+i2) = t2 (i2) - t1 &
        & (i1)
         Call sortx (dt, ind)
!
         idt = 1
         ilag = 1
         Do ! increment bins
            ict = 0
            If (idt > mdt .Or. ilag > nlag) Exit
!
            Do ! increment pairs
               If (idt > mdt) Exit
! check If lag bin finished
               If (ict > nct .And. dt(ind(idt))-dt(ind(Max(1, idt-1))) &
                  > eps .And. mdt-idt > nct/2) Exit
               i1 = (ind(idt)-1) / n2 + 1
               i2 = ind (idt) - (i1-1) * n2
               e = e1 (i1) * e2 (i2)
               If (e > 0._SP) Then
                  Do i = 1, ict ! check dependence
                     If (ib1(i) == i1 .Or. ib2(i) == i2) Exit
                  End Do
                  If (i <= ict) Then ! dependent
                     If (e < e1(ib1(i))*e2(ib2(i))) Then
                        ib1 (i) = i1
                        ib2 (i) = i2
                     End If
                  Else ! independent
                     If (ict >= mct) Exit! this should never happen
                     ict = ict + 1
                     ib1 (ict) = i1
                     ib2 (ict) = i2
                  End If
               End If
               idt = idt + 1
            End Do
!
! evaluate correlation within a bin
            If (ict < 2) Cycle
            Allocate (db1(ict), db2(ict), wb1(ict), wb2(ict), wb(ict), &
               cb(ict), lags(ict))
            db1 = v1 (ib1(:ict))
            wb1 = 0._SP
            Where (e1(ib1(:ict)) > 0._SP) wb1 = 1._SP / e1 (ib1(:ict)) ** 2
            db2 = v2 (ib2(:ict))
            wb2 = 0._SP
            Where (e2(ib2(:ict)) > 0._SP) wb2 = 1._SP / e2 (ib2(:ict)) ** 2
            wb = Sqrt (wb1*wb2)
            sw = sum (wb)
!
            ad1 = sum (db1*wb1) / sum (wb1)
            db1 = db1 - ad1 ! data variance and averages
            ad2 = sum (db2*wb2) / sum (wb2)
            db2 = db2 - ad2
            sd = sum (db1**2*wb1) * sum (db2**2*wb2)
            r = sum (db1*db2*wb)
!
            lags = t2 (ib2(:ict)) - t1 (ib1(:ict))! average & extreme intervals
            lav (ilag) = sum (lags*wb) / sw
            lmx (ilag) = maxval (lags)
            lmi (ilag) = minval (lags)
            Deallocate (db1, db2, wb1, wb2, wb, cb, lags)
!
! correlation and its variance in z-variable
            If (sd <= 0._SP) Cycle
            r = r / Sqrt (sd)
            cc (ilag) = r
            r2 = r * r
            z = (((((r2*3._SP+2._SP)*r2+11._SP)*0.5_SP/(ict-1)+r2+ &
               5._SP)*0.25/(ict-1)+1._SP)*r/(ict-1)+Log(Max(1._SP-r, &
               tiny(z))/(1._SP+r))) * 0.5_SP
            s = Sqrt (((((-r2*3._SP-6._SP)*r2+22._SP)/3._SP/(ict-1)+ &
               4._SP-r2)*0.5/(ict-1)+1._SP)/(ict-1))
            cmi (ilag) = r - Abs (Tanh(z-s)-r)
            cmx (ilag) = r + Abs (Tanh(z+s)+r)
!
            ilag = ilag + 1
         End Do
         ilag = ilag - 1
      End Subroutine covar
!
      Subroutine fgrid (t, fstop, fstep, fr0, nobs)
! Selects an acceptable frequency grid.
! Note: fgrid accounts also for times of invalid observations
!    to maintain the same output for similar batches of data.
         Integer, Intent (In) :: nobs
         Real (TIME), Intent (In) :: t (nobs)
         Real (TIME), Intent (Out) :: fr0, fstep, fstop
!  Copyrights and Distribution:
!  This package is subject to copyrights by its author,
!  (C) A. Schwarzenberg-Czerny 1998-2011,  alex@camk.edu.pl
!
         Integer :: ind (nobs), nobs1, iobs, nfr, nh2
         Real (SP) :: del (nobs-1)
         Real (TIME) :: steps, span
!
         Integer, Parameter :: MAXNFR = 1000000
         Character (Len=*), Parameter :: msg1 (12) = (/ &
'FGRID: *** DANGER *** Data spacing too fine or span too long for good sampling',&
'of periodogramme. You may choose to preset shorter or coarser frequency grid  ',&
'by first loading shorter span data and choosing frequency range and next use  ',& 
'FREEZE & OPEN to load the big data.                                           ',&
'   You may also choose to modify the data: if high frequencies do not matter, ',&
'data could be binned into coarser time grid. Alternatively, you may split data',& 
'into shorter intervals and take average of several periodograms, thus keeping ',&
'the frequency range while degrading resolution.                               ',&
'   This kind of problems are fundamental in time-frequency analysis and depend',&
'on no algorithm. Violation of the uncertainity principle Dt*Df<1, where either',&
'Dt-time span and Df-frequency resolution or Dt-resolution and Df-span yields  ',&
'respectively undersampling or clipping of (effective) Nyquist interval        ' /)
         nh2 = 3
         nobs1 = nobs - 1
         If (nobs1 < 5) Then
            Print *, 'FGRID: Too few observations'
            Return
         End If
!
         Call sortx (Real(t-t(nobs/2+1),SP), ind)
         span = t (ind(nobs)) - t (ind(1))
         If (span <= 0._SP) Then
            Print *, 'FGRID: Input TIMES differ too little'
            Return
         End If
         fstep = 0.3_SP / span / real (nh2/2, kind=SP)
!
         del = t (ind(2:)) - t (ind(:nobs1))
         Call sortx (del, ind(:nobs1))
!
         fstop = Log (real(nobs, kind=SP))
         iobs = Nint (nobs1*(1._SP/(6._SP+0.3_SP*fstop)+0.05_SP)) + 1
         fstop = del (ind(iobs))
         If (fstop <= 0._SP) Then
            Print *, 'FGRID: Observations too finely spaced, ' // &
           'bin them coarsly, if feasible'
            Return
         End If
!
         fstop = 1._SP / fstop * (span/nobs1/fstop) ** (-0.6_SP)
!  fstop=0.5_SP/fstop*(del(ind(nobs1/2))/fstop)**(-0.6_SP)
         steps = fstop / fstep
         If (steps >= MAXNFR) Then
            steps = MAXNFR
            Print '((a))', msg1
         End If
!
         nfr = Nint (steps)
         fstep = fstop / (nfr-1)
         fr0 = 0._SP
         Print *, ' RESULTS OF FREQUENCY BAND EVALUATION:'
         Print '(3(A,1PE10.1))', 'Max. Freq.: ', fstop, '  Resolution: ', &
            fstep, '  Min. Freq.: ', fr0
         Print '(A,i10)', 'No. of points:  ', nfr
      End Subroutine fgrid
!
      Subroutine fitcor (cn, co, cp, chi2, npar1, no)
!
         Integer, Intent (In) :: no, npar1
         Real (SP), Intent (In) :: cn (npar1, no)
         Real (SP), Intent (In) :: co (no+1, no+1)
         Real (SP), Intent (Out) :: cp (npar1, npar1)
         Real (SP), Intent (Out) :: chi2
!
!  Solves linear least sqares problem with weighting by a given full
!  covariance matrix of obserations. This routine applies to observations
!  mutually correlated, with covariance matrix known.
!  Note: obsolete algorithm, kept for testing any upgrade amploying Givens rotations
!
!  Input:
!  co(no,no)         - covariance matrix of observations (to be destroyed)
!  cn(npar1,no)      - condition equations
!  no                - number of observations
!  npar1             - number of parameters+1 (1 for calculation of chi2 only)
!
!  Output:
!  chi2              - chi2
!  cp(npar1,npar1)   - fitted parameters and their covariance matrix
!                        (See CRACOW for details)
!
!  Method:
!  First finds Cholesky triangle root Q of covariance matrix C by a call of CRACOW.
!  Next multiply Q by condition equations. Square the result to get the normal
!  equations. Solve them using  CRACOW and return values of the parameters
!  and their covariance matrix.
!
!  This package is subject to copyrights by its author,
!  (C) A. Schwarzenberg-Czerny 1992-2011,  alex@camk.edu.pl
         Integer :: i, j
         Real (SP) :: cn1 (npar1, no)
         Real (SP) :: co1 (no+1, no+1)
!
         chi2 = 0._SP
!  Find inverse triangle root Q of the observations covariance matrix
         co1 (:no, :no) = co (:no, :no)
         If (cracow(co1, 0) == 0._SP) Then
            Print *, ' Covariance matrix of observations is singular (empty?)'
            Return
         End If
! Multiply condition equations by Q
         Do i = no, 1, - 1
            Do j = 1, npar1
               cn1 (j, i) = dot_product (cn(j, 1:i), co1(1:i, i))
            End Do
         End Do
! Square the result to get normal equations
         Do i = 1, npar1
            Do j = i, npar1
               cp (j, i) = dot_product (cn1(i, 1:no), cn1(j, 1:no))
            End Do
         End Do
!  Fit parameters by LSQ, return them and their covariance matrix and chi2
         If (cracow(cp, no) == 0._SP) Then
            Print *, ' Normal equations are singular, modify trend functions'
            Return
         End If
         chi2 = cp (npar1, npar1) * (no-npar1+1)
      End Subroutine fitcor
!
      Subroutine fouw (t, valin, er, frin, nh2, frout, dfrout, valout, &
         cof, dcof, nobs)
!
! FOUW - Fit Fourier series (cof) and Return vector of residuals
!          (fout), adjusting frequency (fr) by NLSQ
! Note: FOUV ignores invalid observations marked with er<0
!
!  Purpose:
!   Twofold: (i) Find Fourier coefficients and improve frequency
!   (ii) Pre-whiten data with a given frequency
!
!  Method:
!   Fits Fourier series of nh2/2 harmonics with frequency adjustment
!   NLSQ by Newton-Raphson iterations:
!   fin(ph)=cof(1)+sum_{n=1}^{nh2/2}(cof(2n-1)Cos(ph*n)+cof(2n)*Sin(ph*n))
!   where ph=2pi*frout*(t-t0) and t0=(max(t)+min(t))/2
!  Input:
!   t(:),fin(:),er(:) -times, values and weights of observations;
!   frin- initial value of frequency is abs(fr)
!       for frin<0 on Return frout=frin
!       for frin>0 value of frout is adjusted by NLSQ iterations
!   nh2- number of model parameters: nh2/2-no. of harmonics
!  Output:
!   frout-final value of frequency (see frin)
!   cof(:),dcof(:)-Fourier coefficients and errors
!      cof(nh2+2)-approximate t0 (for full accuracy use t0=(max(t)+min(t))/2)
!      dcof(nh2+2)-frequency error (0 for fr<0)
!   fout(:)- residuals from fit/prewhitened data
!
!  Copyrights and Distribution:
!  This package is subject to copyrights by its author,
!  (C) A. Schwarzenberg-Czerny 1998-2011,  alex@camk.edu.pl
!  Please quote A.Schwarzenberg-Czerny, 1995, Astr. & Astroph. Suppl, 110,405 ;
!
         Integer, Parameter :: ITMAX = 100
         Integer, Intent (In) :: nobs, nh2
         Real (SP), Intent (In) :: valin (nobs), er (nobs)
         Real (TIME), Intent (In) :: t (nobs), frin
         Real (TIME), Intent (Out) :: frout, dfrout
         Real (SP), Intent (Out) :: valout (nobs), cof (nh2/2*2+2), &
            dcof (nh2/2*2+2)
!
         Logical :: finish
         Integer :: l, n, nn2, nh, nx, nx1, it, idf, no
         Real (SP), Allocatable :: e (:), a (:, :)
         Real (SP) :: cond, detnth, rhs, sig, avstep
         Real (TIME) :: ph, t0, dfr
!
         nh = Max (1, nh2/2)
         nn2 = nh + nh
!
         If (nobs < nh+nh+2 .Or. size(valin) /= nobs .Or. &
               size(er) /= nobs .Or. size(valout) /= nobs) Then
            Write (*,*) 'FOUW:error: wrong Size of arrays'
            Return
         End If
!
         t0 = (maxval(t)+minval(t)) * 0.5_TIME
         dfr = 1. / (maxval(t)-minval(t))
!
         cof = 0._SP
         dcof = 0._SP
         dfrout = 0._TIME
         frout = Abs (frin)
         Print *, '  it    sigma   step    cond'
         Do it = 1, ITMAX
            nx = nn2 + 2
            nx1 = nx ! adjust only coefficients
            If (it > 1 .And. frin > 0._TIME) nx1 = nx + 1 ! adjust also frequency
            Allocate (e(nx1), a(nx1, nx1))
            a = 0._SP
            no = 0
            Do l = 1, nobs
               ph = frout * (t(l)-t0)
               ph = PI2 * (ph-floor(ph)) * 0.5_SP
               rhs = valin (l) - cof (1)
               e = 0._SP
               e (1) = 1._SP
!
               Do n = 2, nn2, 2
                  e (n) = Cos (ph*n)
                  e (n+1) = Sin (ph*n)
                  e (nx) = e (nx) + (-e(n+1)*cof(n)+e(n)*cof(n+1)) * n
                  rhs = rhs - e (n) * cof (n) - e (n+1) * cof (n+1)
               End Do
!
               e (nx) = e (nx) * PI2 * 0.5_SP * (t(l)-t0)
               e (nx1) = rhs
               valout (l) = rhs
               If (er(l) > 0._SP) Then
                  no = no + 1
                  e = e / er (l)
                  Call givensacc (e, a)
               End If
            End Do
            idf = no - nx1 + 1
            If (idf < 1) Then
               Print *, "Fouw: too few valid observations for requested series"
               Return
            End If
            finish = (frin < 0._TIME .And. it >= 2)
!
            If ( .Not. finish) Then
               cond = givens (a, idf, detnth)
               sig = Sqrt (a(nx1, nx1))
               avstep = Sqrt (sum(a(nx1, :nx1-1)**2)/(nx1-1))
               cof (nx) = 0._SP
               cof (:nx1-1) = cof (:nx1-1) + a (nx1, :nx1-1)
               If (it > 1 .And. frin > 0._SP) frout = frout + cof (nx)
               Do l = 1, nx1 - 1
                  dcof (l) = Sqrt (Max(a(l, l), tiny(a)))
               End Do
               finish = Abs (cof(nx)) > dfr .Or. Abs (cof(nx)) < dcof(nx) * 0.01
            End If
            If (Nint(Exp(Nint(2.5*Log(real(it)))/2.5)) == it .Or. finish) &
               Print '(i5,f10.5,2(1pe8.1))', it, sig, avstep, cond
            If (finish) Exit
            Deallocate (e, a)
         End Do
         cof (nx) = t0
         If (it > 1 .And. frin > 0._TIME) dfrout = dcof (nx)
!
         Print *, 'frout=', frout, '+/-', dfrout
         Do l = 1, nx - 1
            Print '(2(f10.5,a))', cof (l), ' +/-', dcof (l)
         End Do
         Print *, 'Epoch=', t0
      End Subroutine fouw
!
      Subroutine peak (n, fx, xm, fm, dx)
! Scans periodogram for a peak and retutns its parameters
!
! NOTE: peak follows Fortran indexing convention
! For C and Python subtract 1 from xm output
! INPUT:
! n- periodogram length
! fx-periodogram values
! OUTPUT:
! xm-peak location
! fm-peak value
! dx-peak halfwidth
! SPECIAL CONDITIONS:
! fm<0 no peak in valid region
! dx<0 no valid (parabolic) peak
! METHOD:
! in the log(f) space
! fits f(x)=f(0)+f'*x+0.5*f"*x^2 where f'=0.5*(f(+1)-f(-1)) and
! f"=f(+1)+f(-1)-2*f(0)
! x(max)=f'/(-f'')     f(max)=f(0)+0.5*(f')^2/(-f")
! f(max)-b = f(max)+0.5*(-f")(dx)^2 - parabolic error estimate
! dx = sqrt(2*b/(-f")), however in log units:
! b(log) = -log( max( 1-b/exp(f(max)), 0.5), where 0.5 for HWHM
!(C) Alex Schwarzenberg-Czerny, 1999-2005 alex@camk.edu.pl
!
         Integer, Intent (In) :: n
         Real (SP), Intent (In) :: fx (n)
         Real (TIME), Intent (Out) :: xm
         Real (SP), Intent (Out) :: fm
         Real (SP), Intent (Out) :: dx
!
         Integer :: nm, n0
         Real (SP) :: f (3), b, a2, dm
!
! valid maximum?
         nm = maxloc (fx, dim=1)
         fm = fx (nm)
         xm = nm
         nm = Min (n, Max(3, nm+1))
         dx = - 1._SP
         If (fm <= 0._SP .Or. n < 3 .Or. minval(fx(nm-2:nm)) <= 0._SP) Return
! on edge ?
         f = Log (fx(nm-2:nm))
         b = (f(3)-f(1)) * 0.5_SP
         a2 = - (f(1)-2._SP*f(2)+f(3))
         If (a2 <= 0._SP) Return
! parabolic case
         dm = b / a2
         If (Abs(dm) > 1._SP) Return! excess step
         fm = Exp (0.5_SP*b*b/a2+f(2))
         xm = dm + nm - 1._SP
         dx = Sqrt (2._SP*Log(2._SP)/a2)
      End Subroutine peak
!
      Subroutine prew (t, valin, er, frin, nh2, frout, dfrout, valout, nobs)
!
! PREW - Fit trig orthogonal polynomials and Return vector of residuals
!          (fout), adjusted frequency (fr) and its error (dfr)
!
!  Purpose:
!  (i) Pre-whiten data with a given frequency
!  (ii) Fit data with a smooth periodic function
!
!  Method:
!  Finds frequency of the nearby peak of aovmh periodogram and by orthogonal
!  projection expands data into series of trig orthogonal polynomials
!
!  Input:
!   t(:),fin(:),er(:) -times, values and weights of observations;
!   frin- initial value of frequency is abs(fr)
!       for frin<0 on Return frout=frin
!       for frin>0 value of frout is adjusted by NLSQ iterations
!   nh2- number of model parameters: nh2/2-no. of harmonics
!
!  Output:
!   frout-final value of frequency (see frin)
!   cof(:),dcof(:)-Fourier coefficients and errors
!      cof(nh2+2)-approximate t0 (for full accuracy use t0=(max(t)+min(t))/2)
!      dcof(nh2+2)-frequency error (0 for fr<0)
!   fout(:)- residuals from fit/prewhitened data
!
!  Copyrights and Distribution:
!  This package is for free distribution, subject to copyrights by its author,
!  A. Schwarzenberg-Czerny 1998-2011,  alex@camk.edu.pl
!
!  Please quote A.Schwarzenberg-Czerny, 1995, Astr. & Astroph. Suppl, 110,405 ;
!
         Integer, Parameter :: ITMAX = 100, IBAND = 2
         Integer, Intent (In) :: nobs, nh2
         Real (SP), Intent (In) :: valin (nobs), er (nobs)
         Real (TIME), Intent (In) :: t (nobs), frin
         Real (TIME), Intent (Out) :: frout, dfrout
         Real (SP), Intent (Out) :: valout (nobs)
!
         Logical :: finish, mask (nobs)
         Integer :: nn2, nh, it, n, nfr
         Real (SP) :: per (IBAND*2+1), fm, dx, tdf, tfr, rw (nobs), av, sw
         Real (TIME) :: t0, dfr, frmax, fr0
         Complex (cp) :: a (0:nh2/2*2), c (0:nh2/2*2)
!
         nh = nh2 / 2
         nn2 = nh + nh
         nfr = IBAND * 2 + 1
!
         If (nh < 1 .Or. nobs < nh+nh+2 .Or. size(valin) /= nobs .Or. &
               size(er) /= nobs .Or. size(valout) /= nobs) Then
            Write (*,*) 'PREW:error: wrong nh2 or Size of arrays'
            Return
         End If
!
         Call fgrid (t, frout, dfr, fr0, nobs)
         frout = Abs (frin)
!
         t0 = (maxval(t)+minval(t)) * 0.5_TIME
         mask = er > 0._SP
         n = count (mask)
         Where (mask)
            rw = 1._SP / er
         Elsewhere
            rw = 0.
         End Where
         sw = sum (rw**2)
         av = sum (valin*rw**2, mask=mask) / sw
         valout = (valin-av)
!
         If (frin > 0._TIME) Then
            Print *, '  it  dwidth  dfreq'
! Loop to improve frequency
            finish = .True.
            Do it = 1, ITMAX
! Get periodogram peak
               frmax = frout
               fr0 = frout - dfr * IBAND
               Call ortprj (pack(t-t0, mask),pack(valout, mask),pack(rw, mask), &
                  nfr, dfr, fr0, nh2, a, c, per, n)
               Call peak (nfr, per, frout, fm, dx)
               dx = 0.5_SP * dx * dfr
               frout = (frout-1._TIME) * dfr + fr0
!
! Update frequency
               frout = (frout+frmax) * 0.5_TIME ! slowly correct frequency
!
! Check finish
               tfr = real (Max(dfr, tiny(dfr)), kind=SP)
               tdf = real (Abs(dfr-dx), kind=SP) / tfr
               tfr = real (Abs(frmax-frout), kind=SP) / tfr
               finish = (tdf < 0.01 .And. tfr < 0.01) .Or. frin < 0._TIME
               If (Nint(Exp(Nint(2.5*Log(real(it)))/2.5)) == it .Or. finish) &
                  Print '(i5,4(1pe8.1))', it, tdf, tfr
               If (finish) Exit
               dfr = dx
            End Do
!
! Final peak evaluation
            If ( .Not. finish) &
               Print *, "prew: poor convergence of frequency"
            Call peak (IBAND*2+1, per, frout, fm, dx)
            frout = (frout-1._TIME) * dfr + fr0 ! return best frequency  ...
            dfrout = dx * dfr ! ...and error accounting for correlation
         End If ! End of frequency adjustment
!
! Subtract fitted polynomial series
         ! valout as input
         Call ortprj (pack(t-t0, mask), pack(valout, mask), pack(rw, mask), &
            1, dfr, frout, nh2, a, c, per, n)
         Call ortint (t-t0, valout, nobs, frout, nh2, a, c)
         valout = valin - av - valout
!
         Print *, 'frout=', frout, '+/-', dfrout
      End Subroutine prew
!
      Subroutine simul (no, s2n, t, f, er)
!   Simulated light curve
         Integer, Intent (In) :: no
         Real (SP), Intent (In) :: s2n
         Real (TIME), Intent (Out) :: t (no)
         Real (SP), Intent (Out) :: f (no), er (no)
!
         Integer k
         Real (SP) span, omt, r (6)
         Real (TIME) fr, dph
!
         Integer, Allocatable :: s (:)
!
!   start the same random number sequence
         Call random_seed (size=k)
         Allocate (s(k))
         s = 0
         s (1) = 1949
         Call random_seed (put=s)
!
!   simulate data
         span = 40._SP
         fr = 1._TIME / 0.3617869_TIME
         Write (*,*) 'Simulated data, FR0=', fr
         Call random_number (er)
         er = Exp (-er)
         Forall (k=0:no-1) t (k) = span * Sin (k*1._TIME) ** 4.
         t (1) = 0._TIME
         t (no) = span
         Do k = 1, no
            dph = fr * t (k)
            omt = PI2 * (dph-floor(dph))
            Call random_number (r)
            f (k) = er (k) * (sum(r)/6.-3.) / s2n + 11.90_SP + &
           0.434_SP * Cos (omt+4.72_SP) + 0.237_SP * Cos (2._SP*omt+0.741_SP)
         End Do
      End Subroutine simul
!
      Subroutine spectrum (noeff, npar, vr, sptype, th, nfr)
!! Purpose: <br>
!! Converts quadratic norm of fitted model into
!! one of frequency spectrum classic Fisher statistics
         Integer, Intent (In) :: noeff !! number of valid observations
         Integer, Intent (In) :: npar !! number of model parameters
         Integer, Intent (In) :: nfr !! number of frequences
         Real (SP), Intent (In) :: vr !! weighted sum of squares
         Character (Len=*), Intent (In) :: sptype !! type of returned spectrum
         	                             !! 'RAW', 'AOV', 'MODEL', 'RESID'
         Real (SP), Intent (Inout) :: th (nfr)!! periodogram/spectrum
!
         Real (SP) :: fac
         Select Case (trim(sptype))! SP_RAW,SP_AOV,SP_MODEL,SP_RESID
         Case ('RAW')! do nothing
            Continue
         Case ('AOV')! convert to F(npar,no-npar;th)
            fac = (noeff-npar-1) / npar
            th = th / Max (vr-th, epsilon(fac)*vr) * fac
         Case ('MODEL')! convert to ibeta(npar/2,(no-npar)/2;th)
            th = th / vr
         Case ('RESID')! convert to ibeta((no-npar)/2,npar/2;th)
            th = 1._SP - th / vr
         Case Default
            If (VERBOSE) Print *, "Spectrum: wrong SPTYPE"
            Return
         End Select
!
      End Subroutine spectrum
!
      Subroutine totals (xin, e)
! evaluate general properties of x
! note: quantiles computed in a crude way
         Real (SP), Intent (In) :: xin (:), e (:)
!  Copyrights and Distribution:
!  This package is subject to copyrights by its author,
!  (C) A. Schwarzenberg-Czerny 1998-2011,  alex@camk.edu.pl
!
         Real (SP) :: fmax, fmin, fmean, fdev, fq25, fq75, fmed, sw, &
        & wss, wa, wb, a, x (size(xin)), w(size(xin)), ws(size(xin))
         Integer :: ind (size(xin)), signs, n, i
         Character (Len=10) :: tstr
!
         n = count (e > 0._SP)
         If (n <= 0) Then
            Print *, "Totals: No valid observations"
            Return
         End If
         x (:n) = pack (xin, e > 0._SP)
         w (:n) = 1._SP / pack (e, e > 0._SP) ** 2
         signs = 1
         If (n == 1) Then
            fmean = x (1)
            fdev = 0._SP
            fmax = fmean
            fmin = fmean
            fmed = fmean
            fq25 = fmean
            fq75 = fmean
         Else If (n > 1) Then
            sw = sum (w(:n))
            fmean = sum (x(:n)*w(:n)) / sw
            fdev = 0._SP
            x (:n) = x (:n) - fmean
            signs = count (x(:n-1)*x(2:n) <= 0._SP) + 1
            fdev = Max (sum((x(:n))**2*w(:n))*n/(n-1)/sw, 0._SP)
            If (fdev > 0._SP) fdev = Sqrt (fdev)
!
! calculate approx. quantiles by interpolation of
! a roughly binned histogram
            Call sortx (x(:n), ind(:n))
            wss = 0._SP
            Do i = 1, n
               wss = wss + w (ind(i))
               ws (ind(i)) = wss
            End Do
            ws = ws / sw
            fmax = x (ind(n))
            fmin = x (ind(1))
            fmed = fq (0.5_SP, x(ind(:n)), ws(ind(:n)))
            fq25 = fq (0.25_SP, x(ind(:n)), ws(ind(:n)))
            fq75 = fq (0.75_SP, x(ind(:n)), ws(ind(:n)))
         End If
!
         Print '(4(a,i8))', ' length=', n, ' signs=', signs
         Print '(3(a,1pg16.6))', ' min=', fmin, ' max=', fmax, ' mean=', fmean
         Print '(3(a,1pg16.6))', ' q25=', fq25, ' q75=', fq75, ' median=', fmed
         Call date_and_time (TIME=tstr)! standard system routine
         Print '(1(a,1pg16.6),a)', ' std.dev=', fdev, ' time= ' // &
            tstr (1:2) // ':' // tstr (3:4) // ':' // tstr (5:6)
      Contains
         Function fq (q, x, ws)
            Real (SP) :: fq
            Real (SP), Intent (In) :: q, x (:), ws (:)
!
            Integer :: i, n
            n = size (x)
            i = Max (count(ws < q), 1)
            a = (ws(i)-q) / (w(i+1)-ws(i))
            fq = x (i) * (1._SP-a) + x (i+1) * a
         End Function fq
!
      End Subroutine totals
!
End Module aov
!
!
!
!
!
!
