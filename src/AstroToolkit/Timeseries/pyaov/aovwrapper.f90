program aovwrapper
  Use aovconst
  Use aovsub
  use aov
  implicit none
  character (len=16) :: nam
  CHARACTER (LEN=400) :: data
  character (len=100) :: fname
  character (len=100) :: return_fname
  character (len=16) :: sptype 
  character (len=4) :: method 
  real :: mat(10, 20)
  Integer :: nrows,nrows2
  Integer :: ncols 
  Integer :: npar
  Integer :: ncov
  Real (TIME) :: lfr0
  Real (TIME) :: frs,fstop,fstep,fr0,frin,frout, dfrout,xm
  Integer :: lnfr 
  !Integer, Parameter :: n=1113
  !Integer, Parameter :: npar=3
  !Integer, Parameter :: ncov=2
  !Real (TIME), Parameter :: lfr0=2.0
  !Real (TIME), Parameter :: frs=(60.0-2.0)/150000.
  !Integer, Parameter :: lnfr =150000

  Real (SP), allocatable :: lth(:)  
  integer :: ios,expid,ioid,i,lflex,j,len,status,np0,stat,nparameters,narrays,nh2
  real (kind=16) oid,hjd,mjd,mag,magerr,total
  real (sp) :: vr,fm,dx		

  Real (TIME), allocatable  :: t(:),t2(:) !! times of observations
  Real (SP) , allocatable:: v(:),v2(:) !! values of observations
  Real (SP) , allocatable :: w(:),w2(:) !! weights of observations


  Real (SP), allocatable :: valout (:),cof(:),dcof(:)

! Return file format is:
!    Rec1  nparameters, narrays
!    rec2 - rec(2+nparameters-1) values of parameters
!    for each array the first record holds the length of teh array and subsequent records the values in th array

  call get_command_argument(1,fname, len, status)
  Open(15,file=fname)
  Read(15,*,End=1) nam  !  name of routine to invoke
  Read(15,*,End=1) return_fname  !  name of file for output data

  select case (trim(nam))

  Case ('fgrid')
       Read(15,*) nrows
       allocate(t(nrows))
       Do i=1,nrows
             Read(15,*,End=1) t(i)
       End Do
       call fgrid(t, fstop, fstep, fr0, nrows)	
       open(16,file=return_fname,status='replace')
       nparameters=3
       narrays=0
       write (16,*) nparameters
       write (16,*) narrays
       write (16,*) fstop
       write (16,*) fstep
       write (16,*) fr0
       close(16)

  Case ('fouw')
       Read(15,*)frin,nh2
       Read(15,*) nrows
       allocate(t(nrows))
       allocate(v(nrows))
       allocate(w(nrows))
       allocate(valout(nrows))
       allocate(cof(nh2/2*2+2))
       allocate(dcof(nh2/2*2+2))
       Do i=1,nrows
             Read(15,*,End=1) t(i)
       End Do
       Read(15,*) nrows
       Do i=1,nrows
             Read(15,*,End=1) v(i)
       	     v(i)=v(i)-total/nrows
       End do
       Read(15,*) nrows
 
       Do i=1,nrows    
             Read(15,*,End=1) w(i)
       End do
       call fouw(t, v, w, frin, nh2, frout, dfrout, valout,cof, dcof, nrows)


       open(16,file=return_fname,status='replace')
       nparameters=2
       narrays=3
       write (16,*) nparameters
       write (16,*) narrays
       write (16,*) frout
       write (16,*) dfrout
       write (16,*) nrows
       do j=1,nrows
	     write (16,*) valout(j)
       end do
       write (16,*) nh2/2*2+2
       do j=1,nh2/2*2+2
	     write (16,*) cof(j)
       end do
       write (16,*) nh2/2*2+2
       do j=1,nh2/2*2+2
	     write (16,*) dcof(j)
       end do
       close(16)

  Case ('peak')
       Read(15,*) ! no parameters but still need to read parameter line
       Read(15,*) nrows
       allocate(v(nrows))
       Do i=1,nrows
             Read(15,*,End=1) v(i)
       End Do
       call peak(nrows, v, xm, fm, dx)	
       open(16,file=return_fname,status='replace')
       nparameters=3
       narrays=0
       write (16,*) nparameters
       write (16,*) narrays
       write (16,*) xm
       write (16,*) fm
       write (16,*) dx
       close(16)




  Case ('prew')

       Read(15,*) frin,nh2
       Read(15,*) nrows

       allocate(t(nrows))
       allocate(v(nrows))
       allocate(w(nrows))
       allocate(lth(nrows))
       total=0
       Do i=1,nrows
             Read(15,*,End=1) t(i)
   
       End Do
       Read(15,*) nrows

       Do i=1,nrows
             Read(15,*,End=1) v(i)
       	     v(i)=v(i)-total/nrows
       End do
       Read(15,*) nrows
 
       Do i=1,nrows    
             Read(15,*,End=1) w(i)
       End do
       call prew(t, v, w, frin, nh2, frout, dfrout, lth, nrows)

       open(16,file=return_fname,status='replace')
       nparameters=2
       narrays=1
       write (16,*) nparameters
       write (16,*) narrays
       write (16,*) frout
       write (16,*) dfrout
       write (16,*) nrows

       do j=1,nrows
	     write (16,*) lth(j)
       end do
       close(16)


  Case ('aovdrv')

       Read(15,*) method ,npar, ncov, lfr0, frs, lnfr,sptype
       Read(15,*) nrows

       allocate(t(nrows))
       allocate(v(nrows))
       allocate(w(nrows))
       allocate(lth(lnfr))
       total=0
       Do i=1,nrows
             Read(15,*,End=1) t(i)
   
       End Do
       Read(15,*) nrows

       Do i=1,nrows
             Read(15,*,End=1) v(i)
       	     v(i)=v(i)-total/nrows
       End do
       Read(15,*) nrows
 
       Do i=1,nrows    
             Read(15,*,End=1) w(i)
       End do

       Call aovdrv (method,t, v, w,npar, ncov, lfr0, frs, lnfr,sptype, lth,np0,stat,vr,nrows)
       open(16,file=return_fname,status='replace')
       nparameters=3
       narrays=1
       write (16,*) nparameters
       write (16,*) narrays
       write (16,*) np0
       write (16,*) stat
       write (16,*) vr
       write (16,*) lnfr	
       do j=1,lnfr
	     write (16,*) lth(j)
       end do
       close(16)

  End Select

1 Continue
end program aovwrapper
