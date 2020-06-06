      program test
C ----------------------------------------------------------------------------
c   This is a sample program illustrating the use of pebble game routines.
c  It considers a small configuration consisting of 4 atoms with 4 bonds
c  between them, finds the number of floppy modes, decomposes into
c  rigid clusters, and finds if there are any stressed bonds.
c
c  All variables that are not described in comments (tag, btag, block etc.)
c  should be initialized exactly as in this program.
C ----------------------------------------------------------------------------

C ----------------------------------------------------------------------------
c Include the parameter file containing the max number of atoms nsmax=nmax,
c twice the maximum number of bonds nb2max and maximum possible coordination
c number mcoord. Of course, these parameters can be changed.

      include 'par_2dcfrp'

C ----------------------------------------------------------------------------

      integer btag,point(nsmax),linkcf(nb2max),clst(nb2max),tag
      integer pebble(id,0:nmax),block(0:nmax),btrack(0:nmax)
      integer shell(0:nmax)
      integer multcf(nsmax),state(nb2max),stressed,covered
      integer addstress(0:3)
      integer x_dum,y_dum,n_bonds,u
      CHARACTER(len=32) :: arg

      common/rigidity/ pebble,block,tag,shell,nsfail
      common/search/ btrack,btag
      common/triang/ ns,mintri
      common/stressmk/ addstress
      common/network/ n,linkcf,point,multcf,state,covered,stressed

C ------------------ Initialization -------------------------------------------


c     Input Checking

      call getarg(1,arg)


c      n_bonds = 3062
      u=20
      n = 512   ! number of atoms
      ns = n
      nd = 2*n   !  total number of degrees of freedom
      tag = 0
      btag = nd

      iflop = nd   !  initial number of floppy modes; will be updated when
                   !  bonds are inserted


      do s=0,n
       pebble(1,s) = 0   ! 1st pebble of a site; 0 means it is free so far
       pebble(2,s) = 0   ! 2nd pebble
       block(s) = 0
      enddo

      do s=1,n

       multcf(s) = 0            !   initial coordination number for each site

       point(s) = (s-1)*mcoord  !   points where link table for atom s will
                                !   start in linkcf array

      enddo

      stressed = 1
      covered = 2
c      addstress(0) = 1
c      addstress(1) = 0
c      addstress(2) = 1
c      addstress(3) = 0

      addstress(0) = 1
      addstress(1) = 1
      addstress(2) = 3
      addstress(3) = 3



      mintri = 50        ! this is a number of failed pebble searches after
                         ! which triangularization occurs. It can make the
			 ! program faster, but does not affect the results

C -------------- Build the network calling the pebble game ------------------

      open(u,FILE='points.dat')
      do
        read(u,*,iostat=io) x_dum,y_dum
        if (io/=0) then
          exit
        endif
         call build2dCF(x_dum,y_dum,iflop)
      enddo
      close(u)

c  The positions of pebbles, the number of floppy modes iflop, the indicator
c  of stress state(), the array of coordination numbers multcf and the
c  link table linkcf are updated automatically.

C ---------------------------------------------------------------------------

      open(u,FILE='floppymodes',STATUS='OLD')
      write(u,'("Number of floppy modes: ",I5)')iflop ! write resulting number
                                                      !of floppy modes
      close (u)

      call clst2dCF(n,iflop,clst)  ! decompose into clusters;
                                   ! clst(i) indicates the number of
				   ! cluster to which bond i belongs;
				   ! bonds belonging to different clusters
				   ! have different values of clst


      do i=1,n
       ipi = point(i)

       do j=1,multcf(i)    !  loop over all bonds stemming from atom i

        k = linkcf(ipi+j)  !  k is the atom to which the bond points

	if (k .gt. i) then

c	 write(*,'("bond connecting sites ",i1," and ",i1,
c     &   " is in cluster ",i1)')i,k,clst(ipi+j)





c Value of state is even, if the bond is not stressed and odd if it is stressed.
         if (mod(state(ipi+j),2) .eq. 0) then
            write(*, '(I7,A3,I7,A3,I5,I5)') i,' ', k,' ', clst(ipi+j),0

	  else
            write(*, '(I7,A3,I7,A3,I5,I5)') i,' ', k,' ', clst(ipi+j),1
	 endif

c	 write(*,*)

	endif
       enddo
      enddo

      end
