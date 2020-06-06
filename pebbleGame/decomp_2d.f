      subroutine clst2dCF(ns,iflop,clst)
c ---------------------------------------------------------------------
c PROGRAM WRITTEN BY: Donald J. Jacobs                    Sept  7, 1996
c         This program takes the connection information of a 2d generic 
c network as defined by the CF-bond linkage table (ie. linkcf() ) and 
c the pebble index  pebble(,)  containing all rigidity information, and
c decomposes the network into its unique set of rigid clusters. It also 
c calculates the number of floppy modes using the cluster statistics 
c sum rule as an additional self-consistency check. Various moments of
c the cluster statistics is calculated as well. 
c ---------------------------------------------------------------------
      include 'par_2dcfrp'
      integer s,so,s1,sf,stest,btag,rigid,smax,point(nsmax)
      integer linkcf(nb2max),clst(nb2max),tag,tag1
      integer multcf(nsmax),state(nb2max),stressed
      integer pivots(nsmax),mark,covered
      integer   stat(nsmax),gradus
      integer pebble(id,0:nmax),shell(0:nmax),block(0:nmax)
      integer btrack(0:nmax),len,lth
      
      data rigid /2147483644/  mark /32/
c          rigid = 2**31 - 4
c =====================================================================
      common/rigidity/pebble,block,tag,shell,nsfail
      common/network/n,linkcf,point,multcf,state,covered,stressed
      common/search/btrack,btag
c =====================================================================
      smax = 1
      nsvoid = 0
      tag1 = 0
      do s=1,ns
       stat(s) = 0
       pivots(s) = 0
      enddo
c --------------------------------------------- identify isolated sites
         do so=1,ns
          if( multcf(so) .lt. 1 ) stat(1) = stat(1) + 1
         enddo
c ----------------------------------------------------- void correction
      stat(1) = stat(1) - nsvoid
      nscorr = ns - nsvoid
c ----------------------------------------- decompose all clusters with
c                                                    more than one bond
         do so=1,ns
          ipo = point(so)
            do j=1,multcf(so)
               if( state(ipo+j) .lt. mark ) then
               sf = linkcf(ipo+j) 
c --------------------------------------------------------- map cluster
               call free2d2(so)
               call free2d1(sf,so,np1)
      

c +++++++++++++++++++++++++++++++++++++++++++++++ REMOVE
                  if( np1 .eq. 2 ) then
                  write(86,6543)
                  write(*,*)'in decomp',so,sf
 6543             format(5x,'More than three pebbles collected ',
     &                  'across a single bond in mapping cluster') 
                  stop
                  endif
                  if( np1 .ne. 1 ) then
                  write(86,6544)
 6544             format(5x,'Three pebbles can not be collected ',
     &                  'across a single bond in mapping cluster') 
                  stop
                  endif
c +++++++++++++++++++++++++++++++++++++++++++++++ 

c ------------------------------------------------- initialize the burn
c                                      btag marks sites already checked
               btag = btag + 1
               btrack(so) = btag
               btrack(sf) = btag
c -------------------------------------------------- mark rigid cluster
               block(so) = rigid
               block(sf) = rigid
               shell(1) = so
               shell(2) = sf
               kmax = 2
c ----------------------------------------------- burn out rigid region
               k = 0
  100          k = k + 1
                  if( k .gt. kmax ) then
c ------------------------------------------------------ record cluster
                  tag1 = tag1+1
                  stat(kmax) = stat(kmax) + 1
                  if( kmax .gt. smax ) smax = kmax
                     do k=1,kmax
                     s = shell(k)
                     pivots(s) = pivots(s) + 1
                     ips = point(s)
                        do jj=1,multcf(s)
                        stest = linkcf(ips+jj)
                           if( block(stest) .eq. rigid ) then
                           state(ips+jj) = state(ips+jj) + mark
                           clst(ips+jj) = tag1
                           endif
                        enddo
                     enddo
c ------------------------------------------------------- re-initialize
                     do k=1,kmax
                     block( shell(k) ) = 0
                     enddo
                  else
c -------------------------------------------------------- grow cluster
                  s = shell(k) 
                  ips = point(s)
c ------------------------------------------- search and expand about s
                     do jtest=1,multcf(s)
c ------------------------------------------------------------ unmarked
                        if( state(ips+jtest) .lt. mark ) then
                        stest = linkcf(ips+jtest) 
c ----------------------------------------------------------- unchecked
                           if( btrack(stest) .lt. btag ) then
                           btrack(stest) = btag
                           call free2d0(stest,kmax)
                           endif
                        endif
                     enddo
                  goto 100
                  endif
               endif
            enddo
         enddo
c ------------------------------------ calculate number of floppy modes
       ic = 0
          do s=2,smax
          ic = ic + stat(s)* (2*s - 3)
          enddo
       iflop2 = id* nscorr - ic
       npivot = 0 
         do s=1,ns
         if( pivots(s) .gt. 1 ) npivot = npivot+1
         enddo
c --------------------------------------------- self consistency checks
       iflag = 1
       itemp = 0
          do s=1,ns
             if( pivots(s) .gt. mcoord ) then
                if( iflag .eq. 1 ) then
                write(87,6723) 
 6723           format(5x,'PIVOT # exceeds maximum coordination')
                endif
             iflag = 0
             endif
          if( pivots(s) .eq. 0 ) itemp = itemp +  1
          enddo
       itemp2 = stat(1) + nsvoid
          if( itemp .ne. itemp2 ) then
          iflag = 0
          write(87,6724) 
 6724     format(5x,'PIVOT(s) = 0 is inconsistent with ',
     &              ' stat(1) + nsvoid')
          write(87,*) itemp,stat(1),nsvoid 
          endif
       if( iflop2 .ne. iflop ) then
       write(87,6725) iflop,iflop2
 6725  format(5x,'ERROR: F_direct = ',i8,' .NE. F_sumrule = ',i8)
       iflag = 0
       endif
!       if( iflag .ne. 1 ) stop

c ------------------------------------------------------ restore state
         do s=1,ns
          ips = point(s)
            do j=1,multcf(s)
            if( state(ips+j) .ge. mark ) state(ips+j) = 
     &         state(ips+j) - mark
            enddo
         enddo

      return
      end
