      subroutine free2d2(s0)
c ----------------------------------------------------------------
c PROGRAM WRITTEN BY: Donald J. Jacobs                Sept 6, 1996
c     This subroutine collects the MAXIMUM number of pebbles at 
c site s0, without any other free pebbles tied down. Therefore, 
c it is ALWAYS possible to collect two pebbles at site s0.
c ----------------------------------------------------------------
      include 'par_2dcfrp'
      integer s0,tag,btag
      integer pebble(id,0:nmax),shell(0:nmax),block(0:nmax)
      integer btrack(0:nmax) 
c ================================================================
      common/rigidity/ pebble,block,tag,shell,nsfail
      common/search/   btrack,btag
c ================================================================
      nget = 2
      if( pebble(1,s0) .eq. 0 ) nget = nget - 1
      if( pebble(2,s0) .eq. 0 ) nget = nget - 1
      btrack(s0) = 0
         do ip=1,nget
         tag = tag + 1
         block(s0) = tag
         call collect2d(s0)
         enddo
      return
      end

      subroutine free2d1(s0,s1,np)
c --------------------------------------------------------------
c     This subroutine collects the MAXIMUM number of pebbles at 
c site s0, with the (two) free pebbles at site s1 tied down. It 
c is ALWAYS possible to collect at least one pebble. If two 
c pebbles can be collected at site s0, then the CF bond <s0,s1>
c is independent, otherwise it is redundant.
c ----------------------------------------------------------------
      include 'par_2dcfrp'
      integer s,s0,s1,s2,tag,btag
      integer pebble(id,0:nmax),shell(0:nmax),block(0:nmax)
      integer btrack(0:nmax)
c ================================================================
      common/rigidity/ pebble,block,tag,shell,nsfail
      common/search/   btrack,btag
c ================================================================
      np = 0
      if( pebble(1,s0) .eq. 0 ) np = np + 1
      if( pebble(2,s0) .eq. 0 ) np = np + 1
      nget = 2 - np 
      btrack(s0) = 0
         do ip=1,nget
         tag = tag + 1
         block(s0) = tag
         block(s1) = tag
         call collect2d(s0)
c ----------------------------------------------------------------
c                           Note: nsfail = -1 => found free pebble
            if( nsfail .lt. 0 ) then
            np = np + 1
            else
            return
            endif
         enddo
      pebble(1,s0) = s1
      return
      end

      subroutine free2d1t(s0,s1,np)
c ------------------------------------------------------------------------------
c                           Description
c ==============================================================================
c     This subroutine collects the Maximum number of pebbles at site s0, given
c that two pebbles are tied down at site s1. If two free pebbles can be obtained
c at site s0, then the bond <s0,s1> is independent and the state of the network
c is updated to indicate that bond <s0,s1> is covered. If only one free pebble 
c can be found (at least one free pebble can always be found), then bond <s0,s1>
c is redundant and the state of the network is updated to record the entire
c overconstrained region defined by the Laman subgraph (the failed pebble search). 
c The network connectivity is updated within this subrountine. 
c If the number of sites in the failed search is greater than mintri, then the
c region defining the "Laman subgraph" is triangularized using the <s0,s1> bond
c as a base. In this process the topology of the pebble index deviates from that
c of the network. Note that with respect to the original network topology, the 
c <s0,s1> bond is not covered. 
c ------------------------------------------------------------------------------
      include   'par_2dcfrp'
      integer   s,s0,s1,tag,btag
      integer   pebble(id,0:nmax),shell(0:nmax),block(0:nmax)
      integer   linkcf(nb2max),point(nsmax),btrack(0:nmax)
      integer multcf(nsmax),state(nb2max),covered,stressed
      integer addstress(0:3)
c ==============================================================================
      common/network/  n,linkcf,point,multcf,state,covered,stressed
      common/rigidity/ pebble,block,tag,shell,nsfail
      common/search/   btrack,btag
      common/triang/   ns,mintri
      common/stressmk/ addstress
c -------------------------------------------------- update network connectivity 
      multcf(s0) = multcf(s0) + 1
      multcf(s1) = multcf(s1) + 1
      index0 = point(s0) + multcf(s0)
      index1 = point(s1) + multcf(s1)
      linkcf(index0) = s1
      linkcf(index1) = s0
      state(index0) = 0
      state(index1) = 0
c ----------------------------------------------------- apply the 2d pebble game
      np = 0
      if( pebble(1,s0) .eq. 0 ) np = np + 1
      if( pebble(2,s0) .eq. 0 ) np = np + 1
      nget = 2 - np 
      btrack(s0) = 0
         do ip=1,nget
         tag = tag + 1
         block(s0) = tag
         block(s1) = tag
         call collect2d(s0)
c ------------------------------------------------------------------------------
c                           Note: nsfail = -1 => found free pebble
            if( nsfail .lt. 0 ) then
            np = np + 1
            else
            nsfail = nsfail + 1
            shell(nsfail) = s1
c ---------------------------------------------- record over-constrained regions
            do k=1,nsfail
            s = shell(k)
            index = point(s)
               do j=1,multcf(s)
               index = index + 1
                  if( block( linkcf(index) ) .eq. tag ) then 
                  state(index) = addstress( state(index) )
                  endif
               enddo
            enddo

               if( nsfail .ge. mintri ) then
c ------------------------------------------------- triangularize Laman subgraph
               pebble(1,s0) = 0
               pebble(2,s0) = s1
                  do jfail=2,nsfail-1      ! s1 = shell(1) keeps 2 free pebbles
                  s = shell(jfail)
                  pebble(1,s) = s0
                  pebble(2,s) = s1
                  enddo
               endif
            return
            endif
         enddo
c ------------------------------------------------------------------- cover bond
      state(index0) = state(index0) + covered 
      state(index1) = state(index1) + covered 
      pebble(1,s0) = s1
      return
      end

      subroutine free2d0(s0,kmax)
c ----------------------------------------------------------------
c     This subroutine collects a SINGLE free pebble at site s0 if
c possible with the free pebbles at sites s1 and s2 tied down. If
c a free pebble cannot be found, then site s0 is mutually rigid 
c with respect to sites {s1,s2}. If a free pebble is successfully 
c found, then site s0 is floppy with respect to the sites {s1,s2}.
c In this subroutine, we do not collect the maximum number of 
c pebbles possible.
c ----------------------------------------------------------------
      include 'par_2dcfrp'
      integer s0,tag,btag,rigid
      integer pebble(id,0:nmax),shell(0:nmax),block(0:nmax)
      integer btrack(0:nmax) 
      data rigid /2147483644/
c          rigid = 2**31 - 4
c ==============================================================
      common/rigidity/ pebble,block,tag,shell,nsfail
      common/search/   btrack,btag
c ==============================================================
c ______ IMPORTANT!!! ______ IMPORTANT!!! ______ IMPORTANT!!! ____
c ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c %%%%%%%%%%%%%%%%%%%%%%%% Must do this initialization before call
c      block(s1) = rigid
c      block(s2) = rigid
c %%%%%%%%%%%%%%%%%%%%%%%% 
c --------------------------------------------------------------
      btrack(s0) = btag
      if( pebble(1,s0) .eq. 0 ) return
      if( pebble(2,s0) .eq. 0 ) return
      tag = tag + 1
      block(s0) = tag
      kkk = kmax
      call burn2d(s0,kkk)
c --------------------------------------------------------------
         if( nsfail .gt. 0 ) then
            do k=kkk,nsfail
            block( shell(k) ) = rigid
            btrack( shell(k) ) = btag
            enddo
         kmax = nsfail
         return
         endif
      kmax = kkk
      return
      end

      subroutine collect2d(s0)
c ----------------------------------------------------------------
c PROGRAM WRITTEN BY: Donald J. Jacobs                Sept 6, 1996
c This subroutine collects a single pebble at site s0 if possible.
c Free pebbles at other sites {s1,s2,...} that are tied down must 
c be specified within or before using the subroutines free2dN where 
c N=maximum number of free pebbles possible. The subroutine 
c collect2d(s0) is used as the basic pebble retriever for all the
c free2dN subroutines. Before calling collect2d(), tag must be 
c incremented in free2dN as well.
c ----------------------------------------------------------------
      include 'par_2dcfrp'
      integer s,s0,sa,sb,stest,tag,btag
      integer pebble(id,0:nmax),shell(0:nmax),block(0:nmax)
      integer btrack(0:nmax) 
c =========================================++=====================
      common/rigidity/ pebble,block,tag,shell,kmax
      common/search/   btrack,btag
c =======================================++=======================
c ______ IMPORTANT!!! ______ IMPORTANT!!! ______ IMPORTANT!!! ____
c ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c %%%%%%%%%%%%%%%%%%%%%%%% Must do this initialization before call
c         tag = tag + 1
c         block(s0) = tag
c         block(s1) = tag        ...
c         btrack(s0) = 0
c %%%%%%%%%%%%%%%%%%%%%%%%
c ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c ----------------------------------------------------------------
c ------------  variable list  (incomplete) ----------------------
c s0,s,sa,sb,stest   site labels
c kmax = number of sites that have been checked during the current
c      pebble search. Note that kmax = nsfail when there is a 
c      failed pebble search.  nsfail is used as a global variable.
c nsfail = the number of sites within a failed search for a free 
c        pebble where each site is identified by  shell(i) for 
c        ( 1 .le. i .le. kmax ) 
c tag = a reference mark to which site was checked during a search,
c       and it is a running dummy index incremented by one before 
c       each and every new search.
c nsmax = maximum total number of sites
c block(s) = locates all previously tagged sites that were searched
c            block(0) = 0 ALWAYS
c shell(j) = stores the generation of sites in a breath-first search
c btrack( ) = a BACKWARD TRACK that is followed when rearranging pebbles
c pebble(,) = defines the current directed graph of the network
c ----------------------------------------------------------------
c ------------------------------------------ search initialization
         kmax = 1
         shell(kmax) = s0
         k = 0
c -------------------------------------------------- pebble search
  100       k = k + 1
            if( k .gt. kmax ) return
c ----------------------------------------- continue pebble search
            s = shell(k)
            do 150 jj=1,id
            stest = pebble(jj,s)
            if( block(stest) .eq. tag ) goto 150
c ---------------------------------------------- rearrange pebbles
               if( stest .eq. 0 ) then

c ================================================================
               if( s .eq. s0  ) goto 150
c ================================================================

               kmax = -1
                  sa = btrack(s)
                  pebble(jj,s) = sa
  200             sb = btrack(sa)
                     if( pebble(1,sa) .eq. s ) then
                     pebble(1,sa) = sb
                     else
                     pebble(2,sa) = sb
                     endif
                  if( sb .eq. 0 ) return
                  s  = btrack(sb)
                     if( pebble(1,sb) .eq. sa ) then
                     pebble(1,sb) = s
                     else
                     pebble(2,sb) = s
                     endif
                  if( s  .eq. 0 ) return
                  sa = btrack(s)
                     if( pebble(1,s) .eq. sb) then
                     pebble(1,s) = sa
                     else
                     pebble(2,s) = sa
                     endif
                  if( sa .eq. 0 ) return
                  goto 200
c ------------------------------------------------------ grow tree
               else
               kmax = kmax + 1
               shell(kmax) = stest
               btrack(stest) = s
               block(stest) = tag
               endif
  150       continue
            goto 100
      end

      subroutine burn2d(s0,klow)
c ----------------------------------------------------------------
c PROGRAM WRITTEN BY: Donald J. Jacobs               Sept  6, 1996
c This subroutine checks if a single pebble at site so can be made
c free given that the free pebbles on sites {s1,s2} are tied down.
c Note that since a rigid cluster is being burned out, all 
c previously checked sites belonging to the rigid cluster (which 
c are unable to have a free pebble) are marked "rigid" to make
c future searches faster! No rearranging of pebbles is performed.
c ----------------------------------------------------------------
      include 'par_2dcfrp' 
      integer s,s0,sa,sb,stest,tag,btag,nsfail
      integer pebble(id,0:nmax),shell(0:nmax),block(0:nmax)
      integer btrack(0:nmax) 
c ==============================================================
      common/rigidity/ pebble,block,tag,shell,nsfail
      common/search/   btrack,btag
c ==============================================================
c ______ IMPORTANT!!! ______ IMPORTANT!!! ______ IMPORTANT!!! ____
c ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c %%%%%%%%%%%%%%%%%%%%%%%% Must do this initialization before call
c         tag = tag + 1
c         block(s0) = tag
c         block(s_rigid) = rigid
c %%%%%%%%%%%%%%%%%%%%%%%%
c ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c ----------------------------------------------------------------
c ------------  variable list  (incomplete) ----------------------
c s0,s,sa,sb,stest   site labels
c klow = the previous number of sites (used as index for  shell() )
c        that have been found to be mutually rigid.
c kmax = number of sites that have been checked during the current
c      pebble search. Note that kmax = nsfail when there is a 
c      failed pebble search.  nsfail is used as a global variable.
c nsfail = the number of sites within a failed search for a free 
c        pebble where each site is identified by  shell(i) for 
c        ( 1 .le. i .le. kmax ) 
c tag = a reference mark to which site was checked during a search,
c       and it is a running dummy index incremented by one before 
c       each and every new search.
c ns = the total number of sites
c block(s) = locates all previously tagged sites that were searched
c            block(0) = 0 ALWAYS
c shell(j) = stores the generation of sites in a breath-first search
c btrack( ) = a BACKWARD TRACK that is followed when rearranging pebbles
c pebble(,) = defines the current directed graph of the network
c ----------------------------------------------------------------
c ------------------------------------------ search initialization
      kmax = klow + 1
      shell(kmax) = s0
      k = klow
c -------------------------------------------------- pebble search
  100 k = k + 1
         if( k .gt. kmax ) then
         nsfail = kmax
         return
         endif
c ----------------------------------------- continue pebble search
      s = shell(k)
         do jj=1,id
         stest = pebble(jj,s)
            if( block(stest) .lt. tag ) then
               if( stest .eq. 0 ) then
               btrack(stest) = btag
               btrack(s0) = btag
               nsfail = -1
               return
               elseif( btrack(stest) .eq. btag ) then
               btrack(s0) = btag
               nsfail = -1
               return
               else
c ------------------------------------------------------ grow tree
               kmax = kmax + 1
               shell(kmax) = stest
               block(stest) = tag
               endif
            endif
         enddo
      goto 100
      end
