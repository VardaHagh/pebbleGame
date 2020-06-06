c $Id: $
c $Log: $
      subroutine build2dCF(so,sf,iflop)
c ------------------------------------------------------------------------------
c                               Description
c    This subroutine places the CF bond <so,sf> in the network, while checking
c if it is independent or redundant. The network connectivity and the stressed 
c state of all the CF bonds is updated in addition to the pebble index and the 
c number of floppy modes.  When an overconstrained region is larger than
c mintriang,   the topology of the network is changed by the process of
c TRIANGULARIZATION. The triangularization process, and the updating of the 
c network connectivity and stress state is actually contained in the
c subroutine  free2d1t()  which is used in place of  free2d1().
c ------------------------------------------------------------------------------
      integer so,sf
c =====================================================================
c iflop = total number of floppy modes in the network.
c so,sf = site labels defining the incident vertices of a CF bond
c =====================================================================
      call free2d2(so)
      call free2d1t(sf,so,np)
      if( np .eq. 2 ) iflop = iflop - 1 
      return
      end
