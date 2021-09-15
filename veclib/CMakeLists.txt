set (vec_lib_src
  util.c
  memory.c
  xcopy.c   xfill.c   xneg.c    xvneg.c   xsadd.c   xvadd.c   xssub.c
  xvsub.c   xsmul.c   xvmul.c   xsdiv.c   xvrecp.c  xvdiv.c   xzero.c
  xspow.c   xvsgn.c   xznan.c
  xvabs.c   xvamax.c  xvpow.c   xvexp.c   xvlg10.c  xvlog.c   xvatan.c
  xvatn2.c  xvcos.c   xvsin.c   xvsqrt.c  xrand.c   xvhypot.c xvtanh.c
  xsvmvt.c  xsvpvt.c  xsvtsp.c  xsvtvm.c  xsvtvp.c  xsvvmt.c  xsvvpt.c
  xsvvtm.c  xsvvtp.c  xvvmvt.c  xvvpvt.c  xvvtvm.c  xvvtvp.c  xvvvtm.c
  xsvvtt.c  xvvvtt.c
  xseq.c    xsge.c    xsle.c    xslt.c    xsne.c
  icount.c  ifirst.c  ixmax.c   ixmin.c   xsum.c    lany.c    lxsame.c
  vdble.c       xvfloa.c         vsngl.c            xbrev.c
  xscatr.c      xgathr.c         xgathr_scatr.c
  xscatr_sum.c  xgathr_sum.c     xgathr_scatr_sum.c
  xramp.c       xcndst.c         xmask.c
  xvpoly.c      xpolint.c        xclip.c
  xspline.c     xsplquad.c       xiclip.c
  xmxv.c        xmxva.c
  xvvtvvtp.c    xsvvttvp.c       xvvtvvtm.c
  xmxm.f
)

add_library (vec STATIC ${vec_lib_src})
