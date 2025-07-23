#
# This file is a demo file : copy it  and use yoursas an input of adel dataflow !

#
#
#
#Write here R code returning geometric functions in a list named Rfun
#
Rfun = list(
  #azim should return the azimutal deviation of a leaf compared to the previous one (A = axe number [0..n], n = leaf nulmber from the base, ntop = idem from the top
  azim = function(a,n,ntop) {
    180 + 10 * (runif(1)-.5)
  },
  #xysr returns the index of the xysr database associated with a leaf
  xysr = function(a,n,ntop) {
    ntop + 1
  },
  #geomT returns azT,incT and dred.
  #azT is the difference of azimuth between parent axe and tiller
  #incT (deg) is the baseal inclination of tiller
  #dred (mockup length unit,ie cm for wheat) is the distance at flowering between top of parent axe and top of tiller
  geomT = function(a) {
    if (a == 0) {
      azT = 0
      incT = runif(1) * 5
      dred =0
    } else {
      azT = 75 + (runif(1) - .5) * 5
      incT = 82 + (runif(1) - .5) * 5
      dred = runif(1) * 7
    }
    list(azT=azT,incT=incT,dred=dred)
  }
  )
