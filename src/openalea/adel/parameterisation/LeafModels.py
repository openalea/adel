# a collection of equation for simulating leaf shapes

from math import radians,tan,log,sqrt
import numpy as np

from openalea.plantgl.all import NurbsCurve2D


# simple nurbs for training

def simpleNurbs_sr(lw = 0.1):
    """
    construct a simple normalised sr profile, whose parameter is leaf/width ratio
    """

    return NurbsCurve2D([(0,0,1),(0.0508334,lw,1),(0.746667,lw,1),(1,0,1)], width = 2)
    

def simpleNurbs_xy(xins = 0.2,xtop = 0.6, xend = 1, yend = 1):
    """
    construct a simple xy profile, with simplified control points
    """

    return NurbsCurve2D([(0,0,1),(xins,xins,1),(xtop,xtop,1),(xend,yend,1)])


#Geometric model of prevot

def xy_Prevot(phi0,Cp,Pcass,phiM,eps,PsiE):
    """
    """
    b = tan(radians(phi0))
    Cp = max(4./3,min(3.99,Cp))
    xm = b * (4 - 3 * Cp) / (Cp - 4)
    sqxm = sqrt(1 + xm**2)
    sqtan = sqrt(1 + b**2)
    a = 1. / (4 * Pcass) * (  log( (xm + sqxm) / (b + sqtan) ) + xm * sqxm - b * sqtan  )

    x = np.linspace(0,xm)
    
    return x, a * x**2 + b * x,xm

# Geometric model of leaf midrib by Steward And Dwyer(Agricultural & Forest Meteorology, 66 (1993) 247-265.


def xy_StewartDwyer93(insertion_angle, xtop, ytop, xend, yend):
    """
    Construct a (x,y) pair of array defining the midrib

    insertion_angle in the angle between the stem and the leaf (deg)
    xtop,ytop are the coordinates of the highest point of the midrib
    xend, yend are the coordinates of the end of the midrib
    """

    x = np.linspace(0,xend)
    theta = radians(90 - insertion_angle)
    
    Bdiv = ytop**2 * xend**2 + xtop**2 * yend**2 - 2 * xtop * ytop * xend * yend

    if Bdiv == 0 :
	return x, tan(theta) * x

    else :
    
	D = - tan(theta)
	B = ( 2 * xtop * xend * yend + tan(theta) * xtop**2 * xend - (tan(theta) / ytop) * xtop**2 * xend * yend - xtop**2 * yend ) / Bdiv
	A = B * ytop**2 / xtop**2 + ytop / xtop**2
	C = (tan(theta) - 2 * A * xtop) / ytop

	a = B
	b = C * x + 1
	c = A * x**2 + D * x

	delta = b**2 - 4 * a * c

	#if all(delta > 0):
	if (0 > 0):
	    y = ( - b + np.sqrt(delta) ) / (2 * a)
	else: #parabole + parabole
	    b = tan(theta)
	    a =  - b / (2 * xtop)
	    yt = a * xtop**2 + b * xtop
	    y = ytop / yt * (a * x**2 + b * x)
	    a2 = (yend - ytop) / (xend - xtop)**2
	    y[x > xtop] = ytop + a2 *  (x[x > xtop] - xtop)**2
	
	return x,y
