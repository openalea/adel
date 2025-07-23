/*-*-c-*-*/


/*_________________________________________________________L_System mais : version cpfg pour la production de chaine et/ou l'entree de la tortue PlantGL

version adaptee de ADEL SimExp (Gaetan)
*/

#define NOVISUALISATION

/* parametres */

#include "AdelPars.h"
#include "AleaPars.h"


/*_________________________________________Macros
*/

#define MAX(a,b) ((a)>=(b)?a:b)

#define MIN(a,b) ((a)<=(b)?a:b)

/* modele lineaire */
#define yLin(t,tdeb,tfin,ydeb,yfin)				\
  ((t) <= (tdeb) ? ydeb : ((t) >= (tfin) ? yfin : (ydeb) + ((yfin)-(ydeb)) * (t-(tdeb)) / ((tfin)-(tdeb))))

#ifndef TIMESTEP
#define TIMESTEP 1
#endif
#ifndef TEMP
#define TEMP 20
#endif
#ifndef STEPS
#define STEPS 60
#endif

/* ___________________ ADEL SimExp*/

/* REM passer a dd/tip/col en entree */

/* Dates d'apparition des primordias (dd) */
#define TINI(Numf) KinFeu[(Numf)-1][0]
/* Delai initiation - extension rapide des feuilles (dd) */
#define DIE(Numf) KinFeu[(Numf)-1][1]
/* Delai initiation - apparition des feuilles (dd) */
#define DIA(Numf) KinFeu[(Numf)-1][2]
/* Delai initiation - ligulation des feuilles (dd)*/
#define DIL(Numf) KinFeu[(Numf)-1][3]
/* Delai initiation -extension rapide des entrenoeuds */
#define DIE_EN(n) DIL(n)
/* Delai initiation - Fin extension Entrenoeud  : a revoir ?*/
#define DIF_EN(n) DIL(n)+0.6*(DIL(n)-DIE(n))

/* temps thermique, base 9,8 (=dd/d) */
#define TT(Temp) ((Temp>=9.8)?Temp-9.8:0)

/* Profils de longueur de limbe (Fournier and Andrieu, 1998) */
#define L_LIMB(NumPhyto) PhytomerDims[(NumPhyto)-1][0]
/* Profils de longueur de gaine (Fournier and Andrieu, 1998) */
#define L_GAINE(NumPhyto) PhytomerDims[(NumPhyto)-1][1]
/* Profils de longueur d'entrenoeud (Fournier and Andrieu, 1998) */
#define L_ENT(NumPhyto) PhytomerDims[(NumPhyto)-1][2]
/* Largeur max de limbe */
#define Wmax_LI(NumPhyto) PhytomerDims[(NumPhyto)-1][3]



/*Azimuth des phytomere,parametrisation proposee dans Fournier 1998 a partir de Drouet 1997 */
#define AZIM(Numphyto) \
  (((Numphyto)<6) ? 0 :(((Numphyto)<=8) ? (ran(1)-0.5)*20 : (ran(1)-0.5)*60))

/* Fuinctions for whorl*/

/* maj hw par un tube emboité, renvoie distance visible */
#define calcLvisTube(l,lvis)\
  lvis=MAX(0,l-Hw);	    \
  Hw=Hw+lvis

/* maj hw par un cyclindre plein emboité, renvoie distance visible*/
#define calcLvisCyl(l,lvis)\
  lvis=MAX(0,l-Hw);	   \
  Hw=MAX(0,Hw-l)


/*_________________________________________________ Module data */

#define ApexData ad_nump,ad_N,ad_n,ad_age,ad_Tini,ad_Torg
#define ApexRData adr_nump,adr_age,adr_Torg
#define EnData ed_nump,ed_num,ed_l,ed_IsExtended,ed_L,ed_diam,ed_age,ed_Tdeb,ed_Tfin,ed_Torg
#define EnMData emd_nump,emd_num,emd_L,emd_diam
#define LeafData ld_nump,ld_num,ld_Torg,ld_age,ld_IsExtended,ld_l,ld_g,ld_wm,ld_ltube,ld_KinData,ld_HwData,ld_GeomData
#define LeafMData lmd_nump,lmd_num,lmd_L,lmd_G,lmd_Wm,ld_GeomData

#define ld_HwData ld_Tapptip,ld_Tnextip,ld_hwnextip,ld_Tnextcol,ld_hwnextcol
#define ld_KinData ld_L,ld_G,ld_TdebL,ld_TfinL,ld_TfinG
#define ld_GeomData ld_azim,ld_diam,ld_rand

/* _________________________________Macros pour la generation de plante */

#define makeNewApex(nump)			\
  ad_nump=nump;\
  ad_N=NB_PHY;\
  ad_n=0;\
  tt=TT(TEMP);\
  ad_age=tt*TIMESTEP;\
  ad_Tini=TINI(1);\
  ad_Torg=TEMP

/*___________________________________________________ Macros pour l'apex */

#define ageApex(dt) \
tt=TT(ad_Torg);\
ad_age=ad_age+tt*dt

#define ageApexR(dt) \
tt=TT(adr_Torg);\
adr_age=adr_age+tt*dt

#define siNewPhyto ((ad_age >= ad_Tini) && (ad_n < ad_N))

#define makeNewPhyto \
makeNewEn;\
makeNewLeaf;\
ad_n=ad_n+1;\
ad_Tini=TINI(ad_n)

#define siRepro (ad_n >= ad_N)

#define makeRepro \
  adr_nump=ad_nump;\
  adr_age=ad_age;\
  adr_Torg=ad_Torg

#define makeNewEn \
  ed_nump=ad_nump;\
  ed_num=ad_n+1;\
  ed_l=0;\
  ed_IsExtended=0;\
  ed_L=L_ENT(ad_n+1);\
  ed_diam=DIAM;\
  ed_age=ad_age-ad_Tini;\
  ed_Tdeb=DIE_EN(ad_n+1);\
  ed_Tfin=DIF_EN(ad_n+1);\
  ed_Torg=ad_Torg

#define makeNewLeaf \
    ld_nump=ad_nump;\
    ld_num=ad_n+1;\
    ld_Torg=ad_Torg;\
    ld_age=ad_age-ad_Tini;\
    ld_l=0;\
    ld_g=0;\
    ld_wm=Wmax_LI(ad_n+1);\
    ld_ltube=0;\
    ld_IsExtended=0;\
    newKinData(ad_n+1);\
    newGeomData(ad_n+1);\
    numhw=MIN(ad_n+1,ad_N-1);\
    newHwData(numhw)

/* parametres du modèle d'allongement, origine des temps : initiation feuille */
#define newKinData(n)\
    ld_L=L_LIMB(n);\
    ld_G=L_GAINE(n);\
    ld_TdebL=DIE(n);\
    ld_TfinL=DIE(n)+L_LIMB(n)/(L_LIMB(n)+L_GAINE(n))*(DIL(n)-DIE(n));\
    ld_TfinG=DIL(n)

/* modèle de cornet
interpole lineairement entre l'apparition, ou la feuille n'est pas enroullee, l'apparition next tip ou la feuille peute etre enroullee et l'apparition. Les dates sont relative a l'inittiation de la feuille n
*/

#define newHwData(n) \
    ld_Tapptip=DIA(n);\
    ld_Tnextip=TINI(n+1)-TINI(n)+DIA(n+1);\
    ld_hwnextip=(L_LIMB(n+1)+L_GAINE(n+1))*(DIA(n+1)-DIE(n+1))/(DIL(n+1)-DIE(n+1));\
    ld_Tnextcol=TINI(n+1)-TINI(n)+DIL(n+1);\
    ld_hwnextcol=MAX(L_GAINE(n),L_GAINE(n+1))

#define newGeomData(n) \
    ld_azim=AZIM(n);   \
    ld_diam=DIAM;      \
    ld_rand=ran(1)


    /* ___________________________________________________macro pour la feuille */

#define ageLeaf(dt) \
tt=TT(ld_Torg);\
ld_age=ld_age+tt*dt;\
ld_IsExtended=0

#define siExtension (ld_IsExtended <= 0)

#define extendLeaf \
    ld_l=yLin(ld_age,ld_TdebL,ld_TfinL,0,ld_L);		\
    ld_g=yLin(ld_age,ld_TfinL,ld_TfinG,0,ld_G);		\
    ld_ltube=(ld_age <= ld_Tnextip) ? yLin(ld_age,ld_Tapptip,ld_Tnextip,0,ld_hwnextip) : yLin(ld_age,ld_Tnextip,ld_Tnextcol,ld_hwnextip,ld_hwnextcol); \
    ld_IsExtended=1

/* ligulation =  apres next col */
#define siMature (ld_age >= ld_Tnextcol)

#define makeMatureLeaf\
    lmd_nump=ld_nump;\
    lmd_num=ld_num;\
    lmd_L=ld_L;\
    lmd_G=ld_G;\
    lmd_Wm=ld_wm

    /* ___________________________________________________macro pour l'entrenoeud */

#define ageEn(dt)\
tt=TT(ed_Torg);\
ed_age=ed_age+tt*dt;\
ed_IsExtended=0

#define siExtensionEn (ed_IsExtended <= 0)

#define extendEn \
    ed_l=yLin(ed_age,ed_Tdeb,ed_Tfin,0,ed_L);	\
    ed_IsExtended=1

#define siMatureEn (ed_age >= ed_Tfin)

#define makeMatureEn\
    emd_nump=ed_nump;\
    emd_num=ed_num;\
    emd_L=ed_L;\
    emd_diam=ed_diam

/*_________________________________________LSystem */

Lsystem: 0
derivation length: STEPS


/* tableaux parametres : cf macro ci dessus pour indices*/
define: { array PhytomerDims[NB_PHY][4]={DIMVALUE},
    KinFeu[NB_PHY][4]={KINFEUVALUES};}


Start: {
nbplante = NB_PLANTS;	
Hw=0;
alpha =ALPHA_LARG;
beta=-2*(alpha+sqrt(-alpha));
gamma=2*sqrt(-alpha)+alpha;
/* verif input */
 printf("contenu du tableau kinetic:\nNumPhy\tTini\tDie\tDia\tDil\n");
 i=1;
 while (i<=NB_PHY) {printf("%g\t%g\t%g\t%g\t%g\n",i,TINI(i),DIE(i),DIA(i),DIL(i));i=i+1;}
}


#ifndef AXM
Axiom: doplante(0,nbplante)
#else
AXM
#endif

 /*________Productions = mise à jour avec l'environnement et/ou le pas de temps des ages */

Apex(ApexData) : * {ageApex(TIMESTEP);} --> Apex(ApexData)
ApexR(ApexRData) : * {ageApexR(TIMESTEP);} --> ApexR(ApexRData)
Feu(LeafData) : * {ageLeaf(TIMESTEP);} --> Feu(LeafData)
En(EnData) : * {ageEn(TIMESTEP);} --> En(EnData)

/*_________Decomposition = mise a jour de la chaine */

decomposition
maximum depth: 100

doplante(dejnum,reste) : reste > 0 {makeNewApex(dejnum+1);} --> 
	Spacer 
	[newPlant(dejnum+1)Apex(ApexData)]
	doplante(dejnum+1,reste-1)

doplante(dejnum,reste) : reste <=0 --> *

Apex(ApexData) :  siNewPhyto {makeNewPhyto;} -->
	newPhyto 
	En(EnData)
	Feu(LeafData)
	Apex(ApexData)

Apex(ApexData) : siRepro {makeRepro;} --> ApexR(ApexRData)
 
Feu(LeafData) : siExtension {extendLeaf;} --> Feu(LeafData)

Feu(LeafData) : siMature {makeMatureLeaf;} --> FeuM(LeafMData)

En(EnData) : ed_L <= 0 --> *

En(EnData) : siExtensionEn {extendEn;} --> En(EnData)

En(EnData) : siMatureEn {makeMatureEn;} --> EnM(EnMData)
    

homomorphism
maximum depth: 2

  /*_______________________________Calcul visibilités, génération modules pour représentation geom */

newPlant(n) : * {Hw=0;} --> Plant(n)

Feu(LeafData) : * {
	gvis=MAX(0,ld_g-Hw);
	calcLvisTube(ld_ltube,tubevis);
	lvis=MAX(0,ld_l+ld_g-Hw);} --> 
	PStem(ld_nump,ld_num,1,0,tubevis,gvis,ld_diam)
	Leaf(ld_nump,ld_num,1,1,ld_L,lvis,ld_wm,ld_rand,(ld_age-ld_Tapptip)/(ld_Tnextcol-ld_Tapptip),ld_azim)

FeuM(LeafMData) : * {calcLvisTube(lmd_G,gvis);} --> 
	PStem(lmd_nump,lmd_num,1,0,gvis,gvis,ld_diam)
	Leaf(lmd_nump,lmd_num,1,1,lmd_L,lmd_L,lmd_Wm,ld_rand,1,ld_azim)

En(EnData) : {calcLvisCyl(ed_l,lvis);} lvis > 0 --> Stem(ed_nump,ed_num,1,0,lvis,ed_diam)

EnM(EnMData) : {calcLvisCyl(emd_L,lvis);}  lvis > 0 --> Stem(emd_nump,emd_num,1,0,lvis,emd_diam)

/* visualisation si on est sous lstudio */
#ifdef VISUALISATION

Plant(n) : * {Hw=0;} --> ;(9)

newPhyto --> /(180)
PStem(np,n,esp,opak,lg,g,diam) --> F(lg*100)
Leaf(np,n,esp,opak,L,l,wm,random,r,azim) : * {l=l*100;} --> [/(azim)+(35*r)F(l/3)+(25*r)F(l/6)+(25*r)F(l/6)+(45*r)F(l/6)+(35*r)F(l/6)]
Stem(np,n,esp,opak,l,diam) --> F(l*100)
#else
/* sinon on produit une chaine pour alea */

Plant(n) --> newPlant newAxe

newPhyto --> /(180)newMetamer

PStem(np,n,esp,opak,lg,g,diam) : lg > 0  --> StemElement(esp,lg,diam,diam)

Leaf(np,n,esp,opak,L,l,wm,random,r,azim) : l > 0 --> [/(azim)LeafElement(esp,L,l,wm,0,1,n,random)]

Stem(np,n,esp,opak,l,diam) --> StemElement(esp,l,diam,diam)

Spacer --> *
En(EnData) --> *
EnM(EnMData) --> *
#endif

endLsystem
