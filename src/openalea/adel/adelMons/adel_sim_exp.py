from . import IOtable
import os
from os.path import join
from openalea.core.pkgmanager import PackageManager
pm = PackageManager()
pkg = pm.get('openalea.adel.adelMons') 
path = ''
if pkg :
    path = pkg.path


path_ini = os.getcwd()

class adel_sim_exp:
    """  wrapper for adelsimexp2.lsys L-system (graphtal) - return can opt and 8 files """ 

    def __init__(self):
        pass


    def __call__(self, parameters, lsyst, TT, can):
        #creation du AdelPars.h
        os.chdir(path) 
        param = open('AdelPars.h', 'w')
        param.write('#define GEOMFILE "'+parameters+'"\n')
        param.close()

        #lancement d'ADEL et creation du .can
        os.system('graphtal.exe  -DSTEPS='+str(TT)+' -d can '+lsyst+'>'+can)# utiliser popen2

        # recup parametres optique et scene dans parameters
        f = open(parameters, 'r')
        tab = IOtable.table_txt(f) 
        f.close()

        d = {'SOIL_REFLECTANCE':'', 'LEAF_REFLECTANCE':'', 'LEAF_TRANSMITANCE':'', 'SEN_LEAF_REFLECTANCE':'', 'SEN_LEAF_TRANSMITANCE':'', 'STEM_REFLECTANCE':'', 'NB_RANG':'', 'NB_PLANTES':'', 'INTER_RANG':'', 'INTER_PLANTE':''}
        for var in list(d.keys()):
            for i in range(len(tab)):
                if len(tab[i])>1:
                    if tab[i][1]==var:
                        d[var]=tab[i][2]

        #creation du .opt
        f = open('par.opt', 'w')
        self.write_opt(f, d)
        f.close()

        #creation du .8
        lower_edge = [float(d['INTER_RANG'])/2, float(d['INTER_PLANTE'])/2]
        upper_edge = [-(float(d['NB_RANG'])*float(d['INTER_RANG']))+float(d['INTER_RANG'])/2, -(float(d['NB_PLANTES'])*float(d['INTER_PLANTE']))+float(d['INTER_PLANTE'])/2]
        f = open('maize.8', 'w')
        f.write(str(upper_edge[0])+" "+str(upper_edge[1])+"\n")
        f.write(str(lower_edge[0])+" "+str(lower_edge[1])+"\n")
        f.close()

        os.chdir(path_ini)
         
        #renvoie le patn du .can
        if can[0:2]=='.\\':
            return join(path, can[2:]), join(path, 'maize.8'),join(path, 'par.opt')
        else:
            return join(path, can), join(path, 'maize.8'),join(path, 'par.opt')


    @staticmethod
    def write_opt (f, d):
        """ ecriture du .opt a partir du dico contenant parametres optiques """
        f.write("# the number of plant optical species"+"\n")
        f.write("n 2"+"\n")
        f.write("# s stands for soil reflectance, d stands for diffuse and is the only option at this time."+"\n")
        f.write("# This line must be included even if there is no soil in the L-system nor in the simulation."+"\n")
        f.write("s d "+d['SOIL_REFLECTANCE']+"\n")
        f.write("# e stands for species. in this case just 1. therfore one line of data"+"\n")
        f.write("# d has the same meaning as for soil."+"\n")
        f.write("# the first number is stem reflectance"+"\n")
        f.write("# second and third are reflectance and transmittance of upper side"+"\n")
        f.write("# fourth and fifth are reflectance and transmittance of lower side"+"\n")
        f.write("#  this is not singular i.e. this combination of values is ok."+"\n")
        f.write("# optical species index 1 (the opaque part is indexed by using the - sign)"+"\n")
        f.write("e d "+d['STEM_REFLECTANCE']+"   d "+d['LEAF_REFLECTANCE']+" "+d['LEAF_TRANSMITANCE']+" d "+d['LEAF_REFLECTANCE']+" "+d['LEAF_TRANSMITANCE']+"\n")
        f.write("e d "+d['STEM_REFLECTANCE']+"   d "+d['SEN_LEAF_REFLECTANCE']+" "+d['SEN_LEAF_TRANSMITANCE']+" d "+d['SEN_LEAF_REFLECTANCE']+" "+d['SEN_LEAF_TRANSMITANCE']+"\n")
