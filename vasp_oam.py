#!/usr/bin/env python
####################################################################################################
#####                   Calculating the oam from the PROCAR file and plotting the figure       #####
#####                       Qun Yang (Qun.Yang@cpfs.mpg.de)                                    #####
#####                              date: 8 of November 2021                                    #####
####################################################################################################
import os
import numpy as np
import matplotlib.pyplot as plt

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

def main():
    writefile()
    #plot()

def l_matrix():
   # matrix_lx = np.zeros([9,9], dtype=complex)
   # matrix_lx[1,2]=complex(0.0, 1.0)
   # matrix_lx[2,1]=complex(0.0,-1.0)
   # matrix_lx[4,7]=complex(0.0, 1.0)
   # matrix_lx[7,4]=complex(0.0,-1.0)
   # matrix_lx[5,6]=complex(0.0, np.sqrt(3.0))
   # matrix_lx[6,5]=complex(0.0,-np.sqrt(3.0))
   # matrix_lx[5,8]=complex(0.0, 1.0)
   # matrix_lx[8,5]=complex(0.0,-1.0)

   # matrix_ly = np.zeros([9,9], dtype=complex)
   # matrix_ly[2,3]=complex(0.0,-1.0)
   # matrix_ly[3,3]=complex(0.0, 1.0)
   # matrix_ly[4,5]=complex(0.0, 1.0)
   # matrix_ly[5,4]=complex(0.0,-1.0)
   # matrix_ly[6,7]=complex(0.0,-np.sqrt(3.0))
   # matrix_ly[7,6]=complex(0.0, np.sqrt(3.0))
   # matrix_ly[7,8]=complex(0.0, -1.0)
   # matrix_ly[8,7]=complex(0.0,  1.0)

    matrix_lz = np.zeros([9,9], dtype=complex)
    matrix_lz[1,3]=complex(0.0, 1.0)
    matrix_lz[3,1]=complex(0.0,  -1.0)
    matrix_lz[4,8]=complex(0.0, 2.0)
    matrix_lz[8,4]=complex(0.0,  -2.0)
    matrix_lz[5,7]=complex(0.0, 1.0)
    matrix_lz[7,5]=complex(0.0,  -1.0)
    return matrix_lz
#### read OUTCAR and PROCAR ####

def readfiles():

    matrix_lz = l_matrix()

    #######################
    ####   read OUTCAR ####
    #######################

    print('-------- READ E-fermi FROM OUTCAR.scf --------')
    with open("OUTCAR.scf","r") as fp:
        fo = fp.readlines()
    for cnt, line in enumerate(fo):
        ######## fermi energy ########
        if any(line.split()) and line.split()[0] == "E-fermi" :
            Efermi = float(line.split()[2])

    with open('KPOINTS', 'r') as fp:
        fi = fp.readlines()
    ######## k-mesh ########
    mesh  = int(fi[1].strip().split()[0]) # 1001

    print('-------- READ OUTCAR --------')
    with open("OUTCAR","r") as fp:
        fr = fp.readlines()
    pra = []
    prb = []
    prc = []
    for cnt, line in enumerate(fr):
        ######## number of kpoints ########
        if line[22:27] == 'NKPTS':
            NK = int(line.split()[3])
            #print(' NK = ',NK)
            NHSKP = int(NK/mesh)

        ######## number of bands ########
        if line[94:100] == "NBANDS" :
            NBANDS = int(line[101:108])
            #print(' NBANDS = ',NBANDS)

        ######## number of ions ########
        if line[58:63] == "NIONS":
            NATOMS = int( line.split()[11] )
            #print(' NATOMS =', NATOMS)
 
        ######## reciprocal lattice vectors ########
        if line[45:71] == "reciprocal lattice vectors" and not pra:
            cntp = cnt+1
            pra = [ float( fr[cntp].split()[3]), float( fr[cntp].split()[4]), float( fr[cntp].split()[5])]
            cntp = cntp+1
            prb = [ float( fr[cntp].split()[3]), float( fr[cntp].split()[4]), float( fr[cntp].split()[5])]
            cntp = cntp+1
            prc = [ float( fr[cntp].split()[3]), float( fr[cntp].split()[4]), float( fr[cntp].split()[5])]
    
    
        ######## bands ########
        if any(line.split()) and line.split()[0] == "E-fermi" :
            kpoint   = [ [0.0, 0.0, 0.0] for ii in range(NK) ] 
            enk      = [ [0.0]*NK for ii in range(NBANDS) ]

            cntp = cnt + 1
            for ii in range(NK):
                cntp = cntp + 2
                kpoint[ii] = [ float( fr[cntp].split()[3]), float( fr[cntp].split()[4]), float( fr[cntp].split()[5])]    
                cntp = cntp + 1
                for jj in range(NBANDS):
                    cntp = cntp + 1
                    enk[jj][ii] = float( fr[cntp].split()[1] )
    ######## k-path ########
    vkxyz   = [ [0.0, 0.0, 0.0] for ii in range(NK) ] 
    lenk    = [ 0.0 for ii in range(NK) ]
    PR = np.array([pra, prb, prc])
    for cnt, kp in enumerate(kpoint) :
        vkxyz[cnt]   = np.dot( np.array(kp), PR ) 
        if cnt >= 1:
            lenk[cnt]   = np.linalg.norm( vkxyz[cnt] - vkxyz[cnt-1]   ) + lenk[cnt-1]
    
    enk = np.array(enk) - Efermi

            
        #for i in range(0,NK):
        #    fw.write("%14.8f" % lenk[i])
        #    for j in range(NBANDS):
        #        fw.write("%14.8f" % enk[j][i])
        #    fw.write("\n")

    #########################
    ####   read PROCAR   ####
    #########################

    print('-------- READ PROCAR AND DO CALCULATION --------')

    #lx = np.zeros([NK, NBANDS])
    #ly = np.zeros([NK, NBANDS])
    lz = np.zeros([NK, NBANDS])

    with open('PROCAR', 'r') as fp:
        line = fp.readline()
        for ii in range(NK) :
            line = fp.readline() 
            line = fp.readline() 
            line = fp.readline() 
            for jj in range(NBANDS):
                line = fp.readline()
                line = fp.readline()
                line = fp.readline()
                line = fp.readline()
                for kk in range(NATOMS):
                    line = fp.readline()
                if NATOMS > 1 :
                    line = fp.readline()
                # skip x component
                for kk in range(NATOMS):
                    line = fp.readline()
                if NATOMS > 1 :
                    line = fp.readline()
                # skip y component
                for kk in range(NATOMS):
                    line = fp.readline()
                if NATOMS > 1 :
                    line = fp.readline()
                # skip z component
                for kk in range(NATOMS):
                    line = fp.readline()
                if NATOMS > 1 :
                    line = fp.readline()
                line = fp.readline()
            
                for kk in range(NATOMS):
                    fai0 = np.zeros([1,9], dtype=complex)
                    line = fp.readline()
                    line = line.split()
                    for hh in range(9):
                        fai0[0,hh] = complex(float(line[hh+1+hh*1]),float(line[hh+2+hh*1]))
                    #lx[ii,jj] = lx[ii,jj] + np.dot(np.dot(np.conjugate(fai0),matrix_lx),np.transpose(fai0)).real[0][0]
                    #lx[ii,jj] = lx[ii,jj] + np.dot(np.dot(np.conjugate(fai0),matrix_ly),np.transpose(fai0)).real[0][0]
                    lz[ii][jj] = lz[ii][jj] + np.dot(np.dot(np.conjugate(fai0),matrix_lz),np.transpose(fai0)).real[0][0]
                line = fp.readline()
    return NK, NHSKP, NBANDS, lenk, enk, lz

def writefile():
    NK, NHSKP, NBANDS, lenk, enk, lz = readfiles()
    with open("lz.dat", "w") as fw:
        for j in range(NBANDS):
            for i in range(0,NK):
                #fw.write("%14.8f%14.8f%14.8f%14.8f%14.8f" % (lenk[i]*2*np.pi, enk[j][i], lx[i,j].real, ly[i,j].real, lz[i,j]))
                fw.write("%14.8f%14.8f%14.8f" % (lenk[i]*2*np.pi, enk[j][i], lz[i][j]))
                fw.write("\n")
            fw.write("\n")

def plot():
    NK, NHSKP, NBANDS, lenk, enk, lz = readfiles()
    mesh = int(NK/NHSKP)
    ###### read KPOINTS ######
    labels = []
    ind = []
    with open('KPOINTS','r') as fp:
        fileread = fp.readline()
        fileread = fp.readline()
        fileread = fp.readline()
        fileread = fp.readline()

        fileread = fp.readline()
        labels.append('$'+fileread.split('!')[1].strip()+'$')
        fileread = fp.readline()
        kp = str(fileread.split('!')[1].strip())
        labels.append('$'+fileread.split('!')[1].strip()+'$')

        for ii in range(NHSKP-1):
            fileread = fp.readline()
            fileread = fp.readline()
            if str(fileread.split('!')[1].strip()) == kp:
                fileread = fp.readline()
                labels.append('$'+fileread.split('!')[1].strip()+'$')
                kp = str(fileread.split('!')[1].strip())
            else:
                labels.pop()
                labels.append('$'+str(kp)+'|'+fileread.split('!')[1].strip()+'$')
                fileread = fp.readline()
                labels.append('$'+fileread.split('!')[1].strip()+'$')
    ###### read KPOINTS ######

    ######## x index ########
    ind.append(0)
    for ii in range(NHSKP):
        ind.append(lenk[(ii+1)*mesh-1] )
    ######## x index ########
    
    ####### plot pictures ######
    plt.figure(figsize=(8,6))

    miny = -2 
    maxy = 2

    for jj in range(NBANDS):
        p = plt.scatter( lenk, enk[jj], s= (lz[:,jj])*20, marker='o',linewidth='0.5', edgecolors = 'blue', facecolors='None')

    plt.plot( [lenk[0], lenk[-1]], [0.00, 0.00],color='black',linestyle='--',linewidth=0.8 )
    for ii in range(NHSKP):
        plt.plot( [lenk[mesh*(ii+1)-1], lenk[mesh*(ii+1)-1]], [miny, maxy],color='black',linewidth=0.8)

    plt.xlim(lenk[0], lenk[-1])
    plt.ylim(miny, maxy)              
    plt.xticks(ind, labels)      
    plt.ylabel('Energy (eV)', fontsize=12)
    plt.title("OAM_Lz", fontsize=15, fontweight='bold')
    plt.savefig('band.png', dpi=300)

if __name__=='__main__':
    main()
