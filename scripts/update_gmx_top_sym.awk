#!/usr/bin/awk -f
#
{
  l=length(FILENAME);
  extension=substr(FILENAME,l-2,3);
#
# read symmetry file ###################################################
#
  if(extension=="sym") {
#
    if(NF==1) {
      nasym++;
      symid[nasym]=$1;
    }
  }
#
# read topology file ###################################################
#
  if(extension=="top"||extension=="itp") {
#
    if(substr($1,1,1)!=";") {
#
# remove all comments
#
      i=index($0,";");
      if(i>1) {ltmp=substr($0,1,i-1); $0=ltmp;}
#
# where are we in the topology file?
#
      if($2=="atomtypes") { bt=1;bo=0;bb=0;bp=0;ba=0;bd=0;bi=0; }
      if($2=="atoms")     { bt=0;bo=1;bb=0;bp=0;ba=0;bd=0;bi=0; }
      if($2=="bonds")     { bt=0;bo=0;bb=1;bp=0;ba=0;bd=0;bi=0; }
      if($2=="pairs")     { bt=0;bo=0;bb=0;bp=1;ba=0;bd=0;bi=0; }
      if($2=="angles")    { bt=0;bo=0;bb=0;bp=0;ba=1;bd=0;bi=0; }
      if($2=="dihedrals") { ndih2++;                            }
      if(ndih2==1)        { bt=0;bo=0;bb=0;bp=0;ba=0;bd=1;bi=0; }
      if(ndih2==2)        { bt=0;bo=0;bb=0;bp=0;ba=0;bd=0;bi=1; }
      if($2=="system")    { bt=0;bo=0;bb=0;bp=0;ba=0;bd=0;bi=0; ndih2=0 }
#
# collect topology info on atomtypes
#
      if(bt==1 && NF==7 && $2!="atomtypes") {
        natomtypes++;
        atomtype[natomtypes]=$1;
        sigma[natomtypes]=$6;
        epsilon[natomtypes]=$7;
      }
#
# collect topology info on atoms/names/types/charges
#
      if(bo==1 && NF==8 && $2!="atoms") {
        natoms++;
        type[natoms]=$2;
        atomname[natoms]=$5;
        charge0[natoms]=$7;
        mass[natoms]=$8;
      }
#
# collect topology info on pairs
#
      if(bp==1 && NF==3 && $2!="pairs") {
        npairs++;
        pair1[npairs]=$1;
        pair2[npairs]=$2;
      }
#
# collect topology info on bonds
#

      if(bb==1&&NF==5&&$2!="bonds") {
#
        nbonds++;
        bidx1[nbonds]=$1;
        bidx2[nbonds]=$2;
        beq[nbonds]=$4;
        bkb[nbonds]=$5;

        if(type[$2]>=type[$1]) {
          blabelt[nbonds]=type[$1]"-"type[$2];
          blabeln[nbonds]=atomname[$1]"-"atomname[$2];
 	} else {
          blabelt[nbonds]=type[$2]"-"type[$1];
          blabeln[nbonds]=atomname[$2]"-"atomname[$1];
	}

        if(symid[$2]>=symid[$1]) {
          sid=symid[$1]" "symid[$2];
 	} else {
          sid=symid[$2]" "symid[$1];
	}

        blabels[nbonds]=sid;
        nbon[sid]++;
        bindex[sid"-"nbon[sid]]=nbonds;
      }
#
# collect topology info on angles
#
      if(ba==1&&NF==6&&$2!="angles") {
#
        nangles++;
        aidx1[nangles]=$1;
        aidx2[nangles]=$2;
        aidx3[nangles]=$3;
        aeq[nangles]=$5;
        afc[nangles]=$6;

        if(type[$3]>=type[$1]) {
          alabelt[nangles]=type[$1]"-"type[$2]"-"type[$3];
          alabeln[nangles]=atomname[$1]"-"atomname[$2]"-"atomname[$3];
        } else {
          alabelt[nangles]=type[$3]"-"type[$2]"-"type[$1];
          alabeln[nangles]=atomname[$3]"-"atomname[$2]"-"atomname[$1];
        }

        if(symid[$3]>=symid[$1]) {
          sid=symid[$1]"-"symid[$2]"-"symid[$3];
        } else {
          sid=symid[$3]"-"symid[$2]"-"symid[$1];
        }

        alabels[nangles]=sid;
        nang[sid]++;
        aindex[sid"-"nang[sid]]=nangles;
      }
#
# collect topology info on dihedrals
#
      if(bd==1&&NF==8&&$2!="dihedrals") {
#
        ndihedrals++;
        didx1[ndihedrals]=$1;
        didx2[ndihedrals]=$2;
        didx3[ndihedrals]=$3;
        didx4[ndihedrals]=$4;
#
        dphieq[ndihedrals]=$6;
        dfc[ndihedrals]=$7;
        dphase[ndihedrals]=$8;
#
        if(type[$4]==type[$1]) {
          if(type[$3]>=type[$2]) {
            dlabeln[ndihedrals]=atomname[$1]"-"atomname[$2]"-"atomname[$3]"-"atomname[$4];
            dlabelt[ndihedrals]=type[$1]"-"type[$2]"-"type[$3]"-"type[$4]"-"dphase[ndihedrals];
          } else {
            dlabeln[ndihedrals]=atomname[$4]"-"atomname[$3]"-"atomname[$2]"-"atomname[$1];
            dlabelt[ndihedrals]=type[$4]"-"type[$3]"-"type[$2]"-"type[$1]"-"dphase[ndihedrals];
          }
        } else {
          if(type[$4]>=type[$1]) {
            dlabeln[ndihedrals]=atomname[$1]"-"atomname[$2]"-"atomname[$3]"-"atomname[$4];
            dlabelt[ndihedrals]=type[$1]"-"type[$2]"-"type[$3]"-"type[$4]"-"dphase[ndihedrals];
          } else {
            dlabeln[ndihedrals]=atomname[$4]"-"atomname[$3]"-"atomname[$2]"-"atomname[$1];
            dlabelt[ndihedrals]=type[$4]"-"type[$3]"-"type[$2]"-"type[$1]"-"dphase[ndihedrals];
          }
        }

        if(symid[$4]==symid[$1]) {
          if(symid[$3]>=symid[$2]) {
            sid=symid[$1]"-"symid[$2]"-"symid[$3]"-"symid[$4]"-"dphase[ndihedrals];
          } else {
            sid=symid[$4]"-"symid[$3]"-"symid[$2]"-"symid[$1]"-"dphase[ndihedrals];
          }
        } else {
          if(symid[$4]>=symid[$1]) {
            sid=symid[$1]"-"symid[$2]"-"symid[$3]"-"symid[$4]"-"dphase[ndihedrals];
          } else {
            sid=symid[$4]"-"symid[$3]"-"symid[$2]"-"symid[$1]"-"dphase[ndihedrals];
          }
        }

        dlabels[ndihedrals]=sid;
        ndih[sid]++;
        dindex[sid"-"ndih[sid]]=ndihedrals;
      }

#
# collect topology info on impropers
#
      if(bi==1&&NF==8&&$2!="dihedrals") {
#
        nimpropers++;
        iidx1[nimpropers]=$1;
        iidx2[nimpropers]=$2;
        iidx3[nimpropers]=$3;
        iidx4[nimpropers]=$4;
#
        iphieq[nimpropers]=$6;
        ifc[nimpropers]=$7;
        iphase[nimpropers]=$8;
#
#        if(type[$4]==type[$1]) {
#          if(type[$3]>=type[$2]) {
#            ilabeln[nimpropers]=atomname[$1]"-"atomname[$2]"-"atomname[$3]"-"atomname[$4];
#            ilabelt[nimpropers]=type[$1]"-"type[$2]"-"type[$3]"-"type[$4]"-"iphase[nimpropers];
#          } else {
#            ilabeln[nimpropers]=atomname[$4]"-"atomname[$3]"-"atomname[$2]"-"atomname[$1];
#            ilabelt[nimpropers]=type[$4]"-"type[$3]"-"type[$2]"-"type[$1]"-"iphase[nimpropers];
#          }
#        } else {
#          if(type[$4]>=type[$1]) {
#            ilabeln[nimpropers]=atomname[$1]"-"atomname[$2]"-"atomname[$3]"-"atomname[$4];
#            ilabelt[nimpropers]=type[$1]"-"type[$2]"-"type[$3]"-"type[$4]"-"iphase[nimpropers];
#          } else {
#            ilabeln[nimpropers]=atomname[$4]"-"atomname[$3]"-"atomname[$2]"-"atomname[$1];
#            ilabelt[nimpropers]=type[$4]"-"type[$3]"-"type[$2]"-"type[$1]"-"iphase[nimpropers];
#          }
#        }
#
#        if(symid[$4]==symid[$1]) {
#          if(symid[$3]>=symid[$2]) {
#            sid=symid[$1]"-"symid[$2]"-"symid[$3]"-"symid[$4]"-"dphase[ndihedrals];
#          } else {
#            sid=symid[$4]"-"symid[$3]"-"symid[$2]"-"symid[$1]"-"dphase[ndihedrals];
#          }
#        } else {
#          if(symid[$4]>=symid[$1]) {
#            sid=symid[$1]"-"symid[$2]"-"symid[$3]"-"symid[$4]"-"dphase[ndihedrals];
#          } else {
#            sid=symid[$4]"-"symid[$3]"-"symid[$2]"-"symid[$1]"-"dphase[ndihedrals];
#          }
#        }
#
        ilabeln[nimpropers]=atomname[$1]"-"atomname[$2]"-"atomname[$3]"-"atomname[$4];
        ilabelt[nimpropers]=type[$1]"-"type[$2]"-"type[$3]"-"type[$4]"-"iphase[nimpropers];
        sid=symid[$1]"-"symid[$2]"-"symid[$3]"-"symid[$4]"-"dphase[ndihedrals];
#
        ilabels[nimpropers]=sid;
        nimp[sid]++;
        iindex[sid"-"nimp[sid]]=nimpropers;
#
      }
    }
  }
#
# read charge file ###################################################
#
  if(extension=="chg") {
    if(NF==1) {
      ncharge++;
      charge[ncharge]=$1;
    }
  }
#
# read xyz file(s) ###################################################
#
  if(extension=="xyz") {
#
    if(FNR==1) {
      nxyz++;
      natomsxyz=0;
    }
#
    if(FNR>2) {
      natomsxyz++;
      idx=nxyz"-"natomsxyz;
      x[idx]=$2; y[idx]=$3; z[idx]=$4;
    }
  }
#
#
#
#
#
} END {
#
# write topology file
#
  print "[ defaults ]";
  print "1               2               yes             0.5     0.8333";

  printf "\n";

  print "[ atomtypes ]";
  for(i=1;i<=natomtypes;i++) printf "%4s %4s  0.00  0.00  A  %12.6e %12.6e\n",
    atomtype[i],atomtype[i],sigma[i],epsilon[i];

  printf "\n";

  print "[ nonbond_params ]";
  for(i=1;i<=natomtypes;i++) {
    for(j=i;j<=natomtypes;j++) {
      nopt+=2;
      printf "%4s %4s  1   %12.6e %12.6e ; %d %d\n",
        atomtype[i],  atomtype[j], 0.5*(sigma[i]+sigma[j]), sqrt(epsilon[i]*epsilon[j]), nopt-1, nopt;
    }
  }

  printf "\n";

  print "[ moleculetype ]";
  printf " %s   3\n", resname;

  printf "\n";

  print "[ atoms ]";
  for(i=1;i<=natoms;i++) printf "%8d %4s 1  %3s %4s %8d %12.6f %12.4f\n",
    i, type[i], resname, atomname[i], i, charge[i], mass[i];
#
  printf "\n";
#
### write bonds ##############################################################
#
  for(i=1;i<=nbonds;i++) {
    b1=bidx1[i];
    b2=bidx2[i];
    if(symid[b1]>=symid[b2]) bondsym[i]=symid[b2]"-"symid[b1];
    else                     bondsym[i]=symid[b1]"-"symid[b2];
    nbondsym[bondsym[i]]++;
  }
  nsymbondtyes=0;
  for(i in nbondsym) {nsymbondtyes++; bondopt[i]=nsymbondtyes}
#
  for(i=1;i<=nbonds;i++) {
    s=0.0;
    for(j=1;j<=nxyz;j++) {
      m1=j"-"bidx1[i];
      m2=j"-"bidx2[i];
      s+=0.1*sqrt((x[m1]-x[m2])**2+(y[m1]-y[m2])**2+(z[m1]-z[m2])**2);
    }
    avgbeq[bondsym[i]]+=s/nxyz;
  }
#
  print "[ bonds ]";
  for(sc in nbon) {
    for(i=1;i<=nbon[sc];i++) {
      ii=bindex[sc"-"i];
      printf "%6d %6d  1  %12.4e %12.4e ; %d %s %s %s\n", 
      bidx1[ii], bidx2[ii], avgbeq[bondsym[ii]]/nbondsym[bondsym[ii]] ,bkb[ii], nopt+bondopt[bondsym[ii]], 
	blabeln[ii], blabelt[ii],blabels[ii];
    }
  }
#
  nopt+=nsymbondtyes;
#
  printf "\n";
#
### write pairs ##############################################################
#
  print "[ pairs ]";
  for(i=1;i<=npairs;i++) printf "%6d %6d  1\n",pair1[i],pair2[i];
  printf "\n";
#
### write angles ##############################################################
#
  for(i=1;i<=nangles;i++) {
    a1=aidx1[i];
    a2=aidx2[i];
    a3=aidx3[i];
    if(symid[a1]>=symid[a3]) angsym[i]=symid[a3]"-"symid[a2]"-"symid[a1];
    else                     angsym[i]=symid[a1]"-"symid[a2]"-"symid[a3];
    nangsym[angsym[i]]++;
  }
  nsymangtyes=0;
  for(i in nangsym) {nsymangtyes++; angopt[i]=2*nsymangtyes}
#
  print "[ angles ]";
  for(sc in nang) {
    for(i=1;i<=nang[sc];i++) {
      ii=aindex[sc"-"i];
      printf "%6d %6d %6d 1  %12.4e %12.4e ; %d %d  %s %s %s\n", 
      aidx1[ii], aidx2[ii], aidx3[ii], aeq[ii], afc[ii], nopt+angopt[angsym[ii]]-1, nopt+angopt[angsym[ii]],
	alabeln[ii],alabelt[ii],alabels[ii];
    }
  }
#
  nopt+=2*nsymangtyes;
#
### write dihedrals ##############################################################
#
  for(i=1;i<=ndihedrals;i++) {
    d1=didx1[i];
    d2=didx2[i];
    d3=didx3[i];
    d4=didx4[i];
    if(symid[d1]==symid[d4]) {
      if(symid[d2]>=symid[d3]) dihsym[i]=symid[d4]"-"symid[d3]"-"symid[d2]"-"symid[d1]"-"dphase[i];
      else                     dihsym[i]=symid[d1]"-"symid[d2]"-"symid[d3]"-"symid[d4]"-"dphase[i];
    } else {
      if(symid[d1]>=symid[d4]) dihsym[i]=symid[d4]"-"symid[d3]"-"symid[d2]"-"symid[d1]"-"dphase[i];
      else                     dihsym[i]=symid[d1]"-"symid[d2]"-"symid[d3]"-"symid[d4]"-"dphase[i];
    }
    ndihsym[dihsym[i]]++;
  }
  nsymdihtypes=0;
  for(i in ndihsym) {nsymdihtypes++; dihopt[i]=nsymdihtypes}
#
  printf "\n";
  print "[ dihedrals ]";
  for(sc in ndih) {
    for(i=1;i<=ndih[sc];i++) {
      ii=dindex[sc"-"i];
      printf "%6d %6d %6d %6d  9  %12.2f %12.5f %8d ; %d  %s %s %s\n",
        didx1[ii], didx2[ii], didx3[ii], didx4[ii],
        dphieq[ii],dfc[ii],dphase[ii],
        nopt+dihopt[dihsym[ii]],dlabeln[ii],dlabelt[ii],dlabels[ii];
    }
  }
#
  nopt+=nsymdihtypes;
#
### write impropers ##############################################################
#
  for(i=1;i<=nimpropers;i++) {
    i1=iidx1[i];
    i2=iidx2[i];
    i3=iidx3[i];
    i4=iidx4[i];

#    if(symid[i1]==symid[i4]) {
#      if(symid[i2]>=symid[i3]) impsym[i]=symid[i4]"-"symid[i3]"-"symid[i2]"-"symid[i1]"-"iphase[i];
#      else                     impsym[i]=symid[i1]"-"symid[i2]"-"symid[i3]"-"symid[i4]"-"iphase[i];
#    } else {
#      if(symid[i1]>=symid[i4]) impsym[i]=symid[i4]"-"symid[i3]"-"symid[i2]"-"symid[i1]"-"iphase[i];
#      else                     impsym[i]=symid[i1]"-"symid[i2]"-"symid[i3]"-"symid[i4]"-"iphase[i];
#    }

    impsym[i]=symid[i1]"-"symid[i2]"-"symid[i3]"-"symid[i4]"-"iphase[i];
    nimpsym[impsym[i]]++;
  }

  nsymimptypes=0;
  for(i in nimpsym) {nsymimptypes++; impopt[i]=nsymimptypes}
#
  printf "\n";
  print "[ dihedrals ]";
  for(sc in nimp) {
    for(i=1;i<=nimp[sc];i++) {
      ii=iindex[sc"-"i];
      printf "%6d %6d %6d %6d      9  %12.2f %12.5f %8d ; %d  %s %s %s\n",
        iidx1[ii], iidx2[ii], iidx3[ii], iidx4[ii],
        iphieq[ii],ifc[ii],iphase[ii],
        nopt+impopt[impsym[ii]],ilabeln[ii],ilabelt[ii],ilabels[ii];
    }
  }

  printf "\n";

  print "[ system ]";
  print "   system";

  printf "\n";

  print "[ molecules ]";
  printf " %s n_res\n\n", resname;
}
