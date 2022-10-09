#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"
#include "vec.h"


int fitness(const t_parameters* pars, t_genotype* genome, t_topology* top, t_codonlimits* limits, double** coords,
            double* qme, double* w_qme, double* mme, double* mme0) {

  int i,i1,i2,i3,i4,ii,jj,nn,rnn,n_conf,n_mol,n_bond,n_ang,n_dih,n_imp;
  int mmeidx1,mmeidx2,mmeidx3,mmeidx4;

  double sig,eps,sig2,sig6,sig12;
  double r12x,r12y,r12z;
  double r12sq,oor2,oor6,oor12;
  double pot, dpot;
  double f1x,f1y,f1z;
  double f2x,f2y,f2z;
  double f3x,f3y,f3z;
  double kf,a0,diff,c1,c2,r0;
  double dp,fc,dn;
  double r32x,r32y,r32z;
  double br12,br32,cos_theta,sin_theta,angle;
  double rijx,rijy,rijz ,rkjx,rkjy,rkjz,rlkx ,rlky,rlkz ,rklx ,rkly,rklz;
  double brkj,nrkjx,nrkjy,nrkjz,rijnrkj,rlknrkj,Rx,Ry,Rz,Sx,Sy,Sz,bR,bS;
  double nRx,nRy,nRz,nSx,nSy,nSz,rkjrklx,rkjrkly ,rkjrklz,mx,my,mz;
  double nx,ny,nz,mcrossnx,mcrossny,mcrossnz,rkjdotmcrossn,bm2,bn2;
  double fix,fiy,fiz,flx,fly,flz;
  double rijrkjorkj2,rklrkjorkj2,fjx,fjy,fjz,fkx,fky,fkz;
  double phi,sign;
  double sume,avge,f;
  double vdwcosq;
  double acosarg;
  double sum,f1;

  static char dovdw[]="v";
  static char dobon[]="b";
  static char doang[]="a";
  static char dodih[]="d";
  static char doimp[]="i";

/***** set all energies and forces to their initial values (the constant not optimized part) *****/

  for(i=0;i< pars->n_qme+pars->n_qmf;i++) mme[i]=mme0[i];

/***** LJ interactions *****/

  if(strstr(pars->opt,dovdw)!=NULL) {

    vdwcosq = (pars->vdw_cutoff)*(pars->vdw_cutoff);

/***** calc factors for LJ pot *****/

    for(nn=0;nn<top->n_neighbors;nn++) {    
      ii = top->neighbor1[nn]; while(ii>=top->n_atoms) ii -= top->n_atoms;
      jj = top->neighbor2[nn]; while(jj>=top->n_atoms) jj -= top->n_atoms;
      rnn=top->nbidx[nn];
      sig=top->sigma[rnn];
      eps=top->epsilon[rnn];
      if(top->sgrpidx[rnn]!=0 ) sig = genome->codon[limits->codonidx[top->sgrpidx[rnn]]];
      if(top->egrpidx[rnn]!=0 ) eps = genome->codon[limits->codonidx[top->egrpidx[rnn]]];
      sig2=sig*sig;
      sig6=sig2*sig2*sig2;
      sig12=sig6*sig6;
      top->lj12[nn] = top->nb_lj_fac[nn]*4.0*eps*sig12;
      top->lj6[nn]  = top->nb_lj_fac[nn]*4.0*eps*sig6;
    }

/***** loop over all conformers *****/

    for(n_conf=0;n_conf<pars->n_conf;n_conf++) {

      for(nn=0;nn<top->n_neighbors;nn++) {

        i1 = top->neighbor1[nn]+n_conf*top->n_atoms*pars->n_mol;
        i2 = top->neighbor2[nn]+n_conf*top->n_atoms*pars->n_mol;

        r12x = coords[i1][XX]-coords[i2][XX];
        r12y = coords[i1][YY]-coords[i2][YY];
        r12z = coords[i1][ZZ]-coords[i2][ZZ];
        r12sq = r12x*r12x+r12y*r12y+r12z*r12z;

        if(r12sq<vdwcosq) {

          br12 = sqrt(r12sq);
          oor2=1.0/r12sq;
          oor6=oor2*oor2*oor2;
          oor12=oor6*oor6;

          pot  =       top->lj12[nn]*oor12      -     top->lj6[nn]*oor6 ;
          dpot =  12.0*top->lj12[nn]*oor12/br12 - 6.0*top->lj6[nn]*oor6/br12; 

          dpot /= br12;
          f1x =  (dpot)*r12x;
          f1y =  (dpot)*r12y;
          f1z =  (dpot)*r12z;
          f2x = -(dpot)*r12x;
          f2y = -(dpot)*r12y;
          f2z = -(dpot)*r12z;
  
          mmeidx1 = pars->n_qme+3*i1;
          mmeidx2 = pars->n_qme+3*i2;

          if(pars->n_qme>0) mme[n_conf]    += pot;
          if(pars->n_qmf>0) {
            mme[mmeidx1]   += f1x;
            mme[mmeidx1+1] += f1y;
            mme[mmeidx1+2] += f1z;
            mme[mmeidx2]   += f2x;
            mme[mmeidx2+1] += f2y;
            mme[mmeidx2+2] += f2z;
          }
        }
      }
    }
  }

/***** bonds *****/

  if(strstr(pars->opt,dobon)!=NULL) {

    for(n_conf=0;n_conf<pars->n_conf;n_conf++) {

      for(n_bond=0;n_bond<top->n_bonds;n_bond++) {
 
        if(top->isconstrained[n_bond]!=1 && top->bgrpidx[n_bond]!=0) {

          for(n_mol=0;n_mol<pars->n_mol;n_mol++) {

            i1 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->bondi[n_bond] - 1;
            i2 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->bondj[n_bond] - 1;
            r0 = top->bondr[n_bond];
            fc = top->bondk[n_bond];
            if(top->bgrpidx[n_bond]!=0) fc = genome->codon[limits->codonidx[top->bgrpidx[n_bond]]];

            r12x = coords[i1][XX]-coords[i2][XX];
            r12y = coords[i1][YY]-coords[i2][YY];
            r12z = coords[i1][ZZ]-coords[i2][ZZ];
            br12 = sqrt(r12x*r12x+r12y*r12y+r12z*r12z);

            dpot = -fc*(br12-r0);
            pot = 0.5*fc*(br12-r0)*(br12-r0);

            dpot /= br12;
            f1x =  dpot*r12x;
            f1y =  dpot*r12y;
            f1z =  dpot*r12z;
            f2x = -dpot*r12x;
            f2y = -dpot*r12y;
            f2z = -dpot*r12z;

            mmeidx1 = pars->n_qme+3*i1;
            mmeidx2 = pars->n_qme+3*i2;

            if(pars->n_qme>0) mme[n_conf] += pot;
            if(pars->n_qmf>0) {
              mme[mmeidx1]   += f1x;
              mme[mmeidx1+1] += f1y;
              mme[mmeidx1+2] += f1z;
              mme[mmeidx2]   += f2x;
              mme[mmeidx2+1] += f2y;
              mme[mmeidx2+2] += f2z;
            }
          }
        }
      }
    }
  }

/***** angles *****/

  if(strstr(pars->opt,doang)!=NULL) {

    for(n_conf=0;n_conf<pars->n_conf;n_conf++) {

      for(n_ang=0;n_ang<top->n_angles;n_ang++) {

        for(n_mol=0;n_mol<pars->n_mol;n_mol++) {

          i1 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->anglei[n_ang] - 1;
          i2 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->anglej[n_ang] - 1;
          i3 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->anglek[n_ang] - 1;

          a0 = top->angleeq[n_ang];
          kf = top->anglefc[n_ang];
          if(top->cgrpidx[n_ang]!=0) kf = genome->codon[limits->codonidx[top->cgrpidx[n_ang]]];
          if(top->agrpidx[n_ang]!=0) a0 = genome->codon[limits->codonidx[top->agrpidx[n_ang]]];
          a0=a0/180.0*M_PI;

          r12x = coords[i1][XX]-coords[i2][XX];
          r12y = coords[i1][YY]-coords[i2][YY];
          r12z = coords[i1][ZZ]-coords[i2][ZZ];
          r32x = coords[i3][XX]-coords[i2][XX];
          r32y = coords[i3][YY]-coords[i2][YY];
          r32z = coords[i3][ZZ]-coords[i2][ZZ];

          br12 = sqrt(r12x*r12x+r12y*r12y+r12z*r12z);
          br32 = sqrt(r32x*r32x+r32y*r32y+r32z*r32z);
          cos_theta = (r12x*r32x+r12y*r32y+r12z*r32z)/(br12*br32);
          if (cos_theta > 1.0) cos_theta       =  1.0;
          else if (cos_theta < -1.0) cos_theta = -1.0;
          angle = acos(cos_theta);

          diff = angle - a0;
          pot = 0.5*kf*diff*diff;
          sin_theta = sqrt(1.0 - cos_theta*cos_theta);
          diff *= -kf/sin_theta;
          c1 = diff/br12;
          c2 = diff/br32;

          f1x = c1*(r12x*(cos_theta/br12) - r32x/br32);
          f1y = c1*(r12y*(cos_theta/br12) - r32y/br32);
          f1z = c1*(r12z*(cos_theta/br12) - r32z/br32);

          f3x = c2*(r32x*(cos_theta/br32) - r12x/br12);
          f3y = c2*(r32y*(cos_theta/br32) - r12y/br12);
          f3z = c2*(r32z*(cos_theta/br32) - r12z/br12);

          f2x = -f1x-f3x;
          f2y = -f1y-f3y;
          f2z = -f1z-f3z;

          mmeidx1 = pars->n_qme+3*i1;
          mmeidx2 = pars->n_qme+3*i2;
          mmeidx3 = pars->n_qme+3*i3;

          if(pars->n_qme>0) mme[n_conf] += pot;
          if(pars->n_qmf>0) {
            mme[mmeidx1]   += f1x;
            mme[mmeidx1+1] += f1y;
            mme[mmeidx1+2] += f1z;
            mme[mmeidx2]   += f2x;
            mme[mmeidx2+1] += f2y;
            mme[mmeidx2+2] += f2z;
            mme[mmeidx3]   += f3x;
            mme[mmeidx3+1] += f3y;
            mme[mmeidx3+2] += f3z;
          }
        }
      }
    }
  }

/***** dihedrals *****/

  if(strstr(pars->opt,dodih)!=NULL) {

    for(n_conf=0;n_conf<pars->n_conf;n_conf++) {

      for(n_dih=0;n_dih<top->n_dihedrals;n_dih++) {

        for(n_mol=0;n_mol<pars->n_mol;n_mol++) {

          i1 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->dihi[n_dih] - 1;
          i2 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->dihj[n_dih] - 1;
          i3 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->dihk[n_dih] - 1;
          i4 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->dihl[n_dih] - 1;

          dp=top->dihphase[n_dih]/180.0*M_PI;
          fc=top->dihfc[n_dih];
          dn=top->dihpn[n_dih];
          if(top->dgrpidx[n_dih]!=0) fc = genome->codon[limits->codonidx[top->dgrpidx[n_dih]]];

          rijx = coords[i2][XX]-coords[i1][XX];
          rijy = coords[i2][YY]-coords[i1][YY];
          rijz = coords[i2][ZZ]-coords[i1][ZZ];
          rkjx = coords[i2][XX]-coords[i3][XX];
          rkjy = coords[i2][YY]-coords[i3][YY];
          rkjz = coords[i2][ZZ]-coords[i3][ZZ];
          rlkx = coords[i3][XX]-coords[i4][XX];
          rlky = coords[i3][YY]-coords[i4][YY];
          rlkz = coords[i3][ZZ]-coords[i4][ZZ];
          rklx = -rlkx;
          rkly = -rlky;
          rklz = -rlkz;

          brkj = sqrt(rkjx*rkjx+rkjy*rkjy+rkjz*rkjz);
          nrkjx = rkjx/brkj;
          nrkjy = rkjy/brkj;
          nrkjz = rkjz/brkj;
          rijnrkj = rijx*nrkjx+rijy*nrkjy+rijz*nrkjz;
          rlknrkj = rlkx*nrkjx+rlky*nrkjy+rlkz*nrkjz;

          Rx = rijx - rijnrkj*nrkjx;
          Ry = rijy - rijnrkj*nrkjy;
          Rz = rijz - rijnrkj*nrkjz;
          Sx = rlkx - rlknrkj*nrkjx;
          Sy = rlky - rlknrkj*nrkjy;
          Sz = rlkz - rlknrkj*nrkjz;
          bR = sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
          bS = sqrt(Sx*Sx+Sy*Sy+Sz*Sz);

          nRx = Rx/bR;
          nRy = Ry/bR;
          nRz = Rz/bR;
          nSx = Sx/bS;
          nSy = Sy/bS;
          nSz = Sz/bS;

          rkjrklx = rkjy*rklz - rkjz*rkly;
          rkjrkly = rkjx*rklz - rkjz*rklx;
          rkjrklz = rkjx*rklz - rkjz*rklx;
          mx=rijy*rkjz-rijz*rkjy;
          my=rijx*rkjz-rijz*rkjx;
          mz=rijx*rkjy-rijy*rkjx;
          nx=rkjy*rklz-rkjz*rkly;
          ny=rkjx*rklz-rkjz*rklx;
          nz=rkjx*rkly-rkjy*rklx;
          mcrossnx = my*nz-mz*ny;
          mcrossny = mx*nz-mz*nx;
          mcrossnz = mx*ny-my*nx;
          rkjdotmcrossn = rkjx*mcrossnx+rkjy*mcrossny+rkjz*mcrossnz;

          if(rkjdotmcrossn>=0.0) sign= 1.0; else sign=-1.0;

          acosarg=nRx*nSx+nRy*nSy+nRz*nSz;
          if(acosarg>1.0) acosarg=1.0;
          if(acosarg<-1.0) acosarg=-1.0;
          phi = sign*acos(acosarg);

          pot = fc*(1.0+cos(dn*phi - dp));
          dpot = -dn*fc*sin(dn*phi - dp);

          bm2 = mx*mx+my*my+mz*mz;
          bn2 = nx*nx+ny*ny+nz*nz;

          fix = -dpot * brkj/bm2*mx;
          fiy = -dpot * brkj/bm2*my;
          fiz = -dpot * brkj/bm2*mz;
          flx =  dpot * brkj/bn2*nx;
          fly =  dpot * brkj/bn2*ny;
          flz =  dpot * brkj/bn2*nz;

          rijrkjorkj2 = rijnrkj/brkj;
          rklrkjorkj2 = (rklx*rkjx+rkly*rkjy+rklz*rkjz)/brkj/brkj;

          fjx = -fix + rijrkjorkj2 * fix - rklrkjorkj2*flx;
          fjy = -fiy + rijrkjorkj2 * fiy - rklrkjorkj2*fly;
          fjz = -fiz + rijrkjorkj2 * fiz - rklrkjorkj2*flz;
          fkx = -flx - rijrkjorkj2 * fix + rklrkjorkj2*flx;
          fky = -fly - rijrkjorkj2 * fiy + rklrkjorkj2*fly;
          fkz = -flz - rijrkjorkj2 * fiz + rklrkjorkj2*flz;

          mmeidx1 = pars->n_qme+3*i1;
          mmeidx2 = pars->n_qme+3*i2;
          mmeidx3 = pars->n_qme+3*i3;
          mmeidx4 = pars->n_qme+3*i4;

          if(pars->n_qme>0) mme[n_conf]    += pot;
          if(pars->n_qmf>0) {
            mme[mmeidx1]   += fix;
            mme[mmeidx1+1] -= fiy;
            mme[mmeidx1+2] += fiz;
            mme[mmeidx2]   += fjx;
            mme[mmeidx2+1] -= fjy;
            mme[mmeidx2+2] += fjz;
            mme[mmeidx3]   += fkx;
            mme[mmeidx3+1] -= fky;
            mme[mmeidx3+2] += fkz;
            mme[mmeidx4]   += flx;
            mme[mmeidx4+1] -= fly;
            mme[mmeidx4+2] += flz;
          }
        }
      }
    }
  }

/***** impropers *****/

  if(strstr(pars->opt,doimp)!=NULL) {

    for(n_conf=0;n_conf<pars->n_conf;n_conf++) {

      for(n_imp=0;n_imp<top->n_impropers;n_imp++) {

        for(n_mol=0;n_mol<pars->n_mol;n_mol++) {

          i1 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->impi[n_imp] - 1;
          i2 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->impj[n_imp] - 1;
          i3 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->impk[n_imp] - 1;
          i4 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->impl[n_imp] - 1;

          dp=top->impphase[n_imp]/180.0*M_PI;
          fc=top->impfc[n_imp];
          dn=top->imppn[n_imp];
          if(top->igrpidx[n_imp]!=0) fc = genome->codon[limits->codonidx[top->igrpidx[n_imp]]];

          rijx = coords[i2][XX]-coords[i1][XX];
          rijy = coords[i2][YY]-coords[i1][YY];
          rijz = coords[i2][ZZ]-coords[i1][ZZ];
          rkjx = coords[i2][XX]-coords[i3][XX];
          rkjy = coords[i2][YY]-coords[i3][YY];
          rkjz = coords[i2][ZZ]-coords[i3][ZZ];
          rlkx = coords[i3][XX]-coords[i4][XX];
          rlky = coords[i3][YY]-coords[i4][YY];
          rlkz = coords[i3][ZZ]-coords[i4][ZZ];
          rklx = -rlkx;
          rkly = -rlky;
          rklz = -rlkz;

          brkj = sqrt(rkjx*rkjx+rkjy*rkjy+rkjz*rkjz);
          nrkjx = rkjx/brkj;
          nrkjy = rkjy/brkj;
          nrkjz = rkjz/brkj;
          rijnrkj = rijx*nrkjx+rijy*nrkjy+rijz*nrkjz;
          rlknrkj = rlkx*nrkjx+rlky*nrkjy+rlkz*nrkjz;
          Rx = rijx - rijnrkj*nrkjx;
          Ry = rijy - rijnrkj*nrkjy;
          Rz = rijz - rijnrkj*nrkjz;
          Sx = rlkx - rlknrkj*nrkjx;
          Sy = rlky - rlknrkj*nrkjy;
          Sz = rlkz - rlknrkj*nrkjz;
          bR = sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
          bS = sqrt(Sx*Sx+Sy*Sy+Sz*Sz);
          nRx = Rx/bR;
          nRy = Ry/bR;
          nRz = Rz/bR;
          nSx = Sx/bS;
          nSy = Sy/bS;
          nSz = Sz/bS;

          rkjrklx = rkjy*rklz - rkjz*rkly;
          rkjrkly = rkjx*rklz - rkjz*rklx;
          rkjrklz = rkjx*rklz - rkjz*rklx;
          mx=rijy*rkjz-rijz*rkjy;
          my=rijx*rkjz-rijz*rkjx;
          mz=rijx*rkjy-rijy*rkjx;
          nx=rkjy*rklz-rkjz*rkly;
          ny=rkjx*rklz-rkjz*rklx;
          nz=rkjx*rkly-rkjy*rklx;
          mcrossnx = my*nz-mz*ny;
          mcrossny = mx*nz-mz*nx;
          mcrossnz = mx*ny-my*nx;
          rkjdotmcrossn = rkjx*mcrossnx+rkjy*mcrossny+rkjz*mcrossnz;
          if(rkjdotmcrossn>=0.0) sign= 1.0; else sign=-1.0;

          acosarg=nRx*nSx+nRy*nSy+nRz*nSz;
          if(acosarg>1.0) acosarg=1.0;
          if(acosarg<-1.0) acosarg=-1.0;
          phi = sign*acos(acosarg);

          pot = fc*(1.0+cos(dn*phi - dp));
          dpot = -dn*fc*sin(dn*phi - dp);

          bm2 = mx*mx+my*my+mz*mz;
          bn2 = nx*nx+ny*ny+nz*nz;
          fix = -dpot * brkj/bm2*mx;
          fiy = -dpot * brkj/bm2*my;
          fiz = -dpot * brkj/bm2*mz;
          flx =  dpot * brkj/bn2*nx;
          fly =  dpot * brkj/bn2*ny;
          flz =  dpot * brkj/bn2*nz;
          rijrkjorkj2 = rijnrkj/brkj;
          rklrkjorkj2 = (rklx*rkjx+rkly*rkjy+rklz*rkjz)/brkj/brkj;
          fjx = -fix + rijrkjorkj2 * fix - rklrkjorkj2*flx;
          fjy = -fiy + rijrkjorkj2 * fiy - rklrkjorkj2*fly;
          fjz = -fiz + rijrkjorkj2 * fiz - rklrkjorkj2*flz;
          fkx = -flx - rijrkjorkj2 * fix + rklrkjorkj2*flx;
          fky = -fly - rijrkjorkj2 * fiy + rklrkjorkj2*fly;
          fkz = -flz - rijrkjorkj2 * fiz + rklrkjorkj2*flz;

          mmeidx1 = pars->n_qme+3*i1;
          mmeidx2 = pars->n_qme+3*i2;
          mmeidx3 = pars->n_qme+3*i3;
          mmeidx4 = pars->n_qme+3*i4;

          if(pars->n_qme>0) mme[n_conf]    += pot;
          if(pars->n_qmf>0) {
            mme[mmeidx1]   += fix;
            mme[mmeidx1+1] -= fiy;
            mme[mmeidx1+2] += fiz;
            mme[mmeidx2]   += fjx;
            mme[mmeidx2+1] -= fjy;
            mme[mmeidx2+2] += fjz;
            mme[mmeidx3]   += fkx;
            mme[mmeidx3+1] -= fky;
            mme[mmeidx3+2] += fkz;
            mme[mmeidx4]   += flx;
            mme[mmeidx4+1] -= fly;
            mme[mmeidx4+2] += flz;
          }
        }
      }
    }
  }

/***** scale mme *****/

  if(pars->n_qme>0) {
    sume = 0.0;
    for(i=0; i<pars->n_qme; i++) sume += mme[i];
    avge=sume/(double)pars->n_qme;
  }

/***** calculate fitness *****/

  f = 0.0;
  sum = 0.0;
  for(i=0;i<pars->n_qme;i++) {
    f1 = qme[i]-(mme[i]-avge)/pars->stdeve;

    f += f1*f1*w_qme[i];
    sum += w_qme[i];
  }
  for(i=pars->n_qme; i<pars->n_qme+pars->n_qmf;i++) {
    f1 = (qme[i]-mme[i]/pars->stdevf);

    f += f1*f1*w_qme[i];
    sum += w_qme[i];
  }
  genome->fitness = sqrt(f/sum);

  return(1);
}
