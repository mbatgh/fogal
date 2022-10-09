

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"
#include "vec.h"


int setup_mm(const char* filename, const t_parameters* pars, t_genotype* genome, t_topology* top, 
             t_codonlimits* limits, double** coords, double* qme, double* w_qme, double* mme, double* mme0) {

  int i,j,k;
  int ii,jj,nn,rnn,nread;
  int i1, i2, i3, i4,n_conf, n_mol, n_bond, n_ang, n_dih, n_imp;
  int bbond, bangle, bpair;
  int molnr1, molnr2;
  int mmeidx1,mmeidx2,mmeidx3,mmeidx4;

  double r0,kb,debonds;
  double dp,fc,dn;
  double a0,ka;
  double angle;
  double phi;
  double oor2, oor6, oor12, sig, eps, sig2, sig6, sig12;
  double ecoul,evdw,ebond,eimp,edih,eang;
  double oo4pe0 = 138.935458; /* = 1.0/(4.0*M_PI*eps0) in kJ/mol nm e^-2*/
  double vdwcosq;
  double mmestdev;
  double acosarg;
  double sum;
  double f1;

  char dovdw[]="v";
  char dobon[]="b";
  char doang[]="a";
  char dodih[]="d";
  char doimp[]="i";

  int** nbidx;

  double rijx,rijy,rijz,rkjx,rkjy,rkjz,rlkx,rlky,rlkz,rklx,rkly,rklz;
  double brkj;
  double nrkjx,nrkjy,nrkjz,rijnrkj,rlknrkj,Rx,Ry,Rz,Sx,Sy,Sz,bR,bS,nRx,nRy,nRz;
  double nSx,nSy,nSz,dih_pot,dvdphi,mx,my,mz,nx,ny,nz;
  double bm2,bn2,fix,fiy,fiz,flx,fly,flz,rijrkjorkj2,rklrkjorkj2,fjx,fjy,fjz,fkx,fky,fkz;
  double rkjrklx,rkjrkly,rkjrklz;
  double mcrossnx,mcrossny ,mcrossnz,rkjdotmcrossn,sign;
  double imp_pot,ang_pot,bon_pot,cou_pot,vdw_pot;

  double cos_theta,sin_theta,c1,c2,diff,f1x,f1y,f1z,f3x,f3y,f3z,f2x,f2y,f2z;
  double br12,br32;
  double r12sq,dcou_pot,dvdw_pot;
  double r12x,r12y,r12z,r32x,r32y,r32z;
  double sume,avge,f;

  FILE* fef;

  nbidx = i0mat(top->n_atoms*pars->n_mol,top->n_atoms*pars->n_mol);
  for(i=0;i< pars->n_qme+pars->n_qmf;i++) { mme0[i]=0.0; mme[i]=0.0;}

/***** create top->nbidx[i][j] to associate the appropriate LJ parameters with each pair i,j *****/

  for(i=0;i<top->n_atoms*pars->n_mol-1;i++) {
    for(j=i+1;j<top->n_atoms*pars->n_mol;j++) {

      ii=i; jj=j;
      while(ii>=top->n_atoms) ii -= top->n_atoms;
      while(jj>=top->n_atoms) jj -= top->n_atoms;

      for(k=0;k<top->n_nb;k++) {

        if((strcmp(top->nbat1[k],top->type[ii])==0 && strcmp(top->nbat2[k],top->type[jj])==0) ||
           (strcmp(top->nbat2[k],top->type[ii])==0 && strcmp(top->nbat1[k],top->type[jj])==0))
             nbidx[i][j]=k;
      }
    }
  }

/***** set up neighbour list *****/

  top->n_neighbors = 0;

  for(i=0;i<pars->n_mol*top->n_atoms-1;i++) {
    for(j=i+1;j<pars->n_mol*top->n_atoms;j++) {

      molnr1=(int)i/top->n_atoms;
      molnr2=(int)j/top->n_atoms;

      bbond=0;
      bangle=0;
      bpair=0;
      
      ii = i; while(ii>=top->n_atoms) ii -= top->n_atoms;
      jj = j; while(jj>=top->n_atoms) jj -= top->n_atoms;

      for(k=0;k<top->n_bonds;k++) {
        if(((top->ndex[ii]==top->bondi[k] && top->ndex[jj]==top->bondj[k]) || 
            (top->ndex[ii]==top->bondj[k] && top->ndex[jj]==top->bondi[k]) ) && molnr1==molnr2) {
           bbond=1;
        }
      }

      for(k=0;k<top->n_angles;k++) {
        if(((top->ndex[ii]==top->anglei[k] && top->ndex[jj]==top->anglek[k] ) || 
            (top->ndex[ii]==top->anglek[k] && top->ndex[jj]==top->anglei[k] ) ) && molnr1==molnr2) {
           bangle=1;
        }
      }

      for(k=0;k<top->n_pairs;k++) {
        if(((top->ndex[ii]==top->pairi[k] && top->ndex[jj]==top->pairj[k]) || 
            (top->ndex[ii]==top->pairj[k] && top->ndex[jj]==top->pairi[k]) ) && molnr1==molnr2) {
           bpair=1;
        }
      }

      r12x = coords[i][XX]-coords[j][XX];
      r12y = coords[i][YY]-coords[j][YY];
      r12z = coords[i][ZZ]-coords[j][ZZ];
      r12sq = r12x*r12x+r12y*r12y+r12z*r12z;
      br12 = sqrt(r12sq);

      if(bbond==0 && bangle==0) {

        if(bpair==1) {
          top->nb_el_fac[top->n_neighbors]=0.8333;
          top->nb_lj_fac[top->n_neighbors]=0.5;
        } else {
          top->nb_el_fac[top->n_neighbors]=1.0;
          top->nb_lj_fac[top->n_neighbors]=1.0;
        }

        top->neighbor1[top->n_neighbors]=i;
        top->neighbor2[top->n_neighbors]=j;
        top->nbidx[top->n_neighbors]=nbidx[i][j];

        top->n_neighbors++;
      }
    }
  }

/***** calculate energies/forces *****/
/***** loop over all conformations *****/

  for(n_conf=0;n_conf<pars->n_conf;n_conf++) {

    evdw=0.0;
    ecoul=0.0;
    ebond=0.0;
    eimp=0.0;
    edih=0.0;
    eang=0.0;

/***** Coulomb *****/
/***** loop over all neighbours *****/

    for(nn=0;nn<top->n_neighbors;nn++) {

      i1 = top->neighbor1[nn]+n_conf*top->n_atoms*pars->n_mol;
      i2 = top->neighbor2[nn]+n_conf*top->n_atoms*pars->n_mol;
      rnn=top->nbidx[nn];

      r12x = coords[i1][XX]-coords[i2][XX];
      r12y = coords[i1][YY]-coords[i2][YY];
      r12z = coords[i1][ZZ]-coords[i2][ZZ];
      r12sq = r12x*r12x+r12y*r12y+r12z*r12z;
      br12 = sqrt(r12sq);

      ii = top->neighbor1[nn]; while(ii>=top->n_atoms) ii -= top->n_atoms;
      jj = top->neighbor2[nn]; while(jj>=top->n_atoms) jj -= top->n_atoms;

      cou_pot  =  top->nb_el_fac[nn]*oo4pe0*(top->charge[ii])*(top->charge[jj])/br12;
      dcou_pot = -cou_pot/br12;
      ecoul += cou_pot;

      dcou_pot /= br12;
      f1x = -dcou_pot*r12x;
      f1y = -dcou_pot*r12y;
      f1z = -dcou_pot*r12z;
      f2x =  dcou_pot*r12x;
      f2y =  dcou_pot*r12y;
      f2z =  dcou_pot*r12z;

      mmeidx1 = pars->n_qme+3*i1;
      mmeidx2 = pars->n_qme+3*i2;

      if(pars->n_qme>0) mme0[n_conf] += cou_pot;
      if(pars->n_qmf>0) {
        mme0[mmeidx1]   += f1x;
        mme0[mmeidx1+1] += f1y;
        mme0[mmeidx1+2] += f1z;
        mme0[mmeidx2]   += f2x;
        mme0[mmeidx2+1] += f2y;
        mme0[mmeidx2+2] += f2z;
      }
    }

/***** VdW *****/

    vdwcosq = (pars->vdw_cutoff)*(pars->vdw_cutoff);

    for(nn=0;nn<top->n_neighbors;nn++) {

      i1 = top->neighbor1[nn]+n_conf*top->n_atoms*pars->n_mol;
      i2 = top->neighbor2[nn]+n_conf*top->n_atoms*pars->n_mol;
      rnn=top->nbidx[nn];
      sig=top->sigma[rnn];
      eps=top->epsilon[rnn];
      if(strstr(pars->opt,dovdw)!=NULL && top->sgrpidx[rnn]!=0 ) sig = genome->codon[limits->codonidx[top->sgrpidx[rnn]]];
      if(strstr(pars->opt,dovdw)!=NULL && top->egrpidx[rnn]!=0 ) eps = genome->codon[limits->codonidx[top->egrpidx[rnn]]];

      sig2=sig*sig;
      sig6=sig2*sig2*sig2;
      sig12=sig6*sig6;

      r12x = coords[i1][XX]-coords[i2][XX];
      r12y = coords[i1][YY]-coords[i2][YY];
      r12z = coords[i1][ZZ]-coords[i2][ZZ];
      r12sq = r12x*r12x+r12y*r12y+r12z*r12z;

      if(r12sq<=vdwcosq) {

        br12 = sqrt(r12sq);
        oor2=1.0/r12sq;
        oor6=oor2*oor2*oor2;
        oor12=oor6*oor6;

        ii = top->neighbor1[nn]; while(ii>=top->n_atoms) ii -= top->n_atoms;
        jj = top->neighbor2[nn]; while(jj>=top->n_atoms) jj -= top->n_atoms;

        vdw_pot  =  top->nb_lj_fac[nn]*4.0*eps*(      sig12*oor12     -    sig6*oor6      );
        dvdw_pot =  top->nb_lj_fac[nn]*4.0*eps*( 12.0*sig12*oor12/br12-6.0*sig6*oor6/br12 );

        evdw     += vdw_pot;
        dvdw_pot /= br12;

        f1x =  (dvdw_pot)*r12x;
        f1y =  (dvdw_pot)*r12y;
        f1z =  (dvdw_pot)*r12z;
        f2x = -(dvdw_pot)*r12x;
        f2y = -(dvdw_pot)*r12y;
        f2z = -(dvdw_pot)*r12z;

        mmeidx1 = pars->n_qme+3*i1;
        mmeidx2 = pars->n_qme+3*i2;

        if((top->sgrpidx[rnn]==0 && top->egrpidx[rnn]==0) || strstr(pars->opt,dovdw) == NULL ) {
          if(pars->n_qme>0) mme0[n_conf]    += vdw_pot;
          if(pars->n_qmf>0) {
            mme0[mmeidx1]   += f1x;
            mme0[mmeidx1+1] += f1y;
            mme0[mmeidx1+2] += f1z;
            mme0[mmeidx2]   += f2x;
            mme0[mmeidx2+1] += f2y;
            mme0[mmeidx2+2] += f2z;
          }
        }
      }
    }

/***** bonds *****/

    for(n_bond=0;n_bond<top->n_bonds;n_bond++) {
      for(n_mol=0;n_mol<pars->n_mol;n_mol++) {

        i1 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->bondi[n_bond] - 1;
        i2 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->bondj[n_bond] - 1;
        r0 = top->bondr[n_bond];
        kb = top->bondk[n_bond];
        if(strstr(pars->opt,dobon)!=NULL && top->bgrpidx[n_bond]!=0) kb = genome->codon[limits->codonidx[top->bgrpidx[n_bond]]];

        r12x = coords[i1][XX]-coords[i2][XX];
        r12y = coords[i1][YY]-coords[i2][YY];
        r12z = coords[i1][ZZ]-coords[i2][ZZ];
        br12 = sqrt(r12x*r12x+r12y*r12y+r12z*r12z);

        debonds = -kb*(br12-r0);
        bon_pot = 0.5*kb*(br12-r0)*(br12-r0);
        ebond+=bon_pot;

        debonds /= br12;

        f1x =  debonds*r12x;
        f1y =  debonds*r12y;
        f1z =  debonds*r12z;
        f2x = -debonds*r12x;
        f2y = -debonds*r12y;
        f2z = -debonds*r12z;

        mmeidx1 = pars->n_qme+3*i1;
        mmeidx2 = pars->n_qme+3*i2;

        if(top->isconstrained[n_bond]==1) {
          if(top->bgrpidx[n_bond]!=0) {
            printf("bond nr %d (%d-%d) is, both, constrained (X-H) and included in optimization\n", n_bond, i1,i2);
            exit(1);
          }
        } else {
          if(top->bgrpidx[n_bond]==0 || strstr(pars->opt,dobon)==NULL) {
            if(pars->n_qme>0) mme0[n_conf] += bon_pot;
            if(pars->n_qmf>0) {
              mme0[mmeidx1]   += f1x;
              mme0[mmeidx1+1] += f1y;
              mme0[mmeidx1+2] += f1z;
              mme0[mmeidx2]   += f2x;
              mme0[mmeidx2+1] += f2y;
              mme0[mmeidx2+2] += f2z;
            }
          }
        }
      }
    }

/***** angles *****/

    for(n_ang=0;n_ang<top->n_angles;n_ang++) {
      for(n_mol=0;n_mol<pars->n_mol;n_mol++) {

        i1 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->anglei[n_ang] - 1;
        i2 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->anglej[n_ang] - 1;
        i3 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->anglek[n_ang] - 1;
        a0 = top->angleeq[n_ang];
        ka = top->anglefc[n_ang];
        if(top->cgrpidx[n_ang]!=0 && strstr(pars->opt,doang)!=NULL) ka = genome->codon[limits->codonidx[top->cgrpidx[n_ang]]];
        if(top->agrpidx[n_ang]!=0 && strstr(pars->opt,doang)!=NULL) a0 = genome->codon[limits->codonidx[top->agrpidx[n_ang]]];
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
        ang_pot = 0.5*ka*diff*diff;
        eang += ang_pot;

        sin_theta = sqrt(1.0 - cos_theta*cos_theta);
        diff *= -ka/sin_theta;

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

        if((top->cgrpidx[n_ang]==0 && top->agrpidx[n_ang]==0) || strstr(pars->opt,doang)==NULL) {
          if(pars->n_qme>0) mme0[n_conf] += ang_pot;
          if(pars->n_qmf>0) {
            mme0[mmeidx1]   += f1x;
            mme0[mmeidx1+1] += f1y;
            mme0[mmeidx1+2] += f1z;
            mme0[mmeidx2]   += f2x;
            mme0[mmeidx2+1] += f2y;
            mme0[mmeidx2+2] += f2z;
            mme0[mmeidx3]   += f3x;
            mme0[mmeidx3+1] += f3y;
            mme0[mmeidx3+2] += f3z;
          }
        }
      }
    }

/***** dihedrals *****/

    for(n_dih=0;n_dih<top->n_dihedrals;n_dih++) {
      for(n_mol=0;n_mol<pars->n_mol;n_mol++) {

        i1 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->dihi[n_dih] - 1;
        i2 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->dihj[n_dih] - 1;
        i3 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->dihk[n_dih] - 1;
        i4 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->dihl[n_dih] - 1;

        dp=top->dihphase[n_dih]/180.0*M_PI;
        fc=top->dihfc[n_dih];
        dn=top->dihpn[n_dih];
        if(top->dgrpidx[n_dih]!=0 && strstr(pars->opt,dodih)!=NULL) fc = genome->codon[limits->codonidx[top->dgrpidx[n_dih]]];

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
        if(rkjdotmcrossn>=0.0) sign= 1.0;
        if(rkjdotmcrossn< 0.0) sign=-1.0;

        acosarg=nRx*nSx+nRy*nSy+nRz*nSz;
        if(acosarg>1.0) acosarg=1.0;
        if(acosarg<-1.0) acosarg=-1.0;
        phi = sign*acos(acosarg);

        dih_pot = fc*(1.0+cos(dn*phi - dp));
        dvdphi = -dn*fc*sin(dn*phi - dp);
        edih += dih_pot;

        bm2 = mx*mx+my*my+mz*mz;
        bn2 = nx*nx+ny*ny+nz*nz;

        fix = -dvdphi * brkj/bm2*mx;
        fiy = -dvdphi * brkj/bm2*my;
        fiz = -dvdphi * brkj/bm2*mz;

        flx =  dvdphi * brkj/bn2*nx;
        fly =  dvdphi * brkj/bn2*ny;
        flz =  dvdphi * brkj/bn2*nz;

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

        if(top->dgrpidx[n_dih]==0 || strstr(pars->opt,dodih)==NULL) {
          if(pars->n_qme>0) mme0[n_conf]    += dih_pot;
          if(pars->n_qmf>0) {
            mme0[mmeidx1]   += fix;
            mme0[mmeidx1+1] -= fiy;
            mme0[mmeidx1+2] += fiz;
            mme0[mmeidx2]   += fjx;
            mme0[mmeidx2+1] -= fjy;
            mme0[mmeidx2+2] += fjz;
            mme0[mmeidx3]   += fkx;
            mme0[mmeidx3+1] -= fky;
            mme0[mmeidx3+2] += fkz;
            mme0[mmeidx4]   += flx;
            mme0[mmeidx4+1] -= fly;
            mme0[mmeidx4+2] += flz;
          }
        }
      }
    }

/***** impropers *****/

    for(n_imp=0;n_imp<top->n_impropers;n_imp++) {
      for(n_mol=0;n_mol<pars->n_mol;n_mol++) {

        i1 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->impi[n_imp] - 1;
        i2 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->impj[n_imp] - 1;
        i3 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->impk[n_imp] - 1;
        i4 = pars->n_mol*top->n_atoms*n_conf + n_mol*top->n_atoms + top->impl[n_imp] - 1;

        dp=top->impphase[n_imp]/180.0*M_PI;
        fc=top->impfc[n_imp];
        dn=top->imppn[n_imp];
        if(top->igrpidx[n_imp]!=0 && strstr(pars->opt,doimp)!=NULL) fc = genome->codon[limits->codonidx[top->igrpidx[n_imp]]];

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
        if(rkjdotmcrossn>=0.0) sign= 1.0;
        if(rkjdotmcrossn< 0.0) sign=-1.0;

        acosarg=nRx*nSx+nRy*nSy+nRz*nSz;
        if(acosarg>1.0) acosarg=1.0;
        if(acosarg<-1.0) acosarg=-1.0;
        phi = sign*acos(acosarg);

        imp_pot = fc*(1.0+cos(dn*phi - dp));
        dvdphi = -dn*fc*sin(dn*phi - dp);
        eimp += imp_pot;

        bm2 = mx*mx+my*my+mz*mz;
        bn2 = nx*nx+ny*ny+nz*nz;

        fix = -dvdphi * brkj/bm2*mx;
        fiy = -dvdphi * brkj/bm2*my;
        fiz = -dvdphi * brkj/bm2*mz;

        flx =  dvdphi * brkj/bn2*nx;
        fly =  dvdphi * brkj/bn2*ny;
        flz =  dvdphi * brkj/bn2*nz;

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

        if(top->igrpidx[n_imp]==0 || strstr(pars->opt,doimp)==NULL) {
          if(pars->n_qme>0) mme0[n_conf]    += imp_pot;
          if(pars->n_qmf>0) {
            mme0[mmeidx1]   += fix;
            mme0[mmeidx1+1] -= fiy;
            mme0[mmeidx1+2] += fiz;
            mme0[mmeidx2]   += fjx;
            mme0[mmeidx2+1] -= fjy;
            mme0[mmeidx2+2] += fjz;
            mme0[mmeidx3]   += fkx;
            mme0[mmeidx3+1] -= fky;
            mme0[mmeidx3+2] += fkz;
            mme0[mmeidx4]   += flx;
            mme0[mmeidx4+1] -= fly;
            mme0[mmeidx4+2] += flz;
          }
        }
      }
    }
  }

/***** calculate the energies/forces involving variable ff terms, and add to mme0 *****/

  nread = fitness(pars, genome, top, limits, coords, qme, w_qme, mme, mme0);

/***** determine average and normalize energies *****/

  mmestdev = 0.0;

  if(pars->n_qme>0) {
    sume = 0.0;
    for(i=0; i<pars->n_qme; i++) sume += mme[i];
    avge=sume/(double)pars->n_qme;
    for(i=0; i<pars->n_qme; i++) mmestdev += (mme[i]-avge)*(mme[i]-avge);
    for(i=0; i<pars->n_qme; i++) mme[i] = (mme[i]-avge)/pars->stdeve;
  }

  mmestdev = sqrt(mmestdev/(double)pars->n_qme);
  printf("# average total energy and stdev from input FF: <E> = %lf +/- %lf kJ/mol\n", avge, mmestdev);

/***** normalize forces (no need to shift by avg as avg is usually 0) *****/

  for(i=pars->n_qme; i<pars->n_qme+pars->n_qmf; i++) mme[i] = mme[i]/pars->stdevf;

/***** open output file *****/

  if((fef=fopen(filename,"w"))==NULL) error_msg("cannot open ef-init file");

/***** write energies and forces *****/

  for(i=0;i<pars->n_qme;i++) {
    fprintf(fef,"%8d %14.6le %14.6le %14.6le %14.6le %14.6le\n", i, qme[i], mme0[i], mme[i], qme[i]*pars->stdeve, mme[i]*pars->stdeve);
  }
  for(i=pars->n_qme; i<pars->n_qme+pars->n_qmf;i++) {
    fprintf(fef,"%8d %14.6le %14.6le %14.6le %14.6le %14.6le\n", i, qme[i], mme0[i], mme[i], qme[i]*pars->stdevf, mme[i]*pars->stdevf);
  }

/***** calculate fitness *****/

  f = 0.0;
  sum = 0.0;
  for(i=0;i<pars->n_qme;i++) {
    f1 = (qme[i]-mme[i]);
    f += f1*f1*w_qme[i];
    sum += w_qme[i];
  }
  for(i=pars->n_qme; i<pars->n_qme+pars->n_qmf;i++) {
    f1 = (qme[i]-mme[i]);
    f += f1*f1*w_qme[i];;
    sum += w_qme[i];
  }

  genome->fitness = sqrt(f/sum);

  printf("# Fitness at t=0: %lf\n", genome->fitness);
  printf("# saved QM vs GMX energies in %s\n", filename);

  for(i=0;i<top->n_atoms*pars->n_mol;i++) free(nbidx[i]);

  fclose(fef);

  free(nbidx);

  return(0);
}
