// Copyright (c) 2018-2020 Osamu Hirose
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
#include <iostream>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<string.h>
#include<math.h>
#include<unistd.h>
#include<ctype.h>
#include<time.h>
#include<sys/time.h>
#include"../base/util.h"
#include"../base/misc.h"
#include"../base/kdtree.h"
#include"../base/kernel.h"
#include"../base/sampling.h"
#include"bcpd.h"
#include"info.h"
#include"norm.h"

#define SQ(x) ((x)*(x))

void init_genrand64(unsigned long long s);
enum transpose {ASIS=0,TRANSPOSE=1};

void save_variable(const char *prefix, const char *suffix,const double *var, int D, int J, char *fmt, int trans){
  int d,j; char fn[256]; double **buf;
  strcpy(fn,prefix); strcat(fn,suffix);
  if(trans==TRANSPOSE){
    buf=calloc2d(J,D);
    for(j=0;j<J;j++)for(d=0;d<D;d++) buf[j][d]=var[d+D*j];
    write2d(fn,(const double **)buf,J,D,fmt,"NA"); free2d(buf,J,D);
  }
  else {
    buf=calloc2d(D,J);
    for(j=0;j<J;j++)for(d=0;d<D;d++) buf[d][j]=var[d+D*j];
    write2d(fn,(const double **)buf,D,J,fmt,"NA"); free2d(buf,D,J);
  }

  return;
}

void save_corresp(
    const char   *prefix,
    const double *X,
    const double *y,
    const double *a,
    const double *sgm,
    const double  s,
    const double  r,
    pwsz          sz,
    pwpm          pm
  ){
  int i,m,n,D,M,N; int *T,*l,*bi; double *bd; double *p,c,val; char fnP[256],fnc[256],fne[256];
  int S[MAXTREEDEPTH]; int top,ct; double omg,dlt,vol,rad; int si=sizeof(int),sd=sizeof(double);
  FILE *fpP=NULL,*fpc=NULL,*fpe=NULL; int db=pm.opt&PW_OPT_DBIAS; double max,min; int mmax;

  D=sz.D; M=sz.M; N=sz.N; omg=pm.omg; dlt=pm.dlt; rad=dlt*r;
  strcpy(fnP,prefix);strcat(fnP,"P.txt");if(pm.opt&PW_OPT_SAVEP){fpP=fopen(fnP,"w");fprintf(fpP,"[n]\t[m]\t[probability]\n");}
  strcpy(fne,prefix);strcat(fne,"e.txt");if(pm.opt&PW_OPT_SAVEE){fpe=fopen(fne,"w");fprintf(fpe,"[n]\t[m]\t[probability]\n");}
  strcpy(fnc,prefix);strcat(fnc,"c.txt");if(pm.opt&PW_OPT_SAVEC){fpc=fopen(fnc,"w");fprintf(fpc,"[n]\t[1/0]\n");}

  T= (int *)calloc(3*M+1,si); bi= (int *)calloc(6*M,si); bd= (double *)calloc(2*M,sd); p= (double *)calloc(M,sd); l=(int *)calloc(M,si);
  kdtree(T,bi,bd,y,D,M); vol=volume(X,D,N); c=(pow(2.0*M_PI*SQ(r),0.5*D)*omg)/(vol*(1-omg));
  for(n=0;n<N;n++){
    /* compute P, c, e */
    val=c;top=ct=0;do{eballsearch_next(&m,S,&top,X+D*n,rad,y,T,D,M);if(m>=0)l[ct++]=m;}while(top);
    if(!ct){nnsearch(&m,&min,X+D*n,y,T,D,M);l[ct++]=m;}
    for(i=0;i<ct;i++){m=l[i];p[i]=a[m]*gauss(y+D*m,X+D*n,D,rad)*(db?exp(-0.5*D*SQ(sgm[m]*s/r)):1.0);val+=p[i];}
    for(i=0;i<ct;i++){m=l[i];p[i]/=val;}
    max=c/val;mmax=0;for(i=0;i<ct;i++)if(p[i]>max){max=p[i];mmax=l[i]+1;}
    /* print P, c, e */
    if(fpP){for(i=0;i<ct;i++)if(p[i]>1.0f/M){m=l[i];fprintf(fpP,"%d\t%d\t%lf\n",n+1,m+1,p[i]);}}
    if(fpe){fprintf(fpe,"%d\t%d\t%lf\n",n+1,mmax?mmax:l[0],mmax?max:p[0]);}
    if(fpc){fprintf(fpc,"%d\t%d\n",n+1,mmax?1:0);}
  }
  if(fpP){fclose(fpP);} free(l); free(bd);
  if(fpe){fclose(fpc);} free(p); free(bi);
  if(fpc){fclose(fpe);} free(T);
  return;
}

int save_optpath(const char *file, const double *sy, const double *X, pwsz sz, int lp){
  int N=sz.N,M=sz.M,D=sz.D; int si=sizeof(int),sd=sizeof(double);
  FILE *fp=fopen(file,"wb");
  if(!fp){printf("Can't open: %s\n",file);exit(EXIT_FAILURE);}
  fwrite(&N, si,   1,  fp);
  fwrite(&D, si,   1,  fp);
  fwrite(&M, si,   1,  fp);
  fwrite(&lp,si,   1,  fp);
  fwrite(sy, sd,lp*D*M,fp);
  fwrite(X,  sd,   D*N,fp);
  fclose(fp);
  return 0;
}

void scan_beta(double *bet, const char *arg){
  int i=0,ct=0;
  while(arg[i++]!='\0')if(arg[i]==','){ct++;} assert(ct<=1);
  switch(ct){
    case 0: sscanf(arg,"%lf",    bet);        break;
    case 1: sscanf(arg,"%lf,%lf",bet,bet+1);  break;
  }
}

void scan_dwpm(int *dwn, double *dwr, const char *arg){
  char c; int n,m; double r;
  m=sscanf(arg,"%c,%d,%lf",&c,&n,&r);
  if(m!=3) goto err01;
  if(n<=0) goto err03;
  if(r< 0) goto err04;
  if(isupper(c)){r*=-1.0f;} c=tolower(c);
  if(c!='x'&&c!='y'&&c!='b') goto err02;
  if(r<0&&-r<1e-2) r=-1e-2;
  switch(c){
    case 'x': dwn[TARGET]=n; dwr[TARGET]=r; break;
    case 'y': dwn[SOURCE]=n; dwr[SOURCE]=r; break;
    case 'b': dwn[TARGET]=n; dwr[TARGET]=r;
              dwn[SOURCE]=n; dwr[SOURCE]=r; break;
  }
  return;
  err01: printf("ERROR: The argument of '-D' must be 'char,int,real'. \n");         exit(EXIT_FAILURE);
  err02: printf("ERROR: The 1st argument of '-D' must be one of [x,y,b,X,Y,B]. \n");exit(EXIT_FAILURE);
  err03: printf("ERROR: The 2nd argument of '-D' must be positive.     \n");        exit(EXIT_FAILURE);
  err04: printf("ERROR: The 3rd argument of '-D' must be positive or 0.\n");        exit(EXIT_FAILURE);
}

void check_prms(const pwpm pm, const pwsz sz){
  int M=sz.M,N=sz.N,M0=pm.dwn[SOURCE],N0=pm.dwn[TARGET]; M=M0?M0:M; N=N0?N0:N;
  if(pm.nlp<=0){printf("ERROR: -n: Argument must be a positive integer. Abort.\n"); exit(EXIT_FAILURE);}
  if(pm.llp<=0){printf("ERROR: -N: Argument must be a positive integer. Abort.\n"); exit(EXIT_FAILURE);}
  if(pm.omg< 0){printf("ERROR: -w: Argument must be in range [0,1]. Abort.\n");     exit(EXIT_FAILURE);}
  if(pm.omg>=1){printf("ERROR: -w: Argument must be in range [0,1]. Abort.\n");     exit(EXIT_FAILURE);}
  if(pm.lmd<=0){printf("ERROR: -l: Argument must be positive. Abort.\n");           exit(EXIT_FAILURE);}
  if(pm.kpa<=0){printf("ERROR: -k: Argument must be positive. Abort.\n");           exit(EXIT_FAILURE);}
  if(pm.dlt<=0){printf("ERROR: -d: Argument must be positive. Abort.\n");           exit(EXIT_FAILURE);}
  if(pm.lim<=0){printf("ERROR: -e: Argument must be positive. Abort.\n");           exit(EXIT_FAILURE);}
  if(pm.btn<=0){printf("ERROR: -f: Argument must be positive. Abort.\n");           exit(EXIT_FAILURE);}
  if(pm.cnv<=0){printf("ERROR: -c: Argument must be positive. Abort.\n");           exit(EXIT_FAILURE);}
  if(pm.rns< 0){printf("ERROR: -r: Argument must be positive. Abort.\n");           exit(EXIT_FAILURE);}
  if(pm.K<0)   {printf("ERROR: -K: Argument must be a positive integer. Abort.\n"); exit(EXIT_FAILURE);}
  if(pm.K>M)   {printf("ERROR: -K: Argument must be less than M. Abort.\n");        exit(EXIT_FAILURE);}
  if(pm.J<0)   {printf("ERROR: -J: Argument must be a positive integer. Abort.\n"); exit(EXIT_FAILURE);}
  if(pm.J>M+N) {printf("ERROR: -J: Argument must be less than M+N. Abort.\n");      exit(EXIT_FAILURE);}
  if(pm.L<0)   {printf("ERROR: -L: Argument must be a positive integer. Abort.\n"); exit(EXIT_FAILURE);}
  if(pm.L>M)   {printf("ERROR: -L: Argument must be less than M. Abort.\n");        exit(EXIT_FAILURE);}
  if(pm.G>5)   {printf("ERROR: -G: Argument must be one of {0,1,2,3,4}. Abort.\n"); exit(EXIT_FAILURE);}
  if(pm.G<0)   {printf("ERROR: -G: Argument must be one of {0,1,2,3,4}. Abort.\n"); exit(EXIT_FAILURE);}
  if(pm.bet[0]<=0){printf("ERROR: -b: Argument must be positive. Abort.\n");        exit(EXIT_FAILURE);}
  if(pm.bet[1]<=0){printf("ERROR: -b: Argument must be positive. Abort.\n");        exit(EXIT_FAILURE);}
  if(!strchr("exyn",pm.nrm)){printf("\n  ERROR: -u: Argument must be one of 'e', 'x', 'y' and 'n'. Abort.\n\n");exit(EXIT_FAILURE);}
}

void pw_getopt(pwpm *pm, int argc, char **argv){ int opt;
  strcpy(pm->fn[TARGET],"X.txt");   pm->omg=0.0; pm->cnv=1e-4; pm->K=0; pm->opt=0.0; pm->btn=0.20; pm->bet[0]=2.0; pm->L=0;
  strcpy(pm->fn[SOURCE],"Y.txt");   pm->lmd=2.0; pm->nlp= 500; pm->J=0; pm->dlt=7.0; pm->lim=0.15; pm->bet[1]=2.0;
  strcpy(pm->fn[OUTPUT],"output_"); pm->rns=0;   pm->llp=  30; pm->G=0; pm->gma=1.0; pm->kpa=ZERO; pm->nrm='e';
  pm->dwn[SOURCE]=0; pm->dwr[SOURCE]=0.0f;
  pm->dwn[TARGET]=0; pm->dwr[TARGET]=0.0f;
  while((opt=getopt(argc,argv,"D:u:r:w:l:b:k:g:d:e:c:n:N:G:J:K:L:o:x:y:f:s:hpqvaAtW1"))!=-1){
    switch(opt){
      case 'D': scan_dwpm(pm->dwn,pm->dwr,optarg);    break;
      case 'b': scan_beta(pm->bet,optarg);            break;
      case 'w': pm->omg  = atof(optarg);              break;
      case 'l': pm->lmd  = atof(optarg);              break;
      case 'k': pm->kpa  = atof(optarg);              break;
      case 'g': pm->gma  = atof(optarg);              break;
      case 'd': pm->dlt  = atof(optarg);              break;
      case 'e': pm->lim  = atof(optarg);              break;
      case 'f': pm->btn  = atof(optarg);              break;
      case 'c': pm->cnv  = atof(optarg);              break;
      case 'n': pm->nlp  = atoi(optarg);              break;
      case 'N': pm->llp  = atoi(optarg);              break;
      case 'K': pm->K    = atoi(optarg);              break;
      case 'J': pm->J    = atoi(optarg);              break;
      case 'L': pm->L    = atoi(optarg);              break;
      case 'G': pm->G    = atoi(optarg);              break;
      case 'r': pm->rns  = atoi(optarg);              break;
      case 'u': pm->nrm  = *optarg;                   break;
      case 'h': pm->opt |= PW_OPT_HISTO;              break;
      case 'a': pm->opt |= PW_OPT_DBIAS;              break;
      case 'p': pm->opt |= PW_OPT_LOCAL;              break;
      case 'q': pm->opt |= PW_OPT_QUIET;              break;
      case 'A': pm->opt |= PW_OPT_ACCEL;              break;
      case 'W': pm->opt |= PW_OPT_NWARN;              break;
      case '1': pm->opt |= PW_OPT_1NN;                break;
      case 'o': strcpy(pm->fn[OUTPUT],optarg);        break;
      case 'x': strcpy(pm->fn[TARGET],optarg);        break;
      case 'y': strcpy(pm->fn[SOURCE],optarg);        break;
      case 'v': printUsage(); exit(EXIT_SUCCESS);     break;
      case 's':
        if(strchr(optarg,'A')) pm->opt |= PW_OPT_SAVE;
        if(strchr(optarg,'x')) pm->opt |= PW_OPT_SAVEX;
        if(strchr(optarg,'y')) pm->opt |= PW_OPT_SAVEY;
        if(strchr(optarg,'u')) pm->opt |= PW_OPT_SAVEU;
        if(strchr(optarg,'v')) pm->opt |= PW_OPT_SAVEV;
        if(strchr(optarg,'a')) pm->opt |= PW_OPT_SAVEA;
        if(strchr(optarg,'c')) pm->opt |= PW_OPT_SAVEC;
        if(strchr(optarg,'e')) pm->opt |= PW_OPT_SAVEE;
        if(strchr(optarg,'S')) pm->opt |= PW_OPT_SAVES;
        if(strchr(optarg,'P')) pm->opt |= PW_OPT_SAVEP;
        if(strchr(optarg,'T')) pm->opt |= PW_OPT_SAVET;
        if(strchr(optarg,'X')) pm->opt |= PW_OPT_PATHX;
        if(strchr(optarg,'Y')) pm->opt |= PW_OPT_PATHY;
        if(strchr(optarg,'t')) pm->opt |= PW_OPT_PFLOG;
        if(strchr(optarg,'0')) pm->opt |= PW_OPT_VTIME;
      break;
    }
  }
  /* acceleration with default parameters */
  if(pm->opt&PW_OPT_ACCEL) {pm->J=300;pm->K=70;pm->opt|=PW_OPT_LOCAL;}
  /* case: save all */
  if(pm->opt&PW_OPT_SAVE)
    pm->opt|=PW_OPT_SAVEX|PW_OPT_SAVEU|PW_OPT_SAVEC|PW_OPT_SAVEP|PW_OPT_SAVEA|PW_OPT_PFLOG|
             PW_OPT_SAVEY|PW_OPT_SAVEV|PW_OPT_SAVEE|PW_OPT_SAVET|PW_OPT_SAVES|PW_OPT_PATHY;
  /* always save y & info */
  pm->opt |= PW_OPT_SAVEY|PW_OPT_INFO;
  /* for numerical stability */
  pm->omg=pm->omg==0?1e-250:pm->omg;
  /* llp is always less than or equal to nlp */
  if(pm->llp>pm->nlp) pm->llp=pm->nlp;

  return;
}

void memsize(int *dsz, int *isz, pwsz sz, pwpm pm){
  int M=sz.M,N=sz.N,J=sz.J,K=sz.K,D=sz.D; int T=pm.opt&PW_OPT_LOCAL; int L=M>N?M:N,mtd=MAXTREEDEPTH;
  *isz =D;              *dsz =4*M+2*N+D*(5*M+N+13*D+3);  /* common          */
  *isz+=K?M:0;          *dsz+=K?K*(2*M+3*K+D+12):(3*M*M);/* low-rank        */
  *isz+=J?(M+N):0;      *dsz+=J?(D*(M+N)+J+J*J):0;       /* nystrom         */
  *isz+=T?L*6:0;        *dsz+=T?L*2:0;                   /* kdtree (build)  */
  *isz+=T?L*(2+mtd):0;                                   /* kdtree (search) */
  *isz+=T?2*(3*L+1):0;                                   /* kdtree (tree)   */
}

void print_bbox(const double *X, int D, int N){
  int d,n; double max,min; char ch[3]={'x','y','z'};
  for(d=0;d<D;d++){
    max=X[d];for(n=1;n<N;n++) max=fmax(max,X[d+D*n]);
    min=X[d];for(n=1;n<N;n++) min=fmin(min,X[d+D*n]);
    fprintf(stderr,"%c=[%.2f,%.2f]%s",ch[d],min,max,d==D-1?"\n":", ");
  }
}

void print_norm(const double *X, const double *Y, int D, int N, int M, int sw, char type){
    int t=0; char name[4][64]={"for each","using X","using Y","skipped"};
    switch(type){
      case 'e': t=0; break;
      case 'x': t=1; break;
      case 'y': t=2; break;
      case 'n': t=3; break;
    }
    if(sw){
      fprintf(stderr,"  Normalization: [%s]\n",name[t]);
      fprintf(stderr,"    Bounding boxes that cover point sets:\n");
    }
    fprintf(stderr,"    %s:\n",sw?"Before":"After");
    fprintf(stderr,"      Target: "); print_bbox(X,D,N);
    fprintf(stderr,"      Source: "); print_bbox(Y,D,M);
    if(!sw) fprintf(stderr,"\n");
}

double tvcalc(const struct timeval *end, const struct timeval *beg){
  return (end->tv_sec-beg->tv_sec)+(end->tv_usec-beg->tv_usec)/1e6;
}

void fprint_comptime(FILE *fp, const struct timeval *tv, double *tt, int nx, int ny){
  if(fp==stderr) fprintf(fp,"  Computing Time:\n");
  #ifdef MINGW32
  fprintf(fp,"    VB Optimization: %.2lf sec\n",            tvcalc(tv+3,tv+2));
  if(nx||ny) fprintf(fp,"    Downsampling:    %.3lf sec\n", tvcalc(tv+2,tv+1));
  if(ny)     fprintf(fp,"    Interpolation:   %.3lf sec\n", tvcalc(tv+4,tv+3));
  #else
  fprintf(fp,"    VB Optimization: %.2lf sec (real) / %.2lf sec (cpu)\n",           tvcalc(tv+3,tv+2),(tt[3]-tt[2])/CLOCKS_PER_SEC);
  if(nx||ny) fprintf(fp,"    Downsampling:    %.3lf sec (real) / %.3lf sec (cpu)\n",tvcalc(tv+2,tv+1),(tt[2]-tt[1])/CLOCKS_PER_SEC);
  if(ny)     fprintf(fp,"    Interpolation:   %.3lf sec (real) / %.3lf sec (cpu)\n",tvcalc(tv+4,tv+3),(tt[4]-tt[3])/CLOCKS_PER_SEC);
  #endif
  fprintf(fp,"    File reading:    %.3lf sec\n",tvcalc(tv+1,tv+0));
  fprintf(fp,"    File writing:    %.3lf sec\n",tvcalc(tv+5,tv+4));
  if(fp==stderr) fprintf(fp,"\n");
}

void fprint_comptime2(FILE *fp, const struct timeval *tv, double *tt){
  fprintf(fp,"%lf\t%lf\n",tvcalc(tv+3,tv+2),(tt[3]-tt[2])/CLOCKS_PER_SEC);
  fprintf(fp,"%lf\t%lf\n",tvcalc(tv+2,tv+1),(tt[2]-tt[1])/CLOCKS_PER_SEC);
  fprintf(fp,"%lf\t%lf\n",tvcalc(tv+4,tv+3),(tt[4]-tt[3])/CLOCKS_PER_SEC);
  fprintf(fp,"%lf\t%lf\n",tvcalc(tv+1,tv+0),tvcalc(tv+5,tv+4));
}

int main(int argc, char *argv[]){
    std::cout << "hello" << std::endl;
  int d,l,m,n,D,M,N,lp; char mode; double s,r,Np,sgmX,sgmY,*muX,*muY; double *u,*v,*w,*R,*t,*a,*sgm;
  pwpm pm; pwsz sz; double *x,*y,*X,*Y,*wd,**bX,**bY; int *wi; int sd=sizeof(double),si=sizeof(int); FILE *fp; char fn[256];
  int dsz,isz,ysz,xsz; char *ytraj=".optpath.bin",*xtraj=".optpathX.bin"; double tt[6]; struct timeval tv[6];
  int nx,ny,N0,M0=0; double rx,ry,*T,*X0,*Y0=NULL; double sgmT,*muT; double *pf;

  gettimeofday(tv+0,NULL); tt[0]=clock();
  /* read files */
  pw_getopt(&pm,argc,argv);
  bX=read2d(&N,&D,&mode,pm.fn[TARGET],"NA"); X= (double *)calloc(D*N,sd); sz.D=D;
  bY=read2d(&M,&D,&mode,pm.fn[SOURCE],"NA"); Y= (double *)calloc(D*M,sd);
  /* init: random number */
  init_genrand64(pm.rns?pm.rns:clock());
  /* check dimension */
  if(D  != sz.D){printf("ERROR: Dimensions of X and Y are incosistent. dim(X)=%d, dim(Y)=%d\n",sz.D,D);exit(EXIT_FAILURE);}
  if(N<=D||M<=D){printf("ERROR: #points must be greater than dimension\n");exit(EXIT_FAILURE);}
  /* change memory layout */
  for(d=0;d<D;d++)for(n=0;n<N;n++){X[d+D*n]=bX[n][d];}
  for(d=0;d<D;d++)for(m=0;m<M;m++){Y[d+D*m]=bY[m][d];}
  free2d(bX,N,D);
  free2d(bY,M,D);
  /* alias: size */
  sz.M=M; sz.J=pm.J;
  sz.N=N; sz.K=pm.K;
  /* check parameters */
  check_prms(pm,sz);
  /* print: paramters */
  if(!(pm.opt&PW_OPT_QUIET)) printInfo(sz,pm);
  /* normalization */
  muX= (double *)calloc(D,sd); muY= (double *)calloc(D,sd);
  if(!(pm.opt&PW_OPT_QUIET)&&(D==2||D==3)) print_norm(X,Y,D,N,M,1,pm.nrm);
  normalize_batch(X,muX,&sgmX,Y,muY,&sgmY,N,M,D,pm.nrm);
  if(!(pm.opt&PW_OPT_QUIET)&&(D==2||D==3)) print_norm(X,Y,D,N,M,0,pm.nrm);
  /* down-sampling */
  gettimeofday(tv+1,NULL); tt[1]=clock();
  nx=pm.dwn[TARGET]; rx=pm.dwr[TARGET];
  ny=pm.dwn[SOURCE]; ry=pm.dwr[SOURCE];
  if((nx||ny)&&!(pm.opt&PW_OPT_QUIET)) fprintf(stderr,"  Downsampling ...");
  if(nx){X0=X;N0=N;N=sz.N=nx;X=(double *)calloc(D*N,sd);downsample(X,N,X0,D,N0,rx);}
  if(ny){Y0=Y;M0=M;M=sz.M=ny;Y=(double *)calloc(D*M,sd);downsample(Y,M,Y0,D,M0,ry);}
  if((nx||ny)&&!(pm.opt&PW_OPT_QUIET)) fprintf(stderr," done. \n\n");
  gettimeofday(tv+2,NULL); tt[2]=clock();
  /* memory size */
  memsize(&dsz,&isz,sz,pm);
  /* memory size: x, y */
  ysz=D*M; ysz+=D*M*((pm.opt&PW_OPT_PATHY)?pm.nlp:0);
  xsz=D*M; xsz+=D*M*((pm.opt&PW_OPT_PATHX)?pm.nlp:0);
  /* allocaltion */
  wd= (double *)calloc(dsz,sd); x= (double *)calloc(xsz,sd); a= (double *)calloc(M,sd); u= (double *)calloc(D*M,sd); R= (double *)calloc(D*D,sd); sgm= (double *)calloc(M,sd);
  wi= (int *)calloc(isz,si); y= (double *)calloc(ysz,sd); w= (double *)calloc(M,sd); v= (double *)calloc(D*M,sd); t= (double *)calloc(D,  sd); pf = (double *)calloc(3*pm.nlp,sd);
  /* main computation */
  lp=bcpd(x,y,u,v,w,a,sgm,&s,R,t,&r,&Np,pf,wd,wi,X,Y,sz,pm);
  /* interpolation */
  gettimeofday(tv+3,NULL); tt[3]=clock();
  if(ny){
    if(!(pm.opt&PW_OPT_QUIET)) fprintf(stderr,"%s  Interpolating ... ",(pm.opt&PW_OPT_HISTO)?"\n":"");
    T= (double *)calloc(D*M0,sd);
    if(pm.opt&PW_OPT_1NN) interpolate_1nn(T,Y0,M0,v,Y,  &s,R,t,   sz,pm);
    else                  interpolate    (T,Y0,M0,x,Y,w,&s,R,t,&r,sz,pm);
    switch(pm.nrm){
      case 'e': sgmT=sgmX;muT=muX;  break;
      case 'x': sgmT=sgmX;muT=muX;  break;
      case 'y': sgmT=sgmY;muT=muY;  break;
      case 'n': sgmT=1.0f;muT=NULL; break;
      default:  sgmT=sgmX;muT=muX;
    }
    denormlize(T,muT,sgmT,M0,D);
    save_variable(pm.fn[OUTPUT],"y.interpolated.txt",T,D,M0,"%lf",TRANSPOSE); free(T);
    if(!(pm.opt&PW_OPT_QUIET)) fprintf(stderr,"done. \n\n");
  }
  gettimeofday(tv+4,NULL); tt[4]=clock();
  /* save correspondence */
  if((pm.opt&PW_OPT_SAVEP)|(pm.opt&PW_OPT_SAVEC)|(pm.opt&PW_OPT_SAVEE))
    save_corresp(pm.fn[OUTPUT],X,y,a,sgm,s,r,sz,pm);
  /* save trajectory */
  if(pm.opt&PW_OPT_PATHX) save_optpath(xtraj,x+D*M,X,sz,lp);
  if(pm.opt&PW_OPT_PATHY) save_optpath(ytraj,y+D*M,X,sz,lp);
  /* revert normalization */
  denormalize_batch(x,muX,sgmX,y,muY,sgmY,M,M,D,pm.nrm);
  /* save variables */
  if(pm.opt&PW_OPT_SAVEY) save_variable(pm.fn[OUTPUT],"y.txt",y,D,M,"%lf",TRANSPOSE);
  if(pm.opt&PW_OPT_SAVEX) save_variable(pm.fn[OUTPUT],"x.txt",x,D,M,"%lf",TRANSPOSE);
  if(pm.opt&PW_OPT_SAVEU) save_variable(pm.fn[OUTPUT],"u.txt",u,D,M,"%lf",TRANSPOSE);
  if(pm.opt&PW_OPT_SAVEV) save_variable(pm.fn[OUTPUT],"v.txt",v,D,M,"%lf",TRANSPOSE);
  if(pm.opt&PW_OPT_SAVEA) save_variable(pm.fn[OUTPUT],"a.txt",a,M,1,"%e", ASIS);
  if(pm.opt&PW_OPT_SAVET){
    save_variable(pm.fn[OUTPUT],"s.txt",&s,1,1,"%lf",ASIS);
    save_variable(pm.fn[OUTPUT],"R.txt", R,D,D,"%lf",ASIS);
    save_variable(pm.fn[OUTPUT],"t.txt", t,D,1,"%lf",ASIS);
  }
  if((pm.opt&PW_OPT_SAVEU)|(pm.opt&PW_OPT_SAVEV)|(pm.opt&PW_OPT_SAVET)){
    save_variable(pm.fn[OUTPUT],"normX.txt",X,D,N,"%lf",TRANSPOSE);
    save_variable(pm.fn[OUTPUT],"normY.txt",Y,D,M,"%lf",TRANSPOSE);
  }
  if((pm.opt&PW_OPT_SAVES)&&(pm.opt&PW_OPT_DBIAS)){
    for(m=0;m<M;m++) sgm[m]=sqrt(sgm[m]);
    save_variable(pm.fn[OUTPUT],"Sigma.txt",sgm,M,1,"%e",ASIS);
  }
  gettimeofday(tv+5,NULL); tt[5]=clock();
  /* output: computing time */
  if(!(pm.opt&PW_OPT_QUIET)) fprint_comptime(stderr,tv,tt,nx,ny);
  /* save total computing time */
  strcpy(fn,pm.fn[OUTPUT]); strcat(fn,"comptime.txt");
  fp=fopen(fn,"w"); if(pm.opt&PW_OPT_VTIME) fprint_comptime2(fp,tv,tt); else fprint_comptime(fp,tv,tt,nx,ny); fclose(fp);
  /* save computing time for each loop */
  if(pm.opt&PW_OPT_PFLOG){
    strcpy(fn,pm.fn[OUTPUT]); strcat(fn,"profile_time.txt"); fp=fopen(fn,"w");
    #ifdef MINGW32
    for(l=0;l<lp;l++){fprintf(fp,"%f\n",pf[l]);}
    #else
    for(l=0;l<lp;l++){fprintf(fp,"%f\t%f\n",pf[l],pf[l+pm.nlp]);}
    #endif
    fclose(fp);
    strcpy(fn,pm.fn[OUTPUT]); strcat(fn,"profile_sigma.txt"); fp=fopen(fn,"w");
    for(l=0;l<lp;l++){fprintf(fp,"%f\n",pf[l+2*pm.nlp]);}
  }
  /* output info` */
  if(pm.opt&PW_OPT_INFO){
    strcpy(fn,pm.fn[OUTPUT]); strcat(fn,"info.txt");
    fp=fopen(fn,"w");
    fprintf(fp,"loops\t%d\n",  lp);
    fprintf(fp,"sigma\t%lf\n", r);
    fprintf(fp,"N_hat\t%lf\n", Np);
    fclose(fp);
  }
  if((pm.opt&PW_OPT_PATHY)&&!(pm.opt&PW_OPT_QUIET))
    fprintf(stderr,"  ** Search path during optimization was saved to: [%s]\n\n",ytraj);

  free(x);free(X);free(wd);free(muX);free(a);free(v);free(sgm);
  free(y);free(Y);free(wi);free(muY);free(u);free(R);free(t);free(pf);

  return 0;
}

