/* hsmc2d.c -- NVT hard-disk Monte Carlo in 2D, positive control.
 *
 * Purpose: run the SAME pipeline (shell percolation, isoperimetric quotient) on
 * a system whose answer is independently known -- de Graaf reports dynamical
 * arrest at phi_a ~ 0.777 for bidisperse disks over a wide range of size ratio.
 * If the pipeline recovers a feature there, the 3D null result is meaningful.
 * If it does not, the observable is blind and the 3D result says nothing.
 *
 * Composition matches the reference study: 1:1 bidisperse, R^-1 = 1.4.
 *
 *   ./hsmc2d <phi> <n_tri> <mode:0=binary1.4,1=mono,2=binary1.7> <seed> <eq> <prod> <prefix>
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static int N;
static double L,*x,*y,*xu,*yu,*rad,rmax;
static int *head,*lst,ncell; static double cellsize;
static unsigned long long S0,S1v;
static unsigned long long rotl(unsigned long long a,int k){return (a<<k)|(a>>(64-k));}
static unsigned long long nextr(void){unsigned long long s0=S0,s1=S1v,r=s0+s1;
    s1^=s0;S0=rotl(s0,55)^s1^(s1<<14);S1v=rotl(s1,36);return r;}
static double rnd(void){return (nextr()>>11)*(1.0/9007199254740992.0);}

static void build_cells(void){
    ncell=(int)floor(L/(2.0*rmax)); if(ncell<3)ncell=3; cellsize=L/ncell;
    free(head);free(lst); head=malloc(sizeof(int)*ncell*ncell); lst=malloc(sizeof(int)*N);
    for(int i=0;i<ncell*ncell;i++)head[i]=-1;
    for(int i=0;i<N;i++){int cx=(int)(x[i]/cellsize),cy=(int)(y[i]/cellsize);
        if(cx>=ncell)cx=ncell-1; if(cy>=ncell)cy=ncell-1;
        int c=cx*ncell+cy; lst[i]=head[c]; head[c]=i;}
}
static inline int cf(double a){int c=(int)(a/cellsize); if(c>=ncell)c=ncell-1; if(c<0)c=0; return c;}
static void cell_remove(int i,double px,double py){int c=cf(px)*ncell+cf(py);int p=head[c];
    if(p==i){head[c]=lst[i];return;} while(lst[p]!=i)p=lst[p]; lst[p]=lst[i];}
static void cell_add(int i){int c=cf(x[i])*ncell+cf(y[i]); lst[i]=head[c]; head[c]=i;}
static inline double pbc(double d){if(d>0.5*L)d-=L; else if(d<-0.5*L)d+=L; return d;}

static double loc(int i,double px,double py,double ri,int hard){
    double e=0; int cx=cf(px),cy=cf(py);
    for(int a=-1;a<=1;a++)for(int b=-1;b<=1;b++){
        int ix=(cx+a+ncell)%ncell, iy=(cy+b+ncell)%ncell;
        for(int j=head[ix*ncell+iy];j>=0;j=lst[j]){
            if(j==i)continue;
            double dx=pbc(px-x[j]),dy=pbc(py-y[j]),d2=dx*dx+dy*dy,s=ri+rad[j];
            if(d2<s*s){ if(hard)return 1.0; double o=s-sqrt(d2); e+=o*o; }
        }
    }
    return e;
}
static double tot(void){double e=0;for(int i=0;i<N;i++)e+=loc(i,x[i],y[i],rad[i],0);return 0.5*e;}

int main(int argc,char**argv){
    if(argc<8){fprintf(stderr,"usage: hsmc2d phi n_tri mode seed eq prod prefix\n");return 1;}
    double phi=atof(argv[1]); int nt=atoi(argv[2]),mode=atoi(argv[3]);
    unsigned long long seed=strtoull(argv[4],0,10);
    long eq=atol(argv[5]),prod=atol(argv[6]); const char*pre=argv[7];
    S0=seed*6364136223846793005ULL+1442695040888963407ULL; S1v=seed^0x9E3779B97F4A7C15ULL;
    for(int i=0;i<20;i++)nextr();

    int ny=nt, nx=(int)lround(sqrt(3.0)*nt);
    N=2*nx*ny;   /* near-triangular lattice in a square box: 2 sites per rectangle */
    x=malloc(N*8);y=malloc(N*8);xu=malloc(N*8);yu=malloc(N*8);rad=malloc(N*8);
    head=0;lst=0;
    double*rt=malloc(N*8),asum=0;
    for(int i=0;i<N;i++){
        if(mode==0)      rt[i]=(i%2)?0.5:0.5/1.4;
        else if(mode==2) rt[i]=(i%2)?0.5:0.5/1.7;
        else             rt[i]=0.5;
        asum+=M_PI*rt[i]*rt[i];
    }
    L=sqrt(asum/phi);
    /* triangular lattice: a1=(a,0), a2=(a/2, a*sqrt3/2); box is a nt x nt cell */
    double a=L/nt, ay=L/(nt*sqrt(3.0)/2.0*1.0);
    (void)ay;
    /* use an orthogonal cell of nt x nt with 2-site basis (rectangular centred) */
    double bx=L/nx, by=L/ny;
    int k=0;
    for(int i=0;i<nx;i++)for(int j=0;j<ny;j++){
        x[k]=i*bx; y[k]=j*by; k++;
        x[k]=(i+0.5)*bx; y[k]=(j+0.5)*by; k++;
    }
    double r0=sqrt(phi*L*L/(N*M_PI));
    for(int i=0;i<N;i++)rad[i]=r0;
    rmax=r0; build_cells();
    char fn[512]; sprintf(fn,"%s.log",pre); FILE*lg=fopen(fn,"w");
    double nnd=sqrt(bx*bx+by*by)/2.0; double nnd2=bx<by?bx:by;
    double nnmin=nnd<nnd2?nnd:nnd2;
    fprintf(lg,"N=%d phi=%.4f L=%.6f r0=%.6f nn=%.6f\n",N,phi,L,r0,nnmin);
    if(2*r0>nnmin+1e-12){fprintf(lg,"FATAL: init overlaps (phi too high for this cell)\n");
        fclose(lg);return 2;}

    double dmax=0.06*2*r0;
    for(long s=0;s<20000;s++)for(int t=0;t<N;t++){
        int i=(int)(rnd()*N); if(i>=N)i=N-1;
        double px=x[i]+(2*rnd()-1)*dmax,py=y[i]+(2*rnd()-1)*dmax;
        px-=L*floor(px/L); py-=L*floor(py/L);
        if(loc(i,px,py,rad[i],1)==0.0){cell_remove(i,x[i],y[i]);x[i]=px;y[i]=py;cell_add(i);}
    }
    for(int i=0;i<N;i++)rad[i]=rt[i];
    rmax=0; for(int i=0;i<N;i++)if(rad[i]>rmax)rmax=rad[i];
    build_cells();
    double T=1e-2,E=tot(); long an=0;
    while(E>0&&an<6000000){
        for(int t=0;t<N;t++){
            int i=(int)(rnd()*N); if(i>=N)i=N-1;
            double px=x[i]+(2*rnd()-1)*dmax,py=y[i]+(2*rnd()-1)*dmax;
            px-=L*floor(px/L); py-=L*floor(py/L);
            double eo=loc(i,x[i],y[i],rad[i],0),en=loc(i,px,py,rad[i],0);
            if(en<=eo||rnd()<exp(-(en-eo)/T)){cell_remove(i,x[i],y[i]);x[i]=px;y[i]=py;cell_add(i);}
        }
        an++; if(an%50==0){T*=0.9;E=tot();}
    }
    E=tot(); fprintf(lg,"anneal=%ld E=%.3e\n",an,E);
    if(E>0){fprintf(lg,"FATAL: overlaps remain\n");fclose(lg);return 3;}

    long acc=0,att=0;
    for(long s=0;s<eq;s++){
        for(int t=0;t<N;t++){
            int i=(int)(rnd()*N); if(i>=N)i=N-1;
            double px=x[i]+(2*rnd()-1)*dmax,py=y[i]+(2*rnd()-1)*dmax;
            px-=L*floor(px/L); py-=L*floor(py/L); att++;
            if(loc(i,px,py,rad[i],1)==0.0){cell_remove(i,x[i],y[i]);x[i]=px;y[i]=py;cell_add(i);acc++;}
        }
        for(int t=0;t<N/5;t++){
            int i=(int)(rnd()*N),j=(int)(rnd()*N);
            if(i>=N)i=N-1; if(j>=N)j=N-1; if(i==j)continue;
            double ri=rad[i],rj=rad[j]; if(fabs(ri-rj)<1e-15)continue;
            rad[i]=rj;rad[j]=ri;
            if(loc(i,x[i],y[i],rad[i],1)>0.0||loc(j,x[j],y[j],rad[j],1)>0.0){rad[i]=ri;rad[j]=rj;}
        }
        if(s%200==199){double r=(double)acc/att;dmax*=(r>0.3?1.02:0.98);
            if(dmax>0.3*rmax)dmax=0.3*rmax; acc=att=0;}
    }
    for(int i=0;i<N;i++){xu[i]=x[i];yu[i]=y[i];}
    double*u0x=malloc(N*8),*u0y=malloc(N*8);
    for(int i=0;i<N;i++){u0x[i]=xu[i];u0y[i]=yu[i];}
    sprintf(fn,"%s.msd",pre); FILE*fp=fopen(fn,"w");
    sprintf(fn,"%s.cfg",pre); FILE*fc=fopen(fn,"w");
    long every=prod/300; if(every<1)every=1;
    long snap=prod/12; if(snap<1)snap=1;
    for(long s=0;s<prod;s++){
        for(int t=0;t<N;t++){
            int i=(int)(rnd()*N); if(i>=N)i=N-1;
            double ddx=(2*rnd()-1)*dmax,ddy=(2*rnd()-1)*dmax;
            double px=x[i]+ddx,py=y[i]+ddy;
            px-=L*floor(px/L); py-=L*floor(py/L);
            if(loc(i,px,py,rad[i],1)==0.0){cell_remove(i,x[i],y[i]);x[i]=px;y[i]=py;cell_add(i);
                xu[i]+=ddx;yu[i]+=ddy;}
        }
        if(s%every==0){double m=0;
            for(int i=0;i<N;i++){double dx=xu[i]-u0x[i],dy=yu[i]-u0y[i];m+=dx*dx+dy*dy;}
            fprintf(fp,"%ld %.8e\n",s,m/N);}
        if(s%snap==0&&s>0){fprintf(fc,"%d %.10f\n",N,L);
            for(int i=0;i<N;i++)fprintf(fc,"%.8f %.8f %.8f\n",x[i],y[i],rad[i]);}
    }
    fclose(fp);fclose(fc);
    int bad=0; for(int i=0;i<N;i++)if(loc(i,x[i],y[i],rad[i],1)>0.0)bad++;
    fprintf(lg,"final overlap audit: %d\n",bad); fclose(lg);
    return bad?4:0;
}
