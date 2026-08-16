/* hsmc.c -- NVT hard-sphere Monte Carlo in 3D with cell lists.
 *
 * Local single-particle displacement moves give a stochastic dynamics that is a
 * standard proxy for Brownian motion (the experiments this study targets are
 * Brownian).  Swap moves are used ONLY during equilibration; they destroy the
 * physical dynamics and are switched off for the production run that measures
 * the MSD.  This is stated in the spec and must not be changed silently.
 *
 * Initialisation: FCC at the target packing fraction with all radii equal
 * (always overlap-free for eta <= 0.7405), then radii are set to their target
 * polydisperse values and the resulting overlaps are annealed away with a
 * soft potential E = sum (sigma_ij - r_ij)^2 at decreasing temperature until
 * E == 0 exactly.  Production never starts from a configuration with E > 0.
 *
 *   ./hsmc <eta> <n_fcc> <mode:0=poly,1=binary> <seed> <eq_sweeps> <prod_sweeps> <prefix>
 *
 * writes  <prefix>.cfg   x y z r   (final configuration)
 *         <prefix>.msd   sweeps  msd
 *         <prefix>.log   diagnostics
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static int N;
static double L, *x, *y, *z, *xu, *yu, *zu, *rad, rmax;
static int *head, *lst, ncell;
static double cellsize;

/* ---- rng: xoshiro-ish, deterministic and seeded ---- */
static unsigned long long S0, S1v;
static unsigned long long rotl(unsigned long long a, int k){return (a<<k)|(a>>(64-k));}
static unsigned long long nextr(void){
    unsigned long long s0=S0,s1=S1v,r=s0+s1;
    s1^=s0; S0=rotl(s0,55)^s1^(s1<<14); S1v=rotl(s1,36); return r;
}
static double rnd(void){ return (nextr()>>11)*(1.0/9007199254740992.0); }
static double gauss(void){
    double u1=rnd(),u2=rnd();
    if(u1<1e-300)u1=1e-300;
    return sqrt(-2.0*log(u1))*cos(2*M_PI*u2);
}

/* ---- cell list ---- */
static void build_cells(void){
    ncell=(int)floor(L/(2.0*rmax));
    if(ncell<3)ncell=3;
    cellsize=L/ncell;
    free(head); free(lst);
    head=malloc(sizeof(int)*ncell*ncell*ncell);
    lst=malloc(sizeof(int)*N);
    for(int i=0;i<ncell*ncell*ncell;i++)head[i]=-1;
    for(int i=0;i<N;i++){
        int cx=(int)(x[i]/cellsize),cy=(int)(y[i]/cellsize),cz=(int)(z[i]/cellsize);
        if(cx>=ncell)cx=ncell-1; if(cy>=ncell)cy=ncell-1; if(cz>=ncell)cz=ncell-1;
        int c=(cx*ncell+cy)*ncell+cz; lst[i]=head[c]; head[c]=i;
    }
}
static inline int cellof(double a){int c=(int)(a/cellsize); if(c>=ncell)c=ncell-1; if(c<0)c=0; return c;}
static void cell_remove(int i,double px,double py,double pz){
    int c=(cellof(px)*ncell+cellof(py))*ncell+cellof(pz);
    int p=head[c];
    if(p==i){head[c]=lst[i];return;}
    while(lst[p]!=i)p=lst[p];
    lst[p]=lst[i];
}
static void cell_add(int i){
    int c=(cellof(x[i])*ncell+cellof(y[i]))*ncell+cellof(z[i]);
    lst[i]=head[c]; head[c]=i;
}
static inline double pbc(double d){ if(d>0.5*L)d-=L; else if(d<-0.5*L)d+=L; return d; }

/* overlap test / soft overlap energy for particle i placed at (px,py,pz) with radius ri */
static double local_energy(int i,double px,double py,double pz,double ri,int hard){
    double e=0.0;
    int cx=cellof(px),cy=cellof(py),cz=cellof(pz);
    for(int a=-1;a<=1;a++)for(int b=-1;b<=1;b++)for(int c=-1;c<=1;c++){
        int ix=(cx+a+ncell)%ncell, iy=(cy+b+ncell)%ncell, iz=(cz+c+ncell)%ncell;
        for(int j=head[(ix*ncell+iy)*ncell+iz]; j>=0; j=lst[j]){
            if(j==i)continue;
            double dx=pbc(px-x[j]),dy=pbc(py-y[j]),dz=pbc(pz-z[j]);
            double d2=dx*dx+dy*dy+dz*dz, s=ri+rad[j];
            if(d2<s*s){
                if(hard)return 1.0;
                double d=sqrt(d2), o=s-d; e+=o*o;
            }
        }
    }
    return e;
}
static double total_energy(void){
    double e=0;
    for(int i=0;i<N;i++)e+=local_energy(i,x[i],y[i],z[i],rad[i],0);
    return 0.5*e;
}

int main(int argc,char**argv){
    if(argc<8){fprintf(stderr,"usage: hsmc eta n_fcc mode seed eq prod prefix\n");return 1;}
    double eta=atof(argv[1]); int nf=atoi(argv[2]), mode=atoi(argv[3]);
    unsigned long long seed=strtoull(argv[4],0,10);
    long eq=atol(argv[5]), prod=atol(argv[6]); const char*pre=argv[7];
    S0=seed*6364136223846793005ULL+1442695040888963407ULL; S1v=seed^0x9E3779B97F4A7C15ULL;
    for(int i=0;i<20;i++)nextr();

    N=4*nf*nf*nf;
    x=malloc(N*8);y=malloc(N*8);z=malloc(N*8);
    xu=malloc(N*8);yu=malloc(N*8);zu=malloc(N*8);rad=malloc(N*8);
    head=0;lst=0;

    /* target radii */
    double *rt=malloc(N*8), vsum=0;
    for(int i=0;i<N;i++){
        if(mode==0){ double g; do{g=1.0+0.10*gauss();}while(g<0.6||g>1.4); rt[i]=0.5*g; }
        else if(mode==1){ rt[i]=(i%2)?0.5:0.5*0.714; }
        else if(mode==3){ double g; do{g=1.0+0.06*gauss();}while(g<0.7||g>1.3); rt[i]=0.5*g; }
        else            { rt[i]=0.5; }
        vsum+=(4.0/3.0)*M_PI*rt[i]*rt[i]*rt[i];
    }
    L=cbrt(vsum/eta);
    double a=L/nf;                              /* FCC lattice constant */
    double bs[4][3]={{0,0,0},{0,.5,.5},{.5,0,.5},{.5,.5,0}};
    int k=0;
    for(int i=0;i<nf;i++)for(int j=0;j<nf;j++)for(int m=0;m<nf;m++)for(int b=0;b<4;b++){
        x[k]=(i+bs[b][0])*a; y[k]=(j+bs[b][1])*a; z[k]=(m+bs[b][2])*a; k++;
    }
    /* start monodisperse at the same eta -> guaranteed overlap-free on FCC */
    double r0=cbrt(eta*L*L*L/(N*(4.0/3.0)*M_PI));
    for(int i=0;i<N;i++)rad[i]=r0;
    rmax=r0; build_cells();

    FILE*lg; char fn[512]; sprintf(fn,"%s.log",pre); lg=fopen(fn,"w");
    fprintf(lg,"N=%d eta=%.4f L=%.6f a=%.6f r0=%.6f nn=%.6f\n",N,eta,L,a,r0,a/sqrt(2.0));
    if(2*r0 > a/sqrt(2.0)+1e-12){fprintf(lg,"FATAL: monodisperse FCC init overlaps\n");return 2;}

    /* melt the monodisperse crystal */
    double dmax=0.06*2*r0;
    for(long s=0;s<20000;s++)for(int t=0;t<N;t++){
        int i=(int)(rnd()*N); if(i>=N)i=N-1;
        double px=x[i]+(2*rnd()-1)*dmax, py=y[i]+(2*rnd()-1)*dmax, pz=z[i]+(2*rnd()-1)*dmax;
        px-=L*floor(px/L); py-=L*floor(py/L); pz-=L*floor(pz/L);
        if(local_energy(i,px,py,pz,rad[i],1)==0.0){
            cell_remove(i,x[i],y[i],z[i]); x[i]=px;y[i]=py;z[i]=pz; cell_add(i);
        }
    }
    /* switch on the target size distribution and anneal the overlaps to zero */
    for(int i=0;i<N;i++)rad[i]=rt[i];
    rmax=0; for(int i=0;i<N;i++)if(rad[i]>rmax)rmax=rad[i];
    build_cells();
    double T=1e-2, E=total_energy();
    long anneal=0;
    while(E>0 && anneal<4000000){
        for(int t=0;t<N;t++){
            int i=(int)(rnd()*N); if(i>=N)i=N-1;
            double px=x[i]+(2*rnd()-1)*dmax, py=y[i]+(2*rnd()-1)*dmax, pz=z[i]+(2*rnd()-1)*dmax;
            px-=L*floor(px/L); py-=L*floor(py/L); pz-=L*floor(pz/L);
            double eo=local_energy(i,x[i],y[i],z[i],rad[i],0);
            double en=local_energy(i,px,py,pz,rad[i],0);
            if(en<=eo || rnd()<exp(-(en-eo)/T)){
                cell_remove(i,x[i],y[i],z[i]); x[i]=px;y[i]=py;z[i]=pz; cell_add(i);
            }
        }
        anneal++; if(anneal%50==0){T*=0.9; E=total_energy();}
    }
    E=total_energy();
    fprintf(lg,"anneal sweeps=%ld  final soft energy=%.3e\n",anneal,E);
    if(E>0){fprintf(lg,"FATAL: could not remove overlaps at eta=%.4f\n",eta);fclose(lg);return 3;}

    /* ---- equilibration: local moves + swap moves ---- */
    long acc=0,att=0;
    for(long s=0;s<eq;s++){
        for(int t=0;t<N;t++){
            int i=(int)(rnd()*N); if(i>=N)i=N-1;
            double px=x[i]+(2*rnd()-1)*dmax, py=y[i]+(2*rnd()-1)*dmax, pz=z[i]+(2*rnd()-1)*dmax;
            px-=L*floor(px/L); py-=L*floor(py/L); pz-=L*floor(pz/L);
            att++;
            if(local_energy(i,px,py,pz,rad[i],1)==0.0){
                cell_remove(i,x[i],y[i],z[i]); x[i]=px;y[i]=py;z[i]=pz; cell_add(i); acc++;
            }
        }
        for(int t=0;t<N/5;t++){          /* swap moves: equilibration only */
            int i=(int)(rnd()*N), j=(int)(rnd()*N);
            if(i>=N)i=N-1; if(j>=N)j=N-1; if(i==j)continue;
            double ri=rad[i],rj=rad[j];
            if(fabs(ri-rj)<1e-15)continue;
            rad[i]=rj; rad[j]=ri;
            if(local_energy(i,x[i],y[i],z[i],rad[i],1)>0.0 ||
               local_energy(j,x[j],y[j],z[j],rad[j],1)>0.0){ rad[i]=ri; rad[j]=rj; }
        }
        if(s%200==199){ double r=(double)acc/att; dmax*=(r>0.3?1.02:0.98);
                        if(dmax>0.3*rmax)dmax=0.3*rmax; acc=att=0; }
    }
    fprintf(lg,"equilibration done, dmax=%.5f\n",dmax);

    /* ---- production: local moves only, MSD in unfolded coordinates ---- */
    for(int i=0;i<N;i++){xu[i]=x[i];yu[i]=y[i];zu[i]=z[i];}
    double *u0x=malloc(N*8),*u0y=malloc(N*8),*u0z=malloc(N*8);
    for(int i=0;i<N;i++){u0x[i]=xu[i];u0y[i]=yu[i];u0z[i]=zu[i];}
    char fm[512]; sprintf(fm,"%s.msd",pre); FILE*fp=fopen(fm,"w");
    char fcn[512]; sprintf(fcn,"%s.cfg",pre); FILE*fcf=fopen(fcn,"w");
    long every=prod/300; if(every<1)every=1;
    long snap=prod/12; if(snap<1)snap=1;
    acc=att=0;
    for(long s=0;s<prod;s++){
        for(int t=0;t<N;t++){
            int i=(int)(rnd()*N); if(i>=N)i=N-1;
            double ddx=(2*rnd()-1)*dmax, ddy=(2*rnd()-1)*dmax, ddz=(2*rnd()-1)*dmax;
            double px=x[i]+ddx, py=y[i]+ddy, pz=z[i]+ddz;
            px-=L*floor(px/L); py-=L*floor(py/L); pz-=L*floor(pz/L);
            att++;
            if(local_energy(i,px,py,pz,rad[i],1)==0.0){
                cell_remove(i,x[i],y[i],z[i]); x[i]=px;y[i]=py;z[i]=pz; cell_add(i);
                xu[i]+=ddx; yu[i]+=ddy; zu[i]+=ddz; acc++;
            }
        }
        if(s%every==0){
            double m=0;
            for(int i=0;i<N;i++){
                double dx=xu[i]-u0x[i],dy=yu[i]-u0y[i],dz=zu[i]-u0z[i];
                m+=dx*dx+dy*dy+dz*dz;
            }
            fprintf(fp,"%ld %.8e\n",s,m/N);
        }
        if(s%snap==0 && s>0){
            fprintf(fcf,"%d %.10f\n",N,L);
            for(int i=0;i<N;i++)fprintf(fcf,"%.8f %.8f %.8f %.8f\n",x[i],y[i],z[i],rad[i]);
        }
    }
    fclose(fcf);
    fclose(fp);
    fprintf(lg,"production acceptance=%.4f\n",(double)acc/att);

    /* final overlap audit */
    int bad=0;
    for(int i=0;i<N;i++)if(local_energy(i,x[i],y[i],z[i],rad[i],1)>0.0)bad++;
    fprintf(lg,"final overlap audit: %d overlapping particles\n",bad);
    fclose(lg);
    return bad?4:0;
}
