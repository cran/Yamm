
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>

# include <Rinternals.h>
# include <Rmath.h>
# include <R_ext/Linpack.h>



int cmpfunc(const void * a, const void * b) {
    return (*(double*)a > *(double*)b) ? 1 : (*(double*)a < *(double*)b) ? -1:0;
    
}


void CPmedTrapz2D(double *xs, double *ans, int *div, int *ndat)
{
    int i,j,k;
    int ndata = *ndat;
    int d = *div;
    double arg;
    double h=2*M_PI/(double)d;
    double  *a;
    double  *y;

    
    y= malloc(sizeof(double) * ndata);
    a= malloc(sizeof(double) * 2);

    
    
    for(k=0; k<= d; k++){
        
        arg=h*k;
       
        a[0]=cos(arg);
        a[1]=sin(arg);
        
        for (i = 0; i <ndata; i++)
        {
            y[i]=0.0;
            for (j = 0; j <2; j++)
            {
                y[i] += xs[i*2+j]*a[j];
                
            }
        }
        
        
        qsort(y,ndata,sizeof(y[0]),cmpfunc);
        
        
        double projmed;
        if(ndata % 2 ==0){
            int median = ndata/2;
            projmed = (y[median]+y[median-1])/2.0;
        }
        else{
            int median = ((ndata+1)/2)-1;
            projmed = y[median];
        }
        
     
        if(k==0||k==d){
            ans[0] += a[0]*projmed;
            ans[1] += a[1]*projmed;
        }
        else{
            ans[0] += 2*a[0]*projmed;
            ans[1] += 2*a[1]*projmed;
        }
      
    }
    
    
    ans[0] = ans[0]*h/(2.0*M_PI);
    ans[1] = ans[1]*h/(2.0*M_PI);
    
    free((void*) y);
    free((void*) a);
}



void CPmedTrapz3D(double *xs, double *myans, int *divm, int *divl, int *ndat)
{
    int i,j,p,q;
    int ndata = *ndat;
    int m = *divm;
    int l = *divl;
    double phi;
    double theta;
    double projmed;
    double h=M_PI/(double)m;
    double k=2*M_PI/(double)l;
    double  *a;
    double  *y;
    double  *ansx;
    double  *ansy;
    double  *ansz;
    double  *ans;
    
    y= malloc(sizeof(double) * ndata);
    a= malloc(sizeof(double) * 3);
    ansx = malloc(sizeof(double) * (l+1));
    ansy = malloc(sizeof(double) * (l+1));
    ansz = malloc(sizeof(double) * (l+1));
    ans = malloc(sizeof(double) * 3);
    
    
    for(p=0; p<= l; p++){
        theta=k*p;
        ans[0]=0.0;
        ans[1]=0.0;
        ans[2]=0.0;
        for(q=0; q<= m; q++){
            phi=h*q;
            a[0]=cos(phi);
            a[1]=sin(phi)*cos(theta);
            a[2]=sin(phi)*sin(theta);
            
            for (i = 0; i <ndata; i++)
            {
                y[i]=0.0;
                for (j = 0; j <3; j++)
                {
                    y[i] += xs[i*3+j]*a[j];
                    
                }
            }
            
            qsort(y,ndata,sizeof(y[0]),cmpfunc);
            
            if(ndata % 2 ==0){
                int median = ndata/2;
                projmed = (y[median]+y[median-1])/2.0;
            }
            else{
                int median = ((ndata+1)/2)-1;
                projmed = y[median];
            }

            
            projmed=projmed*sin(phi);
            
            if(q==0||q==m){
                ans[0] += a[0]*projmed;
                ans[1] += a[1]*projmed;
                ans[2] += a[2]*projmed;
            }
            else{
                ans[0] += 2*a[0]*projmed;
                ans[1] += 2*a[1]*projmed;
                ans[2] += 2*a[2]*projmed;
                
            }
        }
        
        ansx[p] = ans[0]*h/2.0;
        ansy[p] = ans[1]*h/2.0;
        ansz[p] = ans[2]*h/2.0;
        
        if(p==0||p==l){
            myans[0] += ansx[p];
            myans[1] += ansy[p];
            myans[2] += ansz[p];
        }
        else{
            myans[0] += 2*ansx[p];
            myans[1] += 2*ansy[p];
            myans[2] += 2*ansz[p];
        }
    }
    myans[0] = 3*myans[0]*k/(8.0*M_PI);
    myans[1] = 3*myans[1]*k/(8.0*M_PI);
    myans[2] = 3*myans[2]*k/(8.0*M_PI);
    
    
    free((void*) y);
    free((void*) a);
    free((void*) ansx);
    free((void*) ansy);
    free((void*) ansz);
}



void CPmedMCInt(double *xs, double *ans, int *nproj, int *ndat, int *nd)
{
    int i,j,k;
    int nprojs = *nproj;
    int ndim = *nd;
    int ndata = *ndat;
    double vdp;
    double svdp;
    double  *a;
    double  *y;
    
    
    
    y= malloc(sizeof(double) * ndata);
    a= malloc(sizeof(double) * ndim);
    
    
    
    
    for(k=1; k<= nprojs; k++){

        vdp=0.0;
        svdp=0.0;
        for(i = 0; i <ndim; i++)
        {
            GetRNGstate();
            a[i] = 2*(double)unif_rand()-1.0;
            PutRNGstate();
            vdp += a[i] * a[i];
        }
        
        
        svdp = sqrt(vdp);
        
        for(i = 0; i <ndim; i++)
        {
            a[i] = a[i]/svdp;

        }
        
        
        for (i = 0; i <ndata; i++)
        {
            y[i]=0.0;
            for (j = 0; j <ndim; j++)
            {
                y[i] += xs[i*ndim+j]*a[j];
                
            }
  
        }
        
        
        qsort(y,ndata,sizeof(y[0]),cmpfunc);
        
        
        double projmed;
        if(ndata % 2 ==0){
            int median = ndata/2;
            projmed = (y[median]+y[median-1])/2.0;
        }
        else{
            int median = ((ndata+1)/2)-1;
            projmed = y[median];
        }
      
        
        for(i = 0; i <ndim; i++)
        {
            a[i] = a[i]*projmed;
            ans[i] +=a[i];

        }
        
    }
    
    for(i = 0; i <ndim; i++)
    {
        ans[i] = ndim*ans[i]/(double)nprojs;
    }
    
    free((void*) y);
    free((void*) a);
}



/*
 * Wirth's public domain median function
 */

/*
 * Algorithm from N. Wirth's book, implementation by N. Devillard.
 * This code in public domain.
 */


typedef double elem_type ;

#define ELEM_SWAP(a,b) { register elem_type t=(a);(a)=(b);(b)=t; }


/*---------------------------------------------------------------------------
 Function :   kth_smallest()
 In       :   array of elements, # of elements in the array, rank k
 Out      :   one element
 Job      :   find the kth smallest element in the array
 Notice   :   use the median() macro defined below to get the median.
 
 Reference:
 
 Author: Wirth, Niklaus
 Title: Algorithms + data structures = programs
 Publisher: Englewood Cliffs: Prentice-Hall, 1976
 Physical description: 366 p.
 Series: Prentice-Hall Series in Automatic Computation
 
 ---------------------------------------------------------------------------*/


elem_type kth_smallest(elem_type a[], int n, int k)
{
    register int i,j,l,m ;
    register elem_type x ;
    
    l=0 ; m=n-1 ;
    while (l<m) {
        x=a[k] ;
        i=l ;
        j=m ;
        do {
            while (a[i]<x) i++ ;
            while (x<a[j]) j-- ;
            if (i<=j) {
                ELEM_SWAP(a[i],a[j]) ;
                i++ ; j-- ;
            }
        } while (i<=j) ;
        if (j<k) l=i ;
        if (k<i) m=j ;
    }
    return a[k] ;
}


#define median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))


/*
 * Compute Yamm estimator, via using the actual median
 *
 * x - input data (a "matrix" with *ndat observations on *nd dimensions)
 * mu - suggested "guessed" median value (vector of length *nd)
 * ave - the output value (vector of length *nd)
 * nproj - the number of projections to use to compute it
 * ndat - the number of observations
 * nd - the dimensionality
 * dofabs - whether to use absolute value, or square in the Yamm objective
 * y - workspace (of length number of observations)
 * a - workspace (dimension *nd, holds projection vectors)
 * xs - workspace, stores mu-centred data
 *
 */


void Cyammobj(double *x, double *mu, double *ave, int *nproj, int *ndat, int *nd, int *dofabs, double *y, double *a, double *xs)
{
    int i,j,k;
    int nprojs = *nproj;
    int ndim = *nd;
    int ndata = *ndat;
    double out_sq;
    double vdp;
    double svdp;
    double sum_vec=0.0;
    
    /* Centre the data in x, by mu vector */
    
    for (i = 0; i < ndata; i++)
    {
        for (j = 0; j < ndim; j++)
        {
            xs[i*ndim+j] = x[i*ndim+j]-mu[j];
        }
    }
    
    /* Retrieve the random number generator state */
    
    GetRNGstate();
    
    /* Run through all projections, computing the Yamm objective */
    
    for(k=1; k<= nprojs; k++){
        vdp=0.0;
        svdp=0.0;
        
        /* Generate random projection */
        
        for(i = 0; i <ndim; i++)
        {
            a[i] = 2*(double)unif_rand()-1.0;
            //            printf("%12.6f\n",a[i]);
            vdp += a[i] * a[i];
        }
        
        /* Normalise random projection to have length 1 */
        
        svdp = sqrt(vdp);
        
        for(i = 0; i <ndim; i++)
        {
            a[i] = a[i]/svdp;
        }
        
        /* Project centred data onto projection in a */
        
        for (i = 0; i <ndata; i++)
        {
            y[i]=0.0;
            for (j = 0; j <ndim; j++)
            {
                y[i] += xs[i*ndim+j]*a[j];
                
            }
        }
        
        /*
         * Compute and store the median of y vector in yamm
         */
        
        double yamm;
        
        yamm = median(y, ndata);
        
        if (*dofabs==1)
            out_sq=fabs(yamm);
        else
            out_sq=yamm*yamm;
        sum_vec +=out_sq;
    }
    
    /* Save random number generation state */
    
    PutRNGstate();
    
    double ans=sum_vec/(double)nprojs;
    *ave=ans;
}


