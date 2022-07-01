#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "tools.h"

/* Settings which makes the user do a CTRL-C break out of the loop*/
#if defined(LIBUT) && (defined(_WIN32) || defined(__WIN32__) )

#define STOPMARK utIsInterruptPending()
#define INITBREAK ;
bool utIsInterruptPending(void);

#else

#include <signal.h>
#define INITBREAK   sigint_cont = 1;  (void) signal(SIGINT , ex_sigint);
#define STOPMARK sigint_cont==0
int sigint_cont = 1;
void ex_sigint(int sig) {
	sigint_cont = 0;
}
#endif

void tv_restore_core(double *x,double *y,int *Ii,int *Ic,double gamma,double d,double eps,double L,double mu,int m,int n,int sI,int sIc,int maxiter,double *kf, double *epsilon_kf){

  register double *df,*yk,*wk,*zk,*t1,*t2,*yic,*uij;

  double pobj = 0, dobj = 0, c1, c2;
  int i, j, k, mn = m*n, one = 1;
  double A_kp1=0.5,alpha_kp1,t_k,m1t_k;
  double Ld = L*d;
	register int i1,i2,i3;

  INITBREAK

  df = malloc(m*n*sizeof(double));
  yk = malloc(sI*sizeof(double));
  wk = malloc(sI*sizeof(double));
  zk = malloc(sI*sizeof(double));
  t1  = malloc(sI*sizeof(double));
  t2  = malloc(sIc*sizeof(double));
  yic = malloc(sIc*sizeof(double));
  uij = malloc(2*sizeof(double));

  for (i=0; i<mn; i++){
		x[i] = y[i];
	}

  for (i=0; i<sIc; i++){
		yic[i] = y[Ic[i]];
	}
	
	for (i=0; i<sI; i++){
		wk[i] = 0.0;
	}

  for (k=0; k<maxiter; k++) {

    /* step 1 */
    for (i=0; i<mn; i++) df[i] = 0.0;
    
    pobj = 0.0;
    for (j=0; j<n-1; j++) {
      for (i=0; i<m-1; i++) {
				i1 = (i+1) + j*m;
				i2 = i + (j+1)*m;
				i3= i+j*m;

				uij[0] = x[i1]-x[i3];
				uij[1] = x[i2]-x[i3];
			
				c1 = sqrt(uij[0]*uij[0] + uij[1]*uij[1]);
				pobj += c1; 

				c2 = MAX(mu, c1);
				uij[0] = uij[0]/c2;
				uij[1] = uij[1]/c2;

				df[i1] += uij[0];
				df[i3] -= uij[0];
				df[i2] += uij[1];
				df[i3] -= uij[1];
      }
    }
    
    /* #dobj = -gamma*nrm2(I*D.T*u) + (u.T*D*Ic.T*Ic*b)+ (u.T*D*I.T*d)
			 dobj = -gamma*nrm2(df[I]) + dot(df[Ic],b[Ic]) +dot(df[I],d)*/

		dobj = 0;
		for (i=0; i<sI; i++){			
			t1[i] = df[Ii[i]];
			dobj += t1[i];
		}
		

		dobj=dobj*d;

		for (i=0; i<sIc; i++)
			t2[i] = df[Ic[i]];

    dobj += ddot_(&sIc, t2, &one, yic, &one) - gamma*dnrm2_(&sI, t1, &one);

    if (pobj - dobj < eps || STOPMARK) 
      goto cleanup;
    
    /* step 2 */
		/* t1 = L*x_k[I] - df[I] -L*d
       y_k[I]  = min(1/L, gamma/nrm2(t1))*t1 + d*/

    for (i=0; i<sI; i++) 			
			t1[i] = L*x[Ii[i]] - df[Ii[i]]-Ld;

    c1 = 1/MAX(L, dnrm2_(&sI, t1, &one)/gamma);

		for (i=0; i<sI; i++)
			yk[i] = c1*t1[i]+d;


    /* step 3 */
    /* w_k += (k+1)/2.0*df_I */
    c1 = (1.0+k)/2.0;

		for (i=0; i<sI; i++)
			t1[i] = df[Ii[i]];

    daxpy_(&sI, &c1, t1, &one, wk, &one);
		
		/* t1 = -w_k[I]
       z_k[I]  = min(1/L, gamma/nrm2(t1))*t1+d*/
		
		for (i=0; i<sI; i++) 			
			t1[i] = -wk[i];

    c1 = 1/MAX(L, dnrm2_(&sI, t1, &one)/gamma);

		for (i=0; i<sI; i++) 			
			zk[i] = c1*t1[i]+d;

    /* step 4 */
    /* x_k = z_k*2/(k+3) + y_k*(k+1)/(k+3) */
		alpha_kp1 = (k+2)/2.0;
		A_kp1 += alpha_kp1;

		t_k = alpha_kp1/A_kp1;
		
		for (i=0; i<sI; i++) t1[i] = t_k*zk[i];
    
		m1t_k = 1-t_k;
    daxpy_(&sI, &m1t_k, yk, &one, t1, &one);
		
		/*Update the restored pixels*/
		for (i=0; i<sI; i++) x[Ii[i]] = t1[i];

  }
  
 cleanup:
  free(df);
  free(yk);
  free(wk);
  free(zk);
  free(t1);
  free(t2);
  free(yic);

  kf[0] = (double)(k);
  epsilon_kf[0] = pobj-dobj;

}
