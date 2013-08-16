#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R_ext/Utils.h>
#include <R_ext/Print.h>


struct outputstruct  
    {  
    double logLmax;
	double *chain_warp;
	double *chain_amp;
	double *chain_ker;
	double *chain_shift;
	double *chain_logsigma;
	int indexLmax;
    };  

void R_CheckUserInterrupt(void);
double runif();
void warp3(double *pA,double *pL,double *pR1,double *pR2,double *x, int size,int S, int comps, double ll, double ul);
double normaldens(double t, double sigma, double mu), rtnormal(double mu,double sig,double lower,double upper, double *pnormp);
double denstnormal(double evalt,double mu,double sig,double lower,double upper, double *pnormp);
double kernel(double *kernelpars, double t,int Namp);
double mind(double xx,double yy), maxd(double xx,double yy);
double likely(double *warpp, double *amp,double *kernelpars, int Namp, double *prep, double *X,double *Y, int Q, int S, int N,double *FreqA,double *FreqLl,double *FreqUl,double *FreqL,double *BreakA,double *BreakLl,double *BreakUl,double *BreakL,double ll, double ul,double thresh,double *WX,double *WYinterp,double *R1,double *R2,double *LambdaS, double *Lwarp,double *pnormp);
int imax( int n1, int n2), imin( int n1, int n2);

void amcmc(double *X,double *Y,int Q,int chainlength, int thinning, int burnin, int S, int N, double thresh,double *initvalsigma,double *initvalshift, double *initvalarg,double *initvalpvamp,double *initvalamp,double *initvalker, double *initvalpvker,int Namp,
			  double *FreqA, double *FreqLl,double *FreqUl,double *FreqL,double *BreakA,double *BreakLl,double *BreakUl,double *BreakL,double ll,double ul,double lla,double ula,double *pnormp,double *Chainwarp, double *Chainshift,double *Chainamp, double *Chainker, double *Chainlogsigma, int *iternr,double *WX,double *WYinterp,double *R1,double *R2,double *LambdaS, double *Lwarp,double *sumL);

/*-------------------------MAIN------------------------------ */

void interfaceRC(double *X,double *Y,int *Sp,int *Np, int *chain, int *thin, int *burnin, int *Nampp, double *initvalker, int *Q,double *thresh, double *priorsp, double *pnormp, double *Chainwarp, double *Chainshift,double *Chainamp, double *Chainker, double *Chainlogsigma, int *iternr)
{
double *FreqA, *FreqLl,*FreqUl,*FreqL,*BreakA,*BreakLl,*BreakUl,*BreakL;
int i,k,Namp=*Nampp,N=*Np,S=*Sp;
double *initvalamp, *initvalpvamp,*initvalpvker,*initvalshift, ul,ll,ula,lla;
double *initvaluesw,*initvaluespv; // corresponding pointers
double *WX,*WYinterp,*R1,*R2,*LambdaS,*Lwarp, *sumL;
// struct outputstruct outputmodel;


// memory allocation

initvalamp = (double *) R_alloc(Namp, sizeof(double));
initvalshift = (double *) R_alloc(S-1, sizeof(double));
initvalpvamp = (double *) R_alloc(Namp, sizeof(double));
initvalpvker = (double *) R_alloc(Namp * 3, sizeof(double) );

initvaluesw = (double *) R_alloc(*Q * 3 * (S-1)* *Q +1, sizeof(double));
initvaluespv = (double *) R_alloc(*Q * 3 * (S-1)* *Q +1 , sizeof(double));


WX = (double *) R_alloc(S * N, sizeof(double));
WYinterp = (double *) R_alloc( S * N, sizeof(double));
R1 = (double *) R_alloc(*Q, sizeof(double));
R2 = (double *) R_alloc(*Q, sizeof(double));
LambdaS = (double *) R_alloc(*Q, sizeof(double));
Lwarp = (double *) R_alloc(*Q, sizeof(double));
sumL = (double *) R_alloc(*Q, sizeof(double));


// Respons domain limits
// ***********************************
ula =*Y;
lla =*Y; 
for (i=0;i<(S-1);i++)
	{
	for (k=0; k<N; k++) 
		{
		if (*(Y+i*N + k) > ula){ula = *(Y+i*N + k);}
		if (*(Y+i*N + k) < lla){lla = *(Y+i*N + k);}
		}
	}



// warping domain limits
// ***********************************
ul =*(X+N-1);
ll =*(X); 
for (i=1;i<(S-1);i++)
	{
	if (*(X+i*N + N-1) > ul){ul = *(X+i*S + N-1);}
	if (*(X+i*N ) < ll){ll = *(X+i*S );}
	}

// starting values 
// ***********************************

// for the new warping component (+shift)
*(initvaluesw + *Q-1) = (ul+ll)/2; //a
*(initvaluesw+ 2* *Q -1) = ll+(ul-ll)/4; //w_ll
*(initvaluesw+3* *Q-1) = ul-(ul-ll)/4; // w_ul
for (i=0;i<(S-1);i++)
	{
	*(initvaluesw+3* *Q + (S-1)*(*Q-1)+i)=1/10000; // lambda_i
	*(initvalshift +i)=0.0; //shift
	}

// for proposal density of the new warping component
*(initvaluespv + *Q-1)= (ul-ll)/8.0;
*(initvaluespv+2 * *Q -1)= (ul-ll)/8.0;
*(initvaluespv+3* *Q -1)= (ul-ll)/8.0;
for (i=0;i<(S-1);i++)
	{
	*(initvaluespv+3* *Q + (S-1)*(*Q-1)+i)=0.1;
	}

for (i=0;i<(S-1);i++)
	{
	*(initvaluesw+3* *Q + (S-1)*(*Q-1)+i)=1/10000; // lambda_i
	*(initvalshift +i)=0.0; //shift
	}

if (*Q==1)
	{
	// for proposal density of error variance
	*(initvaluespv+(S-1)+3)=(ula-lla)/10; 
	// for the amplitude part 
	for (i=0; i<Namp; i++)
		{
		*(initvalamp +i)=(ula-lla)/20; //sigma_k
		*(initvalpvamp +i)=(ula-lla)/20; //for proposal density (sigma_k) 
		*(initvalpvker +i*3) = (ul-ll)/20;
		*(initvalpvker +i*3 +1) = (ul-ll)/20;
		*(initvalpvker +i*3 +2) = (ul-ll)/20;
		}
	}
else
	{
	// for proposal density of error variance
	*(initvaluespv+3* *Q + *Q *(S-1))=*(Chainlogsigma+(*Q-1)*3 + (S-1)* (*Q-1)); 
	*(initvaluesw+3* *Q + *Q *(S-1))=*(Chainwarp+(*Q-1)*3 + (S-1)* (*Q-1)); 
	// for the amplitude part
	for (i=0; i<Namp; i++)
		{
		*(initvalamp +i)=*(Chainamp+i); //sigma_k
		*(initvalpvamp +i)=(ula-lla)/20; //for proposal density (sigma_k) 
		*(initvalpvker +i*3) = (ul-ll)/20;
		*(initvalpvker +i*3 +1) = (ul-ll)/20;
		*(initvalpvker +i*3 +2) = (ul-ll)/20;
		}
	// for proposal density of warping parameters & warping parameters
	for (k=0;k<(*Q-1);k++)
		{
		*(initvaluespv + k) = *(Chainlogsigma+k); //a
		*(initvaluespv+ *Q +k) = *(Chainlogsigma+ (*Q-1) +k); //w_ll
		*(initvaluespv+2* *Q +k) = *(Chainlogsigma+ 2* (*Q-1) +k); // w_ul

		*(initvaluesw + k) = *(Chainwarp+k); //a
		*(initvaluesw+ *Q +k) = *(Chainwarp+ (*Q-1)+k); //w_ll
		*(initvaluesw+2* *Q +k) = *(Chainwarp+ 2* (*Q-1) +k); // w_ul

		for (i=0;i<(S-1);i++)
			{
			*(initvaluespv+3* *Q + (S-1)*(k)+i)=*(Chainlogsigma+ 3* (*Q-1) + (S-1)*(k)+i); // lambda_i
			*(initvaluesw+3* *Q + (S-1)*(k)+i)=*(Chainwarp+ 3* (*Q-1) + (S-1)*(k)+i);; // lambda_i
			} 
		}
	for (i=0;i<(S-1);i++)
		{
		*(initvalshift +i)=*(Chainshift+i); // shift
		} 
	}

	for (i=0;i<(*Q);i++)
		{
		*(sumL+i)=0.0;
		}		


	// the priors
	FreqA = priorsp;
	FreqLl = (priorsp+14* *Q);
	FreqUl = (priorsp+14* *Q*2);
	FreqL = (priorsp+14* *Q*3);
	BreakA = (priorsp+14* *Q*(3+S-1));
	BreakLl = (priorsp+14* *Q*(3+S-1) + 15* *Q);
	BreakUl = (priorsp+14* *Q*(3+S-1)+ 15* *Q*2);
	BreakL = (priorsp+14* *Q*(3+S-1)+ 15* *Q*3);

	// model estimation
	amcmc(X,Y,*Q,*chain,*thin, *burnin,S,N,*thresh,initvaluespv, initvalshift,initvaluesw,initvalpvamp,initvalamp,initvalker,initvalpvker, Namp,  FreqA, FreqLl,FreqUl,FreqL,BreakA,BreakLl,BreakUl,BreakL,ll,ul,lla,ula,pnormp,Chainwarp,Chainshift,Chainamp,Chainker, Chainlogsigma,iternr,WX,WYinterp,R1,R2,LambdaS,Lwarp,sumL);

}


/*----------------------------------------------------------- */


void amcmc(double *X,double *Y,int Q,int chainlength, int thinning,int burnin,int S,int N,double thresh, double *initvalsigma, double *initvalshift, double *initvalarg,double *initvalpvamp,double *initvalamp,double *initvalker,double *initvalpvker, int Namp,
    double *FreqA, double *FreqLl,double *FreqUl,double *FreqL,double *BreakA,double *BreakLl,double *BreakUl,double *BreakL, double ll, double ul,double lla, double ula, double *pnormp,double *Chainwarp, double *Chainshift,double *Chainamp, double *Chainker, double *Chainlogsigma, int *iternr,double *WX,double *WYinterp,double *R1,double *R2,double *LambdaS, double *Lwarp,double *sumL)
{
  double kerp= (ul-ll)/100; // minimum distance between kernel parameters
  double *logsigma;
  unsigned short int  *acount, *acountamp, *acountker, thinnum, accepted,acounto=0,u,u2, K,help;
  unsigned int i,ii,jj,jj2,kk,count_c=0;
  double *warpp, *newp, propold, propnew,gammai[2],*gamma=gammai;
  double wshifts,*wshift,*newwshift;
  double logalpha, llw, ulw,lls,uls;
  double oldlikely, newlikely, propold1,propnew1,propold2,propnew2,propold3,propnew3;
  double *amp, *newamp, *ampsigma,*kernelpars, *kersigma,*newkernelpars, logLmax;

// ADDITIONAL SETTINGS
//---------------------
u= 30; //after how many thinnings to update the proposal 
// variances (u*thinning is the thinlength)
u2= 200; 
wshifts = (ul-ll)/100; //proposal variance of horizontal shift
 
llw = ll - 1.0/10.0 *(ul-ll);
ulw = ul + 1.0/10.0 *(ul-ll);
lls = -(ul-ll)/4.0;
uls = -lls;
*gamma=0.01;
*(gamma+1)=0.01;  

K = (S-1) * Q +3* Q +1;
  
    /* INITIALISATIONS. */


    warpp= (double *) R_alloc(K, sizeof(double) );
	newp= (double *) R_alloc(K, sizeof(double) );
	logsigma= (double *) R_alloc(K, sizeof(double) );
	newamp= (double *) R_alloc(Namp, sizeof(double));
	ampsigma=(double *) R_alloc(Namp, sizeof(double));
	amp=(double *) R_alloc(Namp, sizeof(double));
	kernelpars=(double *) R_alloc(Namp*3, sizeof(double) );
	newkernelpars=(double *) R_alloc(Namp*3, sizeof(double) );
	kersigma=(double *) R_alloc((S-1), sizeof(double) );
	wshift = (double *) R_alloc( (S), sizeof(double));
	newwshift = (double *) R_alloc((S), sizeof(double) );

	acount=(unsigned short int *) R_alloc(((S-1) +3), sizeof(unsigned short int));
	acountamp = (unsigned short int *) R_alloc(Namp, sizeof(unsigned short int));
	acountker = (unsigned short int *) R_alloc(Namp, sizeof(unsigned short int));

	for (i=0;i<(S);i++)
		{
		*(wshift+i) = 0;
		*(newwshift+i) = 0;
		}

	for (i=0;i<(3+S);i++)
		{
		acount[i]=0;
		}
	for (i=0;i<K;i++)
		{
		*(warpp+i) = *(initvalarg+i);
		*(logsigma+i) = *(initvalsigma+i);
		*(newp+i) = *(initvalarg+i);
		}
	for (i=0;i<Namp;i++ )
		{
		*(amp+i)=*(initvalamp+i);
		*(ampsigma+i) = *(initvalpvamp+i);
		*(newamp+i) = *(initvalamp+i);

		*(kernelpars+i*3)=*(initvalker+3*i);
		*(kernelpars+i*3+1)=*(initvalker+3*i+1);
		*(kernelpars+i*3+2)=*(initvalker+3*i+2);
		*(newkernelpars+i*3)=*(initvalker+3*i);
		*(newkernelpars+i*3+1)=*(initvalker+3*i+1);
		*(newkernelpars+i*3+2)=*(initvalker+3*i+2);
		*(kersigma+i*3)= *(initvalpvker+3*i);
		*(kersigma+i*3+1)= *(initvalpvker+3*i+1);
		*(kersigma+i*3+2)= *(initvalpvker+3*i+2);		
		acountamp[i]=0;
		acountker[i]=0;
		}
	
	for (jj=0;jj<(S-1);jj++)
			{
			for (i=0;i<Q;i++)
			{
			*(sumL+i)=*(sumL+i)+*(warpp+3*Q+i*(S-1)+jj);
			}
			*(wshift+jj) = *(initvalshift+jj);;
			*(newwshift+jj) = *(initvalshift+jj);;
			}
 

	oldlikely = likely(warpp,amp,kernelpars,Namp,newwshift,X,Y,Q,S,N, FreqA, FreqLl, FreqUl, FreqL, BreakA, BreakLl, BreakUl,BreakL,ll,ul,thresh,WX,WYinterp,R1,R2,LambdaS,Lwarp,pnormp);
    logLmax = oldlikely;
	*iternr = 1;

    /* MAIN ITERATIVE LOOP. */
    for (thinnum=1; thinnum<=chainlength; thinnum++) 
		{
		R_CheckUserInterrupt();

		for (ii=1; ii<=thinning; ii++) 
			{
			
			//shift horizontal
			for (jj=0;jj<(S-1);jj++)
				{
				*(newwshift+jj) = rtnormal(*(wshift+jj),wshifts,lls,uls, pnormp); 
				propold=1;
				propnew=1;
				propnew=propnew*denstnormal(*(wshift+jj),*(newwshift+jj),wshifts,lls,uls, pnormp);
				propold=propold*denstnormal(*(newwshift+jj),*(wshift+jj),wshifts,lls,uls, pnormp);
				newlikely =  likely(newp,amp,kernelpars,Namp,newwshift,X, Y, Q,S,N, FreqA, FreqLl, FreqUl, FreqL, BreakA, BreakLl, BreakUl,BreakL,ll,ul,thresh,WX,WYinterp,R1,R2,LambdaS,Lwarp,pnormp);
				if (newlikely == 0)
					{accepted = 0;}
				else
					{
					logalpha = newlikely - oldlikely+log(propnew)-log(propold);
					accepted = ( log(runif()) < logalpha );
					}
				if (accepted) 
					{
					*(wshift+jj) = *(newwshift+jj);
					} 
				else{
					*(newwshift+jj)= *(wshift+jj);
					}
				}

			/* For the new component. */
			kk=(Q-1);
			/* PROPOSE MOVE OF 'A' COORD. */
				propold=1;
				propnew=1;
				*(newp+kk) = rtnormal(*(warpp+kk),*(logsigma+kk),*(warpp+Q+kk),*(warpp+2*Q+kk), pnormp); 
				while ( *(newp+kk) ==  *(warpp+kk) || *(newp+kk) == *(warpp+Q+kk) || *(newp+kk) == *(warpp+2*Q+kk))
					{*(newp+kk) = rtnormal(*(warpp+kk),*(logsigma+kk),*(warpp+Q+kk),*(warpp+2*Q+kk), pnormp); } 
				propnew= denstnormal(*(warpp+kk),*(newp+kk),*(logsigma+kk),*(warpp+Q+kk),*(warpp+2*Q+kk), pnormp);
				propold= denstnormal(*(newp+kk),*(warpp+kk),*(logsigma+kk),*(warpp+Q+kk),*(warpp+2*Q+kk), pnormp);
				newlikely =  likely(newp,amp,kernelpars,Namp,newwshift,X, Y, Q,S,N, FreqA, FreqLl, FreqUl, FreqL, BreakA, BreakLl, BreakUl,BreakL,ll,ul,thresh,WX,WYinterp,R1,R2,LambdaS,Lwarp,pnormp);
           		if (newlikely == 0)
					{accepted = 0;}
				else
					{
					logalpha = newlikely - oldlikely+log(propnew)-log(propold);
					accepted = ( log(runif()) < logalpha );
					}
				if (accepted) 
					{
					acount[0] = acount[0] + 1;
					oldlikely = newlikely;
					*(warpp+kk) = *(newp +kk);
					} 
				else {*(newp+kk)=*(warpp+kk);}


				/* PROPOSE MOVE OF 'Ll' COORD. */
				propold=1;
				propnew=1;
				*(newp+Q+kk) = rtnormal(*(warpp+Q+kk),*(logsigma+Q+kk),llw,*(warpp+kk), pnormp); 
				while (*(newp+Q+kk) ==  *(warpp+Q+kk)|| *(newp+Q+kk) == llw || *(newp+Q+kk) == *(warpp+kk))
					{*(newp+Q+kk) = rtnormal(*(warpp+Q+kk),*(logsigma+Q+kk),llw,*(warpp+kk), pnormp);}
				propnew= propnew*denstnormal(*(warpp+Q+kk),*(newp+Q+kk),*(logsigma+Q+kk),llw,*(warpp+kk), pnormp);
				propold= propold*denstnormal(*(newp+Q+kk),*(warpp+Q+kk),*(logsigma+Q+kk),llw,*(warpp+kk), pnormp);
				newlikely = likely(newp,amp,kernelpars,Namp,newwshift,X, Y, Q,S,N, FreqA, FreqLl, FreqUl, FreqL, BreakA, BreakLl, BreakUl,BreakL,ll,ul,thresh,WX,WYinterp,R1,R2,LambdaS,Lwarp,pnormp);
           		if (newlikely == 0)
					{accepted = 0;}
				else
					{
					logalpha = newlikely - oldlikely+log(propnew)-log(propold);
					accepted = ( log(runif()) < logalpha );
					}
				if (accepted) 
					{
					acount[1] = acount[1] + 1;
					oldlikely = newlikely;
					*(warpp+Q+kk) = *(newp+Q+kk);
					} 
				else {*(newp+Q+kk)=*(warpp+Q+kk);}


				/* PROPOSE MOVE OF 'Ul' COORD. */
				propold=1;
				propnew=1;
				*(newp+2*Q+kk) = rtnormal(*(warpp+2*Q+kk),*(logsigma+2*Q+kk),*(warpp+kk),ulw,pnormp); 
				while ( *(newp+2*Q+kk) == *(warpp+2*Q+kk)||  *(newp+2*Q+kk) == *(warpp+kk)||  *(newp+2*Q+kk)==ulw)
					{*(newp+2*Q+kk) = rtnormal(*(warpp+2*Q+kk),*(logsigma+2*Q+kk),*(warpp+kk),ulw, pnormp); }
				propnew= propnew*denstnormal(*(warpp+2*Q+kk),*(newp+2*Q+kk),*(logsigma+2*Q+kk),*(warpp+kk),ulw, pnormp);
				propold= propold*denstnormal(*(newp+2*Q+kk),*(warpp+2*Q+kk),*(logsigma+2*Q+kk),*(warpp+kk),ulw, pnormp);
				newlikely =  likely(newp,amp,kernelpars,Namp,newwshift,X, Y, Q,S,N, FreqA, FreqLl, FreqUl, FreqL, BreakA, BreakLl, BreakUl,BreakL,ll,ul,thresh,WX,WYinterp,R1,R2,LambdaS,Lwarp,pnormp);
           		if (newlikely == 0)
					{accepted = 0;}
				else
					{
					logalpha = newlikely - oldlikely+log(propnew)-log(propold);
					accepted = ( log(runif()) < logalpha );
					}
				if (accepted) 
					{
					acount[2] = acount[2] + 1;
					oldlikely = newlikely;
					*(warpp+2*Q+kk) = *(newp+2*Q+kk);
					} 
				else {*(newp+2*Q+kk)= *(warpp+2*Q+kk);}
				/* PROPOSE MOVE OF 'Lambda' COORD. */
				for (jj=0;jj<(S-1);jj++)
					{
					*(sumL+kk)= *(sumL+kk)-*(warpp+3*Q+kk*(S-1)+jj);
					propold=1;
					propnew=1;
					*(newp+3*Q+(kk)*(S-1)+jj) = rtnormal(*(warpp+3*Q+(kk)*(S-1)+jj),*(logsigma+3*Q+(kk)*(S-1)+jj),maxd(-1,-1-*(sumL+kk)),mind(1,1-*(sumL+kk)), pnormp); 
					while (  *(newp+3*Q+(kk)*(S-1)+jj) ==  *(warpp+3*Q+(kk)*(S-1)+jj)  || *(newp+3*Q+(kk)*(S-1)+jj) == mind(1,1-*(sumL+kk)) ||  *(newp+3*Q+(kk)*(S-1)+jj) == maxd(-1,-1-*(sumL+kk)))
						{
						*(newp+3*Q+(kk)*(S-1)+jj) = rtnormal(*(warpp+3*Q+(kk)*(S-1)+jj),*(logsigma+3*Q+(kk)*(S-1)+jj),maxd(-1,-1-*(sumL+kk)),mind(1,1-*(sumL+kk)), pnormp); 
						}
					propnew=propnew*denstnormal(*(warpp+3*Q+(kk)*(S-1)+jj),*(newp+3*Q+(kk)*(S-1)+jj),*(logsigma+3*Q+(kk)*(S-1)+jj),maxd(-1,-1-*(sumL+kk)),mind(1,1-*(sumL+kk)), pnormp);
					propold=propold*denstnormal(*(newp+3*Q+(kk)*(S-1)+jj),*(warpp+3*Q+(kk)*(S-1)+jj),*(logsigma+3*Q+(kk)*(S-1)+jj),maxd(-1,-1-*(sumL+kk)),mind(1,1-*(sumL+kk)), pnormp);
					newlikely =  likely(newp,amp,kernelpars,Namp,newwshift,X, Y, Q,S,N, FreqA, FreqLl, FreqUl, FreqL, BreakA, BreakLl, BreakUl,BreakL,ll,ul,thresh,WX,WYinterp,R1,R2,LambdaS,Lwarp,pnormp);
           			if (newlikely == 0)
						{accepted = 0;}
					else
						{
						logalpha = newlikely - oldlikely+log(propnew)-log(propold);
						accepted = ( log(runif()) < logalpha );
						}
					if (accepted) 
						{
						acount[3+jj] = acount[3+jj] + 1;
						oldlikely = newlikely;
						*(warpp+3*Q+kk*(S-1)+jj) = *(newp+3*Q+kk*(S-1)+jj);
						} 
					else{*(newp+3*Q+kk*(S-1)+jj)= *(warpp+3*Q+kk*(S-1)+jj);}
					*(sumL+kk)=*(sumL+kk)+*(warpp+3*Q+kk*(S-1)+jj);

					}

			/* For the old components. */
			if (Q>1)
			{

            for (kk=0; kk<(Q-1); kk++) 
				{
					/* PROPOSE MOVE OF THE COORDINATES */
				propold=1;
				propnew=1;

				*(newp+kk) = rtnormal(*(warpp+kk),*(logsigma+kk),maxd(*(warpp+Q+kk),*(BreakA+kk)),mind(*(warpp+2*Q+kk),*(BreakA+14*Q+kk)), pnormp); 
				// while ((float) *(newp+kk) == (float) *(warpp+kk) || *(newp+kk) == maxd(*(warpp+Q+kk),*(BreakA+kk)) || *(newp+kk) == mind(*(warpp+2*Q+kk),*(BreakA+14*Q+kk)) )
				//	{*(newp+kk) = rtnormal(*(warpp+kk),*(logsigma+kk),maxd(*(warpp+Q+kk),*(BreakA+kk)),mind(*(warpp+2*Q+kk),*(BreakA+14*Q+kk)), pnormp); }
				propnew= denstnormal(*(warpp+kk),*(newp+kk),*(logsigma+kk),maxd(*(warpp+Q+kk),*(BreakA+kk)),mind(*(warpp+2*Q+kk),*(BreakA+14*Q+kk)), pnormp);
				propold= denstnormal(*(newp+kk),*(warpp+kk),*(logsigma+kk),maxd(*(warpp+Q+kk),*(BreakA+kk)),mind(*(warpp+2*Q+kk),*(BreakA+14*Q+kk)), pnormp);

				*(newp+Q+kk) = rtnormal(*(warpp+Q+kk),*(logsigma+Q+kk),*(BreakLl+kk),mind(*(newp+kk),*(BreakLl+14*Q+kk)), pnormp); 
				//while ((float) *(newp+Q+kk) == (float) *(warpp+Q+kk) || *(newp+Q+kk) == *(BreakLl+kk) || *(newp+Q+kk) == mind(*(newp+kk),*(BreakLl+14*Q+kk)))
				//	{*(newp+Q+kk) = rtnormal(*(warpp+Q+kk),*(logsigma+Q+kk),*(BreakLl+kk),mind(*(newp+kk),*(BreakLl+14*Q+kk)), pnormp); }
				propnew= propnew*denstnormal(*(warpp+Q+kk),*(newp+Q+kk),*(logsigma+Q+kk),*(BreakLl+kk),mind(*(newp+kk),*(BreakLl+14*Q+kk)), pnormp);
				propold= propold*denstnormal(*(newp+Q+kk),*(warpp+Q+kk),*(logsigma+Q+kk),*(BreakLl+kk),mind(*(newp+kk),*(BreakLl+14*Q+kk)), pnormp);
				*(newp+2*Q+kk) = rtnormal(*(warpp+2*Q+kk),*(logsigma+2*Q+kk),maxd(*(newp+kk),*(BreakUl+kk)),*(BreakUl+14*Q+kk), pnormp); 
			
				//while ((float) *(newp+2*Q+kk) == (float) *(warpp+2*Q+kk) || *(newp+2*Q+kk) == maxd(*(newp+kk),*(BreakUl+kk)) || *(newp+2*Q+kk) == *(BreakUl+14*Q+kk))
				//	{*(newp+2*Q+kk) = rtnormal(*(warpp+2*Q+kk),*(logsigma+2*Q+kk),maxd(*(newp+kk),*(BreakUl+kk)),*(BreakUl+14*Q+kk), pnormp); }
				propnew= propnew*denstnormal(*(warpp+2*Q+kk),*(newp+2*Q+kk),*(logsigma+2*Q+kk),maxd(*(newp+kk),*(BreakUl+kk)),*(BreakUl+14*Q+kk), pnormp);
				propold= propold*denstnormal(*(newp+2*Q+kk),*(warpp+2*Q+kk),*(logsigma+2*Q+kk),maxd(*(newp+kk),*(BreakUl+kk)),*(BreakUl+14*Q+kk), pnormp);
				


				for (jj=0;jj<(S-1);jj++)
					{
					*(sumL+kk)=*(sumL+kk)-*(warpp+3*Q+kk*(S-1)+jj);

					propold=1;
					propnew=1;
					*(newp+3*Q+(kk)*(S-1)+jj) = rtnormal(*(warpp+3*Q+(kk)*(S-1)+jj),*(logsigma+3*Q+(kk)*(S-1)+jj),maxd(*(BreakL+kk*(S-1)+jj),-1-*(sumL+kk)),mind(*(BreakL+14*Q*(S-1)+kk*(S-1)+jj),1-*(sumL+kk)), pnormp); 
					while (  *(newp+3*Q+(kk)*(S-1)+jj) ==  *(warpp+3*Q+(kk)*(S-1)+jj)  || *(newp+3*Q+(kk)*(S-1)+jj) == mind(*(BreakL+14*Q*(S-1)+kk*(S-1)+jj),1-*(sumL+kk)) ||  *(newp+3*Q+(kk)*(S-1)+jj) == maxd(*(BreakL+kk*(S-1)+jj),-1-*(sumL+kk)))
						{
						*(newp+3*Q+(kk)*(S-1)+jj) = rtnormal(*(warpp+3*Q+(kk)*(S-1)+jj),*(logsigma+3*Q+(kk)*(S-1)+jj),maxd(*(BreakL+kk*(S-1)+jj),-1-*(sumL+kk)),mind(*(BreakL+14*Q*(S-1)+kk*(S-1)+jj),1-*(sumL+kk)), pnormp); 
						}
					//*(newp+3*Q+(kk)*(S-1)+jj) = rtnormal(*(warpp+3*Q+(kk)*(S-1)+jj),*(logsigma+3*Q+(kk)*(S-1)+jj),*(BreakL+kk*(S-1)+jj),*(BreakL+14*Q*(S-1)+kk*(S-1)+jj), pnormp); 
					//while ( (float) *(newp+3*Q+(kk)*(S-1)+jj) ==  (float) *(warpp+3*Q+(kk)*(S-1)+jj) || *(newp+3*Q+(kk)*(S-1)+jj) ==  *(BreakL+kk*(S-1)+jj) || *(newp+3*Q+(kk)*(S-1)+jj) == *(BreakL+14*Q*(S-1)+kk*(S-1)+jj))
					//	{	*(newp+3*Q+(kk)*(S-1)+jj) = rtnormal(*(warpp+3*Q+(kk)*(S-1)+jj),*(logsigma+3*Q+(kk)*(S-1)+jj),*(BreakL+kk*(S-1)+jj),*(BreakL+14*Q*(S-1)+kk*(S-1)+jj), pnormp); }
					propnew=propnew*denstnormal(*(warpp+3*Q+(kk)*(S-1)+jj),*(newp+3*Q+(kk)*(S-1)+jj),*(logsigma+3*Q+(kk)*(S-1)+jj),maxd(*(BreakL+kk*(S-1)+jj),-1-*(sumL+kk)),mind(*(BreakL+14*Q*(S-1)+kk*(S-1)+jj),1-*(sumL+kk)), pnormp);
					propold=propold*denstnormal(*(newp+3*Q+(kk)*(S-1)+jj),*(warpp+3*Q+(kk)*(S-1)+jj),*(logsigma+3*Q+(kk)*(S-1)+jj),maxd(*(BreakL+kk*(S-1)+jj),-1-*(sumL+kk)),mind(*(BreakL+14*Q*(S-1)+kk*(S-1)+jj),1-*(sumL+kk)), pnormp);
					*(sumL+kk)=*(sumL+kk)+*(newp+3*Q+(kk)*(S-1)+jj);
					}
				newlikely =  likely(newp,amp,kernelpars,Namp,newwshift,X, Y, Q,S,N, FreqA, FreqLl, FreqUl, FreqL, BreakA, BreakLl, BreakUl,BreakL,ll,ul,thresh,WX,WYinterp,R1,R2,LambdaS,Lwarp,pnormp);
           		if (newlikely == 0)
					{accepted = 0;}
				else
					{
					logalpha = newlikely - oldlikely+log(propnew)-log(propold);
					accepted = ( log(runif()) < logalpha );
					}
				if (accepted) 
					{
					acounto = acounto + 1;
					oldlikely = newlikely;
					*(warpp+kk) = *(newp+kk);
					*(warpp+Q+kk) = *(newp+Q+kk);
					*(warpp+2*Q+kk) = *(newp+2*Q+kk);
					for (jj=0;jj<(S-1);jj++)
						{*(warpp+3*Q+kk*(S-1)+jj) = *(newp+3*Q+kk*(S-1)+jj);}
					} 
				else
					{
					for (jj=0;jj<(S-1);jj++)
						{
						*(newp+3*Q+kk*(S-1)+jj)= *(warpp+3*Q+kk*(S-1)+jj);
						*(sumL+kk)=*(sumL+kk)-*(newp+3*Q+(kk)*(S-1)+jj);
						*(sumL+kk)=*(sumL+kk)+*(warpp+3*Q+kk*(S-1)+jj);
						}
					*(newp+kk) = *(warpp+kk);
					*(newp+Q+kk) = *(warpp+Q+kk);
					*(newp+2*Q+kk) = *(warpp+2*Q+kk);
					} 

				} //end of kk loop (components)
			} //end of the if Q>1 loop

				// sigma
				propold=1;
				propnew=1;
				*(newp+K-1) = rtnormal(*(warpp+K-1),*(logsigma+K-1),0,ula-lla, pnormp); 
				while (  *(newp+K-1) == 0 || *(newp+K-1) == ula-lla)
					{*(newp+K-1) = rtnormal(*(warpp+K-1),*(logsigma+K-1),0,ula-lla, pnormp); }
				propnew= propnew*denstnormal(*(warpp+K-1),*(newp+K-1),*(logsigma+K-1),0,ula-lla, pnormp);
				propold= propold*denstnormal(*(newp+K-1),*(warpp+K-1),*(logsigma+K-1),0,ula-lla, pnormp);
				newlikely =  likely(newp,amp,kernelpars,Namp,newwshift,X, Y, Q,S,N, FreqA, FreqLl, FreqUl, FreqL, BreakA, BreakLl, BreakUl,BreakL,ll,ul,thresh,WX,WYinterp,R1,R2,LambdaS,Lwarp,pnormp);
           		if (newlikely == 0)
					{accepted = 0;}
				else
					{
					logalpha = newlikely - oldlikely+log(propnew)-log(propold);
					logalpha = logalpha + log(pow((*(warpp+K-1) / *(newp+K-1)), *(gamma)+1)) - *(gamma+1) * (1/ *(newp+K-1) - 1/ *(warpp+K-1));
					accepted = ( log(runif()) < logalpha );
					}
				if (accepted) 
					{
					acount[3+(S-1)] = acount[3+(S-1)] + 1;
					oldlikely = newlikely;
					*(warpp+K-1) = *(newp +K-1);
					} 
				else {*(newp+K-1)=*(warpp+K-1);}

				// kernel components
				for (kk=0; kk< (Namp);kk++)
					{
					propold1=1;
					propnew1=1;
					if (kk==0)
						{
						*(newkernelpars) = rtnormal(*(kernelpars),*(kersigma),ll-kerp,*(kernelpars+1)-kerp,pnormp); 
						while (*(newkernelpars) == ll-kerp || *(newkernelpars) == *(kernelpars+1)-kerp)
							{*(newkernelpars) = rtnormal(*(kernelpars),*(kersigma),ll-kerp,*(kernelpars+1)-kerp,pnormp); }
						propnew1= denstnormal(*(kernelpars),*(newkernelpars),*(kersigma),ll-kerp,*(kernelpars+1)-kerp, pnormp);
						propold1= denstnormal(*(newkernelpars),*(kernelpars),*(kersigma),ll-kerp,*(kernelpars+1)-kerp, pnormp);
						}
					else 
						{
						*(newkernelpars+kk*3) = rtnormal(*(kernelpars+kk*3),*(kersigma+kk*3),*(kernelpars+(kk-1)*3+2)+kerp,*(kernelpars+(kk)*3+1)-kerp, pnormp);
						while (*(newkernelpars+kk*3) == *(kernelpars+(kk-1)*3+2)+kerp || *(newkernelpars+kk*3) == *(kernelpars+(kk)*3+1)-kerp)
							{*(newkernelpars+kk*3) = rtnormal(*(kernelpars+kk*3),*(kersigma+kk*3),*(kernelpars+(kk-1)*3+2)+kerp,*(kernelpars+(kk)*3+1)-kerp, pnormp); }
						propnew1= denstnormal(*(kernelpars+kk*3),*(newkernelpars+kk*3),*(kersigma+kk*3),*(kernelpars+(kk-1)*3+2)+kerp,*(kernelpars+(kk)*3+1)-kerp, pnormp);
						propold1= denstnormal(*(newkernelpars +kk*3),*(kernelpars+kk*3),*(kersigma+kk*3),*(kernelpars+(kk-1)*3+2)+kerp,*(kernelpars+(kk)*3+1)-kerp, pnormp);
						}
					propold2=1;
					propnew2=1;
					*(newkernelpars+kk*3+1) = rtnormal(*(kernelpars+kk*3+1),*(kersigma+kk*3+1),*(newkernelpars+kk*3)+kerp,*(kernelpars+(kk)*3+2)-kerp, pnormp);
					
					while ( *(newkernelpars+kk*3+1) == *(newkernelpars+kk*3)+kerp || *(newkernelpars+kk*3+1) == *(kernelpars+(kk)*3+2)-kerp)
						{*(newkernelpars+kk*3+1) = rtnormal(*(kernelpars+kk*3+1),*(kersigma+kk*3+1),*(newkernelpars+kk*3)+kerp,*(kernelpars+(kk)*3+2)-kerp, pnormp); }
					propnew2= propnew1*denstnormal(*(kernelpars+kk*3+1),*(newkernelpars+kk*3+1),*(kersigma+kk*3+1),*(newkernelpars+kk*3)+kerp,*(kernelpars+(kk)*3+2)-kerp, pnormp);
					propold2= propold1*denstnormal(*(newkernelpars+kk*3+1),*(kernelpars+kk*3+1),*(kersigma+kk*3+1),*(newkernelpars+kk*3)+kerp,*(kernelpars+(kk)*3+2)-kerp, pnormp);
		
					propold3=1;
					propnew3=1;
					if (kk == (Namp-1))
						{
						*(newkernelpars+kk*3+2) = rtnormal(*(kernelpars+kk*3+2),*(kersigma+kk*3+2),*(newkernelpars+kk*3+1)+kerp,ulw+kerp, pnormp);
						while ( *(newkernelpars+kk*3+2) == *(newkernelpars+kk*3+1)+kerp || *(newkernelpars+kk*3+2) == ulw+kerp )
							{*(newkernelpars+kk*3+2) = rtnormal(*(kernelpars+kk*3+2),*(kersigma+kk*3+2),*(newkernelpars+kk*3+1)+kerp,ulw+kerp, pnormp); }
						propnew3= propnew2*denstnormal(*(kernelpars+kk*3+2),*(newkernelpars+kk*3+2),*(kersigma+kk*3+2),*(newkernelpars+kk*3+1)+kerp,ulw+kerp, pnormp);
						propold3= propold2*denstnormal(*(newkernelpars +kk*3+2),*(kernelpars+kk*3+2),*(kersigma+kk*3+2),*(newkernelpars+kk*3+1)+kerp,ulw+kerp, pnormp);
						}
					else
						{
						*(newkernelpars+kk*3+2) = rtnormal(*(kernelpars+kk*3+2),*(kersigma+kk*3+2),*(newkernelpars+kk*3+1)+kerp,*(kernelpars+(kk+1)*3)-kerp, pnormp);
						while ( *(newkernelpars+kk*3+2) == *(newkernelpars+kk*3+1)+kerp || *(newkernelpars+kk*3+2) == *(kernelpars+(kk+1)*3)-kerp )
							{*(newkernelpars+kk*3+2) = rtnormal(*(kernelpars+kk*3+2),*(kersigma+kk*3+2),*(newkernelpars+kk*3+1)+kerp,*(kernelpars+(kk+1)*3)-kerp, pnormp);}
						propnew3= propnew2*denstnormal(*(kernelpars+kk*3+2),*(newkernelpars+kk*3+2),*(kersigma+kk*3+2),*(newkernelpars+kk*3+1)+kerp,*(kernelpars+(kk+1)*3)-kerp, pnormp);
						propold3= propold2*denstnormal(*(newkernelpars+kk*3+2),*(kernelpars+kk*3+2),*(kersigma+kk*3+2),*(newkernelpars+kk*3+1)+kerp,*(kernelpars+(kk+1)*3)-kerp, pnormp);
						}
					newlikely = likely(warpp,amp,newkernelpars,Namp,newwshift, X, Y, Q,S,N, FreqA, FreqLl, FreqUl, FreqL, BreakA, BreakLl, BreakUl,BreakL,ll,ul,thresh,WX,WYinterp,R1,R2,LambdaS,Lwarp,pnormp);
					
					if (newlikely == 0)
						{accepted = 0;}
					else
						{
						logalpha = newlikely - oldlikely+log(propnew3*propnew2*propnew1)-log(propold3*propold2*propold1);
						accepted = ( log(runif()) < logalpha );
						}
					if (accepted) 
						{
						acountker[kk] = acountker[kk] + 1;
						oldlikely = newlikely;
						*(kernelpars +kk*3) = *(newkernelpars +kk*3);
						*(kernelpars +kk*3+1) = *(newkernelpars +kk*3+1);
						*(kernelpars +kk*3+2) = *(newkernelpars +kk*3+2);
						} 
					else	
						{
						*(newkernelpars+kk)=*(kernelpars+kk);
						*(newkernelpars+kk+1)=*(kernelpars+kk+1);
						*(newkernelpars+kk+2)=*(kernelpars+kk+2);
						}
					}
				// amplitude components
				for (kk=0; kk< Namp;kk++)
					{
					propold=1;
					propnew=1;
					*(newamp+kk) = rtnormal(*(amp+kk),*(ampsigma+kk),0,ula-lla, pnormp); 
					while ( (float) *(newamp+kk) == 0 || (float) *(newamp+kk) == ula-lla)
							{*(newamp+kk) = rtnormal(*(amp+kk),*(ampsigma+kk),0,ula-lla, pnormp);}
					propnew= propnew*denstnormal(*(amp+kk),*(newamp+kk),*(ampsigma+kk),0,ula-lla, pnormp);
					propold= propold*denstnormal(*(newamp +kk),*(amp+kk),*(ampsigma+kk),0,ula-lla, pnormp);
					newlikely = likely(warpp,newamp,kernelpars,Namp,newwshift, X, Y, Q,S,N, FreqA, FreqLl, FreqUl, FreqL, BreakA, BreakLl, BreakUl,BreakL,ll,ul,thresh,WX,WYinterp,R1,R2,LambdaS,Lwarp,pnormp);
           			if (newlikely == 0)
						{accepted = 0;}
					else
						{
						logalpha = newlikely - oldlikely+log(propnew)-log(propold);
						logalpha = logalpha + log(pow((*(amp+kk) / *(newamp+kk)), *(gamma)+1)) - ( *(gamma+1) * (1/ *(newamp+kk) - 1/ *(amp+kk)));
						accepted = ( log(runif()) < logalpha );
						}
					if (accepted) 
						{
						acountamp[kk] = acountamp[kk] + 1;
						oldlikely = newlikely;
						*(amp +kk) = *(newamp +kk);
						} 
					else {*(newamp+kk)=*(amp +kk );}
					}
				} //end of ii (thinloop)
				
				/* DO THE ADAPTING.(batches) */
				help = thinnum/u;
				if (thinnum < burnin && (double) thinnum/ ((double) u) == help )
					{

					if (acount[0] > u*thinning * 0.5 )
							{*(logsigma+Q-1) = mind(*(logsigma+Q-1) + *(logsigma+Q-1)/4,(ul-ll)/2.0);}
					else
						{
						if (acount[0] < u*thinning * 0.4)
								{*(logsigma+Q-1) = maxd((ul-ll)/(4.0*N),*(logsigma+Q-1) - *(logsigma+Q-1)/4);}
						}
					acount[0]=0;

					if (acount[1] > u*thinning * 0.5 )
							{*(logsigma+2*Q-1) = mind((ul-ll)/2.0,*(logsigma+2*Q-1) + *(logsigma+2*Q-1)/4);}
					else
						{
						if (acount[1] < u*thinning * 0.4)
								{*(logsigma+2*Q-1) = maxd((ul-ll)/(4.0*N),*(logsigma+2*Q-1) - *(logsigma+2*Q-1)/4);}
						}
					acount[1]=0;

					if (acount[2] > u*thinning * 0.5)
							{*(logsigma+3*Q-1) = mind((ul-ll)/2.0,*(logsigma+3*Q-1) + *(logsigma+3*Q-1)/4);}
					else
						{
						if (acount[2] < u*thinning * 0.4)
								{*(logsigma+3*Q-1) = maxd((ul-ll)/(4.0*N),*(logsigma+3*Q-1) - *(logsigma+3*Q-1)/4);}
						}
					acount[2]=0;

					for (jj=0; jj<(S-1); jj++) 
						{
						if (acount[3+jj] > u*thinning * 0.5)
							{*(logsigma+3*Q + (Q-1)*(S-1) +jj ) = mind(1,*(logsigma+3*Q + (Q-1)*(S-1) +jj ) + *(logsigma+3*Q + (Q-1)*(S-1) +jj )/4);}
						else
							{
							if (acount[3+jj] < u*thinning * 0.4)
								{*(logsigma+3*Q + (Q-1)*(S-1) +jj ) = maxd(0.0001,*(logsigma+3*Q + (Q-1)*(S-1) +jj ) - *(logsigma+3*Q + (Q-1)*(S-1) +jj )/4);}
							}
						acount[3+jj]=0;
						}

					for (jj=0; jj<Namp; jj++) //amplitude part
						{
						if (acountker[jj] > u*thinning * 0.5)
							{
							*(kersigma+jj*3) = mind((ul-ll)/(2.0),*(kersigma+jj*3) + *(kersigma+jj*3)/4);
							*(kersigma+jj*3+1) = mind((ul-ll)/(2.0),*(kersigma+jj*3+1) + *(kersigma+jj*3+1)/4);
							*(kersigma+jj*3+2) = mind((ul-ll)/(2.0),*(kersigma+jj*3+2) + *(kersigma+jj*3+2)/4);
							}
						else
							{
							if (acountker[jj] < u*thinning * 0.4)
								{
								*(kersigma+jj*3) = maxd((ul-ll)/(4.0*N),*(kersigma+jj*3) - *(kersigma+jj*3)/4);
								*(kersigma+jj*3+1) = maxd((ul-ll)/(4.0*N),*(kersigma+jj*3+1) - *(kersigma+jj*3+1)/4);
								*(kersigma+jj*3+2) = maxd((ul-ll)/(4.0*N),*(kersigma+jj*3+2) - *(kersigma+jj*3+2)/4);
								}
							}
						acountker[jj]=0;

						if (acountamp[jj] > u*thinning * 0.5)
							{*(ampsigma+jj) = mind((ula-lla)/(10.0),*(ampsigma+jj) + *(ampsigma+jj)/4);}
						else
							{
							if (acountamp[jj] < u*thinning * 0.4)
								{*(ampsigma+jj) = maxd((ula-lla)/(1000.0),*(ampsigma+jj) - *(ampsigma+jj)/4);}
							}
						acountamp[jj]=0;
						}
					if (acount[3+(S-1)] > u*thinning * 0.5)
							{*(logsigma+K-1) = mind((ula-lla)/(10.0),*(logsigma+K-1) + *(logsigma+K-1)/4);}
					else
						{
						if (acount[3+(S-1)] < u*thinning * 0.4)
								{*(logsigma+K-1) = maxd((ula-lla)/(1000.0),*(logsigma+K-1) - *(logsigma+K-1)/4);}
						}
					acount[3+(S-1)]=0;
					
					if (Q>1)
						{
						
						if (acounto< u*thinning * 0.4*(Q-1))
							{
							for (jj=0; jj<(Q-1); jj++) 
								{
								*(logsigma+jj) = maxd((ul-ll)/(4.0*N),*(logsigma+jj) - *(logsigma+jj)/4);
								*(logsigma+Q+jj) = maxd((ul-ll)/(4.0*N),*(logsigma+Q+jj) - *(logsigma+Q+jj)/4);
								*(logsigma+2*Q+jj) = maxd((ul-ll)/(4.0*N),*(logsigma+2*Q+jj) - *(logsigma+2*Q+jj)/4);
								for (jj2=0; jj2<(S-1); jj2++) 
									{
									*(logsigma+3*Q+(S-1)*jj+jj2) = maxd(0.0001,*(logsigma+3*Q+(S-1)*jj+jj2) - *(logsigma+3*Q+(S-1)*jj+jj2)/4);
									}
								}
							}
						acounto=0;
						}

					}

if (thinnum > burnin)
	{
	if (count_c == 0)
		{
		for (jj=0; jj<K; jj++) 
			{
			*(Chainlogsigma+jj)= *(logsigma+jj);
			}
		}

		for (jj=0; jj<K; jj++) 
			{
			*(Chainwarp + count_c*K + jj) = *(warpp+jj);
			}

		for (jj=0; jj<Namp; jj++) 
			{
			*(Chainamp + count_c*Namp +jj) = *(amp+jj);
			}

		for (jj=0; jj<(S-1); jj++) 
			{
			*(Chainshift+count_c*(S-1) +jj) = *(wshift+jj);
			}

		for (jj=0; jj<(3*Namp); jj++) 
			{
			*(Chainker+count_c*(3*Namp) +jj) = *(kernelpars+jj);
			}


		/*possibly update maximum loglikelihood value and the iteration of occurance */
		
		if (oldlikely > logLmax) 
			{
			logLmax = oldlikely;
			*iternr= thinnum-burnin;
			}
	count_c++;
	}



help = thinnum/u2;
if ((double) thinnum/ ((double) u2) == help )
{
Rprintf("progress: %d chains.\n", thinnum);
}


} /* End of thinnum for loop. */

 
} // end of amcmc



/*----------------------------------------------------------- */

double rtnormal(double mu,double sig,double lower,double upper, double *pnormp)
{
double pr1,pr2,Prob1,u1,u,v,x,y,rtnorm,r1,r2,p1d,p2d;
int p1,p2;
p1d = (upper-mu)/sig *(double)100.0 +(double)301.0;
p2d = (lower-mu)/sig *(double)100.0 +(double)301.0;
p1 = imin( (int) p1d,600);
p2 = imax( (int) p2d, 0);
pr2 = (double) *(pnormp+p1)-0.5;
pr1 = (double) 0.5-*(pnormp+p2);
Prob1 = pr1/(pr1+pr2);

r1=mu-lower;
r2=upper-mu;

u1 = runif();
if (u1 < Prob1)
	{
	do
		{
		u = runif()*(1-exp(-pow((r1/sig),2))) +exp(-pow((r1/sig),2));
		v = runif();
		x = sqrt(-2*log(u))*cos(2*3.1415926535*v);
		y = sqrt(-2*log(u))*sin(2*3.1415926535*v);
		if (x>0){x=-x;}
		}
	while (x > (r1/sig) || x < -(r1/sig) ||y > (r1/sig) ||y < -(r1/sig) );
	}
else
	{
	do
		{
		u = runif()*(1-exp(-pow((r2/sig),2))) +exp(-pow((r2/sig),2));
		v = runif();
		x = sqrt(-2*log(u))*cos(2*3.1415926535*v);
		y = sqrt(-2*log(u))*sin(2*3.1415926535*v);
		if (x<0){x=-x;}
		}
	while ( x > (r2/sig) || x < -(r2/sig) ||y > (r2/sig) ||y < -(r2/sig));
	}

rtnorm = x*sig+mu;
return(rtnorm);

}
/*----------------------------------------------------------- */

double denstnormal(double evalt,double mu,double sig,double lower,double upper, double *pnormp)
{
int p1,p2;
double pr1,pr2,denstn,p1d,p2d;

denstn=1/sig*(double)1.0/(sig*sqrt(2.0*3.1415926535))*exp(-pow(evalt-mu,2)/(2*pow(sig,2) ) );

p1d= (upper-mu)/sig*(double)100.0 +(double)301.0;
p2d=(lower-mu)/sig*(double)100.0 +(double)301.0;
p1 = imin( (int) p1d,600);
p1 = imax( p1,0);
p2 = imax( (int) p2d, 0);
p2 = imin( p2,600);
pr2 = *(pnormp + p2);
pr1 = *(pnormp + p1);

denstn=denstn/(pr1-pr2);

return(denstn);

}
/*----------------------------------------------------------- */

/* IMIN */
int imin( int n1, int n2)
{
  if (n1<n2)
    return(n1);
  else
    return(n2);
}



/*----------------------------------------------------------- */

/* MIN */
double mind( double xx, double yy)
{
  if (xx<yy)
    return(xx);
  else
    return(yy);
}

/* MAX */
double maxd( double xx, double yy)
{
  if (xx<yy)
    return(yy);
  else
    return(xx);
}


/* IMAX */
int imax( int n1, int n2)
{
  if (n1>n2)
    return(n1);
  else
    return(n2);
}


/*----------------------------------------------------------------------------------------*/
/* pseudo loglikelihood	                                                                 */
/*----------------------------------------------------------------------------------------*/

double likely(double *warpp, double *amp, double *kernelpars, int Namp, double *prep, double *X,double *Y, int Q, int S, int N, double *FreqA,double *FreqLl,double *FreqUl,double *FreqL,double *BreakA,double *BreakLl,double *BreakUl,double *BreakL,double ll, double ul, double thresh,double *WX,double *WYinterp,double *R1,double *R2,double *LambdaS, double *Lwarp,double *pnormp)
{
int r, i,r2,ii,i2;
double sum=0.0, sum2=0.0;
double loglikely, hulp2,prior,result, hulp1, hulp3;
int A_okay ,Ll_okay, Ul_okay, ifast; 
prior = 1;
  


for(i=0; i<(Q); i++) 
	{
	sum=0.0;
	for(ii=0; ii<(S-1); ii++) 
		{
	    sum= sum + *(warpp +3*Q+(i*(S-1))+ii);
		}
	*(LambdaS+i)=-sum;
	}



/* calculate the priors */
 if (Q > 1 && prior != 0)
    {
	for (i=0;i<(Q-1);i++)
	{
	A_okay = 0;
	Ll_okay = 0;
	Ul_okay = 0;
    for(r=0; r < 14; r++) 
       {
		if(A_okay == 1 && Ll_okay == 1 && Ul_okay == 1) {break;}
			
		if(A_okay == 0 && *(warpp + i) >= *(BreakA + (r)*Q+i) && *(warpp + i) <= *(BreakA + (r+1)*Q+i))
			{
			prior = prior* *(FreqA +(r)*Q+i);
			A_okay = 1;
			} 
		if(Ll_okay == 0 && *(warpp +Q +i) >= *(BreakLl +(r)*Q+i) && *(warpp +Q +i) <= *(BreakLl +(r+1)*Q+i) )
			{
			prior = prior* *(FreqLl +(r)*Q+i);
			Ll_okay = 1;
			}   
		if(Ul_okay == 0 && *(warpp +2*Q+i) >= *(BreakUl +(r)*Q+i) && *(warpp +2*Q+i) <= *(BreakUl +(r+1)*Q+i) )
			{
			prior = prior* *(FreqUl +(r)*Q+i);
			Ul_okay = 1;
			} 				 
		}
        for(r2=0; r2 < (S-1); r2++)       
            {  
           	 for(r=0; r < (14); r++) 
            	{
           		if(*(warpp +Q*3 + i*(S-1) +r2) >= *(BreakL + (r)*(Q*(S-1)) + i*(S-1)+r2)  && *(warpp +Q*3 + i*(S-1) +r2) <= *(BreakL + (r+1)*(Q*(S-1))+ i*(S-1)+r2) )
                	{
                	prior = prior* *(FreqL + (r)*(Q*(S-1)) + i*(S-1)+r2);
                	break; 
               	 	}
           		}
			}
		}
	}


for(i=0; i< S-1; i++)
{	
prior = prior * denstnormal( *(prep + i),0,thresh/4.0*(ul-ll),-(ul-ll)/4.0,(ul-ll)/4.0, pnormp);
prior = prior * denstnormal( *(warpp +Q*3 + (Q-1)*(S-1) +i),0,thresh,-1.0,1.0, pnormp);
}

if(prior != 0)  
{

for(i=0; i < Q; i++)
        {
        *(R1+i) = *(warpp + i) - *(warpp +Q +i);
        *(R2+i) = *(warpp +2*Q+i)  - *(warpp + i);
        }    


for(r = 0; r < S; r++)
    {
	if (r < S-1)
		{
		sum2 = sum2 - *(prep+r);
		for(i=Q-1; i >=0; i--)
			{
			*(Lwarp+i) = *(warpp +Q*3 + i*(S-1) +r);
			}
		}
	else 
		{
		for(i=Q-1; i >=0; i--)
			{
			*(Lwarp+i) = *(LambdaS+i);
			}
		*(prep + S-1) = sum2;
		}
	for(i=N-1; i>=0; i--)
		{	
		*(WX + r*N +i) = *(X +r*N +i) + *(prep + r);
		}
	warp3(warpp,Lwarp,R1,R2,(WX +r*N),N,S,Q,ll,ul);
	}


prior = prior * denstnormal( *(prep + S-1),0,thresh/4.0*(ul-ll),-(ul-ll)/4.0,(ul-ll)/4.0, pnormp);
for(i=Q-1; i >=0; i--)
			{
			prior = prior * denstnormal( *(Lwarp+i),0,thresh,-1.0,1.0, pnormp);
			}


/*  data likelihood */

loglikely = 0;


for(r = 0; r < (S); r++)
{
	/* Predicting the warped curves in the warped datapoints of the other curves (linear interpolation)) */
    for(r2 = 0; r2 < S; r2++)
        {
		ifast=0;
        for(i = 0; i < N; i++)
            {
            for(i2 = ifast; i2 < N; i2++)
                {
                if(*(WX + r*N +i) >= *(WX + r2*N +i2) && *(WX + r*N +i) <= *(WX + r2*N +i2+1))
                    {
                    *(WYinterp + r2*N +i)= *(Y + r2*N +i2) + ( (*(Y + r2*N +i2+1)-*(Y + r2*N +i2))/(*(WX + r2*N +i2+1)-*(WX + r2*N +i2)) )*(*(WX + r*N +i) - *(WX + r2*N +i2));
                    ifast = i2;
					break;
                    }
                }
			if (*(WX + r*N +i) > *(WX + r2*N +N-1)) {*(WYinterp + r2*N +i) = *(Y + r2*N +N-1);}
			if (*(WX + r*N +i) < *(WX + r2*N )) {*(WYinterp + r2*N +i) = *(Y + r2*N );}
            }
       }
    if(r > 0)
        {
        for(r2 = 0; r2 < (r); r2++)
            {
            hulp2=0;
			hulp1=0;
            for(i = 0; i < (N); i++)
                {
				hulp3=0;
				for (i2=0;i2<Namp;i2++)
					{
					hulp3=hulp3+( 2.0*pow(kernel(kernelpars,*(WX + r*N +i),Namp)* *(amp+i2),2) );
					}
				hulp2 = hulp2 + (pow( (*(Y + r*N +i)-*(WYinterp + r2*N +i)),2))/(2.0*(hulp3+pow(*(warpp + Q*3 + (S-1)*Q ),2)));
				hulp1= hulp1 + log(  sqrt( 2.0*3.1415926535* (hulp3+pow(*(warpp + Q*3 + (S-1)*Q ),2) ) )  );
				} 
			loglikely = loglikely - (double)1/((S-1)*S)*(hulp2+hulp1);    
            }
        }
    if(r < (S-1))
        {
        for(r2 = (r+1); r2 < (S); r2++)
           {
		   hulp1=0;	
           hulp2=0;
           for(i = 0; i < (N); i++)
                {
				hulp3=0;
 				for (i2=0;i2<Namp;i2++)
					{
					hulp3=hulp3+(2.0*pow(kernel(kernelpars,*(WX + r*N +i),Namp)* *(amp+i2),2) );
					}
				hulp2 = hulp2 + (pow( (*(Y + r*N +i)-*(WYinterp + r2*N +i)),2))/(2.0*(hulp3+pow(*(warpp + Q*3 + (S-1)*Q ),2)));
				hulp1= hulp1 + log(  sqrt( (double) 2.0*3.1415926535* (hulp3+pow(*(warpp + Q*3 + (S-1)*Q ),2) ) )  );
				} 
			loglikely = loglikely - (double)1/((S-1)*S)*(hulp2 + hulp1);      
			}
        }    
}      



result = loglikely+log(prior);


} /* ending of loglikelihood calculation */

else { result=0;} /* for the case of invalid warping pars */

return(result);

}
/*-----------------------------------------------------------------------------------------------------------------------*/



/* NORMALDENS */
double normaldens(double t, double sigma, double mu)
	{
    double dens;
	dens = (double) 1/(sigma* (double) sqrt((double)2*3.1415926535))*exp(-pow(t-mu,2)/(2*pow(sigma,2)));
    return(dens);
	}

/*--------------------------------------------------------------------------------------------------------------------------*/
/* applies a composition of warping components onto timepoints (kernal = Quartic)        	            */
/*--------------------------------------------------------------------------------------------------------------------------*/

void warp3(double *pA, double *pL, double *pR1,double *pR2,double *px, int size, int S, int comps, double ll, double ul)
{
int i,i2;
double r, lambda1, lambda2, c, P, Q, R, U, y,u,lambda,r1,r2,a;

for (i = 0; i < comps; i++)
	{
	a = *(pA+i);
	lambda=*(pL+i);
	r1=*(pR1+i);
	r2=*(pR2+i);

	r = mind(r1,r2);

	if (lambda!=0)
		{
		for (i2 = 0; i2 < size; i2++)
		 {
		 if (*(px+i2)< (a + r2) && *(px+i2) >= (a-3.0*sqrt(3.0)/8.0*lambda*r) )
            {
            *(px+i2) = (*(px+i2)-a)/r2;
            lambda2 = lambda*r/r2;

            c = (float)8.0/(3.0*sqrt(3.0)*lambda2);
            P = -(float)4/3-c* *(px+i2);
            Q = (float)2/27-(float)2/3*(1.0+c* *(px+i2))-(pow(c,2))/8;
            R = Q/2.0 - sqrt( pow((Q/2.0),2)+pow((P/3.0),3));
            if (R==0) 
                {
                U=0;
                y = (float)5/3 - U; 
                }
            else 
                {
                U = (fabs(R)/R)*pow((fabs(R)),(float)1/3 );
                y = (float)5/3 - U;
                y = y + P/(3.0*U);
                }
            *(px +i2) = sqrt(-2.0+2.0*y);
            if (lambda > 0)
                {   u =  (float)1/2.0*(  *(px+i2) - sqrt( -(-6.0+2.0*y +(-2.0*c)/(*(px+i2)))  )  ); }
            else    
                {   u =  (float)1/2.0*(- *(px+i2) + sqrt( -(-6.0+2.0*y -(-2.0*c)/(*(px+i2)))  )  ); }
            *(px+i2) = a+r2*( u+(1.0/c)*pow((1.0-pow(u,2)),2) );
            }
		else{
			if (*(px+i2) < a-3.0*sqrt(3.0)/8.0*lambda*r && *(px+i2) > a-r1)
				{
				*(px+i2) = (*(px+i2)-a)/r1;
				lambda1 = lambda*r/r1;
            
				c = (float)8.0/(3.0*sqrt(3.0)*lambda1);
				P = -(float)4/3-c* *(px+i2);
				Q = (float)2/27-(float)2/3*(1.0+c* *(px+i2))-(pow(c,2))/8;
				R = Q/2.0 - sqrt( pow((Q/2.0),2)+pow((P/3.0),3));
				if (R==0) 
				   {
				   U=0;
				   y = (float)5/3 - U; 
				   }
				else 
				  {
				  U = (fabs(R)/R)*pow((fabs(R)),(float)1/3 );
				  y = (float)5/3 - U;
				  y = y + P/(3.0*U);
				  }
				*(px+i2) = sqrt(-2.0+2.0*y);
				if (lambda > 0)
					{   u =  (float)1/2.0*(  *(px+i2) - sqrt( -(-6.0+2.0*y +(-2.0*c)/(*(px+i2)))  )  ); }
				else    
					{   u =  (float)1/2.0*(- *(px+i2) + sqrt( -(-6.0+2.0*y -(-2.0*c)/(*(px+i2)))  )  ); }
				*(px+i2) = a+r1*( u+(1.0/c)*pow((1.0-pow(u,2)),2) );
				}
			}
		}
	}

    }

}

/*----------------------------------------------------------------------------------------*/
/* Random uniform numbers													    		  */
/*----------------------------------------------------------------------------------------*/

double runif()
{
	double y;
	y= ( (double)rand() )/RAND_MAX;
	return(y);

}

/*-------------------*/
/* Kernel */
/*-------------------*/
double kernel(double *kernelpars, double t,int Namp)
{
int i;
double y=0;

for (i=0; i<Namp; i++)
	{ 
	if (t > *(kernelpars+3*i+1) && t < *(kernelpars+3*i+2))
		{y = y + pow( 1-pow( (t-*(kernelpars+3*i+1))/ (*(kernelpars+3*i+2)-*(kernelpars+3*i+1)),2),2);}
	if (t <= *(kernelpars+3*i+1) && t > *(kernelpars+3*i))
		{y = y + pow(  1-pow( (t-*(kernelpars+3*i+1))/ (*(kernelpars+3*i+1)-*(kernelpars+3*i)),2),2);}
	}
return(y);
}



