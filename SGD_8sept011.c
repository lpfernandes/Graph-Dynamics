#include "include_lib/general.h"

int MEASURE_TRACES=0;
int DEFORMED_MEASURE=0;

int main()
{
FILE *fp,*fopen(),*fp1;
int choice,getint(),answer(),option(),*degrees,*ivector();
int kmax,i,j,N,view;
double **dmatrix(),*Pk,*dvector(),**Prob,RANDOM(),exp(),log(),kav,kvariance,loops3,loops4;
 double DComplexity, PiComplexity;
float getfloat();
 char **c,**cmatrix(),filename[100];//dum;
void test_symmetry(),save_network(),degree_complexity(),calculate_Pi(),calculate_Pi_complexity(),entropy();
void calculate_degree_stats(),run_dynamics();

SETRANDOM();
puts("graph dynamics\nversion of Aug 24, 2011");
puts("A Annibale, Coolen, Fernandes\n");
puts("deformed dynamics to produce graphs with biological Pi or with (dis)assortative Pi");

  printf("Choices: 1 - Literature Edge Switches (No detailed balance); 2 - Target Pi=1; 3 - Target arbitrary Pi\n"); 
  choice=getint("choice=?\n");

 	printf("filename: "); gets(filename); fp=fopen(filename,"r"); 
	if(fp==NULL) { puts("cannot open file"); }
	view=answer("view data while reading");
	/*
	  
	  fscanf(fp,"%d",&N); printf("N=%d",N); puts("");
	  c=cmatrix(1,N,1,N); degrees=ivector(1,N); Pk=dvector(0,N); 
	  
	  for(i=1;i<=N;i++) for(j=1;j<=N;j++) {
	  dum=' '; while((dum==' ')||(dum=='\b')||(dum=='\n')) fscanf(fp,"%c",&dum); 
	  if(view) printf("%c",dum);
	  if((dum!='1')&&(dum!='0')) { 
	  fclose(fp); puts("\nerror in network file data"); 
	  printf("c[%d][%d]=%c",i,j,dum); puts(""); getchar(); 
	  free_cmatrix(c,1,N,1,N);  free_ivector(degrees,1,N); free_dvector(Pk,0,N);   
	  NETWORK_PRESENT=0;
	  }
	  c[i][j]=dum;
	  }
	  fclose(fp);
	*/
	
	fscanf(fp,"%d\n",&N); printf("N=%d\n",N); puts("");
	c=cmatrix(1,N,1,N); degrees=ivector(1,N); Pk=dvector(0,N); 
	
	
	fp1=fopen("InitialConnectivityN.txt","w");
	fprintf(fp1,"%d\n",N);
	for(i=1;i<=N;i++) {
	  for(j=1;j<=N;j++) {
	    fscanf(fp,"%c ",&c[i][j]);
	    if(view) fprintf(fp1,"%c ",c[i][j]);
	  }
	  if(view) fprintf(fp1,"\n");
	}
	fclose(fp1); fclose(fp);
	
	/*
	  for(i=1;i<=N;i++) for(j=1;j<=N;j++) {
	  dum=' '; while((dum==' ')||(dum=='\b')||(dum=='\n')) fscanf(fp,"%c",&dum); 
	  
	  if(view) printf("%c",dum);
	  if((dum!='1')&&(dum!='0')) { 
	  fclose(fp); puts("\nerror in network file data"); 
	  printf("c[%d][%d]=%c",i,j,dum); puts(""); getchar(); 
	  free_cmatrix(c,1,N,1,N);  free_ivector(degrees,1,N); free_dvector(Pk,0,N);   
	  NETWORK_PRESENT=0;
	  }
	  c[i][j]=dum;
	  }
	  fclose(fp);
	  */
	view=0; for(i=1;i<=N;i++) if(c[i][i]=='1') view++; 
	printf("self-interactions:   %d of %d",view,N); puts("");
	if(view>0) if(answer("remove self-interactions")) for(i=1;i<=N;i++) c[i][i]='0';
	view=0; for(i=1;i<=N;i++) for(j=i;j<=N;j++) if(j>i) if(c[i][j]!=c[j][i]) view++; 
	printf("non-symmetric bonds: %d of %d",view,(N*(N-1))/2); puts("");
	puts("");
	
	calculate_degree_stats(Pk,degrees,c,N,&kmax,&kav,&kvariance); 
	degree_complexity(&DComplexity,Pk,degrees,kav,N); 
	
	Prob=dmatrix(1,N,1,N); 
	calculate_Pi("bionet_Pi_raw.txt",c,Prob,Pk,degrees,kmax,kav,N);
	calculate_Pi_complexity(PiComplexity,Prob,degrees,kav,N);
	entropy(DComplexity,PiComplexity,kav,N);
	free_dmatrix(Prob,1,N,1,N); 
	puts("");
	
	switch(choice){
	  
	case 1: puts("edge swap dynamics without mobility correction");
          MEASURE_TRACES=answer("calculate mobility terms");
          run_dynamics(0,c,degrees,Pk,N,&loops3,&loops4,kav,kvariance,kmax);
          if(MEASURE_TRACES) {
            printf("<3-loops/node>=%lf\n",loops3);
            printf("<4-loops/node>=%lf\n",loops4);
          }
          Prob=dmatrix(1,N,1,N);
          calculate_Pi("Pi_false_null.txt",c,Prob,Pk,degrees,kmax,kav,N);
          calculate_Pi_complexity(PiComplexity,Prob,degrees,kav,N);
          entropy(DComplexity,PiComplexity,kav,N);
          free_dmatrix(Prob,1,N,1,N);
          puts(""); break;


	case 2: puts("edge swap dynamics with mobility correction");
	  MEASURE_TRACES=1;
	  run_dynamics(1,c,degrees,Pk,N,&loops3,&loops4,kav,kvariance,kmax);
	  printf("<3-loops/node>=%lf\n",loops3);
	  printf("<4-loops/node>=%lf\n",loops4);
	  Prob=dmatrix(1,N,1,N); 
	  calculate_Pi("Pi_true_null.txt",c,Prob,Pk,degrees,kmax,kav,N);
	  calculate_Pi_complexity(PiComplexity,Prob,degrees,kav,N);
	  entropy(DComplexity,PiComplexity,kav,N);
	  free_dmatrix(Prob,1,N,1,N); 
	  puts(""); break; 
	  
	case 3: puts("canonical edge swap dynamics with mobility correction and degree correlations");
	  MEASURE_TRACES=DEFORMED_MEASURE=1;
	  run_dynamics(1,c,degrees,Pk,N,&loops3,&loops4,kav,kvariance,kmax);
	  printf("<3-loops/node>=%lf\n",loops3);
	  printf("<4-loops/node>=%lf\n",loops4);
	  Prob=dmatrix(1,N,1,N); 
	  calculate_Pi("Pi_true_null.txt",c,Prob,Pk,degrees,kmax,kav,N);
	  calculate_Pi_complexity(PiComplexity,Prob,degrees,kav,N);
	  entropy(DComplexity,PiComplexity,kav,N);
	  free_dmatrix(Prob,1,N,1,N); 
	  puts(""); break; 
	  
	default: break;
	}


	free_cmatrix(c,1,N,1,N); free_ivector(degrees,1,N); free_dvector(Pk,0,N);  
	puts("finished normally"); 
	puts("");
	
	return(0);
}

void save_network(char *filename,char **c,int N) 
{
FILE *fp,*fopen();
int i,j;
	fp=fopen(filename,"w"); fprintf(fp,"%d\n",N);
	for(i=1;i<=N;i++) {
		for(j=1;j<=N;j++) if(c[i][j]=='1') fprintf(fp,"1 "); else fprintf(fp,"0 ");
		fprintf(fp,"\n");
		}
	fclose(fp);
}


void calculate_degree_stats(double *Pk,int *degrees,char **c,int N,int *kmax,double *kav,double *kvar)
{
int i,j,k,Kmax,deg;
double Kav,Kvar,sum;
FILE *fp,*fopen();
	puts("degree distribution ..."); 
	Kmax=0;
	for(i=1;i<=N;i++) { c[i][i]='0'; degrees[i]=0; Pk[i]=0.0; } Pk[0]=0.0; 
	for(i=1;i<=N;i++) { 
		for(j=1;j<=N;j++) if(c[i][j]=='1') degrees[i]++;  
		if(degrees[i]>Kmax) Kmax=degrees[i];
		}
    fp=fopen("degrees.txt","w");
	for(i=1;i<=N;i++) fprintf(fp,"%d %d\n",i,degrees[i]);
	fclose(fp);

	for(k=0;k<=N;k++) {
		for(i=1;i<=N;i++) { if(degrees[i]==k) Pk[k]+=1.0; }
		Pk[k]=Pk[k]/(double)N;
		}
	sum=Pk[0]; Kav=Kvar=0.0; for(i=1;i<=N;i++) { 
		deg=degrees[i]; Kav=Kav+(double)deg; Kvar=Kvar+(double)(deg*deg); sum+=Pk[i];
		}
	Kav=Kav/(double)N; Kvar=Kvar/(double)N;
	printf("kmax=%d, <k>=%lf, <kk>=%lf, sum_k P(k)=%lf",Kmax,Kav,Kvar,sum); puts("");
	*kmax=Kmax; *kav=Kav; *kvar=Kvar;

	puts("saving P(k) ..."); 
    fp=fopen("Pk.txt","w");
	fprintf(fp,"# kmax=%d, <k>=%lf, <kk>=%lf\n",Kmax,Kav,Kvar);
	fprintf(fp,"%lf %lf\n",-0.5,0.0); 
	for(k=0;k<=Kmax;k++) fprintf(fp,"%lf %lf\n%lf %lf\n",(double)k-0.5,Pk[k],(double)k+0.5,Pk[k]);
	fprintf(fp,"%lf %lf\n",(double)Kmax+0.5,0.0);
	fclose(fp);
}


void calculate_Pi(char *filename,char **c,double **Prob,double *Pk,int *degrees,int kcut,double kav,int N)
{
char filename2[100];
int k,l,i,j,KCUT,getint(),filename_extension();
double value,factor;
void test_symmetry();
FILE *fp,*fpl,*fopen();

	puts("calculating Pi(k,k') ..."); 
   
    factor=kav*(1.0-1.0/(double)N)/(double)N;
	for(k=1;k<=kcut;k++) for(l=1;l<=k;l++) Prob[k][l]=Prob[l][k]=0.0;
	for(i=1;i<=N;i++) for(j=1;j<=N;j++) if((i!=j)&&(c[i][j]=='1')) {
		k=degrees[i]; l=degrees[j]; if((k<=kcut)&&(l<=kcut)) Prob[k][l]+=1.0;
		}
	for(k=1;k<=kcut;k++) for(l=1;l<=k;l++) {
		if(Prob[k][l]>0.0) {
			if(k!=l) value=Prob[k][l]*factor/((double)(k*l)*Pk[k]*Pk[l]);
			else value=Prob[k][l]*factor/((double)(k*l)*Pk[k]*(Pk[l]-1.0/(double)N));
			Prob[k][l]=Prob[l][k]=value; 
			}
		else if((Pk[k]*Pk[l])==0.0) Prob[k][l]=Prob[l][k]=1.0; /* this degree combination does not occur */
		}
	test_symmetry(Prob,1,kcut,"Pi");

	KCUT=kcut+1; while((KCUT<1)||(KCUT>kcut)) {
		printf("largest degree in saving Pi(k,k') (max=%d)",kcut);
		KCUT=getint("");
		}

	strcpy(filename2,filename); filename_extension(filename2,"dat");
	fp=fopen(filename,"w"); fpl=fopen(filename2,"w"); 
	for(l=1;l<=KCUT;l++) {
		for(k=1;k<=KCUT;k++) {
			fprintf(fp,"%d %d %lf\n",k,l,Prob[k][l]);
			fprintf(fpl,"%lf ",Prob[k][l]);
			}
		fprintf(fp,"\n"); fprintf(fpl,"\n");
		}
	fclose(fp); fclose(fpl);
}





void degree_complexity(double *DComplexity, double *Pk,int *degrees,double kav,int N) 
{
int i; 
double log(),logfactorial(),dComplexity; 
FILE *fp,*fopen();
	dComplexity=0.0; for(i=1;i<=N;i++) dComplexity=dComplexity+log(Pk[degrees[i]])
		-(double)degrees[i]*log(kav)+kav+logfactorial(degrees[i]);
	dComplexity=dComplexity/(double)N;
	*DComplexity=dComplexity;
	printf("degree complexity: %lf",dComplexity); puts("");
	fp=fopen("complexity.txt","a");
	fprintf(fp,"N=%d, <k>=%lf: C_degrees=%lf\n",N,kav,dComplexity);
	fclose(fp);
}


void calculate_Pi_complexity(double PiComplexity, double **Prob,int *degrees,double kav,int N)
{
double log(),value,Assortativity,Norm,k_i,k_j,Nsqr,variance,average;
int i,j;
FILE *fp,*fopen();
	PiComplexity=Assortativity=Norm=variance=average=0.0; Nsqr=(double)(N*N);
	for(i=1;i<=N;i++) for(j=1;j<=N;j++) if((degrees[i]>0)&&(degrees[j]>0)) {
		value=Prob[degrees[i]][degrees[j]]; k_i=(double)degrees[i]; k_j=(double)degrees[j];
		if(value>0.0) {
			PiComplexity=PiComplexity+k_i*k_j*value*log(value);
			Norm=Norm+(k_i*k_j*value); 
			average=average+0.5*(k_i+k_j)*k_i*k_j*value;
			variance=variance+0.5*(k_i*k_i+k_j*k_j)*k_i*k_j*value;
			k_i=k_i*k_i; k_j=k_j*k_j;
			Assortativity=Assortativity+(k_i*k_j*value);
			}
		}
	PiComplexity=PiComplexity/(2.0*kav*Nsqr);
	Norm=Norm;
	Assortativity=Assortativity/Norm; average=average/Norm; variance=variance/Norm;
	Assortativity=(Assortativity-(average*average))/(variance-(average*average));

	printf("Pi-complexity: %lf",PiComplexity);    puts("");
	printf("Assortativity: %lf",Assortativity); puts("");
	fp=fopen("complexity.txt","a");
	fprintf(fp,"N=%d, <k>=%lf: C_Pi=%lf, R=%lf\n",N,kav,PiComplexity,Assortativity);
	fclose(fp);
}

void entropy(double DComplexity, double PiComplexity, double kav, int N)
{
  double S;
FILE *fp,*fopen();

 S=0.5*kav*(log(N/kav)+1)-DComplexity-PiComplexity;
 printf("Entropy: %lf",S);
 fp=fopen("complexity.txt","a");
 fprintf(fp,"N=%d, <k>=%lf: S=%lf\n",N,kav,S);
 fclose(fp);
}

void test_symmetry(double **matrix,int lower,int upper,char *text)
{
int ok=1,i,j; 
	for(i=lower;i<=upper;i++) for(j=lower;j<=upper;j++) if(matrix[i][j]!=matrix[j][i]) ok=0;
	if(ok) { printf("symmetry test %s ok",text); }
	else   { printf("symmetry test %s failed",text); }
	puts("");
}


double logfactorial(int k)
{
int i;
double log(),value;
	if(k<=1) return(0.0);
	else if(k==2) return(LOGTWO);
	value=LOGTWO; for(i=3;i<=k;i++) value+=log((double)i);
	return(value);
}





void run_dynamics(int correct,char **c,int *degrees,double *Pk,int N,double *loops3,double *loops4,double kav,double kvar,int kmax)
{
int a,b,i,j,ki,lj,k,l,m,n,switchable,do_swap,**imatrix(),**cint,steps1,steps2,max1,max2,MAX2,**neighbours,dum, length;
int *ivector(),*newdegrees,mini,maxi,minj,maxj,NN,mmax,dim,tpic, time;
double RANDOM(),exp(),log(),twoNa,Na,Nb,oldNb,Trc3,Trc4,Loops3,Loops4,dcounter,probability,stopwatch();
double sqrt(),mloops3,mloops4,mcounter,NC,**dmatrix(),dN,factor,value;
double **Prob, **QQ, **L, Hamming;
void count_loops(),start_stopwatch(),save_network(),symmetrize_normalize();
float getfloat();
char **cold, **cmatrix();
FILE *fopen(),*fp,*fp1,*fp2;
double Delta3,Delta4;
 int np, nr, p,r,somma3,somma4,sum3,sum4;
int temp, edges, *T, *nachvec, M;
 int Rl, Ak, Al;
 void interrupt_dynamics();

	NN=N*N; dN=(double)N;
	neighbours=imatrix(1,N,0,kmax); newdegrees=ivector(1,N); cint=imatrix(1,N,1,N);
	for(i=1;i<=N;i++) neighbours[i][0]=0;                        /* we store in [i][0] the degrees k_i of the nodes */
	for(i=1;i<=N;i++) for(j=1;j<=N;j++) if(i!=j) if(c[i][j]=='1') { 
		neighbours[i][0]+=1; neighbours[i][neighbours[i][0]]=j;  /* entries [i][1]->[i][k_i] contain neighbours of i */ 
		}
	///////////
	T=ivector(0,N); T[0]=0;
	temp=0;
	for(i=1;i<=N;i++) 
	  {
	    temp=temp+degrees[i];  
	    T[i]=temp;     //cumulative degrees
	  }
	edges=T[N];
	printf("total edges %d\n",edges);
	nachvec=ivector(1,edges);
	for(i=1;i<=N;i++) for(k=1;k<=(degrees[i]);k++) nachvec[1+T[i-1]+k-1]=neighbours[i][k];

	/////////

	if(MEASURE_TRACES||correct) {
		count_loops(&Trc3,&Trc4,N,c,neighbours,cint); mloops3=Trc3/(double)N; mloops4=Trc4/(double)N;
		Na=0.25*kav*(1.0+(double)N*kav)-0.5*kvar; Na=Na*(double)N; twoNa=2.0*Na;
		Nb=0.25*Trc4+0.5*Trc3;
		for(a=1;a<=N;a++) if(degrees[a]>0) {
			dum=0; for(b=1;b<=degrees[a];b++) dum=dum+degrees[neighbours[a][b]];
			dum=dum*degrees[a];
			Nb=Nb-0.5*(double)dum;
			}
		fp=fopen("loop_stats.txt","a");
		fprintf(fp,"# N=%d, <k>=%lf, <kk>=%lf\n",N,kav,kvar);
		fprintf(fp,"%d %lf %lf %lf\n",0,mloops3,mloops4,(Na+Nb)/(double)(N*N)); fclose(fp);
	}

	k=NN; for(i=1;i<=N;i++) for(j=1;j<=N;j++) if((c[i][j]!='1')&&(c[i][j]!='0')) k--;
	if(k<NN) { puts("problem in matrix c"); getchar(); return; }
	else puts("no problem in matrix c");
	k=0;   for(i=1;i<=N;i++) if(neighbours[i][0]!=degrees[i]) k++;
	if(k>0) { puts("problem in neighbourhood matrix"); getchar(); return; }
	else puts("no problem in neighbourhood matrix");

 ///////////////////////////////////
	if(DEFORMED_MEASURE) {    // run dynamics targeting a certain Pi, here either the one of the network or a (dis)assortative
	  L=dmatrix(1,N,1,N); 	Prob=dmatrix(1,N,1,N);
	  ////////
	  if(answer("Pi biological")){
	    puts("preparing L-matrix ..."); // from the biological Pi and using Q canonical 
	    factor=kav*(1.0-1.0/(double)N)/(double)N;
	    for(k=1;k<=N;k++) for(l=1;l<=k;l++) Prob[k][l]=Prob[l][k]=0.0;
	    for(i=1;i<=N;i++) for(j=1;j<=N;j++) if((i!=j)&&(c[i][j]=='1')) {
	      k=degrees[i]; l=degrees[j]; 
	      Prob[k][l]+=1.0;
	    }
	    for(k=1;k<=kmax;k++) for(l=1;l<=k;l++) {
	      if(Prob[k][l]>0.0) {
		if(k!=l) value=Prob[k][l]*factor/((double)(k*l)*Pk[k]*Pk[l]);
		else value=Prob[k][l]*factor/((double)(k*l)*Pk[k]*(Pk[l]-1.0/(double)N));
		Prob[k][l]=Prob[l][k]=value; 
	      }
	      else if((Pk[k]*Pk[l])==0.0) Prob[k][l]=Prob[l][k]=1.0; /* this degree combination does not occur */
	    } 
	    /////////
	    for(i=1;i<=N;i++) {
	      L[i][i]=0.0;
	      for(j=1;j<i;j++) {
		if((degrees[i]>0)&&(degrees[j]>0)) L[i][j]=L[j][i]=dN*kav/(degrees[i]*degrees[j]*Prob[degrees[i]][degrees[j]])-1.0;
		else L[i][j]=L[j][i]=0.0;
	      }
	    }
	    puts("done");
	  }
	  ////////
	  else{
	    QQ=dmatrix(0,kmax,0,kmax);          // from (dis)assortative Q
	    if(answer("disassortative Pi")){
	      if(answer("strongly")){
		for(k=0;k<=kmax;k++) for(l=0;l<=kmax;l++) QQ[k][l]=(double)(k-l)*(k-l);    // here of the form Q(k,k')=(k-k')^2/C 
	      }
	      else for(k=0;k<=kmax;k++) for(l=0;l<=kmax;l++) QQ[k][l]=sqrt((double)(k-l)*(k-l));    // here of the form Q(k,k')=|k-k'|/C 
	    }
	    else{
	      if(answer("strongly assortative")){
		for(k=0;k<=kmax;k++) for(l=0;l<=kmax;l++) QQ[k][l]=1.0/(1.0+(double)(k-l)*(k-l));    // here of the form Q(k,k')=1/(1+(k-k')^2)/C
	      }
	      else for(k=0;k<=kmax;k++) for(l=0;l<=kmax;l++) QQ[k][l]=1.0/(1.0+sqrt((double)(k-l)*(k-l)));    // here of the form Q(k,k')=1/(1+|k-k'|)/C
	    }
	    printf("symmetrizing Q..\n");
	    symmetrize_normalize(QQ,Pk,kmax);
	    
	    puts("preparing L ...\n"); 
	    
	    for(i=1;i<=N;i++) {
	      L[i][i]=0.0;
	      for(j=1;j<i;j++) {
		if((degrees[i]>0)&&(degrees[j]>0)) L[i][j]=L[j][i]=dN/(kav*QQ[degrees[i]][degrees[j]])-1.0;
		else L[i][j]=L[j][i]=0.0;
	      }
	    }
	    puts("done\n");

	    // useful to know also the Pi one is targeting for plots: use canonical kernel
	    printf("calculating target Pi..\n");
	    for(k=1;k<=kmax;k++) {
	      for(l=1;l<=kmax;l++) {
		Prob[k][l]=(QQ[k][l])*kav*kav/(k*l);
	      }
	    }
	    fp1=fopen("bionet_Pi_target.txt","w"); fp2=fopen("bionet_Pi_target.dat","w"); 
	    for(l=1;l<=kmax;l++) {
	      for(k=1;k<=kmax;k++) {
		fprintf(fp1,"%d %d %lf\n",k,l,Prob[k][l]);
		fprintf(fp2,"%lf ",Prob[k][l]);
	      }
	      fprintf(fp1,"\n"); fprintf(fp2,"\n");
	    }
	    fclose(fp1); fclose(fp2);	
	    

	  }
	  free_dmatrix(Prob,1,N,1,N);
	}
	
	 max1=1+(int)kav; max2=1+N/10;  mmax=100; puts("default: about 10 swaps/bond");   // mmax=was 100 
	if(!answer("OK")) { mmax=-1; while(mmax<1) mmax=10*getint("nr of swaps/bond"); }
	start_stopwatch();
	
	tpic=0;
	Loops3=Loops4=dcounter=0.0; oldNb=INFTY; mini=maxi=minj=maxj=N/2;	

 length=getint("max length of loops");
 interrupt_dynamics(N,c,0,length);
 time=0;
	for(m=0;m<=mmax;m++) { //external loops   
		if(m==tpic) {time++; printf("\n"); interrupt_dynamics(N,c,time,length); tpic=m+mmax/10;}           
		mloops3=mloops4=mcounter=NC=0.0;    /* monitor loops statistics over periods of max1*max2 steps */
		printf("."); fflush(stdout);
		for(steps1=1;steps1<=max1;steps1++){
		  for(steps2=1;steps2<=max2;steps2++) {
/////////////////////////////////////////
		  switchable=0;
		  while(switchable==0){
		    M=edges;
		    Ak=(int)(RANDOM()*M)+1;            
		    k=nachvec[Ak]; 
		    i=1;
		    while((T[i])<Ak) i++;                // got i-k
		    ki=Ak-T[i-1];                       // position of k along the neighbours of i
		    for(a=1;a<=(degrees[i]);a++){                        
		      M=M-degrees[neighbours[i][a]];     // shorten nachvec by chalking out all the segments referring to the neighbours of i..
		    }
		    Al=0; j=1;
		    Rl=(int)(RANDOM()*(M-(degrees[i])))+1;  //.. and to i itself, so seek j!=i,k and disconnected from i  
		    do{                                          
		      if((j==i)||((c[i][j])=='1')) {  // hop segment referring to i or neighbours of i
			Al=Al+degrees[j];   
			j++;
		      }
		      else{
			if(Rl>degrees[j]){      // if in an allowed segment, check if Rl is farther than the length of the segment.  
			  Al=Al+degrees[j];    // if this is the case jump the segment
			  Rl=Rl-degrees[j];    // reduce number of remaining steps along allowed position we have to do
			  j++;
			}
			else{                     // if Al is within k[j] steps then make unitary steps
			  Rl--;
			  Al++;
			}
		      }
		    }                            // so found j
		    while(Rl>=1);
		    l=nachvec[Al];              // l
		    switchable=((l!=i)&&(l!=k)&&((c[l][k])=='0')); // l!=i,k and k,l disconnected  // end if fulfilled
		  }
		  lj=Al-T[j-1];                  // position of l along the neighbour of j
		  switchable=((c[i][k]=='1')&&(c[j][l]=='1')&&(c[i][j]=='0')&&(c[k][l]=='0'));   // just to check
//////////////////////////////////////////

		  if(!switchable) printf("problems: ik=%c, jl=%c, ij=%c, kl=%c\n",c[i][k],c[j][l],c[i][j],c[k][l]);									if(switchable) { 
			  if(i<mini) mini=i; if(i>maxi) maxi=i; if(j<minj) minj=j; if(j>maxj) maxj=j;

			  // Calculate Trace variation

				if(correct||MEASURE_TRACES){	    
					somma3=0;  somma4=0;
					for(np=1;np<=degrees[k];np++){
						p=neighbours[k][np];
						if((p!=i)&&(p!=j)&&(p!=k)&&(p!=l)){
							if(c[l][p]=='1') somma3++;
							if(c[p][i]=='1') somma3--;
							}
						for(nr=1;nr<=degrees[i];nr++){
							r=neighbours[i][nr];
							if((c[r][p])=='1'){
								if((p!=i)&&(p!=k)&&(r!=i)&&(r!=k)&&(r!=j)&&(p!=l)&&((p!=j)||(r!=l))) 
									if(((c[k][p])=='1')&&((c[r][i])=='1')) somma4--;
								}
							}
					}
					for(np=1;np<=degrees[j];np++){
						p=neighbours[j][np];
						if((p!=i)&&(p!=j)&&(p!=k)&&(p!=l)){
							if(c[l][p]=='1') somma3--;
							if(c[p][i]=='1') somma3++;
							}
						for(nr=1;nr<=degrees[l];nr++){
							p=neighbours[j][np];
							r=neighbours[l][nr];
							if((c[r][p])=='1'){
								if((p!=l)&&(p!=j)&&(r!=j)&&(r!=l)&&(p!=i)&&(r!=k)&&((p!=k)||(r!=i))) 
									if(((c[r][l])=='1')&&((c[j][p])=='1')) somma4--;
								}
							}
					}
					for(np=1;np<=degrees[l];np++){
						for(nr=1;nr<=degrees[k];nr++){
							p=neighbours[l][np];
							r=neighbours[k][nr];
							if((c[r][p])=='1'){
								if((p!=k)&&(p!=l)&&(r!=l)&&(r!=k)&&(r!=i)&&(p!=j)&&((p!=i)||(r!=j))) 
									if(((c[r][k])=='1')&&((c[l][p])=='1')) somma4++;
								}
							}
					}
					for(np=1;np<=degrees[i];np++){
						for(nr=1;nr<=degrees[j];nr++){
							p=neighbours[i][np];
							r=neighbours[j][nr];
							if((c[r][p])=='1'){
								if((p!=j)&&(p!=i)&&(p!=k)&&(r!=i)&&(r!=j)&&(r!=l)&&((p!=l)||(r!=k))) 
									if(((c[i][p])=='1')&&((c[r][j])=='1')) somma4++; 
								}
							}
					}
					sum4=8*somma4;
					sum3=6*somma3;

					Delta3=(double)sum3; 
					Delta4=(double)sum4;
				}

				do_swap=1;
				c[i][j]=c[j][i]='1'; c[k][l]=c[l][k]='1'; c[i][k]=c[k][i]='0'; c[j][l]=c[l][j]='0';
                                
				// 
				if(correct) {
					neighbours[i][ki]=j; neighbours[j][lj]=i;
					for(n=1;n<=degrees[k];n++) if(neighbours[k][n]==i) neighbours[k][n]=l;
					for(n=1;n<=degrees[l];n++) if(neighbours[l][n]==j) neighbours[l][n]=k;
					//uncomment for checks start
					//count_loops(&CTrc3,&CTrc4,N,c,neighbours,cint); //Previous counting of loops

					Trc3=Trc3+Delta3;
					Trc4=Trc4+Delta4;

					
					  // printf("%.0lf sec, ",stopwatch()); fflush(stdout); */
					
					Nb=0.25*Trc4+0.5*Trc3;
					for(a=1;a<=N;a++) if(degrees[a]>0) {
						dum=0; for(b=1;b<=degrees[a];b++) dum=dum+degrees[neighbours[a][b]];
						dum=dum*degrees[a];
						Nb=Nb-0.5*(double)dum;
					}
					if(!DEFORMED_MEASURE) probability=(oldNb+Na)/(oldNb+Nb+twoNa);
					else {
						probability=L[i][j]*L[k][l]/(L[i][k]*L[l][j]); 
						probability=(Na+oldNb)/(Na+oldNb+(Na+Nb)*probability);
						}
					do_swap=(RANDOM()<probability);
					if(do_swap) oldNb=Nb;
				}	
				if(!do_swap) { /* reject proposed state and swap back */
					c[i][j]=c[j][i]='0'; c[k][l]=c[l][k]='0'; c[i][k]=c[k][i]='1'; c[j][l]=c[l][j]='1'; 
					
							Trc3=Trc3-Delta3;
							Trc4=Trc4-Delta4;

					neighbours[i][ki]=k; neighbours[j][lj]=l; 
					for(n=1;n<=degrees[k];n++) if(neighbours[k][n]==l) neighbours[k][n]=i;
					for(n=1;n<=degrees[l];n++) if(neighbours[l][n]==k) neighbours[l][n]=j;
				}
				//

				else { // for incorrect and accepted correct 

				  nachvec[Ak]=j; nachvec[Al]=i;    // update bond matrix
				  for(n=1;n<=degrees[k];n++) if(nachvec[n+T[k-1]]==i) nachvec[n+T[k-1]]=l;
				  for(n=1;n<=degrees[l];n++) if(nachvec[n+T[l-1]]==j) nachvec[n+T[l-1]]=k;

					dcounter+=1.0; mcounter+=1.0; 
					if(!correct) {                // for incorrect update neighbours
					  neighbours[i][ki]=j; neighbours[j][lj]=i; 
					  for(n=1;n<=degrees[k];n++) if(neighbours[k][n]==i) neighbours[k][n]=l;
					  for(n=1;n<=degrees[l];n++) if(neighbours[l][n]==j) neighbours[l][n]=k;
					}
					if(MEASURE_TRACES) {
					  if(!correct) {
					    
					    //uncomment for traces checks
					    //		count_loops(&CTrc3,&CTrc4,N,c,neighbours,cint);
					    
					    // Calculate Trace variation
					    somma3=0;  somma4=0;
 for(np=1;np<=degrees[k];np++){
  p=neighbours[k][np];
     if((p!=i)&&(p!=j)&&(p!=k)&&(p!=l)){
       if(c[l][p]=='1') somma3++;
       if(c[p][i]=='1') somma3--;
     }
   for(nr=1;nr<=degrees[i];nr++){
     r=neighbours[i][nr];
     if((c[r][p])=='1'){
       if((p!=i)&&(p!=k)&&(r!=i)&&(r!=k)&&(r!=j)&&(p!=l)&&((p!=j)||(r!=l))) if(((c[k][p])=='1')&&((c[r][i])=='1')) somma4--;
     }
   }
 }
for(np=1;np<=degrees[j];np++){
  p=neighbours[j][np];
     if((p!=i)&&(p!=j)&&(p!=k)&&(p!=l)){
       if(c[l][p]=='1') somma3--;
       if(c[p][i]=='1') somma3++;
     }
   for(nr=1;nr<=degrees[l];nr++){
     p=neighbours[j][np];
     r=neighbours[l][nr];
     if((c[r][p])=='1'){
       if((p!=l)&&(p!=j)&&(r!=j)&&(r!=l)&&(p!=i)&&(r!=k)&&((p!=k)||(r!=i))) if(((c[r][l])=='1')&&((c[j][p])=='1')) somma4--;
     }
   }
 }
for(np=1;np<=degrees[l];np++){
   for(nr=1;nr<=degrees[k];nr++){
     p=neighbours[l][np];
     r=neighbours[k][nr];
     if((c[r][p])=='1'){
       if((p!=k)&&(p!=l)&&(r!=l)&&(r!=k)&&(r!=i)&&(p!=j)&&((p!=i)||(r!=j))) if(((c[r][k])=='1')&&((c[l][p])=='1')) somma4++;
     }
   }
 }
for(np=1;np<=degrees[i];np++){
   for(nr=1;nr<=degrees[j];nr++){
     p=neighbours[i][np];
     r=neighbours[j][nr];
     if((c[r][p])=='1'){
       if((p!=j)&&(p!=i)&&(p!=k)&&(r!=i)&&(r!=j)&&(r!=l)&&((p!=l)||(r!=k))) if(((c[i][p])=='1')&&((c[r][j])=='1')) somma4++; 
     }
   }
 }
 sum4=8*somma4;
 sum3=6*somma3;

 Delta3=(double)sum3; 
 Delta4=(double)sum4;


					Trc3=Trc3+Delta3;
					Trc4=Trc4+Delta4;

      							Nb=0.25*Trc4+0.5*Trc3;
							for(a=1;a<=N;a++) if(degrees[a]>0) {
								dum=0; for(b=1;b<=degrees[a];b++) dum=dum+degrees[neighbours[a][b]];
								dum=dum*degrees[a];
								Nb=Nb-0.5*(double)dum;
								}
					  }
						Loops3=Loops3+Trc3/(double)N;   Loops4=Loops4+Trc4/(double)N; 
						mloops3=mloops3+Trc3/(double)N; mloops4=mloops4+Trc4/(double)N;
						NC=NC+Na+Nb;
					}
				}
		  }
		  }
		}

		if(MEASURE_TRACES) { 
		  fp=fopen("loop_stats.txt","a"); 
		  fprintf(fp,"%d %lf %lf %lf\n",(int)(0.5+dcounter),mloops3/mcounter,mloops4/mcounter,
			  NC/(mcounter*(double)(N*N))); 
		  fclose(fp);
		}
	}

	if(MEASURE_TRACES) {
	  Loops3=Loops3/dcounter; *loops3=Loops3; Loops4=Loops4/dcounter; *loops4=Loops4;
	}
	printf("\nduration: %.0lf seconds",stopwatch());					puts(""); 
	printf("switches executed/bond: %.0lf/%d",dcounter,edges);						puts("");
	printf("node ranges: i=(%d->%d), j=(%d->%d)",mini,maxi,minj,maxj);	puts("");
	
	for(i=1;i<=N;i++) { k=0; for(j=1;j<=N;j++) if(c[i][j]=='1') k++; newdegrees[i]=k; }
	k=N; for(i=1;i<=N;i++) if(degrees[i]!=newdegrees[i]) k--;
	if(k==N) puts("degrees left invariant");
	else printf("%d degrees changed",N-k); puts("");
	
	/*
	if(answer("save graph")) {
	  if(correct)	save_network("true_null.mat",c,N); 
	  else		save_network("false_null.mat",c,N); 
	}
	*/
	
	if(answer("calcualte Hamming distance")){
	  Hamming=0.0;
	  cold=cmatrix(1,N,1,N);
	  fp=fopen("InitialConnectivityN.txt","r");
	  fscanf(fp,"%d\n", &N);
	  for(i=1;i<=N;i++) {
	    for(j=1;j<=N;j++) {
	      fscanf(fp,"%c ",&cold[i][j]);
	      if(cold[i][j]!=c[i][j]) Hamming=Hamming+1.0;
	    }
	  }
	  Hamming=Hamming/(2*edges);
	  printf("Hamming=%lf \n",Hamming);
	  free_cmatrix(cold,1,N,1,N);
	}

	free_imatrix(cint,1,N,1,N); free_imatrix(neighbours,1,N,0,kmax); free_ivector(newdegrees,1,N);
	if(DEFORMED_MEASURE) free_dmatrix(L,1,N,1,N); free_ivector(T,0,N); free_ivector(nachvec,1,edges);
}


void count_loops(double *Trc3,double *Trc4,int N,char **c,int **neighbours,int **cint)
{
int i,j,k,l,nia,nib,nj,k_i,k_j,sum3=0,sum4=0;
	for(i=1;i<=N;i++) {
		k_i=neighbours[i][0]; if(k_i>0) {
			for(nia=1;nia<=k_i;nia++) for(nib=1;nib<=k_i;nib++) {
				j=neighbours[i][nia]; k=neighbours[i][nib];
				if(c[j][k]=='1') sum3++;
				k_j=neighbours[j][0];
				if(k_j>0) for(nj=1;nj<=k_j;nj++) {
					l=neighbours[j][nj]; if(c[k][l]=='1') sum4++;
					}
				}
			}
		}
	*Trc3=(double)sum3; *Trc4=(double)sum4;
}	


void symmetrize_normalize(double **QQ, double *Pk, int kmax)
{
int k,l;
double norm=0.0;
    for(k=0;k<=kmax;k++) for(l=0;l<=kmax;l++) QQ[k][l]=QQ[l][k]=0.5*(QQ[k][l]+QQ[l][k]); 
	for(k=0;k<=kmax;k++) for(l=0;l<=kmax;l++) norm=norm+Pk[k]*Pk[l]*QQ[k][l];
	for(k=0;k<=kmax;k++) for(l=0;l<=kmax;l++) QQ[k][l]=(QQ[k][l])/norm; 
}


void matrix_product(int **A, int **B, int **C, int N)
{
  int i, j, k,sum;

  for(i=1;i<=N;i++) for(j=1;j<=N;j++){
      sum=0;
      for(k=1;k<=N;k++)	sum=sum+(A[i][k])*(B[k][j]);
      C[i][j]=sum;
    }
}

void get_trace(int **C, double *Tr, int N, int clock)
{
  int i;
  double Trace=0.0;
  for(i=1;i<=N;i++) Trace=Trace+(double)(C[i][i]);
  *Tr=(double)Trace/((double)N*2.0*(double)clock);
}

/*
void create_path(char **D, char **PATH, char **NEW_PATH, int N)
{
int i,j,k;

for(i=1;i<=N;i++) for(j>i;j<=N;j++){
if((PATH[i][j])=='1') for(k>j;k<=N;k++){
if((D[j][k])=='1') NEW_PATH[i][k]='1';
else NEW_PATH[i][k]='0';
}
}
}
*/

/*
void create_neigh(char **E, int **neighbours, int N)
{
int i,j;

for(i=1;i<=N;i++){
neighbours[i][0]=0;
for(j=1;j<=N;j++) if(i!=j) if(E[i][j]=='1') { 
neighbours[i][0]+=1; neighbours[i][neighbours[i][0]]=j;   (entries [i][1]->[i][k_i] contain neighbours of i )
}
}
}
*/

void interrupt_dynamics(int N, char **c, int m, int length)
{
  FILE *fp, *fopen();
  int clock,i,j,k,l,loops;
  void create_neigh(), matrix_product(), get_trace(); 
  int **A,**B,**C,*ivector(),**imatrix(), **neighbours;
  char  **cmatrix(), Trace[50], Loop[50], MatrixC[50];
  double Tr;

  A=imatrix(1,N,1,N);
  B=imatrix(1,N,1,N);
  C=imatrix(1,N,1,N);

  // write C[m]
  
  sprintf(MatrixC, "C%d.dat",m);
  fp=fopen(MatrixC,"w");
 
  fprintf(fp,"%d\n",N);
  for(i=1;i<=N;i++){ 
    for(j=1;j<=N;j++){ 
      fprintf(fp, "%c ", c[i][j]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  

 for(i=1;i<=N;i++){
   for(j=1;j<=N;j++){ 
     if((c[i][j])=='1'){ C[i][j]=1; B[i][j]=1;}
     else{ C[i][j]=0; B[i][j]=0;}
   }
 }
 
sprintf(Trace, "Traces%d.dat",m);
fp=fopen(Trace,"w");
for(clock=1; clock<=length; clock++)
  {
    get_trace(C,&Tr,N,clock);
    for(i=1;i<=N;i++) for(j=1;j<=N;j++) A[i][j]=C[i][j];
    matrix_product(A,B,C,N);
    fprintf(fp,"length of closed path=%d, Trace=%lf \n", clock,Tr);
  }
 fclose(fp);


 /*
 PATH=cmatrix(1,N,1,N);  NEWPATH=cmatrix(1,N,1,N); 

 for(i=1;i<=N;i++) for(j=1;j<=N;j++){ 
   PATH[i][j]=D[i][j]=E[i][j]=c[i][j];
 }


length=length-2;


 sprintf(Loop, "Loop%d.dat",m);
fp=fopen(Loop,"w");

 neighbours=imatrix(1,N,0,N);
for(clock=1; clock<=length; clock++)
  {  
    create_neigh(D, neighbours,N); 
    create_path(E,PATH,NEWPATH,N);            // path 1 -> path 2
    for(i=1;i<=N;i++) for(j=1;j<=N;j++) PATH[i][j]=NEWPATH[i][j];
    for(i=1;i<=N;i++) for(j=1;j<=N;j++) D[i][j]=PATH[i][j];
    loops=0;
    for(i=1;i<=N;i++){
      for(l=1;l<=N;l++){
	if(((PATH[i][l])=='1')&&((c[l][i])=='1')) loops++;
      }
    }
    L=(double)loops/N;
    fprintf(fp, "Loops%d=%d\n", length+2,L);
  }
 fclose(fp);
 */

 free_imatrix(A,1,N,1,N); free_imatrix(B,1,N,1,N); free_imatrix(C,1,N,1,N);

 // free_cmatrix(D,1,N,1,N); free_cmatrix(E,1,N,1,N);
 // free_cmatrix(PATH,1,N,1,N); free_cmatrix(NEW_PATH,1,N,1,N);
}
