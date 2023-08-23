#include <iostream>
#include<math.h>
#include<stdlib.h>
#include<conio.h>
#include<string.h>
using namespace std;
int const Gnmax=20000,Gpmax=14000,ptlimit=15;const float pi=3.14159;
float kf=1.0,rho=1.0,dt,dto,volmin=1.0e5,t=0,Norm=1,Rmax,dt2,dtl,emax=0;
int pmax=3,nmax,Nmax,nn,bmax,change=1,pass=0,R=1,test=0,lmax=0,Count=0;
class vec
{
	public:
	float i,j,Tc,Tf,Tx,Ty,Txx,Tyy,Txy,S,Flux;
	bool coarse,deleted,BN;
	int p1,p2,t[ptlimit+1],level;
	
	vec(){deleted=false;BN=false;level=0;coarse=false;p1=0;p2=0;for(int i=0;i<=9;i++){t[i]=0;}}
	
	void pointmarking();
	void nodetemp();
	void nodederivatives();
	void add(int i);
	void del(int i);
	void del();

	vec operator+(vec b);
	vec operator-(vec b);
	float operator*(vec b);
	float operator%(vec b);
	vec operator*(float b);
	vec operator/(float b);
	friend vec operator*(float b,vec a);
};

class triangle
{
	public:
	int N[3],p[3],midpt[3],type,split,j,level;
	bool midptallocation,mod,coarse,deleted,served,dM;

	vec V,C,A[3],e[3],n[3],f[3];
	float S,Vol,Tf,Tc,To,E1,E2,E21,Ec,Ef,E3c,E3f,Sdc,Sdf;
	float Ixx,Iyy,Ixy,Sxx,Syy,Sxy,Tx,Ty,Txx,Tyy,Txy,Txxx,Txxy,Txyy,Tyyy,Ttc,Ttf;
	float de[3],dn[3],dF[3],dc[3],dr[3],ds[3],Flux;
	
	triangle(){type=0;deleted=false;split=0;coarse=false;mod=false;
	Tx=0;Txx=0;Txxx=0; Ty=0;Tyy=0;Tyyy=0; Txxy=0;Txyy=0;Txy=0;};
	void initialize(int p0,int p1,int p2,int j);

//---------triangle functions-------------
	void splitmarking();//iterative
	void colapsebility();
	void coarsemarking();//iterative
	void midpts();
	void cellmodification();
	void setneighbour();
	void parameters();
	void del();

//---------triangle modifications--------
	void halfsplit();
	void fullsplit();
	void boundrysplit();

	void halfcolapse(int i);
	void fullcolapse();
	void boundrycolapse();
	void modified();
//---------diffusion functions------------

	void errors();
	void derivatives1(int A);
	void derivatives2(int A);

	void trucationerrors(int Ar);
	void heatflux(int A);
	
//---------diffusion functions------------

	void errormark();
	void inputcell(FILE *fp,int jj);
	void timestep();
	void gridtaging();

//-----------output functions--------------
	void plot3D(FILE *fp);
	void plot3DN(FILE *fp);
	void map3D(FILE *fp);
	
	void draw2D(FILE *fp);
	void err(FILE *fp,int Ar);
};
triangle cell[Gnmax+1];
vec P[Gpmax+1];

int newpt();
float gauss(vec P,float t,float &S);
float gauss(vec P,float t);

void triangle::modified(){mod=1;cell[N[0]].mod=1;cell[N[1]].mod=1;cell[N[2]].mod=1;}
float d(vec a){float f;f=sqrt(a.i*a.i+a.j*a.j);return f;}
//-----------**************--------------triangle functions--------------*************-----------------
void triangle::initialize(int p0,int p1,int p2,int jj)
{
 	int swap;
	p[0]=p0;	p[1]=p1;	p[2]=p2;
	N[0]=0;		N[1]=0;		N[2]=0;

	C=(P[p[0]]+P[p[1]]+P[p[2]])/3;
	mod=true;	coarse=false;	served=false;	deleted=false;	j=jj;	
	Vol=(P[p1]-P[p0])%(P[p2]-P[p0]);Vol=Vol/2;
	Tc=C.Tc;
	if(Vol<0) {cout<<"-ve volume cell"<<j<<endl;getch();swap=p[1];p[1]=p[2]; p[2]=swap;Vol=-Vol;}  
	P[p0].add(j);	P[p1].add(j);	P[p2].add(j);
}
//-----------**************--------------triangle functions--------------*************-----------------
void triangle::splitmarking()//iterative
{
	int nk,kk,r,s,pa,c1,c2;cell[0].split=0;
	bool Ar,As,Br,Bs;
	if(split>0&&type==0&&0)//full cell addtitional gaurd
	{ 
		kk=split-1;
		if(kk==0) {r=1;s=2;}	if(kk==1) {r=2;s=0;}	if(kk==2) {r=0;s=1;}	
	
		Ar=cell[N[r]].split>0&&cell[N[r]].type==0;
		As=cell[N[s]].split>0&&cell[N[s]].type==0;
		Br=cell[N[r]].type>0&&cell[N[r]].N[0]!=j;
		Bs=cell[N[s]].type>0&&cell[N[s]].N[0]!=j;

		if(Ar||As||Br||Bs) split=-1;
	}
	if(!served&&split<0&&type==0)//full cell
	{ 
		change=1;served=true;
		for(int k=0;k<3;k++)
		{
			if(cell[N[k]].split>0) cell[N[k]].split=-1;
			if(cell[N[k]].split==0){for(kk=0;kk<3;kk++) if(cell[N[k]].N[kk]==j) cell[N[k]].split=kk+1;}
			if(cell[N[k]].type>0&&cell[N[k]].N[0]==j) cell[N[k]].dM=true;
		}	
	}
	if(!served&&split!=0&&type>0)//half cell
	{ 
		change=1;served=true;
		for(int k=0;k<3;k++)
		{
			if(k!=0)
			{

			if(cell[N[k]].split>0) cell[N[k]].split=-1;
			if(cell[N[k]].split==0){for(kk=0;kk<3;kk++) if(cell[N[k]].N[kk]==j) cell[N[k]].split=kk+1;}
			if(cell[N[k]].type>0&&cell[N[k]].N[0]==j&&k!=0) cell[N[k]].dM=true;
			}
		}	
	}cell[0].split=0;
}
//-----------**************--------------triangle functions--------------*************-----------------
void vec::pointmarking()
{
	bool A,B;
	coarse=true;
	if(p1==0||p2==0) coarse=false;
	if(coarse)
	{
	 	for(int j=1;j<=t[0];j++)
		{
			A=cell[t[j]].coarse==false&&cell[t[j]].type==0;
			B=cell[t[j]].split!=0;
			if(A||B){coarse=false;break;}
		}
	}
}
//-----------**************--------------triangle functions--------------*************-----------------
void triangle::colapsebility()
{
	bool A,B,C,D;
	int p0A,p0B,p1A,p1B,p2A,p2B,E;

		p0A=P[p[0]].p1;		p0B=P[p[0]].p2;
		p1A=P[p[1]].p1;		p1B=P[p[1]].p2;
		p2A=P[p[2]].p1;		p2B=P[p[2]].p2;

		D = P[p[0]].coarse||P[p[1]].coarse||P[p[2]].coarse;
		if(!D) coarse=false;

	if(type==0&&D)
	{
		A = p0A==p1A || p0A==p1B || p0B==p1A || p0B==p1B;
		B = p1A==p2A || p1A==p2B || p1B==p2A || p1B==p2B;
		C = p2A==p0A || p2A==p0B || p2B==p0A || p2B==p0B;
		E = p0A*p1A*p2A;
		
		if(A&&B&&C&&E!=0) coarse=true;else coarse=false;

		if(!coarse)
		{
			if(!A&&!B) P[p[1]].coarse=false;
			if(!B&&!C) P[p[2]].coarse=false;
			if(!C&&!A) P[p[0]].coarse=false;
		}
	}
	if(type==1&&D)	
	{
		if(P[p[2]].coarse) coarse=true;
		else coarse=false;
		P[p[0]].coarse=false;
		P[p[1]].coarse=false;
	}
	if(type==2&&D) 
	{
		coarse=false;
		P[p[0]].coarse=false;
		P[p[2]].coarse=false;
	}
}
//-----------**************--------------triangle functions--------------*************-----------------
void triangle::coarsemarking()//iterative
{
	bool A,B,C;
	int jj,E;
	if(coarse&&type==0)
	{
		E=P[p[0]].coarse+P[p[1]].coarse+P[p[2]].coarse;
		if(E<2){P[p[0]].coarse=false;P[p[1]].coarse=false;P[p[2]].coarse=false;coarse=false;change=1;}

		if(E==2&&0)//additional gaurd
		{
			for(int k=0;k<3;k++)
			{
		 	if(P[p[k]].coarse) 
			{
			for(int i=1;i<=P[p[k]].t[0];i++)
			{
			 	jj=P[p[k]].t[i];
				E=P[cell[jj].p[0]].coarse+P[cell[jj].p[1]].coarse+P[cell[jj].p[2]].coarse;
				if(cell[jj].coarse&&cell[jj].type==0&&jj!=j&&E<3)
				{ 
				P[p[0]].coarse=false;
				P[p[1]].coarse=false;
				P[p[2]].coarse=false;
				coarse=false;
				change=1;break;
				}
			}
			}
			}
		}
	}
	if(coarse&&type==1&&!P[p[2]].coarse){coarse=false;change=1;}
}
//-----------**************--------------triangle functions--------------*************-----------------
void triangle::midpts()
{
	int r,s,nk,kk,pa,c1,c2,n0;
	bool A,B,CC;
	float deno;
	if(split!=0)
	{
	midptallocation=true;
	for(int k=0;k<3;k++)
	{
		if(k==0) {r=1;s=2;}		if(k==1) {r=2;s=0;}		if(k==2) {r=0;s=1;}	

		A=(split==k+1||split==-1)&&type==0;
		B=(k!=type)&&k!=0&&type>0;
		CC=dM&&k!=type;

		if(cell[N[k]].midptallocation&&(A||B||CC))//grab
		{for(kk=0;kk<3;kk++) {if(cell[N[k]].N[kk]==j) midpt[k]=cell[N[k]].midpt[kk];}}

		if((A||B||CC)&&!cell[N[k]].midptallocation)//produce
		{
			pa=newpt();if(pmax>Gpmax) {printf("points exeed limit of %i",Gpmax);getch();exit(0);}
			P[pa]=(P[p[r]]+P[p[s]])/2;
			P[pa].Tc=(P[p[r]].Tc+P[p[s]].Tc)/2;

			deno=1/dF[k]+1/dc[k]+4/dn[k];
//			if(N[k]!=0)P[pa].Tc=(4*P[pa].Tc/dn[k]+cell[N[k]].Tc/dF[k]+Tc/dc[k])/deno;
//			if(N[k]!=0) P[pa].Tc==(map(P[pa]-C,Tc,j)/dF[k]+map(P[pa]-cell[N[k]].C,cell[N[k]].Tc,N[k])/dF[k])/deno;
//			if(N[k]!=0) P[pa].Tc==map(P[pa]-C,Tc,j)/2+map(P[pa]-cell[N[k]].C,cell[N[k]].Tc,N[k])/2;		
			P[pa].p1=p[r];P[pa].p2=p[s];
			
			midpt[k]=pa;
			if(P[p[r]].BN&&P[p[s]].BN&&N[k]==0) P[pa].BN=true;	else P[pa].BN=false;

			if(P[p[r]].level>P[p[s]].level)P[pa].level=P[p[r]].level+1;
			else P[pa].level=P[p[s]].level+1;
		}
	}	
	}
}	
//-----------**************--------------triangle functions--------------*************-----------------
void triangle::cellmodification()
{
	bool A,B,C;

	if(split==-1&&type==0) fullsplit();
	if(split>0&&type==0) halfsplit();
	if(split!=0&&type==1) boundrysplit(); 
	if(coarse&&type==0)
	{
		A=P[p[0]].coarse&&P[p[1]].coarse;
		B=P[p[1]].coarse&&P[p[2]].coarse;
		C=P[p[2]].coarse&&P[p[0]].coarse;
		
	 	if(A&&!B) halfcolapse(2);
	 	if(B&&!C) halfcolapse(0);
	 	if(C&&!A) halfcolapse(1);

	 	if(A&&B&&C)fullcolapse();
	}
	if(coarse&&type==1) boundrycolapse();
}

//-----------**************--------------triangle functions--------------*************-----------------
//-----------@@@@@@@@@@@@@@-*************@@@@@@@@@@@@@@@@@@@@***********@@@@@@@@@@@@@@-----------------
float map(vec P,float T,int j)
{
 	float x=P.i,y=P.j,Trec=0;
	Trec=cell[j].Txxx*x*x*x/6+cell[j].Tyyy*y*y*y/6+Trec;
	Trec=cell[j].Txxy*x*x*y/2+cell[j].Txyy*x*y*y/2+Trec;

	Trec=cell[j].Txx*x*x/2+cell[j].Tyy*y*y/2+cell[j].Txy*x*y+Trec;
	Trec=T+cell[j].Tx*x+cell[j].Ty*y+Trec;
 	return Trec;
}
//-----------@@@@@@@@@@@@@@-************triangle modifications***********@@@@@@@@@@@@@-----------------
void triangle::halfsplit()
{
 	modified();vec Co=C;
	R=1;float Trec=Tc,Tr,Ts;
	int k,r,s,pk,pr,ps,mk,Bn=P[p[0]].BN+P[p[1]].BN+P[p[2]].BN;
	if(!dM) k=split-1;
    if(dM) k=0;
	if(k==0) {r=1;s=2;}		if(k==1) {r=2;s=0;}		if(k==2) {r=0;s=1;}	
	
	Tr=(P[p[k]].Tc+2*P[p[r]].Tc)/3;
	Ts=(P[p[k]].Tc+2*P[p[s]].Tc)/3;

	pk=p[k];	pr=p[r];	ps=p[s];	mk=midpt[k];
	P[ps].del(j);

	nmax++;
	if(nmax>Gnmax) {printf("triangles exeed limit");getch();exit(0);}
	cell[nmax].initialize(pk,mk,ps,nmax);
	cell[nmax].type=2;
	cell[nmax].split=0;

//	if(Bn==0) cell[nmax].Tc=map(cell[nmax].C-C,Tc,j);	
	cell[nmax].Tc=Trec/2+Ts/2;

	initialize(pk,pr,mk,j);
//	if(Bn==0) Tc=map(C-Co,Trec,j);	
	Tc=Trec/2+Tr/2;
	type=1;split=0;
}
//-------------@@@@@@@@@@@@-------------triangle modifications-----------@@@@@@@@@@--------------------
void triangle::fullsplit()
{
 	R=1;
 	modified();
	int Bn=P[p[0]].BN+P[p[1]].BN+P[p[2]].BN;
	int k,r,s,m0=midpt[0],m1=midpt[1],m2=midpt[2];float Trec=Tc;
	P[p[0]].del(j);		P[p[1]].del(j);		P[p[2]].del(j);
	for(k=0;k<3;k++)
	{
		if(k==0) {r=1;s=2;}		if(k==1) {r=2;s=0;}		if(k==2) {r=0;s=1;}	

		nmax++;   
		if(nmax>Gnmax) {printf("triangles exeed limit");getch();exit(0);}

		cell[nmax].initialize(p[k],midpt[s],midpt[r],nmax);	
		cell[nmax].split=0;	
		cell[nmax].Tc=Trec/2+P[p[k]].Tc/2;
	}
	initialize(m0,m1,m2,j);Tc=Trec;
	type=0;split=0;
}
//-------------@@@@@@@@@@@@-------------triangle modifications-----------@@@@@@@@@@--------------------
void triangle::boundrysplit()
{	
	R=1;
 	modified();cell[N[1]].modified();
	int Bn=P[p[0]].BN+P[p[1]].BN+P[p[2]].BN;
	int c2=N[1],p0=p[0],p1=p[1],ms=midpt[2],mr=cell[c2].midpt[1],mk=p[2];	
	float Tr=Tc,Ts=cell[c2].Tc;
	vec Cr=C,Cs=cell[c2].C;

	P[p[0]].del(j);			P[p[0]].del(c2);
	initialize(ms,p[1],p[2],j);type=0;
//	Tc=map(C-Cr,Tr,j);	
	Tc=(Tr+Ts+2*P[p[1]].Tc)/4;
	
	cell[c2].initialize(mr,cell[c2].p[1],cell[c2].p[2],c2);
//	if(Bn==0) cell[c2].Tc=map(cell[c2].C-Cs,Ts,c2);	
	cell[c2].Tc=(Tr+Ts+2*P[cell[c2].p[2]].Tc)/4;
	cell[c2].type=0;

	//----------------**********************----------------------------------------	
	nmax++;if(nmax>Gnmax) {printf("triangles exeed limit");getch();exit(0);}
	cell[nmax].initialize(mk,mr,ms,nmax);
	cell[nmax].type=0;
	cell[nmax].split=0;
//	if(Bn==0) cell[nmax].Tc=map(cell[nmax].C-Cr,Tr,j)/2+map(cell[nmax].C-Cs,Ts,c2)/2;
	cell[nmax].Tc=Tr/2+Ts/2;

	nmax++;if(nmax>Gnmax) {printf("triangles exeed limit");getch();exit(0);}
	cell[nmax].initialize(p0,ms,mr,nmax);
	cell[nmax].type=0;
	cell[nmax].split=0;
//	if(Bn==0) cell[nmax].Tc=map(cell[nmax].C-Cr,Tr,j)/2+map(cell[nmax].C-Cs,Ts,c2)/2;
	cell[nmax].Tc=P[p0].Tc/2+cell[nmax-1].Tc/2;

	//----------------**********************----------------------------------------	
	if(dM) 
	{
		halfsplit();
//		if(Bn==0) Tc=map(C-Cr,Tr,j);
//		if(Bn==0) cell[nmax].Tc=map(cell[nmax].C-Cr,Tr,j);
	}
	if(cell[c2].dM) 
	{
		cell[c2].halfsplit();
//		if(Bn==0) cell[c2].Tc=map(cell[c2].C-Cs,Ts,c2);
//		if(Bn==0) cell[nmax].Tc=map(cell[nmax].C-Cs,Ts,c2);
	} 
	split=0;cell[c2].split=0;
}
//-------------@@@@@@@@@@@@-------------triangle modifications-----------@@@@@@@@@@--------------------
void triangle::halfcolapse(int k)
{
	R=1;modified();cell[N[0]].modified();cell[N[1]].modified();cell[N[2]].modified();
	int r,s,pt[3],m,n;
	int Bn=P[p[0]].BN+P[p[1]].BN+P[p[2]].BN;
	if(k==0){r=1;s=2;}	if(k==1){r=2;s=0;}	if(k==2){r=0;s=1;}	

	float Tr=cell[N[r]].Tc,Ts=cell[N[s]].Tc,Tk=cell[N[k]].Tc;
	for(m=0;m<3;m++)
	{for(n=0;n<3;n++){if(cell[N[m]].N[n]==j) pt[m]=cell[N[m]].p[n];}}

	cell[N[r]].initialize(pt[k],pt[r],p[k],N[r]);cell[N[r]].type=1;	
	cell[N[s]].initialize(pt[k],p[k],pt[s],N[s]);cell[N[s]].type=2;	

//	if(Bn==0)cell[N[r]].Tc=map(cell[N[r]].C-C,Tc,j);
//	if(Bn==0)cell[N[s]].Tc=map(cell[N[s]].C-C,Tc,j);

	cell[N[r]].Tc=(2*Tr+Tk)/3;
	cell[N[s]].Tc=(2*Ts+Tk)/3;

	P[pt[k]].del(N[k]);
	P[p[k]].del(j);
	cell[N[k]].del();
	del();
}
//-------------@@@@@@@@@@@@-------------triangle modifications-----------@@@@@@@@@@--------------------
void triangle::fullcolapse()
{
 	R=1;modified();cell[N[0]].modified();cell[N[1]].modified();cell[N[2]].modified();
	int Bn=P[p[0]].BN+P[p[1]].BN+P[p[2]].BN;
	float Trec=Tc;
	int pt[3],k,kk;
	for(k=0;k<3;k++) 
	{
		for(kk=0;kk<3;kk++) if(cell[N[k]].N[kk]==j) pt[k]=cell[N[k]].p[kk];
		P[pt[k]].del(N[k]);
		cell[N[k]].del();
	}
	initialize(pt[0],pt[1],pt[2],j);Tc=Trec;
}
//-------------@@@@@@@@@@@@-------------triangle modifications-----------@@@@@@@@@@--------------------
void triangle::boundrycolapse()
{
	R=1;modified();cell[N[1]].modified();
	int n1=N[1],p2=cell[N[1]].p[2];
	int Bn=P[p[0]].BN+P[p[1]].BN+P[p[2]].BN;
	float Tr=Tc,Ts=cell[n1].Tc;
	vec Cr=C,Cs=cell[n1].C;
	
	P[p2].del(n1);
	P[p[0]].del(n1);
	initialize(p[0],p[1],p2,j);type=0;
//	if(Bn==0)Tc=map(C-Cr,Tr,j)/2+map(C-Cs,Ts,n1)/2;
//	else 
	Tc=Tr/2+Ts/2;
	cell[n1].del();
}
//-------------@@@@@@@@@@@@-------------triangle modifications-----------@@@@@@@@@@@@@-----------------
//-------------@@@@@@@@@@@@@*************@@@@@@@@@@@@@@@@@@@@***********@@@@@@@@@@@@@@-----------------
//-----------**************--------------triangle functions--------------*************-----------------
void triangle::setneighbour()
{
 	int k,kk,r,s,i,l;
 	bool found;
	if(!deleted)
	{
		for(k=0;k<3;k++)
		{
			if(k==0) {r=1;s=2;}		if(k==1) {r=2;s=0;}		if(k==2) {r=0;s=1;}	
	 
			for(i=1;i<=P[p[r]].t[0];i++)
			{for(l=1;l<=P[p[s]].t[0];l++)
			{found=P[p[r]].t[i]==P[p[s]].t[l]&&P[p[s]].t[l]!=j;
			if(found) break;}
			if(found) break;}
			if(found) N[k]=P[p[r]].t[i];
			if(!found&&P[p[r]].BN&&P[p[s]].BN) N[k]=0;
			
			if(!found&&(!P[p[r]].BN||!P[p[s]].BN)) 
			{printf("triangle %i neibour error at N[%i]=%i %i po=%i",j,k,N[k],deleted,P[p[0]].deleted);
			printf(" p1=%i p2=%i.\n",P[p[1]].deleted,P[p[2]].deleted);getch();};
		}
	}
}	
//-----------**************--------------triangle functions--------------*************-----------------
void triangle::parameters()
{
	int r,s,k,nk,l;
	bool Ao,Bo;
	float X,Y,x0,y0,x1,y1,x2,y2,xi,yi,I0,I2;
	vec fk;
	for(k=0;k<3;k++)
	{	
		if(k==0) {r=1;s=2;}    if(k==1) {r=2;s=0;}    if(k==2) {r=0;s=1;}
		Ao=P[p[k]].i>=P[p[r]].i&&P[p[k]].i<=P[p[s]].i;
		Bo=P[p[k]].i<=P[p[r]].i&&P[p[k]].i>=P[p[s]].i;
		if(Ao||Bo) {if(Bo){nk=r;r=s;s=nk;}break;}
	}if(k==3) k--;	

	x0=P[p[r]].i-C.i;		x1=P[p[k]].i-C.i;		x2=P[p[s]].i-C.i;
	y0=P[p[r]].j-C.j;		y1=P[p[k]].j-C.j;		y2=P[p[s]].j-C.j;
	yi=y0*(x2-x1)/(x2-x0)+y2*(x1-x0)/(x2-x0);if(fabs(x2-x0)<1e-8) yi=y0/2+y1/2; 
	
	X=x1-x0;Y=y1-yi;	I0=Y*X*(X*X/4+2*x0*X/3+x0*x0/2);//x1>x0
	X=x1-x2;Y=y1-yi;	I2=Y*X*(X*X/4+2*x2*X/3+x2*x2/2);//x1>x2
	Ixx=fabs(I0-I2);	

	X=x1-x0;Y=y1-yi;	I0=y0*Y*X*(X/3+x0)+Y*(y1+yi)*X*(X/4+x0)/2;
	X=x1-x2;Y=y1-yi;	I2=y2*Y*X*(X/3+x2)+Y*(y1+yi)*X*(X/4+x2)/2;
	Ixy=I0-I2;if(Y<0) Ixy=-Ixy;

	for(k=0;k<3;k++)
	{	
		if(k==0) {r=1;s=2;}    if(k==1) {r=2;s=0;}    if(k==2) {r=0;s=1;}
		Ao=P[p[k]].j>=P[p[r]].j&&P[p[k]].j<=P[p[s]].j;
		Bo=P[p[k]].j<=P[p[r]].j&&P[p[k]].j>=P[p[s]].j;
		if(Ao||Bo) {if(Bo){nk=r;r=s;s=nk;}break;}
	}if(k==3) k--;	
 
	x0=P[p[r]].i-C.i;		x1=P[p[k]].i-C.i;		x2=P[p[s]].i-C.i;
	y0=P[p[r]].j-C.j;		y1=P[p[k]].j-C.j;		y2=P[p[s]].j-C.j;
	xi=x0*(y2-y1)/(y2-y0)+x2*(y1-y0)/(y2-y0);if(fabs(y2-y0)<1e-8) xi=x0/2+x1/2; 

	X=x1-xi;Y=y1-y0;	I0=Y*X*(Y*Y/4+2*y0*Y/3+y0*y0/2);//x1>x0
	X=x1-xi;Y=y1-y2;	I2=Y*X*(Y*Y/4+2*y2*Y/3+y2*y2/2);//x1>x2
	Iyy=fabs(I0-I2);

//	I0=(y1-yi)*((x1*x1+x0*x0)*(x1+x0)/4-x0*(x1*x1+x0*x0+x1*x0)/3);
//	printf("Ixx=%G Iyy=%G Ixy=%G\n",Ixx,Iyy,Ixy);getch();

	vec a,c,b,B,Xo,Yo;float tt;
	for(int k=0;k<=2;k++)
	{
	    if(k==0){r=1;s=2;}	if(k==1) {r=2;s=0;}	if(k==2){r=0;s=1;}
		nk=N[k];if(nk==0) cell[nk].C = (P[p[s]]+P[p[r]])/2;
		
		A[k].i = P[p[s]].j - P[p[r]].j;
		A[k].j = P[p[r]].i - P[p[s]].i;
	
		de[k]= d(cell[nk].C-C);
		dn[k]= d(P[p[s]]-P[p[r]]);

		e[k] = (cell[nk].C - C)/de[k];
		n[k] = (P[p[s]]-P[p[r]])/dn[k];
	
		tt=((P[p[s]]-C)%n[k])/(e[k]%n[k]);
		f[k]=C+tt*e[k];

		dc[k]=d(C-f[k]);		dF[k]=de[k]-dc[k];
		dr[k]=d(P[p[r]]-f[k]);	ds[k]=dn[k]-dr[k];

	}
}

//-----------**************--------------triangle functions--------------*************-----------------
//-----------**************---------------Tempreture functions---------------*************-----------------
void triangle::heatflux(int Ar)
{		
	int r,s,nk,kk,Bn;
	float Te,Tn,ap,w,Ri,tx,ty;if(Ar<3)Flux=0;if(Ar==3)E2=0;
	Bn = P[p[0]].BN + P[p[1]].BN + P[p[2]].BN;
	for(kk=0;kk<3;kk++) {if(P[p[kk]].BN&&Bn==1) break;if(!P[p[kk]].BN&&Bn==2) break;if(!N[kk]!=0&&Bn==3) break;}

	for(int k=0;k<=2;k++)//Jf.Af calculations.
	{
		if(k==0) {r=1;s=2;}		if(k==1) {r=2;s=0;}		if(k==2) {r=0;s=1;}
		w=e[k]*n[k];		nk=N[k];

		if(Ar==1)
		{
			if(nk==0) cell[nk].Tc=(P[p[s]].Tc+P[p[r]].Tc)/2;
			Te=(cell[nk].Tc-Tc)/de[k];
			Tn=(P[p[s]].Tc-P[p[r]].Tc)/dn[k];

		}
		
		if(Ar==2)
		{
			if(nk==0) cell[nk].Tf=(P[p[s]].Tf+P[p[r]].Tf)/2;
			Te=(cell[nk].Tf-Tf)/de[k];
			Tn=(P[p[s]].Tf-P[p[r]].Tf)/dn[k];
		}

		Flux=Flux+(Te-w*Tn)*e[k]*A[k]/(1-w*w);
	}Flux=kf*Flux/Vol;//Flux=Flux-E3c;
	Ri=(Tf-Tc)/dt-Flux-S;

	if(Rmax<Ri)Rmax=Ri;
	if(Ar==1){Tf=Tc+(Flux+S)*dt;//if(Tf>150) printf("Tf=%f\n",Tf);Ttc=Flux+S;
	}
	if(Ar==2){Ttf=Flux+S; E1=(Ttf-Ttc)*.5/dt;}
//	if(Bn>0)E2=0;
}

/*void triangle::errors(int Ar)
{
	float Ac[3][5],err,w,xe,ye,xn,yn,Te,Tn,Te2,Tn2,axc[3],ayc[3];
	int nk,r,s;
	for(int k=0;k<3;k++)
	{
		if(k==0) {r=1;s=2;}		if(k==1) {r=2;s=0;}		if(k==2) {r=0;s=1;}	
		w=e[k]*n[k];		nk=N[k];

		Ac[k][0]=1;		Ac[k][1]=f[k].i-C.i;		Ac[k][2]=f[k].j-C.j;	

		if(Ar==1)
		{
			if(nk==0) cell[nk].Tc=(P[p[s]].Tc+P[p[r]].Tc)/2;
			Te=(cell[nk].Tc-Tc)/de[k];
			Tn=(P[p[s]].Tc-P[p[r]].Tc)/dn[k];
	
			Ac[k][3]=((Te-w*Tn)*e[k].i+(Tn-w*Te)*n[k].i)/(1-w*w);
			Ac[k][4]=((Te-w*Tn)*e[k].j+(Tn-w*Te)*n[k].j)/(1-w*w);
		}
		if(Ar==2)
		{
			if(nk==0) cell[nk].Tf=(P[p[s]].Tf+P[p[r]].Tf)/2;
			Te=(cell[nk].Tf-Tf)/de[k];
			Tn=(P[p[s]].Tf-P[p[r]].Tf)/dn[k];
	
			Ac[k][3]=((Te-w*Tn)*e[k].i+(Tn-w*Te)*n[k].i)/(1-w*w);
			Ac[k][4]=((Te-w*Tn)*e[k].j+(Tn-w*Te)*n[k].j)/(1-w*w);
		}
	}solution(Ac,axc,ayc);

	Tx=axc[0];Ty=ayc[0];Txx=axc[1];Tyy=ayc[2];Txy=(axc[2]+ayc[1])/2;
	if(Ar==1){Ec = Ixx*Txx/2+ Iyy*Tyy/2+Ixy*Txy;	Ec = Ec/Vol;}
	if(Ar==2){Ef = Ixx*Txx/2+ Iyy*Tyy/2+Ixy*Txy;	Ef = Ef/Vol;	E2=(Ec-Ef)/dt;}
}
*/
//-----------**************--------------error functions-------------*************-----------------
void triangle::trucationerrors(int Ar)
{
	int nk,r,s;
	float xe,ye,xn,yn,deno,Cxx,Cyy,Cxy,Fxx,Fyy,Fxy,Pxx,Pyy,Pxy,Tee,Tnn,Ten,w,Ce,Cn,E3=0;//E3c=0;E3f=0;
	for(int k=0;k<3;k++)
	{
		if(k==0) {r=1;s=2;}		if(k==1) {r=2;s=0;}		if(k==2) {r=0;s=1;}	

		xe=e[k].i;	ye=e[k].j;	xn=n[k].i;	yn=n[k].j;

		Tee=(cell[N[k]].Tx-Tx)*xe	+	(cell[N[k]].Ty-Ty)*ye;			Tee=Tee/de[k];
		Tnn=(P[p[s]].Tx-P[p[r]].Tx)*xn+	(P[p[s]].Ty-P[p[r]].Ty)*yn;	Tnn=Tnn/dn[k];

/*		deno=(1/dc[k]+1/dF[k]);
		Cxx=Txx;			Cyy=Tyy;			Cxy=Txy;
		Fxx=cell[N[k]].Txx;	Fyy=cell[N[k]].Tyy;	Fxy=cell[N[k]].Txy;
		
		if(N[k]!=0)
		{
			Pxx=(Cxx/dc[k]+Fxx/dF[k])/deno;	
			Pyy=(Cyy/dc[k]+Fyy/dF[k])/deno;
			Pxy=(Cxy/dc[k]+Fxy/dF[k])/deno;
		}	
		else {Pxx=Txx;	Pyy=Tyy;	Pxy=Txy;}

		//Te=xeTx+yeTy;
	 	Tee=xe*xe*Pxx+ye*ye*Pyy+2*xe*ye*Pxy;
		Tnn=xn*xn*Pxx+yn*yn*Pyy+2*xn*yn*Pxy;
		Ten=xe*xn*Pxx+ye*yn*Pyy+(xe*yn+ye*xn)*Pxy;
*/
		w=e[k]*n[k];Ce=e[k]*A[k]/(1-w*w);Cn=-w*Ce;
		if(N[k]!=0)
		{
	  		E3=E3+Ce*(dF[k]-dc[k])*Tee/2+Cn*(ds[k]-dr[k])*Tnn/2;
			E3=E3+(ds[k]-dr[k])*(Ce*Ten+Cn*Tnn)/2;
		}
	}
	E3 = kf*E3/Vol;
	if(Ar==1)E3c=E3;if(Ar==2)E3f=E3;
	
/*	vec a,dv;float EE,tt,Sd,Sxf,Sxb,Syf,Syb,Sff,Sbb,Sfb,Sbf;
	if(Ar==1)tt=t;if(Ar==2)tt=t+dt;

	dv.i=.01;dv.j=0;a=C+dv;gauss(a,tt,Sxf);		a=C-dv;gauss(a,tt,Sxb);
	dv.i=0;dv.j=.01;a=C+dv;gauss(a,tt,Syf);		a=C-dv;gauss(a,tt,Syb);

	dv.i=.01;dv.j=.01;a=C+dv;gauss(a,tt,Sff);	a=C-dv;gauss(a,tt,Sbb);
	dv.i=.01;dv.j=-.01;a=C+dv;gauss(a,tt,Sfb);	a=C-dv;gauss(a,tt,Sbf);

	Sxx=(Sxf+Sxb-2*S)/1e-4;	
	Syy=(Syf+Syb-2*S)/1e-4;	
	Sxy=(Sff+Sbb-Sbf-Sfb)/4e-4;

	Sd = Ixx*Sxx/2 + Iyy*Syy/2 + Ixy*Sxy;		Sd = Sd/Vol;
//	EE = Ixx*Txx/2 + Iyy*Tyy/2 + Ixy*Txy;		EE = EE/Vol;

//	Bn=P[p[0]].BN+P[p[1]].BN+P[p[2]].BN;
	if(Ar==1) Sdc=Sd;
//	if(Ar==2){Ef=EE;E2=Ef-Ec;Sdf=Sd;	if(Bn>0)E2=0;}*/
}
/*void triangle::derivatives1(int A)
{
	float x,y,xx=0,xy=0,yy=0,Tcx1=0,Tcy1=0,A1[2][3];

	for(int k=0;k<3;k++)
	{
	 	x=cell[N[k]].C.i-C.i;y=cell[N[k]].C.j-C.j;
		xx=xx+x*x;	xy=xy+x*y;	yy=yy+y*y;
		Tcx1=Tcx1+cell[N[k]].Tc*x;
		Tcy1=Tcy1+cell[N[k]].Tc*y;
	}
	A1[0][0]=xx;	A1[0][1]=xy;	A1[0][2]=Tcx1;
	A1[1][0]=xy;	A1[1][1]=yy;	A1[1][2]=Tcy1;
	Tx=(Tcx1*yy-Tcy1*xy)/(xx*yy-xy*xy);		
	Ty=(Tcy1*xx-Tcx1*xy)/(xx*yy-xy*xy);
}
*/
void vec::nodetemp()
{
	float S=0,Sf=0,Sx=0,Bn;int pt;
	for(int j=1;j<=t[0];j++) 
	{
		pt=t[j];Bn=P[cell[pt].p[0]].BN+P[cell[pt].p[1]].BN+P[cell[pt].p[2]].BN;

	  	Sf=Sf+cell[pt].Tf/d(*this-cell[pt].C);Sx=Sx+cell[pt].Flux/d(*this-cell[pt].C);S=S+1/d(*this-cell[pt].C);

// 		Sf=Sf+cell[pt].Tf*cell[pt].Vol;	Sx=Sx+cell[pt].Flux*cell[pt].Vol;S=S+cell[pt].Vol;
	}Tf=Sf/S;Flux=Sx/S;
}
void triangle::derivatives1(int A)
{
	vec n,a,b;
	a=P[p[1]]-P[p[0]]; a.Flux=P[p[1]].Flux-P[p[0]].Flux; 			
	b=P[p[2]]-P[p[0]]; b.Flux=P[p[2]].Flux-P[p[0]].Flux; 			

	if(A==1){n.i=a.j*b.Tc-a.Tc*b.j;		n.j=a.Tc*b.i-a.i*b.Tc;}
	if(A==2){n.i=a.j*b.Tf-a.Tf*b.j;		n.j=a.Tf*b.i-a.i*b.Tf;}	
	if(A==3){n.i=a.j*b.Flux-a.Flux*b.j;	n.j=a.Flux*b.i-a.i*b.Flux;}	

	n.Tc=a.i*b.j-a.j*b.i;		

	Tx=-n.i/n.Tc;
	Ty=-n.j/n.Tc;
}

void triangle::errors()
{
	float Trec;int nk;E2=0;
	for(int k=0;k<3;k++)
	{
		nk=N[k];
//		Trec=(P[cell[nk].p[0]].Flux+P[cell[nk].p[1]].Flux+P[cell[nk].p[2]].Flux)/3;
		Trec=(P[cell[nk].p[0]].Tc+P[cell[nk].p[1]].Tc+P[cell[nk].p[2]].Tc)/3;

		Trec = Trec + (P[p[k]].i-cell[nk].C.i)*cell[nk].Tx;
		Trec =          Trec + (P[p[k]].j-cell[nk].C.j)*cell[nk].Ty;

		if(fabs(P[p[k]].Tc-Trec)>E2&&nk!=0) E2=fabs(P[p[k]].Tc-Trec);
//		if(fabs(P[p[k]].Flux-Trec)>E2&&nk!=0) E2=fabs(P[p[k]].Flux-Trec);

		if(N[k]==0)cell[N[k]].Vol=0;
	}
//	E21=(Flux*Vol+cell[N[0]].Flux*cell[N[0]].Vol+cell[N[1]].Flux*cell[N[1]].Vol+cell[N[2]].Flux*cell[N[2]].Vol);
//	E21=E21/(Vol+cell[N[0]].Vol+cell[N[1]].Vol+cell[N[2]].Vol);
//	E21=(Flux-E21);
	E21=Flux-(P[p[0]].Flux+P[p[1]].Flux+P[p[2]].Flux)/3;E21=E21*3/4;
}

void vec::nodederivatives()
{
	float S=0,Sx=0,Sy=0;int pt;
	for(int jj=1;jj<=t[0];jj++) 
	{
		pt=t[jj];
		Sx=Sx+cell[pt].Tx/d(*this-cell[pt].C);
		Sy=Sy+cell[pt].Ty/d(*this-cell[pt].C);
		S=S+1/d(*this-cell[pt].C);
	}Tx=Sx/S;Ty=Sy/S;

	if(i==0||i==1)Ty=0;
	if(j==0||j==1)Tx=0;
}

void solution(float A[][4],float x[])
{
	float D0=0,D1=0,D2=0,D=0;
	
	D=D + A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]);
	D=D + A[0][1]*(A[1][2]*A[2][0]-A[1][0]*A[2][2]);
	D=D + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
	
	D0=D0+A[0][3]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]);
	D0=D0+A[0][1]*(A[1][2]*A[2][3]-A[1][3]*A[2][2]);
	D0=D0+A[0][2]*(A[1][3]*A[2][1]-A[1][1]*A[2][3]);
	
	D1=D1+A[0][0]*(A[1][3]*A[2][2]-A[1][2]*A[2][3]);
	D1=D1+A[0][3]*(A[1][2]*A[2][0]-A[1][0]*A[2][2]);
	D1=D1+A[0][2]*(A[1][0]*A[2][3]-A[1][3]*A[2][0]);
	
	D2=D2+A[0][0]*(A[1][1]*A[2][3]-A[1][3]*A[2][1]);
	D2=D2+A[0][1]*(A[1][3]*A[2][0]-A[1][0]*A[2][3]);
	D2=D2+A[0][3]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
	x[0]=D0/D;x[1]=D1/D;x[2]=D2/D;
}
void triangle::derivatives2(int A)
{
	vec n,no,a,b,dv;float Tyx,EE,tt,Sd;
	a=P[p[1]]-P[p[0]];	b=P[p[2]]-P[p[0]]; 			

	n.i=a.j*b.Tx-a.Tx*b.j;		
	n.j=a.Tx*b.i-a.i*b.Tx;		
	n.Tx=a.i*b.j-a.j*b.i;		

	Txx=-n.i/n.Tx;
	Txy=-n.j/n.Tx;

	n.i=a.j*b.Ty-a.Ty*b.j;		
	n.j=a.Ty*b.i-a.i*b.Ty;		
	n.Ty=a.i*b.j-a.j*b.i;		

	Tyx=-n.i/n.Ty;
	Tyy=-n.j/n.Ty;
	Txy=(Txy+Tyx)/2;

	if(A==1)tt=t;if(A==2)tt=t+dt;
	float Sxf,Sxb,Syf,Syb,Sff,Sbb,Sfb,Sbf;

	dv.i=.01;dv.j=0;a=C+dv;gauss(a,tt,Sxf);		a=C-dv;gauss(a,tt,Sxb);
	dv.i=0;dv.j=.01;a=C+dv;gauss(a,tt,Syf);		a=C-dv;gauss(a,tt,Syb);

	dv.i=.01;dv.j=.01;a=C+dv;gauss(a,tt,Sff);	a=C-dv;gauss(a,tt,Sbb);
	dv.i=.01;dv.j=-.01;a=C+dv;gauss(a,tt,Sfb);	a=C-dv;gauss(a,tt,Sbf);

	Sxx=(Sxf+Sxb-2*S)/1e-4;	
	Syy=(Syf+Syb-2*S)/1e-4;	
	Sxy=(Sff+Sbb-Sbf-Sfb)/4e-4;

	Sd = Ixx*Sxx/2 + Iyy*Syy/2 + Ixy*Sxy;		Sd = Sd/Vol;
	EE = Ixx*Txx/2 + Iyy*Tyy/2 + Ixy*Txy;		EE = EE/Vol;

	if(A==1) {Ec=EE;Sdc=Sd;}
	if(A==2){Ef=EE;E2=Ef-Ec;Sdf=Sd;}
}


/*
void triangle::derivatives2(int A)
{
	float x,y,xx=0,xy=0,yy=0,Txx1=0,Txy1=0,Tyx1=0,Tyy1=0,A1[3][4],x1[3];
	for(int k=0;k<3;k++)
	{
	 	x=cell[N[k]].C.i-C.i;y=cell[N[k]].C.j-C.j;
		xx=xx+x*x;	xy=xy+x*y;	yy=yy+y*y;
		Txx1=Txx1+cell[N[k]].Tx*x;
		Tyx1=Tyx1+cell[N[k]].Ty*x;
		Txy1=Txy1+cell[N[k]].Tx*y;
		Tyy1=Tyy1+cell[N[k]].Ty*y;
	}
	
	A1[0][0]=xx;	A1[0][1]=xy;		A1[0][2]=0;		A1[0][3]=Txx1;
	A1[1][0]=xy;	A1[1][1]=xx+yy;		A1[1][2]=xy;	A1[1][3]=Txy1+Tyx1;
	A1[2][0]=0;		A1[2][1]=xy;		A1[2][2]=yy;	A1[2][3]=Tyy1;
	solution(A1,x1);Txx=x1[0];			Txy=x1[1];		Tyy=x1[2];
}
void triangle::derivatives2(int A)
{
	float A1[3][4],x1[3];
	for(int k=0;k<3;k++)
	{
	 	A1[k][0]=e[k].i*e[k].i;A1[k][1]=e[k].j*e[k].j;A1[k][2]=2*e[k].i*e[k].j;
		A1[k][3]=(cell[N[k]].Tx-Tx)*e[k].i/de[k]+(cell[N[k]].Ty-Ty)*e[k].j/de[k];
		printf("Txx=%f Tyy=%f\t",(cell[N[k]].Tx-Tx)*e[k].i,(cell[N[k]].Ty-Ty)*e[k].j);
	}getch();printf("\n");
	solution(A1,x1);Txx=x1[0];Tyy=x1[1];Txy=x1[2];
if(N[0]==0||N[1]==0||N[2]==0){Txx=0;Tyy=0;Txy=0;}
}
*/

void errors(int A)
{
 	for(int j=1;j<=nmax;j++) cell[j].derivatives1(A);
	for(int i=1;i<=pmax;i++) P[i].nodederivatives();
	for(int j=1;j<=nmax;j++) cell[j].errors();
}
void derivatives(int A)
{
	for(int j=1;j<=nmax;j++) cell[j].derivatives2(A);
}

//-----------**************--------------diffusion functions-------------*************---
//-------------@@@@@@@@@@@@@*************@@@@@@@@@@@@@@@@@@@@***********@@@@@@@@@@@@@@---
//-----------**************---------------general functions--------------*************-----------------
void reinitialize()
{
	for(int j=0;j<=Nmax;j++) 
	{
		cell[j].split=0;
		cell[j].served=false;
		cell[j].coarse=false;
		cell[j].mod=false;
		cell[j].midptallocation=false;
		cell[j].dM=false;
	}Nmax=nmax;
}
//-----------**************---------------general functions--------------*************-----------------
void renumber()
{
	int p,nk,r,s;
	while(cell[nmax].deleted) {nmax--;}
	for(int j=1;j<nmax;j++) 
	{
		if(cell[j].deleted)
		{
			cell[j]=cell[nmax];
			cell[j].modified();
			cell[j].j=j;
			cell[j].deleted=false;

			for(int k=0;k<3;k++)
			{
				p=cell[j].p[k];
				P[p].del(nmax);
				P[p].add(j);
			}cell[nmax].del();
			while(cell[nmax].deleted) {nmax--;}
		}
	}
	/*while(P[pmax].coarse) {P[pmax].del();pmax--;}
	for(int i=1;i<=pmax;i++) 
	{
		if(P[i].coarse)
		{
			P[i]=P[pmax];
			for(int j=1;j<=P[i].t[0];j++)
			{
		 	nk=P[i].t[j];
 			for(int k=0;k<3;k++)
			{
				if(cell[nk].p[k]==pmax) 
				{
					if(k==0){r=1;s=2;}	if(k==1){r=2;s=0;}	if(k==2){r=0;s=1;}	
		    		cell[nk].p[k]=i;
					if(P[cell[nk].p[r]].p1==pmax)P[cell[nk].p[r]].p1=i;
					if(P[cell[nk].p[r]].p2==pmax)P[cell[nk].p[r]].p2=i;

					if(P[cell[nk].p[s]].p1==pmax)P[cell[nk].p[s]].p1=i;
					if(P[cell[nk].p[s]].p2==pmax)P[cell[nk].p[s]].p2=i;

				}
			}
			}
			P[pmax].coarse=true;
			while(P[pmax].coarse) {P[pmax].del();pmax--;}
			P[i].del();
		}
	}*/
}


//-----------**************---------------general functions--------------*************-----------------
//-------------@@@@@@@@@@@@@*************@@@@@@@@@@@@@@@@@@@@***********@@@@@@@@@@@@@@-----------------
void vec::add(int i)
{
	bool found=false;
	if(t[0]==ptlimit)
	{
		printf("error .Triangle %i exceded point limit of 12.\n",i);
		getch();exit(1);
	}
	for(int j=1;j<=t[0];j++)if(i==t[j]) found=true;
	if(i==0) found=true;
	if(!found){t[0]++;t[t[0]]=i;}

}
//-----------**************---------------vector functions---------------*************-----------------
void vec::del(int i)
{
	bool found=false;
	for(int j=1;j<=t[0];j++)
	{
		if(t[j]==i) found=true;
		if(found&&j!=t[0]) t[j]=t[j+1];
	}
	if(!found)
	{
	for(int j=1;j<=t[0];j++){printf("\n%i ",t[j]);}

		printf("\nerror .Triangle %i not found in triangle removal from point\n",i);
		printf("\nmax=%i\tpmax=%i",nmax,pmax);
		
		getch();
		exit(1);
	}
	t[t[0]]=0;t[0]--;
}
void triangle::del()
{
	deleted=true;coarse=false;
	split=0;type=0;
	mod=false;dM=false;
	midptallocation=false;
	served=true;
	modified();
	
//	P[p[0]].del(j);P[p[1]].del(j);P[p[2]].del(j);

	N[0]=0;N[1]=0;N[2]=0;
	p[0]=0;p[1]=0;p[2]=0;
}
	
//-----------**************---------------vector functions---------------*************-----------------
void vec::del()
{
	deleted=true;
	coarse=false;
	BN=true;
	p1=0;p2=0;i=0;j=0;level=0;
	for(int i=0;i<=9;i++){t[i]=0;}
}
//-----------**************---------------general functions--------------*************-----------------
int newtri()
{
	int j;
	for(j=nmax;j>0;j--) if(cell[j].deleted) break;

	if(j==0) {nmax++;return nmax;}
	else{cell[j].deleted=false;return j;}
}
//-----------**************---------------general functions--------------*************-----------------
int newpt()
{
	int i;
	for(i=pmax;i>0;i--) if(P[i].deleted) break;

	if(i==0) {pmax++;return pmax;}
	else{P[i].deleted=false;return i;}
}

//-----------**************---------------general functions--------------*************-----------------
//-----------**************---------------general functions--------------*************-----------------
//-----------**************---------------vector functions---------------*************-----------------
//-------------@@@@@@@@@@@@@*************@@@@@@@@@@@@@@@@@@@@***********@@@@@@@@@@@@@@-----------------
