#include <stdio.h>
#include <stdlib.h>
#include "myfunction_ppm.c"
#include "rccount.c"
#include <string.h>
#include <math.h>

/*--------------Data structures-----------*/
int *mat;       /* input image */
int *cdata,*visited,*taken;     
#include "line-draw.c"
typedef struct
{
        int i,j;
        int yes;
} vertx;

//vertx *cv;
vertx five[5],*v;
float *vd,*rd,ep=0;
int row,col,vdcount=0,rdcount=0,vcount=0;
int mark,count0=0,*direc,relax=1,x;
float vdev=0.0,omega=0.0;
float beta,alpha,totvd=0.0,totrd=0.0;

/*---------------Functions-------------*/
void get_the_mat(char * ,int ,int );    /* to read the input image and fill the 'mat' matrix*/
void build_image(int,int,char*);              /*construct the image to show curvature*/

void curve(int, int);
void start(int,int);
void getDist(vertx,vertx);
void aidp(vertx *,int,int);
void calcDist();
void calcPerp(float *,float *,vertx *, int, int, int);
float findPath(vertx *,vertx *,int *,int *,vertx, vertx);
void findNextPt(int,int,int *,int *, int);
void findn(vertx *,int *,int,int,vertx);
void poly(int i,int j);
int checkR1R2(int,int);

/*---------------------------------------*/
main(int argc,char* argv[])
                        {
                        int *x,flag=0;
                        int i,j,vc=0;
			FILE *fm=fopen("mat.dat","w");
			FILE *fvs=fopen("vis.dat","w");
			FILE *fv=fopen("ver.dat","w");

                        x=count_rc(argv[1]);   /* the pgm file is given as input ; returns #row, #col*/
                        row=*(x+1);
                        col=*x;
			printf("%d\t%d\n",row,col);
                          
			get_the_mat(argv[1],row,col); /* input image read -- fill mat*/
                        cdata   =(int *)malloc(row*col*sizeof(int));
                        visited =(int *)malloc(row*col*sizeof(int));
                        taken   =(int *)malloc(row*col*sizeof(int));

			printf("Enter alpha: ");
			scanf("%f",&alpha);
			printf("Enter beta:  ");
			scanf("%f",&beta);

			for(i=0;i<row;i++)
			{
				for(j=0;j<col;j++)
				{
					if(*(mat+i*col+j)==0 || *(mat+i*col+j)==1 || *(mat+i*col+j)==2)
						{*(mat+i*col+j)=0;
						*(cdata+i*col+j)=0;
						fprintf(fv,"%d\t%d\n",i,j);
						count0++;
						}
					fprintf(fm,"%d\t",*(mat+i*col+j));
					*(cdata+i*col+j)=255;
					*(visited+i*col+j)=0; 
					*(taken+i*col+j)=0;
				}
				fprintf(fm,"\n");
			}
			fprintf(fm,"\n");
			fclose(fm);
			fclose(fv);			

			v=(vertx *)malloc(count0*sizeof(int));
			direc=(int *)malloc(count0*sizeof(int));

			for(i=0;i<row;i++)
			   for(j=0;j<col;j++)
				if(*(mat+i*col+j)==0)
					*(cdata+i*col+j)=0;
			
			for(i=0;i<row;i++)
			{
				for(j=0;j<col;j++)
				{
					if(*(mat+i*col+j)==0 && *(taken+i*col+j)==0)
						{omega=0.0;ep=0.0;vdev=0.0;vcount=0;totvd=0.0;totrd=0.0;mark=0;vdcount=0;rdcount=0;x=0;
						printf("Calling for %d row, %d col\n",i,j);
						poly(i,j);
						//flag=1;
						}//if(flag) break;
				}//if(flag) break;
			}
			//poly(i,j);

			for(i=0;i<row;i++)
			   for(j=0;j<col;j++)
			      if(*(visited+i*col+j)==1)
				   {fprintf(fvs,"%d\t%d\n",i,j);
				   *(cdata+i*col+j)=20;
				   }
			fclose(fvs);

                        build_image(row,col, argv[1]);/*build image -- uses cdata*/   
                        }

void poly(int r,int c)
{
int i=r,j=c,k,chng,vc,ni,nj;
int start1=0,dir,pol=0;
mark=0;
FILE *ft=fopen("taken.dat","w");
//cv=(vertx *)malloc(row*col*sizeof(vertx));

v[mark].i=i;
v[mark].j=j;
v[mark].yes=0;
mark++;
fprintf(ft,"%d\t%d\n",i,j);
do
	{
	
	//cv[vcount].i=i;
	//cv[vcount].j=j;
	//cv[vcount].yes=0;

	chng=0;
	*(taken+i*col+j)=1;  //fprintf(ft,"%d\t%d\n",i,j);
	//printf("poly:: ver: row=%d col=%d\t",i,j);
	ni=i; nj=j;

	if(j>0 && *(mat+i*col+j-1)==0 && *(taken+i*col+j-1)==0){nj=j-1; chng++; dir=4;}
   
	if(*(mat+i*col+j+1)==0 && *(taken+i*col+j+1)==0){nj=j+1; chng++; dir=0;}

	if(i>0 && *(mat+(i-1)*col+j+1)==0 && *(taken+(i-1)*col+j+1)==0){ni=i-1; nj=j+1; chng++; dir=1;}

	if(i>0 && j>0 && *(mat+(i-1)*col+j-1)==0 && *(taken+(i-1)*col+j-1)==0){ni=i-1; nj=j-1; chng++; dir=3;}

	if(i>0 && *(mat+(i-1)*col+j)==0 && *(taken+(i-1)*col+j)==0){ni=i-1; chng++; dir=2;}

	if(j>0 && *(mat+(i+1)*col+j-1)==0 && *(taken+(i+1)*col+j-1)==0){ni=i+1; nj=j-1; chng++; dir=5;}

	if(*(mat+(i+1)*col+j)==0 && *(taken+(i+1)*col+j)==0){ni=i+1; chng++; dir=6;}

	if(*(mat+(i+1)*col+j+1)==0 && *(taken+(i+1)*col+j+1)==0){ni=i+1; nj=j+1; chng++; dir=7;}

	printf("chng=%d\n",chng);
	if(chng>1 && vcount>0) break;
	i=ni; j=nj;

	*(direc+vcount)=dir;
	vcount++;
	
	k=checkR1R2(start1,vcount);
	if(k)
		{start1=vcount;
		v[mark].i=i; v[mark].j=j; v[mark].yes=0; mark++;
		fprintf(ft,"%d\t%d\n",i,j);}

	if(!chng) {printf("Breaking..\n");break;}

	}while((i!=r || j!=c) && i<row && j<col);

if((i<=r+1 && i>=r-1) && (j<=c+1 && j>=c-1) && vcount>2) pol=1;

if(chng>1) pol=0;

if(!k && (chng==1 || (!chng && !pol))){v[mark].i=i; v[mark].j=j; v[mark].yes=0; mark++; fprintf(ft,"%d\t%d\n",i,j);}
fclose(ft);
printf("vcount= %d mark=%d\n",vcount,mark);

if(mark==1)
	{v[mark-1].yes=1;pol=2; bresenham_line(i,j,v[mark-1].i,v[mark-1].j,col);}
else if(mark==2)
	{v[mark-2].yes=v[mark-1].yes=1;
	bresenham_line(v[mark-2].i,v[mark-2].j,v[mark-1].i,v[mark-1].j,col);
	pol=2;}

	vc=vcount;
	vcount>>1;

	vd=(int *)malloc(vcount*(sizeof(int)));
	rd=(int *)malloc(vcount*(sizeof(int)));
	vcount=vc;

//	printf("------------------------Vertices------\n");
//	for(vc=0;vc<mark;vc++)
//		printf("Taken = i: %d  j: %d\n",v[vc].i,v[vc].j);

if(pol==1)
        start(0,mark); 
else if(!pol)
	curve(0,mark);

	for(vc=0;vc<mark;vc++)
	{
		if(v[vc].yes==1)
		{//fprintf(ft,"%d\t%d\n",v[vc].i,v[vc].j);
		*(cdata+(v[vc].i)*col+(v[vc].j))=30;
		}
	}  //fclose(ft);

//	free(v);
//	free(cv);
	free(vd);
	free(rd);

}

//--------------------------------------------DSS Check----------------------------------------------------------------------


int checkR1R2(int start,int ct)
{
   
    if(ct==0 || ct==1)return 0;
    int second,edge_break,first_run,prev,second_run,nonsing_run,prev_sec,first,third;
        second=300;edge_break=0;
        first_run=prev=0;second_run=0;nonsing_run=0;
for(x=start;x<ct;x++)
printf(" %d ",*(direc+x));
printf("\n");
        for(x=start;x<ct;x++)// checking rule 1
        {

                              if(x==0)
                                      first=*(direc+x);
                              if(*(direc+x)!=first && second==300)
                                      {second=*(direc+x);}
                              else if(*(direc+x)!=first && *(direc+x)!=second)
                                      {edge_break=1;/*count=x;*/break;}
       
         }

    if(second==300)return 0;
             if(edge_break==1)
             {
             return edge_break;
             }
            prev_sec=-1;prev=0;nonsing_run=0;x=start;first_run=0;
        while(*(direc+x)==first && x<ct)
            {first_run++;x++;}
        while(*(direc+x)==second && x<ct)
            {nonsing_run++;x++;}
        if(first_run<nonsing_run)
        {
            prev=nonsing_run;//3
            nonsing_run=0;
            while(*(direc+x)==first && x<ct)   x++;
            if(x==ct)return 0;
            while(*(direc+x)==second && x<ct){nonsing_run++;x++;}//3

            if(nonsing_run!=0)
            {
                if(x==ct)return 0;
                if(nonsing_run>=(prev-relax) && nonsing_run<=(prev+relax))
                    prev_sec=nonsing_run;
                else return 1;
            }
            else return 0;

            nonsing_run=0;
            while(x<ct)
            {
                nonsing_run=0;
                if(*(direc+x+1)==first && x!=ct-1)
                {
                    return 1;
                }
                 else x++;
            if(x==ct)return 0;
            while(x<ct && *(direc+x)==second)
            {
                nonsing_run++;
                x++;
            }
            if(!(nonsing_run>=(prev-relax) && nonsing_run<=(prev+relax)))
            {
                return 1;
            }
            nonsing_run=0;
        }
        return 0;
        }
        else if(nonsing_run==1)
        {
            prev=first_run;//2
            first_run=0;
            while(*(direc+x)==second && x<ct)
            {
                x++;
            }
        if(x==ct)return 0;
        while(*(direc+x)==first && x<ct)
        {
            first_run++;
            x++;
        }
        first_run=0;
        while(x<ct)
        {
            if(*(direc+x)==second && *(direc+x+1)==second && x+1<ct)
                return 1;
            else
                 x++;
            while(x<ct && *(direc+x)==first)
                {
                first_run++;
                x++;
                }
            if(x==ct )return 0;
            if(!(first_run>=(prev-relax) && first_run<=(prev+relax)))return 1;
            first_run=0;
        }
    return 0;
        }
       
             return 1;
}
//--------------------------------------------Adaptively Improved Douglas Peucker Algorithm-----------------------------------

void start(int start, int end)
{
vertx center, p1, p2;
int cx=0,cy=0,q,l1max=0,l2max=0,l=0,a1,a2;

for(q=start;q<end;q++)
  {
     cx+=v[q].i;
     cy+=v[q].j;//printf("Sum: %d , %d\n",x,y);
  }
if(end==start) {center.i=cx; center.j=cy;}
else{
cx=cx/(end-start);
center.i=cx;
cy=cy/(end-start);
center.j=cy;}

for(q=start;q<end;q++)
   {
      l=sqrt(pow((cx-v[q].i),2)+ pow((cy-v[q].j),2));
      if(l>l1max)
                {
                l1max=l;
                p1=v[q];a1=q;
                }
   }
   
for(q=start;q<end;q++)
   {
      l=sqrt(pow((p1.i-v[q].i),2)+ pow((p1.j-v[q].j),2));
      if(l>l1max)
                {
                l1max=l;
                p2=v[q];a2=q;
                }
   }
   
*(visited+p1.i*col+p1.j)=1;
*(visited+p2.i*col+p2.j)=1;
printf("x=%d\ty=%d\na1=%d  p1-i=%d  j=%d\na2=%d  p2-i=%d  j=%d\n",cx,cy,a1,p1.i,p1.j,a2,p2.i,p2.j);
getDist(p1,p2);

}

void getDist(vertx p1, vertx p2)//Finding the start point
{
    int i,j,index=0,a=0,k,l;
    vertx *cvA, *cvB;
    int cntA=0,cntB=0;
    vertx vprev;int prev=0;

    cvA=(vertx *)malloc((vcount/2)*sizeof(vertx));
    cvB=(vertx *)malloc((vcount/2)*sizeof(vertx));

    ep=findPath(cvA,cvB,&cntA,&cntB,p1,p2);
 
    if(ep==-1)
		printf("error\n");
    omega+=ep*beta;                 //Radial dist. threshold Omega
    
    printf("Parameters Omega: %f\tEpsilon:%f\n",omega,ep);
    aidp(cvA,0,cntA);
    cvA[cntA-1].yes=1;
    
    aidp(cvB,0,cntB);
    cvB[cntB-1].yes=1;

    for(k=0;k<cntA;k++)
    {
    for(l=0;l<mark;l++)
	if(v[l].i==cvA[k].i && v[l].j==cvA[k].j && cvA[k].yes==1)
		{prev++;v[l].yes=1;
		if(prev>1)
		bresenham_line(vprev.i,vprev.j,cvA[k].i,cvA[k].j,col);
		else
		bresenham_line(p1.i,p1.j,cvA[k].i,cvA[k].j,col);
		vprev=cvA[k];}
    }

bresenham_line(p2.i,p2.j,vprev.i,vprev.j,col);

prev=0;
    for(k=0;k<cntB;k++)
    {
     for(l=0;l<mark;l++)
	if(v[l].i==cvB[k].i && v[l].j==cvB[k].j && cvB[k].yes==1)
		{prev++;v[l].yes=1;
		if(prev>1)
		bresenham_line(vprev.i,vprev.j,cvB[k].i,cvB[k].j,col);
		else
		bresenham_line(p1.i,p1.j,cvB[k].i,cvB[k].j,col);
		vprev=cvB[k];}
    }
bresenham_line(p2.i,p2.j,vprev.i,vprev.j,col);

//free(cvA);
//free(cvB);
}

void aidp(vertx *points,int start,int end)
{
    float vd=0.0,vdist=0.0,rd=0.0;
    int index,i,chck=0;
    printf("end:%d\tstart:%d\n",end,start);
    if(end-start<3)
       return;

    for(i=start;i<end;i++)
       {
          calcPerp(&vd,&rd,points,start,i,end);
          if(vd>vdist)
                         {
                         vdist=vd;
                         index=i;
                         }
       }
       printf("Max dist:%f\tIndex:%d\n",vdist,index);
    
    if(vdist>=ep || (vdist<ep && rd>=omega))
       {
	   points[start].yes=1;
           points[index].yes=1;
	   if(index){
           aidp(points,start,index-1);printf("BRANCH\n");
           aidp(points,index,end);}
                
       }
    else
        {
        points[start].yes=1;
        //points[end].yes=1;
        }
     
}

//----------------------------------------------------------------------------------------------------------------------------
float findPath(vertx *cvA,vertx *cvB,int *acnt,int *bcnt,vertx p1, vertx p2)//Traverse the curves A and B
{
int i,j;
int a=0,err=0,next1row,next1col,next2row,next2col,dir1,dir2;
float epA=0, epB=0, avA=0, avB=0;
int cntA,cntB;

i=p1.i;j=p1.j;

cntA=*(acnt);
cntB=*(bcnt);

for(a=0;a<mark && (v[a].i!=i || v[a].j!=j);a++);

printf("fP: a=%d row=%d col=%d\n",a,i,j);

next1row=v[(a+1)%mark].i;
next1col=v[(a+1)%mark].j;

next2row=v[(mark+a-1)%mark].i;
next2col=v[(mark+a-1)%mark].j;

printf("next1row:%d\tnext1col:%d\nnext2row:%d\tnext2col:%d\n",next1row,next1col,next2row,next2col);


//------------------Traverse curve A
if(p2.i!=next1row || p2.j!=next1col)
   {
   *(visited+next1row*col+next1col)=1;
   printf("entered 1st call\n");
   findn(cvA,&cntA,next1row,next1col,p2);}
else printf("couldn't enter");

if(cntA)
   {if(cntA!=1)
	   avA=(totvd/(cntA-1));                                //Get the Avg. Vertical Dist.
   else avA=totvd;

if(cntA<10)
	{
	for(i=0;i<cntA;i++)
		{cvA[i].yes=1;}	
	}
   for(i=0;i<vdcount;i++)                           //Calculate the mean deviation
      vdev+=abs(vd[i]-avA);
   printf("vdev before division=%f and vdcount=%d\n",vdev,vdcount);
   if(cntA!=1) vdev/=(cntA-1);

   epA = avA + alpha*vdev;                           //Vertical dist. threshold Epsilon for A
   omega=totrd;                             
   printf("cntA=%d\tavA=%f\nvdev=%f\tepA=%f\ntotvd=%f\ttotrd=%f\n",cntA,avA,vdev,epA,totvd,totrd);
   }
else
   {epA=0;omega=0;}

//--------------------Traverse curve B
totvd=0;totrd=0;
if(p2.i!=next2row || p2.j!=next2col)
   {
   *(visited+next2row*col+next2col)=1;
   findn(cvB,&cntB,next2row,next2col,p2);
   }

if(cntB)
   {if(cntB!=1)
	   avB=(totvd/(cntB-1));                                 //Get the Avg. Vertical Dist.
    else
	avB=totvd;

if(cntB<10)
	{
	for(i=0;i<cntB;i++)
		{cvB[i].yes=1;}	
	}
   vdev=0.0;
   for(i=0;i<vdcount;i++)                             //Calculate the mean deviation
      vdev+=abs(vd[i]-avB);
   printf("vdev before division=%f and vdcount=%d\n",vdev,vdcount);
   if(cntB!=1) vdev/=(cntB-1);

   epB = avB + alpha*vdev;                            //Vertical dist. threshold Epsilon for B
   omega+=totrd;
   omega/=(cntA+cntB-1);                               
   printf("cntB=%d\tavB=%f\nvdev=%f\tepB=%f\tomega=%f\ntotvd=%f\ttotrd=%f\n",cntB,avB,vdev,epB,omega,totvd,totrd);
   }

else
   {epB=0;omega=0;}

*(acnt)=cntA; *(bcnt)=cntB;

if(epA>epB) return epA;
return epB;

}

void findn(vertx *c,int *cnt,int i,int j,vertx end)//------------------Searches the neighbourhood of a point to get the next
{
int ch=0;
int cnt1,k,a,b;
*(cnt)=0;
cnt1=*(cnt);

for(a=0;a<mark && (v[a].i!=i || v[a].j!=j);a++);

printf("ver: a=%d row=%d col=%d\n",a,i,j);

while(v[a].i!=end.i || v[a].j!=end.j)
	{

	c[cnt1].i=i;
	c[cnt1].j=j;
	c[cnt1].yes=0;

	for(k=0;(k<cnt1 && cnt1<4) || k<4; k++)
		five[k]=five[k+1];
	five[k].i=i;   five[k].j=j;
	printf("ver in loop: a=%d row=%d col=%d\n",a,i,j);
	*(visited+i*col+j)=1;

	if(cnt1%2==0 && cnt1>=4)  calcDist();
	cnt1++;

	a=(a+1)%mark;//printf("a1=%d\n",a);
	i=v[a].i; j=v[a].j;
	if(*(visited+i*col+j)==1) 
		{a=(mark+a-2)%mark;/*printf("a2=%d\n",a);*/
		i=v[a].i; j=v[a].j;
		if(*(visited+i*col+j)==1) 
			{printf("Breaking, no choice\n");break;}
		}
	}


*(cnt)=cnt1;
}

void calcDist()//-----------------------------------calculate the vertical and radial distances------------------------
{
    float vdist,rdist,rd1,rd2;
    float m=0.0,c=0.0;
    vertx p1,p2,p;
    p1=five[0];
    p2=five[4];
    p=five[2];
    
    if(p1.j==p2.j) {/*printf("mark\n");*/vdist=abs(p.j-p1.j);}
    else{
    m=(p1.i-p2.i)/(p1.j-p2.j);
    c=p1.i-p1.j*m;
    
    vdist=abs((p.i-(p.j*m)-c))/sqrt(1+m*m);}
    rd1=sqrt(pow((p.i-p1.i),2)+pow((p.j-p1.j),2));
    rd2=sqrt(pow((p.i-p2.i),2)+pow((p.j-p2.j),2));
    //printf("P1: r=%d c=%d P: r=%d c=%d P2: r=%d c=%d\n",p1.i,p1.j,p.i,p.j,p2.i,p2.j);
    //printf("slope=%f  intercept=%f  vert dist=%f rad 1=%f rad2=%f\n",m,c,vdist,rd1,rd2);
    //printf("vert dist=%f rad 1=%f rad2=%f\n",vdist,rd1,rd2);
    vd[vdcount++]=vdist;
    totvd+=vdist;
    if(rd1>=rd2) rdist=rd1; else rdist=rd2;
    rd[rdcount++]=rdist;
    totrd+=rdist;
}

void calcPerp(float *vd ,float *rd ,vertx *points, int start, int x, int end)
{
    float ve,r,rd1,rd2,m,c;
    vertx p1,p2,p;
    p1=points[start];
    p2=points[end];
    p=points[x];
    
    if(p1.j==p2.j) 
	ve=abs(p.j-p1.j);
    else
    {
    m=(p1.i-p2.i)/(p1.j-p2.j);
    c=p1.i-p1.j*m;
    ve=abs((p.i-(p.j*m)-c)/sqrt(1+m*m));
    }

    rd1=sqrt(pow((p.i-p1.i),2)+pow((p.j-p1.j),2));
    rd2=sqrt(pow((p.i-p2.i),2)+pow((p.j-p2.j),2));

    //printf("P1: r=%d c=%d P: r=%d c=%d P2: r=%d c=%d\n",p1.i,p1.j,p.i,p.j,p2.i,p2.j);
    //printf("slope=%f  intercept=%f  vert dist=%f rad 1=%f rad2=%f\n",m,c,ve,rd1,rd2);
    //printf("CalcPerp vert dist=%d  rad dist1=%d rad dist2=%d\n",ve,rd1,rd2);

    if(rd1>rd2)   r=rd1; else r=rd2;
    
    *(vd)=ve;
    *(rd)=r;
}

//---------------------------------------AIDP on curves----------------------------------------------------------------------

void curve(int start, int end)
{
vertx st,en,vprev;
vertx *c;
int i,j,k,l,cnt=0,prev=0;
float avg;

st=v[start];
en=v[end];

c=(vertx *)malloc(vcount*sizeof(vertx));

findn(c,&cnt,st.i,st.j,en);

if(cnt)
   {if(cnt!=1)
	   avg=(totvd/(cnt-1));                                //Get the Avg. Vertical Dist.
   else avg=totvd;
if(cnt<10)
	{
	for(i=0;i<cnt;i++)
		{c[i].yes=1;}	
	}

   for(i=0;i<vdcount;i++)                           //Calculate the mean deviation
      vdev+=abs(vd[i]-avg);
   printf("vdev before division=%f and vdcount=%d\n",vdev,vdcount);
   if(cnt!=1) vdev/=(cnt-1);

   ep = avg + alpha*vdev;                           //Vertical dist. threshold Epsilon for A
   omega=(totrd/cnt)+beta*vdev;                             
   printf("cnt=%d\tavg=%f\nvdev=%f\tep=%f\ntotvd=%f\ttotrd=%f\n",cnt,avg,vdev,ep,totvd,totrd);
   }
else
   {ep=0;omega=0;}

aidp(c,0,cnt);
c[cnt-1].yes=1;

    for(k=0;k<cnt;k++)
    {
    for(l=0;l<mark;l++)
	if(v[l].i==c[k].i && v[l].j==c[k].j && c[k].yes==1)
		{prev++;v[l].yes=1;
		if(prev>1)
		bresenham_line(vprev.i,vprev.j,c[k].i,c[k].j,col);
		vprev=c[k];}
    }

free(c);
}

//---------------------------------------Image input and output--------------------------------------------------------------
void get_the_mat(char *s,int r,int c)
                        {
                          FILE *fp;
			  unsigned char ch;
                          int line_count=0;
                          int r1=0,c1=0;
                          //int ch;
                          int i,j;
                         /*----------open the image file------------*/
                          fp=fopen(s,"rb");
                          mat=(int*)malloc(r*c*sizeof(int));
                          while(fread(&ch,sizeof(ch),1,fp)==1) /* read all values and place in mat*/
                             {
                               if(ch=='\n') line_count++;
                                if(line_count==4) break;
                              }

                           while(fread(&ch,sizeof(ch),1,fp)==1)
                                {
                                  *(mat+r1*c+c1)=ch;
                                  c1++;
                                  if(c1>(c-1)) /* matrix building*/
                                      {c1=c1%c;
                                       r1=(r1+1)%r;}
                                 }

                           fclose(fp);

                         } /*func. end*/
/*---------------------------------------*/
void build_image (int r,int c, char* s) /*construct the image to show curvature*/
                     {
                          FILE *fp1,*fp2;
                          int red,green,blue;
                          int i,j,mydata;
                          unsigned char ch,color[3];
                          char t[200];
                          int d,m;
                          char mystring[30]; /* output image name*/
                          build_my_file(r,c);/*to make myfile.dat */
                          /*--------------------------------*/
                 for(i=strlen(s)-1;i>=0;i--)
                      {
                       if (s[i]=='.' ) break; /*to reflect the input image name to output*/
                       }
                  for(j=0;j<i;j++)
                       t[j]=s[j];
                    t[j]='\0';
               /*-------------------------------*/
			  sprintf(mystring,"-%.1f-%.2f.ppm",alpha,beta);
                          strcat(t,mystring);

                          fp2=fopen("myfile.dat","rb"); /* contains the header lines of the image*/
                          fp1=fopen(t,"wb+"); /*update t if needed- to change name*/
                          while(fread(&ch,sizeof(ch),1,fp2)==1)
                          fwrite(&ch,sizeof(ch),1,fp1);
                          for(i=0;i<r;i++)
                          for(j=0;j<c;j++)
                           {
                            //ch=*(cdata+i*c+j);
				mydata=*(cdata+i*c+j);//printf("%d\n",mydata);
				if(mydata==0 ){color[0]=127;color[1]=127;color[2]=127;/*printf("0\n");*/fwrite(color,1,3,fp1);}
				else if(mydata==255 ){color[0]=255;color[1]=255;color[2]=255;/*printf("255\n");*/fwrite(color,1,3,fp1);}
				else if(mydata==9){color[0]=0;color[1]=0;color[2]=255;/*printf("9\n");*/fwrite(color,1,3,fp1);}
				else if(mydata==20){color[0]=255;color[1]=0;color[2]=0;/*printf("20\n");*/fwrite(color,1,3,fp1);}
				else if(mydata==30){color[0]=0;color[1]=255;color[2]=0;/*printf("30\n");*/fwrite(color,1,3,fp1);}
                            //fwrite(&ch,sizeof(ch),1,fp1);
                           }
                          

                         fclose(fp2);
                         fclose(fp1);
                       } /*func. end*/
/*---------------------------------------*/

