
#include <stdio.h>
#include <stdlib.h>
#include "myfunction_pgm.c"
#include "rccount.c"
#include <string.h>
#include <math.h>

//--------------Data structures-----------
typedef struct {
	int i; int j;
	} Vertex;
int eresponse[1000000],finval;
int *mat,*direc,*tmat;       /* input image */
int *cdata,*visited,*newimage;
int st,count,start_i,start_j,vcount,start,ecount,ethresh,thresh,new_i,new_j,new_dir,relax;
int x=0,*dirs,r=0,expval=0; vcount=0;
ethresh=10;  
int row,col;  
Vertex *v;


//---------------Functions-------------
void get_the_mat(char * ,int ,int );    /* to read the input image and fill the 'mat' matrix*/
void build_image(int,int,char*);              /*construct the image to show curvature*/
void find_8N(int,int,int,int);
void find3N(int, int, int);
int check_R1(int);
void find_max(int,int);
void myfilter(int r,int c);

int mysqrt(int x)
   {
     int i;
     int r;
      for(i=0;i<x;i++)
           {
            if(i*i>=x)
          {
             r=i;
              break;
            }
             }
      return i;
  }

int prewitt(int i,int j,int r,int c)
   {
    // int i,j;
     
      int dvx,dvy;
       int t;
           if(i==0||i==r-1)
                        {
                        if(j==c-1)
                                dvx=255-*(mat+i*c+j);
                        else dvx=*(mat+i*c+j+1)-*(mat+i*c+j);
                        }
                     else
              dvx= *(mat+(i-1)*c+j-1)*(-1)+*(mat+(i-1)*c+j)*(-1)+*(mat+(i-1)*c+j+1)*(-1)
                              +*(mat+(i)*c+j-1)*0+ *(mat+(i)*c+j)*0+*(mat+(i)*c+j+1)*0+
                              *(mat+(i+1)*c+j-1)*(1)+*(mat+(i+1)*c+j)*(1)+*(mat+(i+1)*c+j+1)*(1);
              
                                                                                                                            
      
                                                                                                                            
        if(j==0||j==c-1)
               {
                if(i==r-1)
                        dvy=255-*(mat+i*c+j);
                else dvy=*(mat+(i+1)*c+j)-*(mat+i*c+j);
               }
            else
                dvy= *(mat+(i-1)*c+j-1)*(-1)+*(mat+(i-1)*c+j)*(0)+*(mat+(i-1)*c+j+1)*(1)
                              +*(mat+(i)*c+j-1)*(-1)+ *(mat+(i)*c+j)*0+*(mat+(i)*c+j+1)*1+
                              *(mat+(i+1)*c+j-1)*(-1)+*(mat+(i+1)*c+j)*(0)+*(mat+(i+1)*c+j+1)*(1);
      t= testpixel(dvx,dvy);
           return (t);
}
int testpixel(int x,int y)
                 {
                    int s;
                    s=mysqrt(abs(x)*abs(x)+abs(y)*abs(y));
                   return s;   
                   }
//Main: Computes Prewitt response matrix and builds the new image



int main(int argc,char* argv[])
                        {
                          int *x;

                          int i,j;int result;
			  int fx,fy;
                          x=count_rc(argv[1]);   /* the pgm file is given as input ; returns #row, #col*/
                          row=*(x+1);
                          col=*x;
//row=row/2;col=col/2;
                          cdata   =(int*)malloc(row*col*sizeof(int));
			  newimage=(int*)malloc(row*col*sizeof(int));
			  direc=(int*)malloc(row*col*sizeof(int));
			  mat=(int*)malloc(row*col*sizeof(int));
			  v=(Vertex*)malloc(row*col*sizeof(Vertex));
 get_the_mat(argv[1],row,col); /* input image read -- fill mat*/
                          myfilter(row,col);
			 
printf("enter the relaxation\n");
scanf("%d",&relax);
			for(i=0;i<row;i++)
			{
				for(j=0;j<col;j++)
				{
					*(newimage+i*col+j)=255;
				}
			}
                           for(i=0;i<row;i++)
                           for(j=0;j<col;j++)
                            {

/*			    fy= (*(mat+(i-1)*col+(j-1))-*(mat+(i+1)*col+(j-1)))	+	(*(mat+(i-1)*col+(j))-*(mat+(i+1)*col+(j)))	+	(*(mat+(i-1)*col+(j+1))-*(mat+(i+1)*col+(j+1)));
			    fy=fy/6;

			    fx=	(*(mat+(i-1)*col+(j+1))-*(mat+(i-1)*col+(j-1))) +	(*(mat+(i)*col+(j+1))-*(mat+(i)*col+(j-1)))	+	(*(mat+(i+1)*col+(j+1))-*(mat+(i+1)*col+(j-1)));
			    fx=fx/6;

			    result=sqrt(pow(fx,2)+pow(fy,2));
*(cdata+i*col+j)=result;  */
*(cdata+i*col+j)=prewitt(i,j,row,col);
                          } 
                           


find_max(row,col);
                           build_image(row,col, argv[1]);/*build image -- uses cdata*/   
                           return 0;
                         } /* main end*/  

void myfilter(int r,int c)
              {
                  int i,j;
                  
                  for(i=0;i<=(r-1);i++)
                  for(j=0;j<=(c-1);j++)
                    {
                       if(i==0||j==0||i==r-1||j==c-1 )
                        *(mat+i*c+j)=*(tmat+i*c+j);
                     else
                         *(mat+i*c+j)=(*(tmat+i*c+j)*4+*(tmat+(i-1)*c+j-1)
                                     +*(tmat+(i-1)*c+j)*2+*(tmat+(i-1)*c+j+1)
                                     +*(tmat+i*c+j-1)*2+*(tmat+i*c+j+1)*2
                                     +*(tmat+(i+1)*c+j-1)+*(tmat+(i+1)*c+j)*2
                                     +*(tmat+(i+1)*c+j+1))/16;

                         /*  *(mat+i*c+j)=(*(tmat+i*c+j)+*(tmat+(i-1)*c+j-1)
                                     +*(tmat+(i-1)*c+j)+*(tmat+(i-1)*c+j+1)
                                     +*(tmat+i*c+j-1)+*(tmat+i*c+j+1)
                                     +*(tmat+(i+1)*c+j-1)+*(tmat+(i+1)*c+j)
                                     +*(tmat+(i+1)*c+j+1))/9;*/
                    }
               }


int expavg(int r,int c) //-------------Calculating exponential average
{
     int pre,cnt,i,p,j,val,ct,change;
  //   printf("in expavg and start is %d and vcount is %d. The count of direction is %d.ecount is %d\n",start,vcount,count,ecount);
   val=0;p=1;ct=change=0;
if(r<0 || c<0 || r>row|| c>col)return 0;
if(ecount<7)
{
//printf("PREWITT is %d\n",*(cdata+r*col+c));
finval=*(cdata+r*col+c);
return *(cdata+r*col+c);
}
for(cnt=0;cnt<=ecount-1;cnt++)
{
if(eresponse[cnt]<=0 || eresponse[cnt]==eresponse[cnt+1])
{

pre=cnt;
while(pre<ecount)
	{eresponse[pre]=eresponse[pre+1];pre++;}
	change++;
}
}
ecount-=change;
if(ecount<7){finval=*(cdata+r*col+c);return *(cdata+r*col+c);}
for(cnt=ecount-1;cnt>=0;cnt--)
           {
	   i=v[cnt].i;
           j=v[cnt].j;
	  // printf("in expavg %d  %d",eresponse[cnt],eresponse[cnt]>>p);
           val+=(eresponse[cnt]>>p);p++;
	   //printf("..\n");ct++;if(ct==8)break;
           }

finval=(val+*(cdata+r*col+c));//finval=finval>>1;
//printf("%d\n",*(cdata+r*col+c));
//printf("finval is %d\n",finval);
return finval;
}
int check(int r,int c) //---------------Checks whether a vertex obeys the threshold
{

     if(ethresh<=expavg(r,c))
       { return 1;}
     if(*(cdata+r*col+c)>=thresh)
        return 1;
     else
         {return 0;}
}

//Finds the point having local maximum over its 8 neighbours (Start Point) ----------------------------------------

void find_max(int row,int col)
{

	int i,j,item;	
	start=0;
	printf("Enter the threshold:\n");
	scanf("%d",&thresh);
x=0;
	for(i=1;i<row-1;i++)
	{
		for(j=1;j<col-1;j++)
		{
			if(*(cdata+(i*col)+j)>thresh)
			{
				item=*(cdata+(i*col)+j);
				if( *(cdata+(i*col)+j-1)<=item && *(cdata+(i*col)+j+1)<=item && *(cdata+((i-1)*col)+j-1)<=item && *(cdata+((i-1)*col)+j)<=item && *(cdata+((i-1)*col)+j+1)<=item  && *(cdata+((i+1)*col)+j-1)<=item && *(cdata+((i+1)*col)+j)<=item && *(cdata+(i*col)+j+1)<=item)
					{ //printf("calling find_8N\n");
				//	v[vcount].i=i;v[vcount].j=j;vcount++;
				//	*(newimage+i*col+j)=0;
					//eresponse[ecount++]=expavg(i,j);
					start_i=i;start_j=j;ecount=0;
					find_8N(i,j,row,col);
					find_8N(start_i,start_j,row,col);
					}
			}
		}
	}
}

//Finds direction from the Start Point---------------------------------------------

void find_8N(int i,int j,int row,int col)
{
//*(newimage+i*col+j)=0;
//start=vcount;
	v[vcount].i=i;v[vcount].j=j;vcount++;
eresponse[ecount++]=expavg(i,j);
if(i<0 || i>row ||j<0 ||j > col-1)return;
	int max,max_i,max_j,dir;
	max=-1;max_i=0;max_j=0;dir=0;
	if(i!=0 && j!=0 && *(cdata+(i-1)*col+j-1)>max && *(newimage+(i-1)*col+j-1)!=0)
		{
		max=*(cdata+(i-1)*col+j-1);max_i=i-1;
		max_j=j-1;
		dir=3;
		}

	if(i!=0 && *(cdata+(i-1)*col+j)>max && *(newimage+(i-1)*col+j)!=0)
		{
		max=*(cdata+(i-1)*col+j);
		max_i=i-1;
		max_j=j;
		dir=2;
		}

	if(i!=0  && *(cdata+(i-1)*col+j+1)>max && *(newimage+(i-1)*col+j+1)!=0)
		{
		max=*(cdata+(i-1)*col+j+1);
		max_i=i-1;
		max_j=j+1;
		dir=1;
		}

	if(j!=0 && *(cdata+(i)*col+j-1)>max && *(newimage+(i)*col+j-1)!=0)
		{
		max=*(cdata+(i)*col+j-1);
		max_i=i;
		max_j=j-1;
		dir=4;
		}

	if(*(cdata+(i)*col+j+1)>max && *(newimage+(i)*col+j+1)!=0)
		{
		max=*(cdata+(i)*col+j-1);
		max_i=i;
		max_j=j+1;
		dir=0;
		}

	if(j!=0 && *(cdata+(i+1)*col+j-1)>max && *(newimage+(i+1)*col+j-1)!=0)
		{
		max=*(cdata+(i+1)*col+j-1);
		max_i=i+1;
		max_j=j-1;
		dir=5;
		}

	if(*(cdata+(i+1)*col+j)>max && *(newimage+(i+1)*col+j)!=0)
		{
		max=*(cdata+(i+1)*col+j);
		max_i=i+1;
		max_j=j;
		dir=6;
		}

	if(*(cdata+(i+1)*col+j+1)>max && *(newimage+(i+1)*col+j+1)!=0)
		{
		max=*(cdata+(i+1)*col+j+1);
		max_i=i+1;
		max_j=j+1;
		dir=7;
		}
	if(max==-1)return;
 //find3N(max_i,max_j,dir);
if(i!=0 && j!=0 && *(cdata+(i-1)*col+j-1)==max && *(newimage+(i-1)*col+j-1)!=0 && check(i-1,j-1)==1)
{
	v[vcount].i=i-1;
	v[vcount].j=j-1;
	eresponse[ecount++]=expavg(i-1,j-1);
	vcount++;
	*(direc+count)=3;count++;//start_i=i-1;start_j=j-1; printf("calling 3N\n");
	if(j-1<col && j-1>=0)
	find3N(i-1,j-1,3);
}
//count=0;
if(i!=0 && *(cdata+(i-1)*col+j)==max && *(newimage+(i-1)*col+j)!=0 && check(i-1,j)==1)
{
	v[vcount].i=i-1;
	v[vcount].j=j;
	vcount++;eresponse[ecount++]=expavg(i-1,j);
	*(direc+count)=2;
	count++;//start_i=i-1;start_j=j;//printf("calling 3N\n");
	find3N(i-1,j,2);}
//count=0;
if(i!=0 && *(cdata+(i-1)*col+j+1)==max && *(newimage+(i-1)*col+j+1)!=0 && check(i-1,j+1)==1)
{
	v[vcount].i=i-1;
	v[vcount].j=j+1;
	vcount++;eresponse[ecount++]=expavg(i-1,j+1);
	*(direc+count)=1;
	count++;//start_i=i-1;start_j=j+1; printf("calling 3N\n");
	find3N(i-1,j+1,1);}
//count=0;
if(*(cdata+(i)*col+j-1)==max && check(i,j-1)==1 && *(newimage+(i)*col+j-1)!=0)
{
	v[vcount].i=i;
	v[vcount].j=j-1;eresponse[ecount++]=expavg(i,j-1);
	vcount++;*(direc+count)=4;
	count++;//start_i=i;start_j=j-1; printf("calling 3N\n");
	if(j-1<col)find3N(i,j-1,4);}
//count=0;
if(*(cdata+(i)*col+j+1)==max && check(i,j+1)==1 && *(newimage+(i)*col+j+1)!=0)
{
	v[vcount].i=i;
	v[vcount].j=j+1;eresponse[ecount++]=expavg(i,j+1);
	vcount++;*(direc+count)=0;
	count++;//start_i=i;start_j=j+1; printf("calling 3N\n");
	find3N(i,j+1,0);}
//count=0;
if(j!=0 && *(cdata+(i+1)*col+j-1)==max && check(i+1,j-1)==1 && *(newimage+(i+1)*col+j-1)!=0)
{
	v[vcount].i=i+1;
	v[vcount].j=j-1;
	vcount++;*(direc+count)=5;eresponse[ecount++]=expavg(i+1,j-1);
	count++;//start_i=i+1;start_j=j-1; printf("calling 3N\n");
	if(j-1<col)find3N(i+1,j-1,5);}
//count=0;
if(*(cdata+(i+1)*col+j)==max && check(i+1,j)==1 && *(newimage+(i+1)*col+j)!=0)
{	v[vcount].i=i+1;
	v[vcount].j=j;eresponse[ecount++]=expavg(i+1,j);
	vcount++;*(direc+count)=6;
	count++;//start_i=i+1;start_j=j; printf("calling 3N\n");	
	find3N(i+1,j,6);}
//count=0;
if(*(cdata+(i+1)*col+j+1) ==max && check(i+1,j+1)==1 && *(newimage+(i+1)*col+j+1)!=0)
{
	v[vcount].i=i+1;
	v[vcount].j=j+1;eresponse[ecount++]=expavg(i+1,j+1);
	vcount++;*(direc+count)=7;
	count++;//start_i=i+1;start_j=j+1; printf("calling 3N\n");
	find3N(i+1,j+1,7);}
}

int findNextPt(int i,int j,int *nexti,int *nextj, int dir) // -----------Gets coordinates of next point based on direction
{
 //printf("in findnext\n");

if(dir==0)
	{*(nexti)=i; *(nextj)=j+1;}
if(dir==1)
	{*(nexti)=i-1; *(nextj)=j+1;}
if(dir==2)
	{*(nexti)=i-1; *(nextj)=j;}
if(dir==3)
	{*(nexti)=i-1; *(nextj)=j-1;}
if(dir==4)
	{*(nexti)=i; *(nextj)=j-1;}
if(dir==5)
	{*(nexti)=i+1; *(nextj)=j-1;}
if(dir==6)
	{*(nexti)=i+1; *(nextj)=j;}
if(dir==7)
	{*(nexti)=i+1; *(nextj)=j+1;}
 //printf("next_i %d and next_j %d\n",*nexti,*nextj);
if(*nexti<0 || *nextj<0)return 0;
if(*nexti>=row || *nextj>=col)return 0;
//new_i=*nexti;new_j=*nextj;
if(check(*nexti,*nextj)==0)return 0;
else 
return 1;
}

void edge_creation(int ct)
{
//count--;
//FILE *fp;
 //printf("entered edge creation\n");
int i,prev_i,prev_j,now_i,in,m,de,inj,dej,now_j;
unsigned char ch;
start=vcount;
if(count<3){st=1;*(direc)=*(direc+count-1);count=1;return;}
*(newimage+start_i*col+start_j)=1;
*(cdata+start_i*col+start_j)=-1;
prev_i=start_i;prev_j=start_j;

for(i=0;i<count;i++)
{
//printf("direction is %d\n",*(direc+i));
//printf("in loop\n");
switch (*(direc+i))
{
case 0:now_i=prev_i;now_j=prev_j+1;break;
case 1:now_i=prev_i-1;now_j=prev_j+1;break;
case 2:now_i=prev_i-1;now_j=prev_j;break;
case 3:now_i=prev_i-1;now_j=prev_j-1;break;
case 4:now_i=prev_i;now_j=prev_j-1;break;
case 5:now_i=prev_i+1;now_j=prev_j-1;break;
case 6:now_i=prev_i+1;now_j=prev_j;break;
case 7:now_i=prev_i+1;now_j=prev_j+1;break;
}
//printf("now_i %d now_j %d\n",now_i,now_j);
if(now_i>=0 && now_i<row && now_j>=0 && now_j<col)
{
*(newimage+col*now_i+now_j)=0;
*(cdata+col*now_i+now_j)=-1;
prev_i=now_i;prev_j=now_j;}
else return;
if(*(direc+i)==0)
{
	//printf("direction 0\n");
//printf("now_i is %d\n",now_i);
	in=now_i+1;de=now_i-1;
	for(m=1;m<=2 && de>=0 && in<row;m++)
	{
//		*(newimage+in*col+now_j)=255;*(newimage+de*col+now_j)=255;	
		*(cdata+in*col+now_j)=0;*(cdata+de*col+now_j)=0;in++;de--;
	}
}

if(*(direc+i)==1)
{
//	printf("direction 1\n");
	in=now_i+1;inj=now_j+1;de=now_i-1;dej=now_j-1;
	for(m=1;m<=2 && de>=0 && dej>=0 && in<row && inj<col;m++)
	{
//		*(newimage+in*col+inj)=255;*(newimage+de*col+dej)=255;	
		*(cdata+in*col+inj)=0;*(cdata+de*col+dej)=0;in++;de--;inj++;dej--;
	}
}

if(*(direc+i)==2)
{
//	printf("direction 2\n");
	in=now_j+1;de=now_j-1;
//printf("in %d de %d\n",in,de);
	for(m=1;m<=2 && de>=0 && in<col;m++)
	{	//printf("hi\n");printf("ppp%d\t%d\n",*(cdata+now_i*col+in),*(cdata+now_i*col+de));
//		*(newimage+now_i*col+in)=255;
		*(cdata+now_i*col+in)=0;printf("in %d\n",in);
//		*(newimage+now_i*col+de)=255;
		*(cdata+now_i*col+de)=0;printf("de %d\n",de);in++;de--;
	}
//printf("over\n");
}

if(*(direc+i)==3)
{
//	printf("direction 3\n");
	in=now_i-1;inj=now_j+1;de=now_i+1;dej=now_j-1;
	for(m=1;m<=2 && in>=0 && dej>=0 && de<row && inj<col;m++)
	{
//		*(newimage+in*col+inj)=255;*(newimage+de*col+dej)=255;	
		*(cdata+in*col+inj)=0;*(cdata+de*col+dej)=0;in--;de++;inj++;dej--;
	}
}

if(*(direc+i)==4)
{
//	printf("direction 4\n");
	in=now_i+1;de=now_i-1;
	for(m=1;m<=2 && in<row && de>=0;m++)
	{
//		*(newimage+in*col+now_j)=255;*(newimage+de*col+now_j)=255;	
		*(cdata+in*col+now_j)=0;*(cdata+de*col+now_j)=0;in++;de--;
	}
}

if(*(direc+i)==5)
{
//	printf("direction 5\n");
	in=now_i+1;inj=now_j+1;de=now_i-1;dej=now_j-1;
	for(m=1;m<=2 && in<row && de>=0 && inj<col && dej>=0;m++)
	{
//		*(newimage+in*col+inj)=255;*(newimage+de*col+dej)=255;	
		*(cdata+in*col+inj)=0;*(cdata+de*col+dej)=0;in++;de--;inj++;dej--;
	}
}

if(*(direc+i)==6)
{
//	printf("direction 6\n");
	in=now_j+1;de=now_j-1;
	for(m=1;m<=2 && in<col && de>=0;m++)
	{
//		*(newimage+now_i*col+in)=255;*(newimage+now_i*col+de)=255;	
		*(cdata+now_i*col+in)=0;*(cdata+now_i*col+de)=0;in++;de--;
	}
}

if(*(direc+i)==7)
{
//	printf("direction 7\n");
	in=now_i-1;inj=now_j+1;de=now_i+1;dej=now_j-1;
	for(m=1;m<=2 && in>=0 && dej>=0 && de<row && inj<col;m++)
	{
//		*(newimage+in*col+inj)=255;*(newimage+de*col+dej)=255;	
		*(cdata+in*col+inj)=0;*(cdata+de*col+dej)=0;in--;de++;inj++;dej--;
	}
}
}

st=1;
*(newimage+col*now_i+now_j)=2;
*(direc)=*(direc+count-1);count=0;
start_i=now_i;start_j=now_j;
 //printf("done edge creation\n");
}


void find3N(int i,int j,int dir)//  ---------------------------------Finds maximum among 3 neighbours
{
int k;
int p1dir,p2dir,p3dir;
int p1,p2,p3,maxval;
p1=p2=p3=-1;
int nexti1,nextj1,nexti2,nextj2,nexti3,nextj3;
if(check(i,j)==0)return;

 if(count==0)return;
if(i<0 || j<0 ||i >row||j>col)return;
p1dir=dir;
p2dir=(dir+1)%8;
p3dir=(dir+7)%8;
//printf("p1dir %d p2dir %d p3dir %d\n",p1dir,p2dir,p3dir);
 //printf("calling findnext\n");

if(findNextPt(i,j,&nexti1,&nextj1,p1dir)==1)
{if(nexti1<0 || nextj1<0 || nexti1>row || nextj1>col)return;

p1=*(cdata+nexti1*col+nextj1);
}

 //printf("calling findnext\n");
if(findNextPt(i,j,&nexti2,&nextj2,p2dir)==1)
{if(nexti2<0 || nextj2<0 || nexti2>row || nextj2>col)return;
p2=*(cdata+nexti2*col+nextj2);
}

 //printf("calling findnext\n");
if(findNextPt(i,j,&nexti3,&nextj3,p3dir)==1)
{if(nexti3<0 || nextj3<0 || nexti3>row || nextj3>col)return;
p3=*(cdata+nexti3*col+nextj3);
}

if(p1==-1 || p2==-1 || p3==-1)
{//printf("no max\n");
edge_creation(count);
//*(direc)=dir;count=1;
//start_i=i;start_j=j;
count=0;
st=0;ecount=0;return;}
if(p1>=p2 && p1>=p3)
	maxval=p1;
else if(p2>=p1 && p2>=p3)
	maxval=p2;
else if(p3>=p1 && p3>=p2)
	maxval=p3;

//*(newimage+i*col+j)=0;
if(p1==maxval)
	{
              *(direc+count)=p1dir;count++;
              k=check_R1(count);// printf("returned %d\n",k);
	         if(k==1){edge_creation(count);*(direc)=p1dir;count=1;
//		printf("new_i %d new_j %d\n",new_i,new_j);
			{start_i=i;start_j=j;}	
		if(check(nexti1,nextj1)==1)
		{
		eresponse[ecount++]=expavg(nexti1,nextj1);	
		find3N(nexti1,nextj1,p1dir);
		}
		return;}
              else{
 //	      printf("calling 3N recursively\n");
		if(nexti1<row && nextj1<col)
              {v[vcount].i=i;v[vcount].j=j;vcount++;
		eresponse[ecount++]=expavg(nexti1,nextj1);			
		find3N(nexti1,nextj1,p1dir);return;}}
    }
if(p2==maxval)
{
              *(direc+count)=p2dir;count++;
              k=check_R1(count); //printf("returned %d\n",k);
              if(k==1){edge_creation(count);*(direc)=p2dir;count=1;
//			printf("new_i %d new_j %d\n",new_i,new_j);				
		{start_i=i;start_j=j;}
		if(check(nexti2,nextj2)==1){eresponse[ecount++]=expavg(nexti2,nextj2);	
		find3N(nexti2,nextj2,p2dir);}return;}
	      else
		{
 //	      printf("calling 3N recursively\n");
		if(nexti2<row && nextj2<col){
		v[vcount].i=i;v[vcount].j=j;vcount++;
		eresponse[ecount++]=expavg(nexti2,nextj2);
	      find3N(nexti2,nextj2,p2dir);return;}}
}
if(p3==maxval)
{
    *(direc+count)=p3dir;count++;
    k=check_R1(count);// printf("returned %d\n",k);
    if(k==1){edge_creation(count);*(direc)=p3dir;count=1;
//	printf("new_i %d new_j %d\n",new_i,new_j);		
if(check(nexti3,nextj3)==1){eresponse[ecount++]=expavg(nexti3,nextj3);	
	{start_i=i;start_j=j;}
	find3N(nexti3,nextj3,p3dir);}return;}
else {// printf("calling 3N recursively\n");
	if(nexti3<row && nextj3<col){
v[vcount].i=i;v[vcount].j=j;vcount++;
eresponse[ecount++]=expavg(nexti3,nextj3);
find3N(nexti3,nextj3,p3dir);return;}}
}

}

int check_R1(int ct)
{
// printf("entering check_R1\n");
if(ct==0 || ct==1)return 0;
 int second,edge_break,first_run,prev,second_run,nonsing_run,prev_sec,first,third;
             second=300;edge_break=0;
             first_run=prev=0;second_run=0;nonsing_run=0;
             for(x=0;x<ct;x++)// checking rule 1
             {
 //				printf("%d ",*(direc+x));
                              if(x==0)first=*(direc+x);
                              if(*(direc+x)!=first && second==300){second=*(direc+x);}
                              else if(*(direc+x)!=first && *(direc+x)!=second){edge_break=1;count=x;break;}
		
             }
 //printf("\n");
if(second==300)return 0;
             if(edge_break==1)
             {

 //			printf("more than 3.returning *%d\n",edge_break);                             
			 return edge_break;
             }
 //printf(".\n");//return 0;
             prev_sec=-1;prev=0;nonsing_run=0;x=0;first_run=0;
		while(*(direc+x)==first && x<ct){first_run++;x++;}//1
		while(*(direc+x)==second && x<ct){nonsing_run++;x++;}//2
		//if(first_run==1 || nonsing_run==1){;}
		//else return 1;
		if(first_run<nonsing_run)
		{
		prev=nonsing_run;//3
		nonsing_run=0;
		while(*(direc+x)==first && x<ct){//if(*(direc+x+1)==first)return 1;
		x++;}
		if(x==ct)return 0;
		while(*(direc+x)==second && x<ct){nonsing_run++;x++;}//3
		if(nonsing_run!=0)
		{
		if(x==ct)return 0;
		if(nonsing_run>=(prev-relax) && nonsing_run<=(prev+relax))prev_sec=nonsing_run;
		else return 1;
		}
		else
		{return 0;}
//printf("___");
		nonsing_run=0;
		while(x<ct)
		{
		nonsing_run=0;
			if(*(direc+x+1)==first && x!=ct-1){//printf("next is same\n");
				return 1;} else x++;
			//while(*(direc+x)==first && x<ct)x++;			
			if(x==ct)return 0;
			while(x<ct && *(direc+x)==second){nonsing_run++;x++;}
if(!(nonsing_run>=(prev-relax) && nonsing_run<=(prev+relax)))
{//printf("nonsing is %d. prev is %d.beyond bounds.1 \n",nonsing_run,prev);
return 1;}
			nonsing_run=0;
		}return 0;
		}
		else if(nonsing_run==1)
		{
		prev=first_run;//2
		first_run=0;
		while(*(direc+x)==second && x<ct){//if(*(direc+x+1)==second && x!=ct-1)return 1;
		x++;}
		if(x==ct)return 0;
		while(*(direc+x)==first && x<ct){first_run++;x++;}//3
//		printf("firstrun is %d and prev is %d\n",first_run,prev);
		/*if(first_run!=0)
		{
		if(x==ct || x==ct-1)return 0;
		if(first_run>=(prev-relax) && first_run<=(prev+relax))prev_sec=first_run;
		else return 1;
                }
		else
		{return 0;}//printf("prev_sec is %d\n",prev_sec);*/
		first_run=0;
		while(x<ct)
		{
			if(*(direc+x)==second && *(direc+x+1)==second && x+1<ct)return 1; else x++;
			while(x<ct && *(direc+x)==first){first_run++;x++;}//printf("firstrun is %d\n",first_run);
			if(x==ct )return 0;
			if(!(first_run>=(prev-relax) && first_run<=(prev+relax)))return 1;
			first_run=0;
		}
	return 0;
		}
		
// 	printf("returning %d\n",edge_break);                             
             return 1;

}

/*---------------------------------------*/
void get_the_mat(char *s,int r,int c)
                        {
                          FILE *fp;
                          int line_count=0;
                          int row=0,col=0;
                          unsigned char ch;
                          int i,j;
                         /*----------open the image file------------*/
                          fp=fopen(s,"rb");
                          tmat=(int*)malloc(r*c*sizeof(int));
                          while(fread(&ch,sizeof(ch),1,fp)==1) /* read all values and place in mat*/
                             {
                               if(ch=='\n') line_count++;
                                if(line_count==4) break;
                              }
                           while( fread(&ch,sizeof(ch),1,fp)==1)
                                {
                                  *(tmat+row*c+col)=ch;
                                     col++;
                                    if(col>(c-1)) /* matrix building*/
                                      {col=col%c;
                                       row=(row+1)%r;}
                                 }
				
                           fclose(fp);

                         } /*func. end*/
/*---------------------------------------*/
void build_image (int r,int c, char* s) /*construct the image to show curvature*/
                     {
                          FILE *fp1,*fp2,*fp;
                          unsigned char color[3];
                          int i,j;char mydata;
                          unsigned char ch;
                          char t[200];
                          int d,m;
			fp=fopen("prewittmatrix.txt","w");
                          char *mystring ="-output.pgm"; /* output image name*/
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
                          strcat(t,mystring);
                          fp2=fopen("myfile.dat","rb"); /* contains the header lines of the image*/
                          fp1=fopen(t,"wb+"); /*update t if needed- to change name*/
                          while(fread(&ch,sizeof(ch),1,fp2)==1)
                          fwrite(&ch,sizeof(ch),1,fp1);
                          for(i=0;i<r;i++)
                          for(j=0;j<c;j++)
                           {
                            mydata=*(newimage+i*c+j);
//255 255 255 is white; 0 0 0 is black; 255 0 0 is red; 0 255 0 is green; 0 0 255 is blue. 
				/*if(mydata==0 ){color[0]=0;color[1]=0;color[2]=0;}
				else if(mydata==10 ){color[0]=255;color[1]=255;color[2]=255;}
				else if(mydata==1){color[0]=255;color[1]=0;color[2]=0;}
				else if(mydata==2){color[0]=0;color[1]=255;color[2]=0;}*/
				                    printf("%d\n",mydata); 
				                     fprintf(fp1,"%c",mydata);                                                        
							
			   //fwrite(color,1,3,fp1);                      
				 }
                          

                         fclose(fp2);
                         fclose(fp1);
                       } /*func. end*/
/*---------------------------------------*/

