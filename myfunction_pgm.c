int  build_my_file(int,int);
char mychar(int);
//---------------------------------

   int build_my_file(int r, int c)
{         
	  char ch[50];
	  int i,d;
	  char mc;
	  FILE *fp;
	  char rs[25],cs[25];
	  int rs_limit;
	  int cs_limit;
	  fp=fopen("myfile.dat","wb+");
	  ch[0]='P';
	  ch[1]='5';// for ppm file this character will be 6.
	  ch[2]='\n';
	  ch[3]='#';
	  ch[4]='t';
	  ch[5]='\n';
	  ch[6]='\0';
	   //-------------------------------
	  
	   i=0;
	   while(ch[i]!='\0')
	   {      mc=ch[i];
		 fwrite(&mc,sizeof(mc),1,fp);
		 i++;
	   }	 
	//------------------------------------------
	i=0;
	while(c>=10)
	{
		d=c%10;
		c=c/10;
		cs[i]=mychar(d);
		i++;
	}
	cs[i]=mychar(c);
	cs_limit=i;
	for(i=cs_limit;i>=0;i--)
	{
		mc=cs[i];
		fwrite(&mc,sizeof(mc),1,fp);
	}
	mc =' ';
	fwrite(&mc,sizeof(mc),1,fp);

	//-------------------------------------
	     i=0;
          while(r>=10)
	  {
          d=r%10;
	  r=r/10;
	  rs[i]=mychar(d);
	  i++;
	  }
	  rs[i]=mychar(r);
	  rs_limit=i;

	  for(i=rs_limit;i>=0;i--)
	  { mc=rs[i];
		  fwrite(&mc,sizeof(mc),1,fp);
	  }
	  mc='\n';
	  fwrite(&mc,sizeof(mc),1,fp);
//----------------------------------------------
     ch[1]='2';
     ch[2]='5';
     ch[3]='5';
     for(i=1;i<=3;i++){mc=ch[i];
	     fwrite(&mc,sizeof(mc),1,fp);}
	     mc='\n';
	     fwrite(&mc,sizeof(mc),1,fp);
	     
	  fclose(fp);


	   
	
return 0;
}
      
	   

 char mychar(int x)
{
	 char c;
	 if (x==1) c='1';
	 if(x==2) c='2';
         if(x==3) c='3';
	 if(x==4) c='4';
	 if(x==5) c='5';
	 if(x==6) c='6';
	 if(x==7) c='7';
	 if(x==8) c='8';
	 if(x==9) c='9';
	 if(x==0) c='0';
	  return c;

}
