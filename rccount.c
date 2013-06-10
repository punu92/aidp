int * count_rc(char*);
//------------------------------------

int* count_rc(char* filename)
   {
     FILE *fp;
     int line_count=0;
     char ch;
     int i=0,j=0;
     char str[50];
     int c;
     int a[2];
     char r[25],s[25];
     fp=fopen(filename,"rb");
     while(fread(&ch,sizeof(ch),1,fp)==1)
     {
        if(ch=='\n') line_count++;
        if (line_count==2) break;
     }
     while(fread(&ch,sizeof(ch),1,fp)==1)
     {
       if(ch!='\n'&&ch!=' ')
       { str[i]=ch;i++;}
       else if(ch==' '&&str[i-1]!='*') {str[i]='*';i++;}
       else if(ch==' '&&str[i-1]=='*'){}
       else if (ch=='\n') {str[i]='\0'; break;}
     }
     fclose(fp);
     c=0;i=0;
     while(str[c]!='*')
     {
         r[i]=str[c];i++;c++;
      }
     r[i]='\0';
     c++;
     while(str[c]!='\0')
     {
             s[j]=str[c];j++;c++;
      }
     s[j]='\0';
     a[0]=atoi(r);
     a[1]=atoi(s);   /// return array a[0]= row , a[1]= col.
  return (&a[0]);
}




//98360 65336
