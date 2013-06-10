
/*--------------- Bresenham Line ------------------*/ 
/*cdata: matrix where the line will be drawn */ 

void bresenham_line(int x1,int y1,int x2,int y2,int c)
{
	int dx,dy,x,y,d,dE,dNE;
	dx=(x1-x2<0)?x2-x1:x1-x2;
	dy=(y1-y2<0)?y2-y1:y1-y2;

	if(dx>=dy)
	{
		if(y1<y2)
		{
			x1=x1+x2;
			x2=x1-x2;
			x1=x1-x2;

			y1=y1+y2;
			y2=y1-y2;
			y1=y1-y2;
		}

		if(x1<x2)	/* quad 1 */
		{
			/*printf("quad 1.\n");*/

			d=2*dy+dx;
			dx=x2-x1;
			dy=y2-y1;
			dE=2*dy;
			dNE=2*(dy+dx);
			x=x1,y=y1;
			
			//if (*(cdata+x*c+y)==0)  // set according to your requirement, 
                           *(cdata+x*c+y)=9;    // my background : *(cdata+x*c+y)==0, you may have diff. val.
                                                      // to avoid overwriting
			         
			while(x<x2)
			{
				if(d>=0)
				{
					x++;
					d+=dE;
				}
				else
				{
					x++;
					y--;
					d+=dNE;
				}
				
				//if (*(cdata+x*c+y)==0) 
			           *(cdata+x*c+y)=9;
			}

		}
		else		/* quad 4*/
		{
			/*printf("quad 4.\n");*/

			dx=x2-x1;
			dy=y2-y1;
			d=dx-2*dy;
			dE=(2)*(-dy);
			dNE=2*(dx-dy);
			x=x1,y=y1;
			
                //    if (*(cdata+x*c+y)==0) 
			  *(cdata+x*c+y)=9;
			while(x>x2)
			{
				if(d<=0)	/* W */
				{
					x--;
					d+=dE;
				}
				else		/* NW*/
				{
					x--;
					y--;
					d+=dNE;
				}
				
		//		if (*(cdata+x*c+y)==0) 
			  *(cdata+x*c+y)=9; 
			}
		}

	}

	else
	{
		if(y1<y2)
		{
			x1=x1+x2;
			x2=x1-x2;
			x1=x1-x2;

			y1=y1+y2;
			y2=y1-y2;
			y1=y1-y2;
		
		}
		
		if(x1<=x2)
		{
			/*printf("quad 2.\n");*/
			dx=x2-x1;
			dy=y2-y1;
			d=dy+2*dx;
			dE=2*dx;
			dNE=2*(dy+dx);
			x=x1,y=y1;
			
                  // if (*(cdata+x*c+y)==0) 
			  *(cdata+x*c+y)=9;
			while(y>y2)
			{
				if(d<=0)
				{
					y--;
					d+=dE;
				}

				else
				{
					x++;
					y--;
					d+=dNE;
				}

				
			//	if (*(cdata+x*c+y)==0) 
			  *(cdata+x*c+y)=9;
			}

		}
		else
		{
			/*printf("quad 3.\n");*/

			dx=x2-x1;
			dy=y2-y1;
			d=2*dx-dy;
			dE=2*dx;
			dNE=2*(dx-dy);
			x=x1,y=y1;
			
			
			//if (*(cdata+x*c+y)==0) 
			  *(cdata+x*c+y)=9;

			while(y>y2)
			{
				if(d>=0)
				{
					y--;
					d+=dE;
				}
				else
				{
					x--;
					y--;
					d+=dNE;
				}
				
				
			//	if (*(cdata+x*c+y)==0) 
			  *(cdata+x*c+y)=9;
			}

		}
	}

       }
