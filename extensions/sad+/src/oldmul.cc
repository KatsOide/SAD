
map_da oldmul(const map_da& x,const map_da& y)
{
   map_da z;
   da z1,z2,z3,z4,z5,z6;
   
   z=0.;
   //   Nvr2=Nvar/2;
   int Ncv2=N_cv/2;
   int ic1,ic2,Iaddr,i;

   z1=1.;
   for(int ipz=0;ipz<=N_ord;ipz++) {
      *(iv2+2)=ipz;
      //   z1=pz^ipz
      if(ipz>0) z1*=y.m[Ncv2+2];

      z2=z1;
      for(int ipy=0;ipy<=N_ord-ipz;ipy++) {
	 *(iv2+1)=ipy;
      //   z2=pz^ipz*py^ipy
	 if(ipy>0) z2*=y.m[Ncv2+1];
	 
	 z3=z2;
         for(int ipx=0;ipx<=N_ord-ipz-ipy;ipx++) {
	    *(iv2)=ipx;
      //  z3=pz^ipz*py^ipy*px^ipx
	    if(ipx>0) z3*=y.m[Ncv2];
            ic2=ccal(Nord,Nvr2,iv2);

	    z4=z3;
	    for(int iz=0;iz<=N_ord-ipz-ipy-ipx;iz++) {
	       *(iv1+2)=iz;
	       if(iz>0) z4*=y.m[2];
	       
	       z5=z4;
	       for(int iy=0;iy<=N_ord-ipz-ipy-ipx-iz;iy++) {
	          *(iv1+1)=iy;
	          if(iy>0) z5*=y.m[1];
		  
		  z6=z5;
		  for(int ix=0;ix<=N_ord-ipz-ipy-ipx-iz-iy;ix++) {
		     *iv1=ix;
		     if(ix>0) z6*=y.m[0];
		     ic1=ccal(Nord,Nvr2,iv1);
		     Iaddr=*(d1+ic1)+*(d2+ic2);
		     for(i=0;i<N_cv;i++) z.m[i]+=x.m[i].co(Iaddr)*z6;
		  }
	       }
	    }
	 }
      }
   }
   return z;
}


map_da oldmul(const map_da& x,const map_da& y,int No)
{
   map_da z;
   da z1,z2,z3,z4,z5,z6;
   
   z=0.;
   //   Nvr2=Nvar/2;
   int Ncv2=N_cv/2;
   int ic1,ic2,Iaddr,i;

   z1=1.;
   for(int ipz=0;ipz<=N_ord;ipz++) {
      *(iv2+2)=ipz;
      //   z1=pz^ipz
      if(ipz>0) z1=muln(z1,y.m[Ncv2+2],No);

      z2=z1;
      for(int ipy=0;ipy<=N_ord-ipz;ipy++) {
	 *(iv2+1)=ipy;
      //   z2=pz^ipz*py^ipy
	 if(ipy>0) z2=muln(z2,y.m[Ncv2+1],No);
	 
	 z3=z2;
         for(int ipx=0;ipx<=N_ord-ipz-ipy;ipx++) {
	    *(iv2)=ipx;
      //  z3=pz^ipz*py^ipy*px^ipx
	    if(ipx>0) z3=muln(z3,y.m[Ncv2],No);
            ic2=ccal(Nord,Nvr2,iv2);

	    z4=z3;
	    for(int iz=0;iz<=N_ord-ipz-ipy-ipx;iz++) {
	       *(iv1+2)=iz;
	       if(iz>0) z4=muln(z4,y.m[2],No);
	       
	       z5=z4;
	       for(int iy=0;iy<=N_ord-ipz-ipy-ipx-iz;iy++) {
	          *(iv1+1)=iy;
	          if(iy>0) z5=muln(z5,y.m[1],No);
		  
		  z6=z5;
		  for(int ix=0;ix<=N_ord-ipz-ipy-ipx-iz-iy;ix++) {
		     *iv1=ix;
		     if(ix>0) z6=muln(z6,y.m[0],No);
		     ic1=ccal(Nord,Nvr2,iv1);
		     Iaddr=*(d1+ic1)+*(d2+ic2);
		     for(i=0;i<N_cv;i++) z.m[i]+=x.m[i].co(Iaddr)*z6;
		  }
	       }
	    }
	 }
      }
   }
   return z;
}
