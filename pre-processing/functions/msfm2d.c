#include "mex.h"
#include "math.h"
#include "common.c"

/*
 * This function MSFM2D calculates the shortest distance from a list of
 * points to all other pixels in an image, using the
 *Multistencil Fast Marching Method (MSFM). This method gives more accurate
 *distances by using second order derivatives and cross neighbours.
 *
 *T=msfm2d(F, SourcePoints, UseSecond, UseCross)
 *
 *inputs,
 *  F: The speed image
 *  SourcePoints : A list of starting points [2 x N] (distance zero)
 *  UseSecond : Boolean Set to true if not only first but also second
 *               order derivatives are used (default)
 *  UseCross: Boolean Set to true if also cross neighbours
 *               are used (default)
 *outputs,
 *  T : Image with distance from SourcePoints to all pixels
 *
 *Function is written by D.Kroon University of Twente (June 2009)
 */

double CalculateDistance(double *T, double Fij, int *dims, int i, int j, bool usesecond, bool usecross, bool *Frozen) {
    /* Derivatives */
    double Tm[4]={0, 0, 0, 0};
    double Tm2[4]={0, 0, 0, 0};
    double Coeff[3];
    
    /* local derivatives in distance image */
    double Tpatch_2_3, Txm2, Tpatch_4_3, Txp2;
    double Tpatch_3_2, Tym2, Tpatch_3_4, Typ2;
    double Tpatch_2_2, Tr1m2, Tpatch_4_4, Tr1p2;
    double Tpatch_2_4, Tr2m2, Tpatch_4_2, Tr2p2;
    
    /* Return values root of polynomial */
    double ansroot[2]={0, 0};
    
    /* Loop variables  */
    int q, t;
    
    /* Derivative checks */
    bool ch1, ch2;
    
    /* Order derivatives in a certain direction */
    int Order[4]={0, 0, 0, 0};
    
    /* Current location */
    int in, jn;
    
    /* Constant cross term */
    const double c1=0.5;
    
    double Tt, Tt2;
    /*Get First order derivatives (only use frozen pixel)  */
    in=i-1; jn=j+0; if(isfrozen2d(in, jn, dims, Frozen)) { Tpatch_2_3=T[mindex2(in, jn, dims[0])]; } else { Tpatch_2_3=INF; }
    in=i+0; jn=j-1; if(isfrozen2d(in, jn, dims, Frozen)) { Tpatch_3_2=T[mindex2(in, jn, dims[0])]; } else { Tpatch_3_2=INF; }
    in=i+0; jn=j+1; if(isfrozen2d(in, jn, dims, Frozen)) { Tpatch_3_4=T[mindex2(in, jn, dims[0])]; } else { Tpatch_3_4=INF; }
    in=i+1; jn=j+0; if(isfrozen2d(in, jn, dims, Frozen)) { Tpatch_4_3=T[mindex2(in, jn, dims[0])]; } else { Tpatch_4_3=INF; }
    if(usecross) {
        in=i-1; jn=j-1; if(isfrozen2d(in, jn, dims, Frozen)) { Tpatch_2_2=T[mindex2(in, jn, dims[0])]; } else { Tpatch_2_2=INF; }
        in=i-1; jn=j+1; if(isfrozen2d(in, jn, dims, Frozen)) { Tpatch_2_4=T[mindex2(in, jn, dims[0])]; } else { Tpatch_2_4=INF; }
        in=i+1; jn=j-1; if(isfrozen2d(in, jn, dims, Frozen)) { Tpatch_4_2=T[mindex2(in, jn, dims[0])]; } else { Tpatch_4_2=INF; }
        in=i+1; jn=j+1; if(isfrozen2d(in, jn, dims, Frozen)) { Tpatch_4_4=T[mindex2(in, jn, dims[0])]; } else { Tpatch_4_4=INF; }
    }
    /*The values in order is 0 if no neighbours in that direction  */
    /*1 if 1e order derivatives is used and 2 if second order  */
    /*derivatives are used  */
    Order[0]=0; Order[1]=0; Order[2]=0; Order[3]=0;
    /*Make 1e order derivatives in x and y direction  */
    Tm[0] = min( Tpatch_2_3 , Tpatch_4_3); if(IsFinite(Tm[0])){ Order[0]=1; }
    Tm[1] = min( Tpatch_3_2 , Tpatch_3_4); if(IsFinite(Tm[1])){ Order[1]=1; }
    /*Make 1e order derivatives in cross directions  */
    if(usecross) {
        Tm[2] = min( Tpatch_2_2 , Tpatch_4_4); if(IsFinite(Tm[2])){ Order[2]=1; }
        Tm[3] = min( Tpatch_2_4 , Tpatch_4_2); if(IsFinite(Tm[3])){ Order[3]=1; }
    }
	
    /*Make 2e order derivatives  */
    if(usesecond) {
        /*Get Second order derivatives (only use frozen pixel) */
        in=i-2; jn=j+0; if(isfrozen2d(in, jn, dims, Frozen)) { Txm2=T[mindex2(in, jn, dims[0])]; } else { Txm2=INF; }
        in=i+2; jn=j+0; if(isfrozen2d(in, jn, dims, Frozen)) { Txp2=T[mindex2(in, jn, dims[0])]; } else { Txp2=INF; }
        in=i+0; jn=j-2; if(isfrozen2d(in, jn, dims, Frozen)) { Tym2=T[mindex2(in, jn, dims[0])]; } else { Tym2=INF; }
        in=i+0; jn=j+2; if(isfrozen2d(in, jn, dims, Frozen)) { Typ2=T[mindex2(in, jn, dims[0])]; } else { Typ2=INF; }
        if(usecross) {
            in=i-2; jn=j-2; if(isfrozen2d(in, jn, dims, Frozen)) { Tr1m2=T[mindex2(in, jn, dims[0])]; } else { Tr1m2=INF; }
            in=i-2; jn=j+2; if(isfrozen2d(in, jn, dims, Frozen)) { Tr2m2=T[mindex2(in, jn, dims[0])]; } else { Tr2m2=INF; }
            in=i+2; jn=j-2; if(isfrozen2d(in, jn, dims, Frozen)) { Tr2p2=T[mindex2(in, jn, dims[0])]; } else { Tr2p2=INF; }
            in=i+2; jn=j+2; if(isfrozen2d(in, jn, dims, Frozen)) { Tr1p2=T[mindex2(in, jn, dims[0])]; } else { Tr1p2=INF; }
        }
        
        Tm2[0]=0; Tm2[1]=0;Tm2[2]=0; Tm2[3]=0;
        /*pixels with a pixeldistance 2 from the center must be */
        /*lower in value otherwise use other side or first order */
        ch1=(Txm2<Tpatch_2_3)&&IsFinite(Tpatch_2_3); ch2=(Txp2<Tpatch_4_3)&&IsFinite(Tpatch_4_3);
        if(ch1&&ch2) {
            Tm2[0] =min( (4.0*Tpatch_2_3-Txm2)/3.0 , (4.0*Tpatch_4_3-Txp2)/3.0);  Order[0]=2;
        }
        else if (ch1) {
            Tm2[0]=(4.0*Tpatch_2_3-Txm2)/3.0; Order[0]=2;
        }
        else if(ch2) {
            Tm2[0] =(4.0*Tpatch_4_3-Txp2)/3.0; Order[0]=2;
        }
        
        ch1=(Tym2<Tpatch_3_2)&&IsFinite(Tpatch_3_2); ch2=(Typ2<Tpatch_3_4)&&IsFinite(Tpatch_3_4);
        
        if(ch1&&ch2) {
            Tm2[1] =min( (4.0*Tpatch_3_2-Tym2)/3.0 , (4.0*Tpatch_3_4-Typ2)/3.0); Order[1]=2;
        }
        else if(ch1) {
            Tm2[1]=(4.0*Tpatch_3_2-Tym2)/3.0; Order[1]=2;
        }
        else if(ch2) {
            Tm2[1]=(4.0*Tpatch_3_4-Typ2)/3.0; Order[1]=2;
        }
        if(usecross) {
            ch1=(Tr1m2<Tpatch_2_2)&&IsFinite(Tpatch_2_2); ch2=(Tr1p2<Tpatch_4_4)&&IsFinite(Tpatch_4_4);
            if(ch1&&ch2) {
                Tm2[2] =min( (4.0*Tpatch_2_2-Tr1m2)/3.0 , (4.0*Tpatch_4_4-Tr1p2)/3.0); Order[2]=2;
            }
            else if(ch1) {
                Tm2[2]=(4.0*Tpatch_2_2-Tr1m2)/3.0; Order[2]=2;
            }
            else if(ch2){
                Tm2[2]=(4.0*Tpatch_4_4-Tr1p2)/3.0; Order[2]=2;
            }
            
            ch1=(Tr2m2<Tpatch_2_4)&&IsFinite(Tpatch_2_4); ch2=(Tr2p2<Tpatch_4_2)&&IsFinite(Tpatch_4_2);
            if(ch1&&ch2){
                Tm2[3] =min( (4.0*Tpatch_2_4-Tr2m2)/3.0 , (4.0*Tpatch_4_2-Tr2p2)/3.0); Order[3]=2;
            }
            else if(ch1) {
                Tm2[3]=(4.0*Tpatch_2_4-Tr2m2)/3.0; Order[3]=2;
            }
            else if(ch2) {
                Tm2[3]=(4.0*Tpatch_4_2-Tr2p2)/3.0; Order[3]=2;
            }
        }
    }
    /*Calculate the distance using x and y direction */
    Coeff[0]=0; Coeff[1]=0; Coeff[2]=-1/(max(pow2(Fij),eps));
    
    for (t=0; t<2; t++) {
        switch(Order[t]) {
            case 1:
                Coeff[0]+=1; Coeff[1]+=-2*Tm[t]; Coeff[2]+=pow2(Tm[t]);
                break;
            case 2:
                Coeff[0]+=(2.2500); Coeff[1]+=-2.0*Tm2[t]*(2.2500); Coeff[2]+=pow2(Tm2[t])*(2.2500);
                break;
        }
    }
    roots(Coeff, ansroot);
    Tt=max(ansroot[0], ansroot[1]);
    /*Calculate the distance using the cross directions */
    if(usecross) {
        /* Original Equation */
        /*    Coeff[0]=0; Coeff[1]=0; Coeff[2]=-1/(max(pow2(Fij),eps)) */
        Coeff[0]+=0; Coeff[1]+=0; Coeff[2]+=-1/(max(pow2(Fij),eps));
        for (t=2; t<4; t++) {
            switch(Order[t]) {
                case 1:
                    Coeff[0]+=c1; Coeff[1]+=-2.0*c1*Tm[t]; Coeff[2]+=c1*pow2(Tm[t]);
                    break;
                case 2:
                    Coeff[0]+=c1*2.25; Coeff[1]+=-2*c1*Tm2[t]*(2.25); Coeff[2]+=pow2(Tm2[t])*c1*2.25;
                    break;
            }
        }
        if(Coeff[0]>0) {
            roots(Coeff, ansroot);
            Tt2=max(ansroot[0], ansroot[1]);
            /*Select minimum distance value of both stensils */
            Tt=min(Tt, Tt2);
        }
    }
    /*Upwind condition check, current distance must be larger */
    /*then direct neighbours used in solution */
    /*(Will this ever happen?) */
    if(usecross) {
        for(q=0; q<4; q++) { 
            if(IsFinite(Tm[q])&&(Tt<Tm[q])) 
            { 
                Tt=Tm[minarray(Tm, 4)]+(1/(max(Fij,eps)));
            }
        }
    }
    else {
        for(q=0; q<2; q++)
        { 
            if(IsFinite(Tm[q])&&(Tt<Tm[q])) {
                Tt=Tm[minarray(Tm, 2)]+(1/(max(Fij,eps)));}
        }
    }
    return Tt;
}

/* The matlab mex function  */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] ) {
    /* The input variables */
    double *F, *SourcePoints;
    bool *useseconda, *usecrossa;
    bool usesecond=true;
    bool usecross=true;
    
    /* The output distance image */
    double *T;
    
    /* Euclidian distance image */
    double *Y;
    
    /* Current distance values */
    double Tt, Ty;
    
    /* Matrix containing the Frozen Pixels" */
    bool *Frozen;
    
    /* Augmented Fast Marching (For skeletonize) */
    bool Ed;
    
    /* Size of input image */
    const mwSize *dims_c;
    mwSize dims[2];
    
    /* Size of  SourcePoints array */
    const mwSize *dims_sp_c;
    mwSize dims_sp[2];
    
    /* Number of pixels in image */
    int npixels;
    
    /* Neighbour list */
    int neg_free;
    int neg_pos;
    double *neg_listv;
    double *neg_listx;
    double *neg_listy;
    double *neg_listo;
    
    int *listprop;
    double **listval;
    
    /* Neighbours 4x2 */
    int ne[8]={-1, 1, 0, 0, 0, 0, -1, 1};

    /* Loop variables  */
    int z, k, itt, q;
    
    /* Current location */
    int x, y, i, j;
    
    /* Index */
    int IJ_index, XY_index, index;
    
    /* Check for proper number of input and output arguments. */
    if(nrhs<3) {
        mexErrMsgTxt("2 to 4 inputs are required.");
    }
    if(nlhs==1) { Ed=0; }
    else if (nlhs==2) { Ed=1; }
    else {
        mexErrMsgTxt("One or Two outputs required");
    }
    
    /* Check data input types  */
    if(mxGetClassID(prhs[0])!=mxDOUBLE_CLASS) {
        mexErrMsgTxt("Speed image must be of class double");
    }
    if(mxGetClassID(prhs[1])!=mxDOUBLE_CLASS) {
        mexErrMsgTxt("SourcePoints must be of class double");
    }
    
    if((nrhs>2)&&(mxGetClassID(prhs[2])!= mxLOGICAL_CLASS)) {
        mexErrMsgTxt("UseSecond must be of class boolean / logical");
    }
    if((nrhs>3)&&(mxGetClassID(prhs[3])!= mxLOGICAL_CLASS)) {
        mexErrMsgTxt("UseCross must be of class boolean / logical");
    }
        
    /* Get the sizes of the input image */
    if(mxGetNumberOfDimensions(prhs[0])==2) {
        dims_c = mxGetDimensions(prhs[0]);
        dims[0]=dims_c[0]; dims[1]=dims_c[1];
        npixels=dims[0]*dims[1];
    }
    else {
        mexErrMsgTxt("Speed image must be 2d.");
    }
    
    /* Get the sizes of the  SourcePoints */
    dims_sp_c = mxGetDimensions(prhs[1]);
    
    if(dims_sp_c[0]!=2) {
        mexErrMsgTxt("SourcePoints must be a 2xn matrix.");
    }
    dims_sp[0]=dims_sp_c [0]; dims_sp[1]=dims_sp_c[1];
    
    
    /* Get pointers/data from  to each input. */
    F=(double*)mxGetPr(prhs[0]);
    SourcePoints=(double*)mxGetPr(prhs[1]);
    if(nrhs>2){ useseconda = (bool*)mxGetPr(prhs[2]); usesecond=useseconda[0];}
    if(nrhs>3){ usecrossa = (bool*)mxGetPr(prhs[3]); usecross=usecrossa[0];}
    
    /* Create the distance output array */
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    /* Assign pointer to output.  */
    /*Distance image, also used to store the index of narrowband pixels  */
    /*during marching process  */
    T= mxGetPr(plhs[0]);
    if(Ed) { 
        plhs[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        Y= mxGetPr(plhs[1]);
    }
    
    /* Pixels which are processed and have a final distance are frozen */
    Frozen = (bool*)malloc( npixels* sizeof(bool) );
    for(q=0;q<npixels;q++){Frozen[q]=0; T[q]=-1;}
    if(Ed)
    {
    for(q=0;q<npixels;q++){Y[q]=-1;}
    }
    
    /*Free memory to store neighbours of the (segmented) region */
    neg_free = 100000;
    neg_pos=0;
    neg_listx = (double *)malloc( neg_free*sizeof(double) );
    neg_listy = (double *)malloc( neg_free*sizeof(double) );
    if(Ed) {
        neg_listo = (double *)malloc( neg_free*sizeof(double) );
        for(q=0;q<neg_free;q++) { neg_listo[q]=0; }
    }
    
    /* List parameters array */
    listprop=(int*)malloc(3* sizeof(int));
    /* Make jagged list to store a maximum of 2^64 values */
    listval= (double **)malloc( 64* sizeof(double *) );
    /* Initialize parameter list */
    initialize_list(listval, listprop);
    neg_listv=listval[listprop[1]-1];
    
    /*(There are 3 pixel classes:
     *  - frozen (processed)
     *  - narrow band (boundary) (in list to check for the next pixel with smallest distance)
     *  - far (not yet used)
     */
    
    /* set all starting points to distance zero and frozen  */
    /* and add all neighbours of the starting points to narrow list  */
    for (z=0; z<dims_sp[1]; z++) {
        /*starting point  */
        x= (int)SourcePoints[0+z*2]-1;
        y= (int)SourcePoints[1+z*2]-1;
        XY_index=x+y*dims[0];
        
        /*Set starting point to frozen and distance to zero  */
        Frozen[XY_index]=1;
        T[XY_index]=0;
        if(Ed) { Y[XY_index]=0; }
    }
    
    for (z=0; z<dims_sp[1]; z++) {
        /*starting point  */
        x= (int)SourcePoints[0+z*2]-1;
        y= (int)SourcePoints[1+z*2]-1;
        XY_index=x+y*dims[0];
        
        /* Add neigbours of starting points  */
        for (k=0; k<4; k++) {
            /*Location of neighbour  */
            i=x+ne[k]; j=y+ne[k+4];
            IJ_index=i+j*dims[0];
            
            /*Check if current neighbour is not yet frozen and inside the
             *picture  */
            if(isntfrozen2d(i, j, dims, Frozen)) {
                Tt=(1/(max(F[IJ_index],eps)));
					
                Ty=1;
                /*Update distance in neigbour list or add to neigbour list */
                if(T[IJ_index]>0) {
                    if(neg_listv[(int)T[IJ_index]]>Tt) {
                        listupdate(listval, listprop, (int)T[IJ_index], Tt);
                    }
                    if(Ed)
                    {
                        neg_listo[(int)T[IJ_index]]=min(neg_listo[(int)T[IJ_index]],Ty);
                    }
                }
                else {
                    /*If running out of memory at a new block  */
                    if(neg_pos>=neg_free) {
                        neg_free+=100000;
                        neg_listx = (double *)realloc(neg_listx, neg_free*sizeof(double) );
                        neg_listy = (double *)realloc(neg_listy, neg_free*sizeof(double) );
                        if(Ed) {
                            neg_listo = (double *)realloc(neg_listo, neg_free*sizeof(double) );
                        }
                    }
                    list_add(listval, listprop, Tt);
                    neg_listv=listval[listprop[1]-1];
                    neg_listx[neg_pos]=i;
                    neg_listy[neg_pos]=j;
                    if(Ed){
                        neg_listo[neg_pos]=Ty;
                    }
                    T[IJ_index]=neg_pos;
                    neg_pos++;
                }
            }
        }
    }
    /*Loop through all pixels of the image  */
    for (itt=0; itt<npixels; itt++) {
        /*Get the pixel from narrow list (boundary list) with smallest
         *distance value and set it to current pixel location  */
        index=list_minimum(listval, listprop);
        neg_listv=listval[listprop[1]-1];
 		
        /* Stop if pixel distance is infinite (all pixels are processed)  */
        if(IsInf(neg_listv[index])) {  break; }
        x=(int)neg_listx[index]; y=(int)neg_listy[index];
        
        XY_index=x+y*dims[0];
        Frozen[XY_index]=1;
        T[XY_index]=neg_listv[index];
        if(Ed) { Y[XY_index]=neg_listo[index]; }
      
     
        /*Remove min value by replacing it with the last value in the array  */
        list_remove_replace(listval, listprop, index) ;
        neg_listv=listval[listprop[1]-1];
        if(index<(neg_pos-1)) {
            neg_listx[index]=neg_listx[neg_pos-1];
            neg_listy[index]=neg_listy[neg_pos-1];
            if(Ed){
                neg_listo[index]=neg_listo[neg_pos-1];
            }
            T[(int)(neg_listx[index]+neg_listy[index]*dims[0])]=index;
        }
        neg_pos =neg_pos-1;
       
    
        /*Loop through all 4 neighbours of current pixel  */
        for (k=0;k<4;k++) {
            
            /*Location of neighbour  */
            i=x+ne[k]; j=y+ne[k+4];
            IJ_index=i+j*dims[0];

            /*Check if current neighbour is not yet frozen and inside the  */
            /*picture  */
            if(isntfrozen2d(i, j, dims, Frozen)) {
				
                Tt=CalculateDistance(T, F[IJ_index], dims, i, j, usesecond, usecross, Frozen);
				        
				if(Ed) {
                    Ty=CalculateDistance(Y, 1, dims, i, j, usesecond, usecross, Frozen);
                }

                /*Update distance in neigbour list or add to neigbour list */
                IJ_index=i+j*dims[0];
                if((T[IJ_index]>-1)&&T[IJ_index]<=listprop[0]) {
                    if(neg_listv[(int)T[IJ_index]]>Tt) {
                        listupdate(listval, listprop,    (int)T[IJ_index], Tt);
                    }
                    if(Ed)
                    {
                        neg_listo[neg_pos]=min(neg_listo[neg_pos],Ty);
                    }
                }
                else {
                    /*If running out of memory at a new block */
                    if(neg_pos>=neg_free) {
                        neg_free+=100000;
                        neg_listx = (double *)realloc(neg_listx, neg_free*sizeof(double) );
                        neg_listy = (double *)realloc(neg_listy, neg_free*sizeof(double) );
                        if(Ed) {
                            neg_listo = (double *)realloc(neg_listo, neg_free*sizeof(double) );
                        }
                    }
                    list_add(listval, listprop, Tt);
                    neg_listv=listval[listprop[1]-1];
                    neg_listx[neg_pos]=i; neg_listy[neg_pos]=j;
                    if(Ed) {
                        neg_listo[neg_pos]=Ty;
                    }
                    T[IJ_index]=neg_pos;
                    neg_pos++;
                }
            }
        }
        
    }
    /* Free memory */
    /* Destroy parameter list */
    destroy_list(listval, listprop);
    free(neg_listx);
    free(neg_listy);
    if(Ed) {
        free(neg_listo);
    }
    free(Frozen);
}




