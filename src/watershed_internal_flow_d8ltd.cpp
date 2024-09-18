// C/C++ code
// Author: Emanuele Cordano
// Date: February 2021
// Flow direction D8-LTD algorithm implementation (under development)
//


//See this former code: 
//line 217 qui: https://github.com/TheHortonMachine/hortonmachine/blob/master/hmachine/src/main/java/org/hortonmachine/hmachine/modules/geomorphology/draindir/OmsDrainDir.java per la prescisione 🙂

// TO BE IMPLEMENTED
#include "spatRaster.h"
#include "NA.h"
#include <cmath> // per la funzione pow
#include "watershed_internal_flow_d8ltd.h"
#include "watershed_internal.h"
#define NODATA_2 -999999

// ESRI terra covension:
//   //     # x-1 x 	x+1
//   // y-1 ## 32	64	128
//   // y   ## 16	0	1
//   // y+1 ## 8	4	2
//   

//  Li at al, 2022's convection 
//   //     # x-1 x 	x+1
//   // y-1 ## 1	2	3
//   // y   ## 8	0	4
//   // y+1 ## 7	6	5

/*
 * 
 * nx,ny number of cols,rows 
 * x,y    row,col
 * pv flow direction value 
 * conv_type convention type used , see code
 * 
 
 * 
 */
int nextcell_point_conv1(int nx, int ny,int x,int y,double pv,int conv_type) {
  
  
  
   
    std::vector<double> idir={0,1,2,4,8,16,32,64,128}; // ESRI conv: 0 cell, ... (clockwise from right) 
    int iout=offset(nx,ny,x, y);
    if (conv_type==2) {
     // Li at al 2021 (from Orlandini at al 20)
     idir={0,4,5,6,7,8,1,2,3}; // Li-Orlandini conv 0 cell, ... (clockwise from right) 
    }
    
    if (inRaster(nx, ny, x + 1, y) && pv==idir[1]) {
        iout=offset(nx, ny, x + 1, y);
    } else if (inRaster(nx, ny, x+1,y+1) && pv==idir[2]) {
        
      iout=offset(nx, ny, x+1, y+1);
    } else if (inRaster(nx, ny, x,y+1) && pv==idir[3]) {
      iout=offset(nx, ny, x, y+1);
    } else if (inRaster(nx, ny, x-1,y+1) && pv==idir[4]) {
      iout=offset(nx, ny, x-1, y+1);
    } else if (inRaster(nx, ny, x-1,y) && pv==idir[5]) {
      iout=offset(nx, ny, x-1, y);
    } else if (inRaster(nx, ny, x-1,y-1) && pv==idir[6]) {
      iout=offset(nx, ny, x-1, y-1);
    } else if (inRaster(nx, ny, x,y-1) && pv==idir[7]) {
      iout=offset(nx, ny, x, y-1);
    } else if (inRaster(nx, ny, x+1,y-1) && pv==idir[8]) {
      iout=offset(nx, ny, x+1, y-1);
    } else {
      iout=offset(nx, ny, x, y);
    }
    //printf("iout=%d\n",iout);
    return(iout);
    }
  
  /*
   * 
   * e elevation
   * nx,ny number of cols,rows
   * sr SR in Li et, 2022 (https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021WR030948)
   * sm SM in Li et al,2022 (https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021WR030948)
   * sfacet drainage facet 
   * tdc TD_c in Li et, 2022 (https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021WR030948)
   * tdd TD_c in Li et, 2022 (https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021WR030948)
   *
   * 
   */

void slope_direction(double* e, int nx, int ny, double *sr,double *sm,int *sfacet,
                     double* tdc, double *tdd,double L,
                     std::vector<double> ddp1,std::vector<double> ddp2,int nncell,int conv_type) {
  
  int i;
  ///int L=1; // dx=dy=L=1
  int ncell=nx*ny;
  double x,y;
 // int efacet=0;
  double e0,e1,e2;
  int facet=0;
//  std::vector<double> ddp1 = {2,2,4,4,6,6,8,8}; // this routine uses Orlandini-Li et, 2022's  convention
//  std::vector<double> ddp2 = {1,3,3,5,5,7,7,1};

 // int nncell=8;
  //std::vector<double> vOut(nx*ny,0);
  //std::vector<double> slope_mgn(ncell,0);
  
 // double slope=
 
 
  for (int i = 0; i < nx*ny; i++) {
    
   
    e0=*(e+i);
  
    
    std::vector<int> offset1(nncell,i);
    std::vector<int> offset2(nncell,i);
    x = getCol(nx, ny, i);   // ATTENTION: base 0 or 1?
    y = getRow(nx, ny, i);
    double slope_mgn_temp=0;
    double slope_mgn=0;
    int facet=0;
 
    
    
    for (int j=0; j<nncell; j++){
      
      //if (ddp1[j]==1) {
   /// printf("ba2 j=%d facet=%d i=%d e0=%f e1=%f e2=%f\n",j,facet,i,e0,e1,e2);
    e1=*(e+nextcell_point_conv1(nx,ny,x,y,ddp1[j],conv_type));
    e2=*(e+nextcell_point_conv1(nx,ny,x,y,ddp2[j],conv_type));
    //if (std::isnan(e1) | std::isnan(e0) | e1==e0) e1=1.05*e0; 
    //if (std::isnan(e2) | std::isnan(e0) | e2==e0) e2=1.05*e0;
       
    /// facet=j;
    double mean_e=(e0+e1+e2)/3; //EXPERIMENTAL    /////pow(pow(e0-e1,2)+pow(e1-e2,2),0.5)/L;
    slope_mgn_temp=pow(pow(e0-e1,2)+pow(e1-e2,2),0.5)/L;
    if (e0<mean_e) {
      
      slope_mgn_temp=0;
      
    }
    double flow_angle_tan=(e1-e2)/(e0-e1);
    if (flow_angle_tan>1)  {
      slope_mgn_temp=0;
    } else if (flow_angle_tan<0) {
      
      slope_mgn_temp=0;
    } // 20240913
    
    //if (slope_mgn_temp<0) slope_mgn_temp=0;
   
    if (slope_mgn_temp>slope_mgn) {
        slope_mgn=slope_mgn_temp;
        facet=j;
      
       }
    }
    
    double e_temp=e0;
    if (facet==0)  for (int j=0; j<nncell; j++){
      
      //if (ddp1[j]==1) {
      /// printf("ba2 j=%d facet=%d i=%d e0=%f e1=%f e2=%f\n",j,facet,i,e0,e1,e2);
      e1=*(e+nextcell_point_conv1(nx,ny,x,y,ddp1[j],conv_type));
      e2=*(e+nextcell_point_conv1(nx,ny,x,y,ddp2[j],conv_type));
      if (!(std::isnan(e1) | std::isnan(e2))) {
        if (e1<e_temp) {
          e_temp=e1;
          facet=j;
        }
        if (e2<e_temp) {
          e_temp=e2;
          facet=j;
        } // added on 2024 09 16
      }
      
    }
    
    
    
    
    
    e1=*(e+nextcell_point_conv1(nx,ny,x,y,ddp1[facet],conv_type));
    e2=*(e+nextcell_point_conv1(nx,ny,x,y,ddp2[facet],conv_type)); 
    //if (std::isnan(e1) | std::isnan(e0) | e1==e0) e1=1.05*e0;
    //if (std::isnan(e2) | std::isnan(e0) | e2==e0) e2=1.05*e0;
    *(sfacet+i)=facet;
    *(sr+i)=(std::atan((e1-e2)/(e0-e1))); // RAD 20240912
    *(sm+i)=pow(pow(e0-e1,2)+pow(e1-e2,2),0.5)/L;
    if (*(sr+i)<=0) {
      
      *(sr+i)=0;
      *(sm+i)=(e0-e1)/L;
      
    } else if (*(sr+i)>=M_PI/4) {
      
        *(sr+i)=M_PI/4;
        *(sm+i)=(e0-e2)/(L*sqrt(2));
    }
    
    
    *(tdc+i)=L*std::sin(*(sr+i));
    *(tdd+i)=sqrt(2)*L*std::sin(M_PI/4-*(sr+i)); // see qq 5 on Li at al,2022
  //    printf("ba33  facet=%d i=%d e0=%f e1=%f e2=%f slope_mgn=%f slope_mgn_temp=%f d1=%f,d2=%f tdc=%f tdd=%f,\n",facet,i,e0,e1,e2,slope_mgn,slope_mgn_temp,ddp1[facet],ddp2[facet],*(tdc+i),*(tdd+i));
   // }
    
    
    
    
      //// to do 
}
      
      
      
      
      
}
 
 
void transverse_deviation(double *e, double *tdc, double *tdd,double *sr,double *sm, int *sfacet,int nx, int ny, double L,
                          
                          double *atdc, double *atdd, double *atdplus,double *pflow,int *kupdate,double lambda,
                          std::vector<double> ddp1,std::vector<double> ddp2,std::vector<double> sigma,int nncell,int conv_type)
                          
                          
                          
                           {   
   
  // TO DO
  
 // std::vector<double> ddp1 = {2,2,4,4,6,6,8,8}; // this routine uses Orlandini-Li et, 2022's  convention
 // std::vector<double> ddp2 = {1,3,3,5,5,7,7,1};
//  std::vector<double> sigma = {1,-1,1,-1,1,-1,1,-1};
  int x,y;
//  int nncell=8;

    //     
    //
    //
  int k=1;
  int niter=nx*ny;
  int exit_cond=0;
  int facet,nextc,nextd,nextp;
  double atdplus_temp;
  double e0,e1,e2;
  double pflow_estimate=ddp1[0];
  
  
  for (int i = 0; i < nx*ny; i++) {  
   facet=*(sfacet+i);
 //  *(kcoeff+i)=1;
   *(atdc+i)=*(tdc+i)*sigma[facet];
   *(atdplus+i)=0;
   *(atdd+i)=*(tdd+i)*sigma[facet]*(-1); // 20240912
   *(pflow+i)=ddp1[0];
   *(kupdate+i)=0;
  }
  int cnt=0;
 // for (int iter=0;iter<niter,iter++) {
 for (int j = 0; j < nx*ny; j++) {
  int i=j; 
  cnt++;
  if (*(kupdate+j)==0) do {
    exit_cond=0;
    x = getCol(nx, ny, i);   // ATTENTION: base 0 or 1?
    y = getRow(nx, ny, i);
    nextp=i;
    pflow_estimate=*(pflow+i);
    e0=*(e+i);
    facet=*(sfacet+i);
    nextc=nextcell_point_conv1(nx,ny,x,y,ddp1[facet],conv_type);
    nextd=nextcell_point_conv1(nx,ny,x,y,ddp2[facet],conv_type);
    e1=*(e+nextc);
    e2=*(e+nextd); 
    //if (std::isnan(e1) | std::isnan(e0)  | e1==e0) e1=1.05*e0;
    //if (std::isnan(e2) | std::isnan(e0)  | e2==e0) e2=1.05*e0;
    if (e2<e1) {
      pflow_estimate=ddp2[facet];
     // nextp=nextd;
      //*(atdplus+nextp)
      //atdplus_temp=*(atdd+nextp);
      
    } else  {
      pflow_estimate=ddp1[facet];
      
    } 
    nextp=nextcell_point_conv1(nx,ny,x,y,pflow_estimate,conv_type);
    int facetc=*(sfacet+nextc);
    int facetd=*(sfacet+nextd);
    *(atdc+nextc)=*(tdc+nextc)*sigma[facetc]+*(atdplus+i)*lambda;
    *(atdd+nextd)=*(tdd+nextd)*sigma[facetd]*(-1)+*(atdplus+i)*lambda; // 20240924
      
 

  
   
    if ((e0<e1) & (e0>e2)) { //// pit ??? 
        //*(atdplus+nextd)=*(atdplus+nextd)+*(atdd+i); //DA VEDERE BENE
      pflow_estimate=ddp2[facet];
      nextp=nextd;
        //*(atdplus+nextp)
      atdplus_temp=*(atdd+nextp);
        
        //*(pflow2+i)=ddp1[facet]; //tocontinue...
    } else if ((e0<e2) & (e0>e1)) {
      pflow_estimate=ddp1[facet];
      nextp=nextc;
      atdplus_temp=*(atdc+nextp);
    } else if (abs(*(atdc+i))<=abs(*(atdd+i))) {
        
      pflow_estimate=ddp1[facet];
      nextp=nextc;
      atdplus_temp=*(atdc+nextp);
   // } else if (abs(*(atdd+i))<abs(*(atdc+i))) {
    } else  {
      pflow_estimate=ddp2[facet];
      nextp=nextd;
      atdplus_temp=*(atdd+nextp);
        
    }
   //21 printf("zzx i=%d pflow=%f %f\n",i,pflow_estimate,*(pflow+i));
  //  if (pflow_estimate==*(pflow+i)) {
    //  exit_cond=1;
  //21    printf("zzz pflow=%f pfow=%f\n",pflow_estimate,*(pflow+i));
  //21    printf("zzz i=%d nextp=%d su j=%d facet=%d \n",i,nextp,j,facet);
      
 //   }  else {
      
 //     exit_cond=0;
 //   }///else if (nextp!=i) {
    
    
    
    if (abs(atdplus_temp)>abs(*(atdplus+nextp))){
        
      *(atdplus+nextp)=atdplus_temp;
        
    }
    ///exit_cond=0;
    *(pflow+i)=pflow_estimate;
    *(kupdate+i)=cnt;
    printf("x=%d y=%d atdd=%f atdc=%f tdd=%f tdc=%f sm=%f sr=%f pflow=%f kup=%d i=%d j=%d nextp=%d\n facet=%d e0=%f e1=%f e2=%f \n",x,y,*(atdd+i),*(atdc+i),*(tdd+i),*(tdc+i),*(sm+i),*(sr+i),*(pflow+i),*(kupdate+i),i,j,nextp,facet,e0,e1,e2);
    if (nextp!=i) {
      
      i=nextp;
      
      exit_cond=0;
    } else {
      
      exit_cond=2;
    }
      
      
 // }
    
    //}  else {
    //  printf("L2 i=%d nextp=%d su j=%d facet=%f pflow=%f\n",i,nextp,j,facet,*(pflow+i));
  //    exit_cond=1;
  //  }
  ///  printf("j=%d su i=%d exit_cond=%d lambda=%f facet=%d\n",j,i,exit_cond,lambda,facet);
    
  } while (exit_cond==0); 
  }
  // NOVALUE
  for (int i=0;i<nx*ny;i++) {
  //  
  
    double e0=*(e+i);
  
    if (std::isnan(e0)) {
      *(pflow+i)=e0;
   }
  

 
  
  
  
}
 


    
    
    
    

  
  
}
  
  
SpatRaster  SpatRaster::d8ltd(double lambda,SpatOptions &opt) {
    // DA TESTARE
    SpatRaster out=geometry();
    //std::vector<std::string> oname="watershed";
    //out.setNames(oname);
    int nx=ncol();
    int ny=nrow();
    double Lx=xres();
    double Ly=yres();
    double L=Lx; // to check 
    //printf("nx=%d ny=%d\n",nx,ny);
    //Rcpp::IntegerVector pOut(nx*ny);
    // https://www.codeguru.com/cpp/cpp/cpp_mfc/stl/article.php/c4027/C-Tutorial-A-Beginners-Guide-to-stdvector-Part-1.htm 
    std::vector<double> e=getValues(0,opt); //EC 20211203 //see https://www.delftstack.com/howto/cpp/how-to-convert-vector-to-array-in-cpp/
    ///std::vector<double> sr=weight.getValues(0,opt);
    
    std::vector<double> pOutv(nx*ny,0);
    std::vector<int> sfacet(nx*ny,0);
    std::vector<double> tdc(nx*ny,0);
    std::vector<double> tdd(nx*ny,0);
    std::vector<double> atdc(nx*ny,0);
    std::vector<double> atdd(nx*ny,0);
    std::vector<double> atdplus(nx*ny,0);
   
    std::vector<double> sr(nx*ny,0);
    std::vector<double> sm(nx*ny,0);
    
    std::vector<double> pflow(nx*ny,0);
    std::vector<int> kupdate(nx*ny,0);
    
    /* FLOW Direction Convention */ 
    
    //std::vector<double> ddp1 = {0,2,2,4,4,6,6,8,8}; // this routine uses Orlandini-Li et, 2022's  convention
    //std::vector<double> ddp2 = {0,1,3,3,5,5,7,7,1};
    //int conv_type=2;
    // ESRI terra covension:
    //   //     # x-1 x 	x+1
    //   // y-1 ## 32	64	128
    //   // y   ## 16	0	1
    //   // y+1 ## 8	4	2
    //   
    
    
    
    int conv_type=1;
    std::vector<double> ddp1 = {0,64,64,1,1,4,4,16,16};// this routine ESRI  convention
    std::vector<double> ddp2 = {0,32,128,128,2,2,8,8,32};
    
    
    //  Li at al, 2022's convection 
    //   //     # x-1 x 	x+1
    //   // y-1 ## 1	2	3
    //   // y   ## 8	0	4
    //   // y+1 ## 7	6	5
    
   
    std::vector<double> sigma = {0,1,-1,1,-1,1,-1,1,-1};
    int nncell=9;
    
    
    
    slope_direction(&e[0],nx,ny,&sr[0],&sm[0],&sfacet[0],
                    &tdc[0],&tdd[0],L,ddp1,ddp2,nncell,conv_type);
    
    
    
    
    
    transverse_deviation(&e[0],&tdc[0],&tdd[0],&sr[0],&sm[0],&sfacet[0],nx,ny,L,
                         
                         &atdc[0], &atdd[0], &atdplus[0], &pflow[0],&kupdate[0],lambda,ddp1,ddp2,sigma,nncell,conv_type);
    
    
    // NOVALUE
    //for (int i=0;i<nx*ny;i++) {
    //  
      
    //  double e0=*(e+i);
      
    //  if (e0==0) {
    //    *(pflow+i)=*(e+i);
     // }
      
 //   }
    
    
    
    if (!out.writeStart(opt,filenames())) {
      readStop();
      return out;
    }
    
    // convert pflow ....
    // out.writeValues(pOutv,0,ny,0,nx); UNTIL 20220725
    out.writeValues(pflow,0,ny); //,0,nx); // LOOK AT writeValuesGDAL
    out.writeStop();
    
    
    // DEALLOC 
    
    
    
    return out;
    
    
  }  
  

