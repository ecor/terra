// C/C++ code
// Author: Emanuele Cordano
// Date: February 2021
// Flow dirrection D8-LTD algorithm implementation (under development)
//

// TO BE IMPLEMENTED

//#define NODATA0 -9999



int nextcell_point_conv1(int nx, int ny,int x,int y,int pv,int conv_type);
//void slope_direction(double* e, int nx, int ny, double *sr,double *sm,int *sfacet,double L);


void slope_direction(double* e, int nx, int ny, double *sr,double *sm,int *sfacet,
                     double* tdc, double *tdd,double L,
                     std::vector<double> ddp1,std::vector<double> ddp2,int nncell,int conv_type);


void transverse_deviation(double *e, double *tdc, double *tdd, int *sfacet,int nx, int ny, double L,
                          
                          double *atdc, double *atdd, double *atdplus,double *pflow,int *kupdate,double lambda,
                          std::vector<double> ddp1,std::vector<double> ddp2,std::vector<double> sigma,int nncell,int conv_type);