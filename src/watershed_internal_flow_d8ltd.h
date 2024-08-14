// C/C++ code
// Author: Emanuele Cordano
// Date: February 2021
// Flow dirrection D8-LTD algorithm implementation (under development)
//

// TO BE IMPLEMENTED

//#define NODATA0 -9999



int nextcell_point_conv1(int nx, int ny,int x,int y,int pv,int conv_type);
void slope_direction(double* e, int nx, int ny, double *sr,double *sm,int *sfacet,double L);

////void slope_magnitude(double* v, int nx, int ny, double* vOut,double dx,double dy);


// 
// 
// void NextCell(double* p, int nx, int ny,int* pnext) {
//   
//   //int* q;           // A pointer to a queue of raster cells (offset in memory) to be processed
//   //  int qSize = 50;   // Starting queue size, that can be dinamically incremented if needed
//   //int delta;        // Offset in memory from base queue address of a raster cell
//   //int n = 0;        // Number of raster cells to be processed in queue
//   //  int nLoop = 0;    // Counter for loops over cells
//   int i;
//   
//   // ## 32	64	128
//   // ## 16	x	1
//   // ## 8	4	2
//   
//   for (i=0;i<nx*ny;i++) {
//     *(pnext+i)=i; //-9999;
//   }
//   for (int x=0;x<nx;x++) {
//     for (int y=0;y<ny;y++) {
//       
//       i = offset(nx,ny,x,y);
//       if (inRaster(nx, ny, x + 1, y) && *(p+i)==1) {
//         *(pnext+i)=offset(nx, ny, x + 1, y);
//       } else if (inRaster(nx, ny, x+1,y+1) && *(p+i)==2) {
//         
//         *(pnext+i)=offset(nx, ny, x+1, y+1);
//       } else if (inRaster(nx, ny, x,y+1) && *(p+i)==4) {
//         *(pnext+i)=offset(nx, ny, x, y+1);
//       } else if (inRaster(nx, ny, x-1,y+1) && *(p+i)==8) {
//         *(pnext+i)=offset(nx, ny, x-1, y+1);
//       } else if (inRaster(nx, ny, x-1,y) && *(p+i)==16) {
//         *(pnext+i)=offset(nx, ny, x-1, y);
//       } else if (inRaster(nx, ny, x-1,y-1) && *(p+i)==32) {
//         *(pnext+i)=offset(nx, ny, x-1, y-1);
//       } else if (inRaster(nx, ny, x,y-1) && *(p+i)==64) {
//         *(pnext+i)=offset(nx, ny, x, y-1);
//       } else if (inRaster(nx, ny, x+1,y-1) && *(p+i)==128) {
//         *(pnext+i)=offset(nx, ny, x+1, y-1);
//       } 
//     }
//     
//   }  
//   
// }  
// 
