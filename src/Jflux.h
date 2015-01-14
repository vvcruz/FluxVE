
//Entire J vector calculation routines
double prob_current1D(double *J, double *Re, double *Im, int nx, double stepx, int dyn, double m);
double prob_current2D(double **Jx, double **Jy, double **Re, double **Im, int nx, int ny, double stepx, double stepy, double m1, double m2, int order);
//Flux along a line without splines.
double integflux_x(double **Re,double **Im,int fp0,int fp1,int fp2,double stepx, double stepy,double mx,double yi,int order);
double integflux_y(double **Re,double **Im,int fp0,int fp1,int fp2,double stepx, double stepy,double my,double xi,int order);
//Flux along a lines applying splines.
double integflux_x_spl(double x0, double yi,int fp1, int fp2,double *XKNOT, double *YKNOT,int nx,int ny,double *BCOEFRe, double *BCOEFIm, double stepx, double stepy,double mx,int order);
double integflux_y_spl(double y0, double xi,int fp1, int fp2,double *XKNOT, double *YKNOT,int nx,int ny,double *BCOEFRe, double *BCOEFIm, double stepx, double stepy,double my,int order);
//Flux across a line applying splines and Jacobi coordinates transformation
double integflux_x_jac(double R0, double ri,int fp1, int fp2,double *XKNOT, double *YKNOT,int nx,int ny,double *BCOEFRe, double *BCOEFIm, double stepR, double stepr,double mx,int order,double ma,double mb,double mc);
//Numerical derivative routines
double deriv(double y1, double y2, double step);
double derivO4(double y1, double y2, double y3, double y4, double step);

// Energy resolved flux routine
double spl_integflux_x(double *bcoefre,double *bcoefim, double *xknot, double *yknot, double *eknot, int kx,double x0, double *Y, double E,int nxr,int nyr,int NYG,int nE, double stepx,double stepy, int *fp,double mx);

// Energy resolved flux routine + jacobi coordinates
double spl_jac_integflux_x(double *bcoefre,double *bcoefim, double *xknot, double *yknot, double *eknot, int kx,double R0, double ri, double E,int nxr,int nyr, int NYG, int nE,double stepR,double stepr,double ma,double mb, double mc, double mx);

//Energy resolved flux routine 2
double spl_integflux_x2(double *WPERe,double *WPEIm,double *XX,double x0, double *Y, double E, int kE,int nxr,int nyr, int NYG, int nE,double stepx,double stepy, int *fp,double mx);

//numerical integration routines
double comp_simpson38(double *fx,double stepx,int n);
double simpson38(double *y,double stepx);
double simpson13(double *y,double stepx);
