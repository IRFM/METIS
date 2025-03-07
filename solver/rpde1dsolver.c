static char mc_version[] = "MATLAB Compiler 1.2.1 Jan 21 1999 infun";
/*
 *  MATLAB Compiler: 1.2.1
 *  Date: Jan 21 1999
 *  Arguments: -Z -r rpde1dsolver 
 */
#ifndef ARRAY_ACCESS_INLINING
#error You must use the -inline option when compiling MATLAB compiler generated code with MEX or MBUILD
#endif
#ifndef MATLAB_COMPILER_GENERATED_CODE
#define MATLAB_COMPILER_GENERATED_CODE
#endif

#include <math.h>
#include "mex.h"
#include "mcc.h"

/* static array S0_ (1 x 3) int, line 231 */
static int S0r_[] =
{
	3, 1, 2, 
};
static mxArray S0_ = mccCINIT( mccINT, 1, 3, S0r_, 0 );
/* static array S1_ (1 x 3) int, line 232 */
static int S1r_[] =
{
	3, 1, 2, 
};
static mxArray S1_ = mccCINIT( mccINT, 1, 3, S1r_, 0 );
/* static array S2_ (1 x 3) int, line 234 */
static int S2r_[] =
{
	3, 1, 
	2, 
};
static mxArray S2_ = mccCINIT( mccINT, 1, 3, S2r_, 0 );
/* static array S3_ (1 x 3) int, line 235 */
static int S3r_[] =
{
	3, 1, 2, 
};
static mxArray S3_ = mccCINIT( mccINT, 1, 3, S3r_, 0 );
/* static array S4_ (1 x 3) int, line 237 */
static int S4r_[] =
{
	3, 1, 2, 
};
static mxArray S4_ = mccCINIT( mccINT, 1, 3, S4r_, 0 );
/* static array S5_ (1 x 3) int, line 238 */
static int S5r_[] =
{
	3, 
	1, 2, 
};
static mxArray S5_ = mccCINIT( mccINT, 1, 3, S5r_, 0 );
/* static array S6_ (1 x 3) int, line 248 */
static int S6r_[] =
{
	3, 1, 2, 
};
static mxArray S6_ = mccCINIT( mccINT, 1, 3, S6r_, 0 );
/* static array S7_ (1 x 3) int, line 249 */
static int S7r_[] =
{
	3, 1, 2, 
	
};
static mxArray S7_ = mccCINIT( mccINT, 1, 3, S7r_, 0 );
/* static array S8_ (1 x 3) int, line 251 */
static int S8r_[] =
{
	3, 1, 2, 
};
static mxArray S8_ = mccCINIT( mccINT, 1, 3, S8r_, 0 );
/* static array S9_ (1 x 3) int, line 252 */
static int S9r_[] =
{
	3, 1, 2, 
};
static mxArray S9_ = mccCINIT( mccINT, 1, 3, S9r_, 0 );
/* static array S10_ (1 x 3) int, line 254 */
static int S10r_[] =
{
	3, 1, 
	2, 
};
static mxArray S10_ = mccCINIT( mccINT, 1, 3, S10r_, 0 );
/* static array S11_ (1 x 3) int, line 255 */
static int S11r_[] =
{
	3, 1, 2, 
};
static mxArray S11_ = mccCINIT( mccINT, 1, 3, S11r_, 0 );
/* static array S12_ (1 x 3) int, line 257 */
static int S12r_[] =
{
	3, 1, 2, 
};
static mxArray S12_ = mccCINIT( mccINT, 1, 3, S12r_, 0 );
/* static array S13_ (1 x 3) int, line 258 */
static int S13r_[] =
{
	3, 
	1, 2, 
};
static mxArray S13_ = mccCINIT( mccINT, 1, 3, S13r_, 0 );
/* static array S14_ (1 x 3) int, line 273 */
static int S14r_[] =
{
	3, 1, 2, 
};
static mxArray S14_ = mccCINIT( mccINT, 1, 3, S14r_, 0 );
/* static array S15_ (1 x 3) int, line 274 */
static int S15r_[] =
{
	3, 1, 2, 
	
};
static mxArray S15_ = mccCINIT( mccINT, 1, 3, S15r_, 0 );
/* static array S16_ (1 x 3) int, line 276 */
static int S16r_[] =
{
	3, 1, 2, 
};
static mxArray S16_ = mccCINIT( mccINT, 1, 3, S16r_, 0 );
/* static array S17_ (1 x 3) int, line 277 */
static int S17r_[] =
{
	3, 1, 2, 
};
static mxArray S17_ = mccCINIT( mccINT, 1, 3, S17r_, 0 );
/* static array S18_ (1 x 3) int, line 279 */
static int S18r_[] =
{
	3, 1, 
	2, 
};
static mxArray S18_ = mccCINIT( mccINT, 1, 3, S18r_, 0 );
/* static array S19_ (1 x 3) int, line 280 */
static int S19r_[] =
{
	3, 1, 2, 
};
static mxArray S19_ = mccCINIT( mccINT, 1, 3, S19r_, 0 );
/* static array S20_ (1 x 3) int, line 314 */
static int S20r_[] =
{
	3, 1, 2, 
};
static mxArray S20_ = mccCINIT( mccINT, 1, 3, S20r_, 0 );
/* static array S21_ (1 x 3) int, line 315 */
static int S21r_[] =
{
	3, 
	1, 2, 
};
static mxArray S21_ = mccCINIT( mccINT, 1, 3, S21r_, 0 );
/* static array S22_ (1 x 3) int, line 317 */
static int S22r_[] =
{
	3, 1, 2, 
};
static mxArray S22_ = mccCINIT( mccINT, 1, 3, S22r_, 0 );
/* static array S23_ (1 x 3) int, line 318 */
static int S23r_[] =
{
	3, 1, 2, 
	
};
static mxArray S23_ = mccCINIT( mccINT, 1, 3, S23r_, 0 );
/* static array S24_ (1 x 3) int, line 320 */
static int S24r_[] =
{
	3, 1, 2, 
};
static mxArray S24_ = mccCINIT( mccINT, 1, 3, S24r_, 0 );
/* static array S25_ (1 x 3) int, line 321 */
static int S25r_[] =
{
	3, 1, 2, 
};
static mxArray S25_ = mccCINIT( mccINT, 1, 3, S25r_, 0 );
/* static array S26_ (1 x 3) int, line 330 */
static int S26r_[] =
{
	3, 1, 
	2, 
};
static mxArray S26_ = mccCINIT( mccINT, 1, 3, S26r_, 0 );
/* static array S27_ (1 x 3) int, line 331 */
static int S27r_[] =
{
	3, 1, 2, 
};
static mxArray S27_ = mccCINIT( mccINT, 1, 3, S27r_, 0 );
/* static array S28_ (1 x 3) int, line 334 */
static int S28r_[] =
{
	3, 1, 2, 
};
static mxArray S28_ = mccCINIT( mccINT, 1, 3, S28r_, 0 );
/* static array S29_ (1 x 3) int, line 335 */
static int S29r_[] =
{
	3, 
	1, 2, 
};
static mxArray S29_ = mccCINIT( mccINT, 1, 3, S29r_, 0 );
/* static array S30_ (1 x 3) int, line 337 */
static int S30r_[] =
{
	3, 1, 2, 
};
static mxArray S30_ = mccCINIT( mccINT, 1, 3, S30r_, 0 );
/* static array S31_ (1 x 3) int, line 338 */
static int S31r_[] =
{
	3, 1, 2, 
	
};
static mxArray S31_ = mccCINIT( mccINT, 1, 3, S31r_, 0 );
/* static array S32_ (1 x 3) int, line 340 */
static int S32r_[] =
{
	3, 1, 2, 
};
static mxArray S32_ = mccCINIT( mccINT, 1, 3, S32r_, 0 );
/* static array S33_ (1 x 3) int, line 341 */
static int S33r_[] =
{
	3, 1, 2, 
};
static mxArray S33_ = mccCINIT( mccINT, 1, 3, S33r_, 0 );
/* static array S34_ (1 x 3) int, line 355 */
static int S34r_[] =
{
	3, 1, 
	2, 
};
static mxArray S34_ = mccCINIT( mccINT, 1, 3, S34r_, 0 );
/* static array S35_ (1 x 3) int, line 356 */
static int S35r_[] =
{
	3, 1, 2, 
};
static mxArray S35_ = mccCINIT( mccINT, 1, 3, S35r_, 0 );
/* static array S36_ (1 x 3) int, line 358 */
static int S36r_[] =
{
	3, 1, 2, 
};
static mxArray S36_ = mccCINIT( mccINT, 1, 3, S36r_, 0 );
/* static array S37_ (1 x 3) int, line 359 */
static int S37r_[] =
{
	3, 
	1, 2, 
};
static mxArray S37_ = mccCINIT( mccINT, 1, 3, S37r_, 0 );
/* static array S38_ (1 x 3) int, line 361 */
static int S38r_[] =
{
	3, 1, 2, 
};
static mxArray S38_ = mccCINIT( mccINT, 1, 3, S38r_, 0 );
/* static array S39_ (1 x 3) int, line 362 */
static int S39r_[] =
{
	3, 1, 2, 
	
};
static mxArray S39_ = mccCINIT( mccINT, 1, 3, S39r_, 0 );
/* static array S40_ (1 x 47) text, line 183: '=>>>> NaN ou Inf dans les...' */
static unsigned short S40__r_[] =
{
         61,   62,   62,   62,   62,   32,   78,   97,
         78,   32,  111,  117,   32,   73,  110,  102,
         32,  100,   97,  110,  115,   32,  108,  101,
        115,   32,  100,  111,  110,  110,  101,  101,
        115,   32,  100,  101,  115,   32,   80,   68,
         69,   32,   33,   33,   33,   92,  110,
};
static mxArray S40_ = mccCINIT( mccTEXT,  1, 47, S40__r_, 0);
/* static array S41_ (1 x 47) text, line 188: '=>>>> NaN ou Inf dans les...' */
static unsigned short S41__r_[] =
{
         61,   62,   62,   62,   62,   32,   78,   97,
         78,   32,  111,  117,   32,   73,  110,  102,
         32,  100,   97,  110,  115,   32,  108,  101,
        115,   32,  100,  111,  110,  110,  101,  101,
        115,   32,  100,  101,  115,   32,   80,   68,
         69,   32,   33,   33,   33,   92,  110,
};
static mxArray S41_ = mccCINIT( mccTEXT,  1, 47, S41__r_, 0);
/* static array S42_ (1 x 3) text, line 216: 'off' */
static unsigned short S42__r_[] =
{
        111,  102,  102,
};
static mxArray S42_ = mccCINIT( mccTEXT,  1, 3, S42__r_, 0);
/* static array S43_ (1 x 2) text, line 223: 'on' */
static unsigned short S43__r_[] =
{
        111,  110,
};
static mxArray S43_ = mccCINIT( mccTEXT,  1, 2, S43__r_, 0);
/* static array S44_ (1 x 3) text, line 299: 'off' */
static unsigned short S44__r_[] =
{
        111,  102,  102,
};
static mxArray S44_ = mccCINIT( mccTEXT,  1, 3, S44__r_, 0);
/* static array S45_ (1 x 2) text, line 306: 'on' */
static unsigned short S45__r_[] =
{
        111,  110,
};
static mxArray S45_ = mccCINIT( mccTEXT,  1, 2, S45__r_, 0);
/***************** Compiler Assumptions ****************
 * M-File: /tmp_mnt/usr/drfc/cgc/matlab5/zineb/v1.5/solver/rpde1dsolver.m
 *
 *       A           	real vector/matrix
 *       ALPHA       	real vector/matrix
 *       ALPHAP      	real vector/matrix
 *       AP          	real vector/matrix
 *       B           	real vector/matrix
 *       B0_         	boolean scalar temporary
 *       B1_         	boolean scalar temporary
 *       B2_         	boolean scalar temporary
 *       BM0_        	boolean vector/matrix temporary
 *       BM1_        	boolean vector/matrix temporary
 *       BP          	real vector/matrix
 *       C           	real vector/matrix
 *       CP          	real vector/matrix
 *       D           	real vector/matrix
 *       DP          	real vector/matrix
 *       F           	real vector/matrix
 *       FP          	real vector/matrix
 *       FPa         	real vector/matrix
 *       FPb         	real vector/matrix
 *       FPc         	real vector/matrix
 *       FS          	real vector/matrix
 *       Fa          	real vector/matrix
 *       Fb          	real vector/matrix
 *       Fc          	real vector/matrix
 *       I0_         	integer scalar temporary
 *       I1_         	integer scalar temporary
 *       I2_         	integer scalar temporary
 *       I3_         	integer scalar temporary
 *       IDENTITE    	real vector/matrix
 *       IM0_        	integer vector/matrix temporary
 *       IM10_       	integer vector/matrix temporary
 *       IM11_       	integer vector/matrix temporary
 *       IM1_        	integer vector/matrix temporary
 *       IM2_        	integer vector/matrix temporary
 *       IM3_        	integer vector/matrix temporary
 *       IM4_        	integer vector/matrix temporary
 *       IM5_        	integer vector/matrix temporary
 *       IM6_        	integer vector/matrix temporary
 *       IM7_        	integer vector/matrix temporary
 *       IM8_        	integer vector/matrix temporary
 *       IM9_        	integer vector/matrix temporary
 *       K           	integer scalar
 *       M           	integer scalar
 *       R0_         	real scalar temporary
 *       RM0_        	real vector/matrix temporary
 *       RM10_       	real vector/matrix temporary
 *       RM11_       	real vector/matrix temporary
 *       RM12_       	real vector/matrix temporary
 *       RM13_       	real vector/matrix temporary
 *       RM14_       	real vector/matrix temporary
 *       RM1_        	real vector/matrix temporary
 *       RM2_        	real vector/matrix temporary
 *       RM3_        	real vector/matrix temporary
 *       RM4_        	real vector/matrix temporary
 *       RM5_        	real vector/matrix temporary
 *       RM6_        	real vector/matrix temporary
 *       RM7_        	real vector/matrix temporary
 *       RM8_        	real vector/matrix temporary
 *       RM9_        	real vector/matrix temporary
 *       S           	real vector/matrix
 *       S0_         	<constant>
 *       S10_        	<constant>
 *       S11_        	<constant>
 *       S12_        	<constant>
 *       S13_        	<constant>
 *       S14_        	<constant>
 *       S15_        	<constant>
 *       S16_        	<constant>
 *       S17_        	<constant>
 *       S18_        	<constant>
 *       S19_        	<constant>
 *       S1_         	<constant>
 *       S20_        	<constant>
 *       S21_        	<constant>
 *       S22_        	<constant>
 *       S23_        	<constant>
 *       S24_        	<constant>
 *       S25_        	<constant>
 *       S26_        	<constant>
 *       S27_        	<constant>
 *       S28_        	<constant>
 *       S29_        	<constant>
 *       S2_         	<constant>
 *       S30_        	<constant>
 *       S31_        	<constant>
 *       S32_        	<constant>
 *       S33_        	<constant>
 *       S34_        	<constant>
 *       S35_        	<constant>
 *       S36_        	<constant>
 *       S37_        	<constant>
 *       S38_        	<constant>
 *       S39_        	<constant>
 *       S3_         	<constant>
 *       S4_         	<constant>
 *       S5_         	<constant>
 *       S6_         	<constant>
 *       S7_         	<constant>
 *       S8_         	<constant>
 *       S9_         	<constant>
 *       T0          	real vector/matrix
 *       T1          	real vector/matrix
 *       U           	real vector/matrix
 *       V           	real vector/matrix
 *       V0          	real vector/matrix
 *       V1          	real vector/matrix
 *       W           	real vector/matrix
 *       X           	real vector/matrix
 *       Y           	real vector/matrix
 *       Z           	real vector/matrix
 *       a           	real vector/matrix
 *       alpha       	real vector/matrix
 *       alphap      	real vector/matrix
 *       ap          	real vector/matrix
 *       b           	real vector/matrix
 *       bp          	real vector/matrix
 *       c           	real vector/matrix
 *       cat         	<function>
 *       comp        	integer vector/matrix
 *       cp          	real vector/matrix
 *       d           	real vector/matrix
 *       d2fpdx2     	real vector/matrix
 *       dfpdx       	real vector/matrix
 *       dp          	real vector/matrix
 *       dt          	real vector/matrix
 *       dt2x        	real vector/matrix
 *       dtx2        	real vector/matrix
 *       dx          	real vector/matrix
 *       f           	real vector/matrix
 *       find        	<function>
 *       fp          	real vector/matrix
 *       fp0         	real vector/matrix
 *       fpkp1       	real vector/matrix
 *       fpp         	real vector/matrix
 *       fpv         	real vector/matrix
 *       i0          	integer vector/matrix
 *       i1          	integer vector/matrix
 *       i2          	integer vector/matrix
 *       ia          	integer vector/matrix
 *       ib          	integer vector/matrix
 *       ic          	integer vector/matrix
 *       im          	integer vector/matrix
 *       ind         	integer vector/matrix
 *       ip          	integer vector/matrix
 *       isempty     	<function>
 *       isfinite    	<function>
 *       iz          	integer vector/matrix
 *       izp         	integer vector/matrix
 *       km          	real scalar
 *       length      	<function>
 *       m           	integer scalar
 *       mode        	real vector/matrix
 *       n           	integer vector/matrix
 *       n           	integer scalar  => n_1
 *       nargout     	<function>
 *       ones        	<function>
 *       p           	real vector/matrix
 *       permute     	<function>
 *       pp          	real vector/matrix
 *       q           	real vector/matrix
 *       qp          	real vector/matrix
 *       reshape     	<function>
 *       rpde1dsolver	<function being defined>
 *       rpde1dsolver_creux	<function>
 *       s           	real vector/matrix
 *       size        	<function>
 *       squeeze     	<function>
 *       sum         	<function>
 *       u           	real vector/matrix
 *       v           	real vector/matrix
 *       w           	real vector/matrix
 *       warning     	<function>
 *       zeros       	<function>
 *       zverbose    	<function>
 *******************************************************/

void
mexFunction(
    int nlhs_,
    mxArray *plhs_[],
    int nrhs_,
    const mxArray *prhs_[]
)
{
   mxArray *Mplhs_[4];
   mxArray *Mprhs_[11];
   

   if (nrhs_ > 18 )
   {
      mexErrMsgTxt( "Too many input arguments." );
   }

   if (nlhs_ > 8 )
   {
      mexErrMsgTxt( "Too many output arguments." );
   }

   mcmSetLineNumber(0);
   {
      mxArray fp;
      mxArray dfpdx;
      mxArray d2fpdx2;
      mxArray FS;
      mxArray S;
      mxArray ALPHA;
      mxArray ALPHAP;
      mxArray IDENTITE;
      mxArray A;
      mxArray B;
      mxArray C;
      mxArray D;
      mxArray AP;
      mxArray BP;
      mxArray CP;
      mxArray DP;
      mxArray F;
      mxArray FP;
      mxArray V0;
      mxArray T0;
      mxArray V1;
      mxArray T1;
      mxArray mode;
      mxArray f;
      mxArray dx;
      mxArray dt;
      int K = 0;
      int M = 0;
      mxArray n;
      mxArray ind;
      mxArray dtx2;
      mxArray dt2x;
      mxArray a;
      mxArray b;
      mxArray c;
      mxArray d;
      mxArray ap;
      mxArray bp;
      mxArray cp;
      mxArray dp;
      mxArray s;
      mxArray u;
      mxArray v;
      mxArray w;
      mxArray U;
      mxArray V;
      mxArray W;
      mxArray X;
      mxArray Y;
      mxArray Z;
      mxArray M_comp;
      mxArray i1;
      mxArray i2;
      mxArray Fa;
      mxArray FPa;
      mxArray Fb;
      mxArray FPb;
      mxArray Fc;
      mxArray FPc;
      mxArray ia;
      mxArray ib;
      mxArray ic;
      mxArray alpha;
      mxArray alphap;
      mxArray p;
      mxArray q;
      int m = 0;
      int n_1 = 0;
      mxArray pp;
      mxArray qp;
      mxArray iz;
      mxArray izp;
      double km = 0.0;
      mxArray fpv;
      mxArray fp0;
      mxArray fpkp1;
      mxArray fpp;
      mxArray im;
      mxArray i0;
      mxArray ip;
      mxArray BM0_;
      unsigned short B0_ = 0;
      unsigned short B1_ = 0;
      mxArray IM0_;
      int I0_ = 0;
      mxArray IM1_;
      mxArray BM1_;
      unsigned short B2_ = 0;
      mxArray IM2_;
      mxArray IM3_;
      mxArray IM4_;
      mxArray IM5_;
      mxArray RM0_;
      mxArray RM1_;
      mxArray RM2_;
      mxArray RM3_;
      mxArray RM4_;
      mxArray RM5_;
      mxArray IM6_;
      mxArray IM7_;
      mxArray IM8_;
      mxArray IM9_;
      mxArray IM10_;
      mxArray IM11_;
      mxArray RM6_;
      mxArray RM7_;
      mxArray RM8_;
      mxArray RM9_;
      mxArray RM10_;
      mxArray RM11_;
      mxArray RM12_;
      mxArray RM13_;
      mxArray RM14_;
      int I1_ = 0;
      int I2_ = 0;
      int I3_ = 0;
      double R0_ = 0.0;
      
      mccRealInit(A);
      mccImportCopy(&A, ((nrhs_>0) ? prhs_[0] : 0), 0, 0);
      mccRealInit(B);
      mccImportCopy(&B, ((nrhs_>1) ? prhs_[1] : 0), 0, 0);
      mccRealInit(C);
      mccImportCopy(&C, ((nrhs_>2) ? prhs_[2] : 0), 0, 0);
      mccRealInit(D);
      mccImportCopy(&D, ((nrhs_>3) ? prhs_[3] : 0), 0, 0);
      mccRealInit(AP);
      mccImportCopy(&AP, ((nrhs_>4) ? prhs_[4] : 0), 0, 0);
      mccRealInit(BP);
      mccImportCopy(&BP, ((nrhs_>5) ? prhs_[5] : 0), 0, 0);
      mccRealInit(CP);
      mccImportCopy(&CP, ((nrhs_>6) ? prhs_[6] : 0), 0, 0);
      mccRealInit(DP);
      mccImportCopy(&DP, ((nrhs_>7) ? prhs_[7] : 0), 0, 0);
      mccRealInit(F);
      mccImportCopy(&F, ((nrhs_>8) ? prhs_[8] : 0), 0, 0);
      mccRealInit(FP);
      mccImportCopy(&FP, ((nrhs_>9) ? prhs_[9] : 0), 0, 0);
      mccRealInit(V0);
      mccImport(&V0, ((nrhs_>10) ? prhs_[10] : 0), 0, 0);
      mccRealInit(T0);
      mccImport(&T0, ((nrhs_>11) ? prhs_[11] : 0), 0, 0);
      mccRealInit(V1);
      mccImport(&V1, ((nrhs_>12) ? prhs_[12] : 0), 0, 0);
      mccRealInit(T1);
      mccImport(&T1, ((nrhs_>13) ? prhs_[13] : 0), 0, 0);
      mccRealInit(mode);
      mccImport(&mode, ((nrhs_>14) ? prhs_[14] : 0), 0, 0);
      mccRealInit(f);
      mccImport(&f, ((nrhs_>15) ? prhs_[15] : 0), 0, 0);
      mccRealInit(dx);
      mccImport(&dx, ((nrhs_>16) ? prhs_[16] : 0), 0, 0);
      mccRealInit(dt);
      mccImport(&dt, ((nrhs_>17) ? prhs_[17] : 0), 0, 0);
      mccRealInit(fp);
      mccRealInit(dfpdx);
      mccRealInit(d2fpdx2);
      mccRealInit(FS);
      mccRealInit(S);
      mccRealInit(ALPHA);
      mccRealInit(ALPHAP);
      mccRealInit(IDENTITE);
      mccIntInit(n);
      mccIntInit(ind);
      mccRealInit(dtx2);
      mccRealInit(dt2x);
      mccRealInit(a);
      mccRealInit(b);
      mccRealInit(c);
      mccRealInit(d);
      mccRealInit(ap);
      mccRealInit(bp);
      mccRealInit(cp);
      mccRealInit(dp);
      mccRealInit(s);
      mccRealInit(u);
      mccRealInit(v);
      mccRealInit(w);
      mccRealInit(U);
      mccRealInit(V);
      mccRealInit(W);
      mccRealInit(X);
      mccRealInit(Y);
      mccRealInit(Z);
      mccIntInit(M_comp);
      mccIntInit(i1);
      mccIntInit(i2);
      mccRealInit(Fa);
      mccRealInit(FPa);
      mccRealInit(Fb);
      mccRealInit(FPb);
      mccRealInit(Fc);
      mccRealInit(FPc);
      mccIntInit(ia);
      mccIntInit(ib);
      mccIntInit(ic);
      mccRealInit(alpha);
      mccRealInit(alphap);
      mccRealInit(p);
      mccRealInit(q);
      mccRealInit(pp);
      mccRealInit(qp);
      mccIntInit(iz);
      mccIntInit(izp);
      mccRealInit(fpv);
      mccRealInit(fp0);
      mccRealInit(fpkp1);
      mccRealInit(fpp);
      mccIntInit(im);
      mccIntInit(i0);
      mccIntInit(ip);
      mccBoolInit(BM0_);
      mccIntInit(IM0_);
      mccIntInit(IM1_);
      mccBoolInit(BM1_);
      mccIntInit(IM2_);
      mccIntInit(IM3_);
      mccIntInit(IM4_);
      mccIntInit(IM5_);
      mccRealInit(RM0_);
      mccRealInit(RM1_);
      mccRealInit(RM2_);
      mccRealInit(RM3_);
      mccRealInit(RM4_);
      mccRealInit(RM5_);
      mccIntInit(IM6_);
      mccIntInit(IM7_);
      mccIntInit(IM8_);
      mccIntInit(IM9_);
      mccIntInit(IM10_);
      mccIntInit(IM11_);
      mccRealInit(RM6_);
      mccRealInit(RM7_);
      mccRealInit(RM8_);
      mccRealInit(RM9_);
      mccRealInit(RM10_);
      mccRealInit(RM11_);
      mccRealInit(RM12_);
      mccRealInit(RM13_);
      mccRealInit(RM14_);
      
      
      /* % pas de verifications des entrees - gain de temps */
      
      /* % les dimensions */
      /* K=size(F,1); */
      if(mccNOTSET(&F))
      {
         mexErrMsgTxt( "variable F undefined, line 153" );
      }
      K = mccGetDimensionSize(&F, 1);
      /* M=size(F,2); */
      M = mccGetDimensionSize(&F, 2);
      
      
      /* % mise en forme de FP (nul pour le mode predictif) */
      /* % modification du 07/08/2000 -> securite */
      /* n=find(~mode); */
      if(mccNOTSET(&mode))
      {
         mexErrMsgTxt( "variable mode undefined, line 159" );
      }
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         unsigned short *p_BM0_;
         int I_BM0_=1;
         double *p_mode;
         int I_mode=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&mode), mccN(&mode));
         mccAllocateMatrix(&BM0_, m_, n_);
         I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
         p_BM0_ = mccSPR(&BM0_);
         I_mode = (mccM(&mode) != 1 || mccN(&mode) != 1);
         p_mode = mccPR(&mode);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_BM0_+=I_BM0_, p_mode+=I_mode)
               {
                  *p_BM0_ = (!*p_mode);
               }
            }
         }
      }
      mccFind(&n, &BM0_);
      /* if ~isempty(n) */
      B0_ = mccIsEmpty(&n);
      B1_ = (!B0_);
      if ((double)B1_)
      {
         /* FP(:,n) = zeros(K,length(n)); */
         I0_ = mccGetLength(&n);
         mccZerosMN(&IM1_, K, I0_);
         mccFindIndex(&IM0_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_FP;
            int I_FP=1;
            int *p_IM0_;
            int I_IM0_=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM1_), mccN(&IM1_));
            if (!mccNOTSET(&FP)) m_ = mcmCalcResultSize(m_, &n_, mccM(&FP), (mccM(&IM0_) * mccN(&IM0_)));
            mccGrowMatrix(&FP, m_, mccGetMaxIndex(&IM0_ ,mccN(&FP)));
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM0_ += I_IM0_)
               {
                  p_FP = mccPR(&FP) + mccM(&FP) * ((int)(*p_IM0_ - .5)) + 0;
                  for (i_=0; i_<m_; ++i_, p_FP+=I_FP, p_IM1_+=I_IM1_)
                  {
                     *p_FP = ((int)*p_IM1_);
                  }
               }
            }
         }
         /* end */
      }
      /* % fin modification du 07/08/2000 */
      
      /* % modification du 23/08/2000 */
      /* % gain de temps */
      /* n=find(mode); */
      mccFind(&n, &mode);
      /* if ~isempty(n) */
      B1_ = mccIsEmpty(&n);
      B0_ = (!B1_);
      if ((double)B0_)
      {
         /* A(:,:,n)   = zeros(K,M,length(n)); */
         I0_ = mccGetLength(&n);
         Mprhs_[0] = mccTempMatrix(K, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(M, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mplhs_[0] = &IM0_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 169);
         mccFindIndex(&IM1_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_A;
            int I_A=1;
            int *p_IM1_;
            int I_IM1_=1;
            int *p_IM0_;
            int I_IM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&A), (mccM(&IM1_) * mccN(&IM1_)));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM0_), mccN(&IM0_));
            mccGrowMatrix(&A, m_, mccGetMaxIndex(&IM1_ ,mccN(&A)));
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_A = mccPR(&A) + mccM(&A) * ((int)(*p_IM1_ - .5)) + 0;
                  for (i_=0; i_<m_; ++i_, p_A+=I_A, p_IM0_+=I_IM0_)
                  {
                     *p_A = ((int)*p_IM0_);
                  }
               }
            }
         }
         /* B(:,:,n)   = zeros(K,M,length(n)); */
         I0_ = mccGetLength(&n);
         Mprhs_[0] = mccTempMatrix(K, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(M, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mplhs_[0] = &IM1_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 170);
         mccFindIndex(&IM0_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_B;
            int I_B=1;
            int *p_IM0_;
            int I_IM0_=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&B), (mccM(&IM0_) * mccN(&IM0_)));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM1_), mccN(&IM1_));
            mccGrowMatrix(&B, m_, mccGetMaxIndex(&IM0_ ,mccN(&B)));
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM0_ += I_IM0_)
               {
                  p_B = mccPR(&B) + mccM(&B) * ((int)(*p_IM0_ - .5)) + 0;
                  for (i_=0; i_<m_; ++i_, p_B+=I_B, p_IM1_+=I_IM1_)
                  {
                     *p_B = ((int)*p_IM1_);
                  }
               }
            }
         }
         /* C(:,:,n)   = zeros(K,M,length(n)); */
         I0_ = mccGetLength(&n);
         Mprhs_[0] = mccTempMatrix(K, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(M, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mplhs_[0] = &IM0_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 171);
         mccFindIndex(&IM1_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_C;
            int I_C=1;
            int *p_IM1_;
            int I_IM1_=1;
            int *p_IM0_;
            int I_IM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&C), (mccM(&IM1_) * mccN(&IM1_)));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM0_), mccN(&IM0_));
            mccGrowMatrix(&C, m_, mccGetMaxIndex(&IM1_ ,mccN(&C)));
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_C = mccPR(&C) + mccM(&C) * ((int)(*p_IM1_ - .5)) + 0;
                  for (i_=0; i_<m_; ++i_, p_C+=I_C, p_IM0_+=I_IM0_)
                  {
                     *p_C = ((int)*p_IM0_);
                  }
               }
            }
         }
         /* D(:,n)     = zeros(K,length(n)); */
         I0_ = mccGetLength(&n);
         mccZerosMN(&IM1_, K, I0_);
         mccFindIndex(&IM0_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_D;
            int I_D=1;
            int *p_IM0_;
            int I_IM0_=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM1_), mccN(&IM1_));
            if (!mccNOTSET(&D)) m_ = mcmCalcResultSize(m_, &n_, mccM(&D), (mccM(&IM0_) * mccN(&IM0_)));
            mccGrowMatrix(&D, m_, mccGetMaxIndex(&IM0_ ,mccN(&D)));
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM0_ += I_IM0_)
               {
                  p_D = mccPR(&D) + mccM(&D) * ((int)(*p_IM0_ - .5)) + 0;
                  for (i_=0; i_<m_; ++i_, p_D+=I_D, p_IM1_+=I_IM1_)
                  {
                     *p_D = ((int)*p_IM1_);
                  }
               }
            }
         }
         /* AP(:,:,n)  = zeros(K,M,length(n)); */
         I0_ = mccGetLength(&n);
         Mprhs_[0] = mccTempMatrix(K, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(M, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mplhs_[0] = &IM0_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 173);
         mccFindIndex(&IM1_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_AP;
            int I_AP=1;
            int *p_IM1_;
            int I_IM1_=1;
            int *p_IM0_;
            int I_IM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&AP), (mccM(&IM1_) * mccN(&IM1_)));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM0_), mccN(&IM0_));
            mccGrowMatrix(&AP, m_, mccGetMaxIndex(&IM1_ ,mccN(&AP)));
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_AP = mccPR(&AP) + mccM(&AP) * ((int)(*p_IM1_ - .5)) + 0;
                  for (i_=0; i_<m_; ++i_, p_AP+=I_AP, p_IM0_+=I_IM0_)
                  {
                     *p_AP = ((int)*p_IM0_);
                  }
               }
            }
         }
         /* BP(:,:,n)  = zeros(K,M,length(n)); */
         I0_ = mccGetLength(&n);
         Mprhs_[0] = mccTempMatrix(K, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(M, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mplhs_[0] = &IM1_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 174);
         mccFindIndex(&IM0_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_BP;
            int I_BP=1;
            int *p_IM0_;
            int I_IM0_=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&BP), (mccM(&IM0_) * mccN(&IM0_)));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM1_), mccN(&IM1_));
            mccGrowMatrix(&BP, m_, mccGetMaxIndex(&IM0_ ,mccN(&BP)));
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM0_ += I_IM0_)
               {
                  p_BP = mccPR(&BP) + mccM(&BP) * ((int)(*p_IM0_ - .5)) + 0;
                  for (i_=0; i_<m_; ++i_, p_BP+=I_BP, p_IM1_+=I_IM1_)
                  {
                     *p_BP = ((int)*p_IM1_);
                  }
               }
            }
         }
         /* CP(:,:,n)  = zeros(K,M,length(n)); */
         I0_ = mccGetLength(&n);
         Mprhs_[0] = mccTempMatrix(K, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(M, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mplhs_[0] = &IM0_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 175);
         mccFindIndex(&IM1_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_CP;
            int I_CP=1;
            int *p_IM1_;
            int I_IM1_=1;
            int *p_IM0_;
            int I_IM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&CP), (mccM(&IM1_) * mccN(&IM1_)));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM0_), mccN(&IM0_));
            mccGrowMatrix(&CP, m_, mccGetMaxIndex(&IM1_ ,mccN(&CP)));
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_CP = mccPR(&CP) + mccM(&CP) * ((int)(*p_IM1_ - .5)) + 0;
                  for (i_=0; i_<m_; ++i_, p_CP+=I_CP, p_IM0_+=I_IM0_)
                  {
                     *p_CP = ((int)*p_IM0_);
                  }
               }
            }
         }
         /* DP(:,n)    = zeros(K,length(n)); */
         I0_ = mccGetLength(&n);
         mccZerosMN(&IM1_, K, I0_);
         mccFindIndex(&IM0_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_DP;
            int I_DP=1;
            int *p_IM0_;
            int I_IM0_=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM1_), mccN(&IM1_));
            if (!mccNOTSET(&DP)) m_ = mcmCalcResultSize(m_, &n_, mccM(&DP), (mccM(&IM0_) * mccN(&IM0_)));
            mccGrowMatrix(&DP, m_, mccGetMaxIndex(&IM0_ ,mccN(&DP)));
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM0_ += I_IM0_)
               {
                  p_DP = mccPR(&DP) + mccM(&DP) * ((int)(*p_IM0_ - .5)) + 0;
                  for (i_=0; i_<m_; ++i_, p_DP+=I_DP, p_IM1_+=I_IM1_)
                  {
                     *p_DP = ((int)*p_IM1_);
                  }
               }
            }
         }
         /* end */
      }
      
      /* % retrait des NaN et des inf dans les donnees d'entrees -> securite */
      /* ind =find(~isfinite(F)); */
      Mprhs_[0] = &F;
      Mplhs_[0] = &BM0_;
      mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "isfinite", 180);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         unsigned short *p_BM1_;
         int I_BM1_=1;
         unsigned short *p_BM0_;
         int I_BM0_=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&BM0_), mccN(&BM0_));
         mccAllocateMatrix(&BM1_, m_, n_);
         I_BM1_ = (mccM(&BM1_) != 1 || mccN(&BM1_) != 1);
         p_BM1_ = mccSPR(&BM1_);
         I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
         p_BM0_ = mccSPR(&BM0_);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_BM1_+=I_BM1_, p_BM0_+=I_BM0_)
               {
                  *p_BM1_ = (!*p_BM0_);
               }
            }
         }
      }
      mccFind(&ind, &BM1_);
      /* if ~isempty(ind) */
      B0_ = mccIsEmpty(&ind);
      B1_ = (!B0_);
      if ((double)B1_)
      {
         /* F(ind) = zeros(1,length(ind)); */
         I0_ = mccGetLength(&ind);
         mccZerosMN(&IM0_, 1, I0_);
         mccFindIndex(&IM1_, &ind);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_F;
            int *p_IM1_;
            int I_IM1_=1;
            int *p_IM0_;
            int I_IM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM0_), mccN(&IM0_));
            if (m_ == 1 && n_ == 1) m_ = mccCalcSubscriptDimensions(m_, &n_, mccM(&IM1_), mccN(&IM1_), &F);
            mccGrowVector(&F, mccGetMaxIndex(&IM1_ ,mccM(&F)*mccN(&F)));
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_IM1_+=I_IM1_, p_IM0_+=I_IM0_)
                  {
                     mccPR(&F)[((int)(*p_IM1_ - .5))] = ((int)*p_IM0_);
                  }
               }
            }
         }
         /* zverbose('=>>>> NaN ou Inf dans les donnees des PDE !!!\n'); */
         Mprhs_[0] = &S40_;
         Mplhs_[0] = 0;
         mccCallMATLAB(0, Mplhs_, 1, Mprhs_, "zverbose", 183);
         /* end */
      }
      /* ind =find(~isfinite(FP)); */
      if(mccNOTSET(&FP))
      {
         mexErrMsgTxt( "variable FP undefined, line 185" );
      }
      Mprhs_[0] = &FP;
      Mplhs_[0] = &BM1_;
      mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "isfinite", 185);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         unsigned short *p_BM0_;
         int I_BM0_=1;
         unsigned short *p_BM1_;
         int I_BM1_=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&BM1_), mccN(&BM1_));
         mccAllocateMatrix(&BM0_, m_, n_);
         I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
         p_BM0_ = mccSPR(&BM0_);
         I_BM1_ = (mccM(&BM1_) != 1 || mccN(&BM1_) != 1);
         p_BM1_ = mccSPR(&BM1_);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_BM0_+=I_BM0_, p_BM1_+=I_BM1_)
               {
                  *p_BM0_ = (!*p_BM1_);
               }
            }
         }
      }
      mccFind(&ind, &BM0_);
      /* if ~isempty(ind) */
      B1_ = mccIsEmpty(&ind);
      B0_ = (!B1_);
      if ((double)B0_)
      {
         /* FP(ind) = zeros(1,length(ind)); */
         I0_ = mccGetLength(&ind);
         mccZerosMN(&IM1_, 1, I0_);
         mccFindIndex(&IM0_, &ind);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_FP;
            int *p_IM0_;
            int I_IM0_=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM1_), mccN(&IM1_));
            if (m_ == 1 && n_ == 1) m_ = mccCalcSubscriptDimensions(m_, &n_, mccM(&IM0_), mccN(&IM0_), &FP);
            mccGrowVector(&FP, mccGetMaxIndex(&IM0_ ,mccM(&FP)*mccN(&FP)));
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_IM0_+=I_IM0_, p_IM1_+=I_IM1_)
                  {
                     mccPR(&FP)[((int)(*p_IM0_ - .5))] = ((int)*p_IM1_);
                  }
               }
            }
         }
         /* zverbose('=>>>> NaN ou Inf dans les donnees des PDE !!!\n'); */
         Mprhs_[0] = &S41_;
         Mplhs_[0] = 0;
         mccCallMATLAB(0, Mplhs_, 1, Mprhs_, "zverbose", 188);
         /* end */
      }
      /* % fin modification du 23/08/2000  */
      
      
      /* % constante dimentionnelle */
      /* dtx2 = dt ./ dx .^ 2; */
      if(mccNOTSET(&dt))
      {
         mexErrMsgTxt( "variable dt undefined, line 194" );
      }
      if(mccNOTSET(&dx))
      {
         mexErrMsgTxt( "variable dx undefined, line 194" );
      }
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_dtx2;
         int I_dtx2=1;
         double *p_dt;
         int I_dt=1;
         double *p_dx;
         int I_dx=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dt), mccN(&dt));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dx), mccN(&dx));
         mccAllocateMatrix(&dtx2, m_, n_);
         I_dtx2 = (mccM(&dtx2) != 1 || mccN(&dtx2) != 1);
         p_dtx2 = mccPR(&dtx2);
         I_dt = (mccM(&dt) != 1 || mccN(&dt) != 1);
         p_dt = mccPR(&dt);
         I_dx = (mccM(&dx) != 1 || mccN(&dx) != 1);
         p_dx = mccPR(&dx);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_dtx2+=I_dtx2, p_dt+=I_dt, p_dx+=I_dx)
               {
                  *p_dtx2 = (*p_dt / (double) mcmRealPowerInt(*p_dx, 2));
               }
            }
         }
      }
      /* dt2x = dt ./ 2 ./ dx; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_dt2x;
         int I_dt2x=1;
         double *p_dt;
         int I_dt=1;
         double *p_dx;
         int I_dx=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dt), mccN(&dt));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dx), mccN(&dx));
         mccAllocateMatrix(&dt2x, m_, n_);
         I_dt2x = (mccM(&dt2x) != 1 || mccN(&dt2x) != 1);
         p_dt2x = mccPR(&dt2x);
         I_dt = (mccM(&dt) != 1 || mccN(&dt) != 1);
         p_dt = mccPR(&dt);
         I_dx = (mccM(&dx) != 1 || mccN(&dx) != 1);
         p_dx = mccPR(&dx);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_dt2x+=I_dt2x, p_dt+=I_dt, p_dx+=I_dx)
               {
                  *p_dt2x = ((*p_dt / (double) 2) / (double) *p_dx);
               }
            }
         }
      }
      
      /* % 1 - combinaison des coeficients */
      /* a = dtx2 .* A - dt2x .* B; */
      if(mccNOTSET(&A))
      {
         mexErrMsgTxt( "variable A undefined, line 198" );
      }
      if(mccNOTSET(&B))
      {
         mexErrMsgTxt( "variable B undefined, line 198" );
      }
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_a;
         int I_a=1;
         double *p_dtx2;
         int I_dtx2=1;
         double *p_A;
         int I_A=1;
         double *p_dt2x;
         int I_dt2x=1;
         double *p_B;
         int I_B=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dtx2), mccN(&dtx2));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&A), mccN(&A));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dt2x), mccN(&dt2x));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&B), mccN(&B));
         mccAllocateMatrix(&a, m_, n_);
         I_a = (mccM(&a) != 1 || mccN(&a) != 1);
         p_a = mccPR(&a);
         I_dtx2 = (mccM(&dtx2) != 1 || mccN(&dtx2) != 1);
         p_dtx2 = mccPR(&dtx2);
         I_A = (mccM(&A) != 1 || mccN(&A) != 1);
         p_A = mccPR(&A);
         I_dt2x = (mccM(&dt2x) != 1 || mccN(&dt2x) != 1);
         p_dt2x = mccPR(&dt2x);
         I_B = (mccM(&B) != 1 || mccN(&B) != 1);
         p_B = mccPR(&B);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_a+=I_a, p_dtx2+=I_dtx2, p_A+=I_A, p_dt2x+=I_dt2x, p_B+=I_B)
               {
                  *p_a = ((*p_dtx2 * (double) *p_A) - (*p_dt2x * (double) *p_B));
               }
            }
         }
      }
      /* b = - 2 .* dtx2 .* A  + dt .* C; */
      if(mccNOTSET(&C))
      {
         mexErrMsgTxt( "variable C undefined, line 199" );
      }
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_b;
         int I_b=1;
         double *p_dtx2;
         int I_dtx2=1;
         double *p_A;
         int I_A=1;
         double *p_dt;
         int I_dt=1;
         double *p_C;
         int I_C=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dtx2), mccN(&dtx2));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&A), mccN(&A));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dt), mccN(&dt));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&C), mccN(&C));
         mccAllocateMatrix(&b, m_, n_);
         I_b = (mccM(&b) != 1 || mccN(&b) != 1);
         p_b = mccPR(&b);
         I_dtx2 = (mccM(&dtx2) != 1 || mccN(&dtx2) != 1);
         p_dtx2 = mccPR(&dtx2);
         I_A = (mccM(&A) != 1 || mccN(&A) != 1);
         p_A = mccPR(&A);
         I_dt = (mccM(&dt) != 1 || mccN(&dt) != 1);
         p_dt = mccPR(&dt);
         I_C = (mccM(&C) != 1 || mccN(&C) != 1);
         p_C = mccPR(&C);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_b+=I_b, p_dtx2+=I_dtx2, p_A+=I_A, p_dt+=I_dt, p_C+=I_C)
               {
                  *p_b = ((( -2 * (double) *p_dtx2) * (double) *p_A) + (*p_dt * (double) *p_C));
               }
            }
         }
      }
      /* c = dtx2 .* A + dt2x .* B; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_c;
         int I_c=1;
         double *p_dtx2;
         int I_dtx2=1;
         double *p_A;
         int I_A=1;
         double *p_dt2x;
         int I_dt2x=1;
         double *p_B;
         int I_B=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dtx2), mccN(&dtx2));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&A), mccN(&A));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dt2x), mccN(&dt2x));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&B), mccN(&B));
         mccAllocateMatrix(&c, m_, n_);
         I_c = (mccM(&c) != 1 || mccN(&c) != 1);
         p_c = mccPR(&c);
         I_dtx2 = (mccM(&dtx2) != 1 || mccN(&dtx2) != 1);
         p_dtx2 = mccPR(&dtx2);
         I_A = (mccM(&A) != 1 || mccN(&A) != 1);
         p_A = mccPR(&A);
         I_dt2x = (mccM(&dt2x) != 1 || mccN(&dt2x) != 1);
         p_dt2x = mccPR(&dt2x);
         I_B = (mccM(&B) != 1 || mccN(&B) != 1);
         p_B = mccPR(&B);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_c+=I_c, p_dtx2+=I_dtx2, p_A+=I_A, p_dt2x+=I_dt2x, p_B+=I_B)
               {
                  *p_c = ((*p_dtx2 * (double) *p_A) + (*p_dt2x * (double) *p_B));
               }
            }
         }
      }
      /* d = dt .* D; */
      if(mccNOTSET(&D))
      {
         mexErrMsgTxt( "variable D undefined, line 201" );
      }
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_d;
         int I_d=1;
         double *p_dt;
         int I_dt=1;
         double *p_D;
         int I_D=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dt), mccN(&dt));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&D), mccN(&D));
         mccAllocateMatrix(&d, m_, n_);
         I_d = (mccM(&d) != 1 || mccN(&d) != 1);
         p_d = mccPR(&d);
         I_dt = (mccM(&dt) != 1 || mccN(&dt) != 1);
         p_dt = mccPR(&dt);
         I_D = (mccM(&D) != 1 || mccN(&D) != 1);
         p_D = mccPR(&D);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_d+=I_d, p_dt+=I_dt, p_D+=I_D)
               {
                  *p_d = (*p_dt * (double) *p_D);
               }
            }
         }
      }
      
      /* ap = dtx2 .* AP - dt2x .* BP; */
      if(mccNOTSET(&AP))
      {
         mexErrMsgTxt( "variable AP undefined, line 203" );
      }
      if(mccNOTSET(&BP))
      {
         mexErrMsgTxt( "variable BP undefined, line 203" );
      }
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_ap;
         int I_ap=1;
         double *p_dtx2;
         int I_dtx2=1;
         double *p_AP;
         int I_AP=1;
         double *p_dt2x;
         int I_dt2x=1;
         double *p_BP;
         int I_BP=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dtx2), mccN(&dtx2));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&AP), mccN(&AP));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dt2x), mccN(&dt2x));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&BP), mccN(&BP));
         mccAllocateMatrix(&ap, m_, n_);
         I_ap = (mccM(&ap) != 1 || mccN(&ap) != 1);
         p_ap = mccPR(&ap);
         I_dtx2 = (mccM(&dtx2) != 1 || mccN(&dtx2) != 1);
         p_dtx2 = mccPR(&dtx2);
         I_AP = (mccM(&AP) != 1 || mccN(&AP) != 1);
         p_AP = mccPR(&AP);
         I_dt2x = (mccM(&dt2x) != 1 || mccN(&dt2x) != 1);
         p_dt2x = mccPR(&dt2x);
         I_BP = (mccM(&BP) != 1 || mccN(&BP) != 1);
         p_BP = mccPR(&BP);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_ap+=I_ap, p_dtx2+=I_dtx2, p_AP+=I_AP, p_dt2x+=I_dt2x, p_BP+=I_BP)
               {
                  *p_ap = ((*p_dtx2 * (double) *p_AP) - (*p_dt2x * (double) *p_BP));
               }
            }
         }
      }
      /* bp = - 2 .* dtx2 .* AP  + dt .* CP; */
      if(mccNOTSET(&CP))
      {
         mexErrMsgTxt( "variable CP undefined, line 204" );
      }
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_bp;
         int I_bp=1;
         double *p_dtx2;
         int I_dtx2=1;
         double *p_AP;
         int I_AP=1;
         double *p_dt;
         int I_dt=1;
         double *p_CP;
         int I_CP=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dtx2), mccN(&dtx2));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&AP), mccN(&AP));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dt), mccN(&dt));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&CP), mccN(&CP));
         mccAllocateMatrix(&bp, m_, n_);
         I_bp = (mccM(&bp) != 1 || mccN(&bp) != 1);
         p_bp = mccPR(&bp);
         I_dtx2 = (mccM(&dtx2) != 1 || mccN(&dtx2) != 1);
         p_dtx2 = mccPR(&dtx2);
         I_AP = (mccM(&AP) != 1 || mccN(&AP) != 1);
         p_AP = mccPR(&AP);
         I_dt = (mccM(&dt) != 1 || mccN(&dt) != 1);
         p_dt = mccPR(&dt);
         I_CP = (mccM(&CP) != 1 || mccN(&CP) != 1);
         p_CP = mccPR(&CP);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_bp+=I_bp, p_dtx2+=I_dtx2, p_AP+=I_AP, p_dt+=I_dt, p_CP+=I_CP)
               {
                  *p_bp = ((( -2 * (double) *p_dtx2) * (double) *p_AP) + (*p_dt * (double) *p_CP));
               }
            }
         }
      }
      /* cp = dtx2 .* AP + dt2x .* BP; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_cp;
         int I_cp=1;
         double *p_dtx2;
         int I_dtx2=1;
         double *p_AP;
         int I_AP=1;
         double *p_dt2x;
         int I_dt2x=1;
         double *p_BP;
         int I_BP=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dtx2), mccN(&dtx2));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&AP), mccN(&AP));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dt2x), mccN(&dt2x));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&BP), mccN(&BP));
         mccAllocateMatrix(&cp, m_, n_);
         I_cp = (mccM(&cp) != 1 || mccN(&cp) != 1);
         p_cp = mccPR(&cp);
         I_dtx2 = (mccM(&dtx2) != 1 || mccN(&dtx2) != 1);
         p_dtx2 = mccPR(&dtx2);
         I_AP = (mccM(&AP) != 1 || mccN(&AP) != 1);
         p_AP = mccPR(&AP);
         I_dt2x = (mccM(&dt2x) != 1 || mccN(&dt2x) != 1);
         p_dt2x = mccPR(&dt2x);
         I_BP = (mccM(&BP) != 1 || mccN(&BP) != 1);
         p_BP = mccPR(&BP);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_cp+=I_cp, p_dtx2+=I_dtx2, p_AP+=I_AP, p_dt2x+=I_dt2x, p_BP+=I_BP)
               {
                  *p_cp = ((*p_dtx2 * (double) *p_AP) + (*p_dt2x * (double) *p_BP));
               }
            }
         }
      }
      /* dp = dt .* DP; */
      if(mccNOTSET(&DP))
      {
         mexErrMsgTxt( "variable DP undefined, line 206" );
      }
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_dp;
         int I_dp=1;
         double *p_dt;
         int I_dt=1;
         double *p_DP;
         int I_DP=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dt), mccN(&dt));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&DP), mccN(&DP));
         mccAllocateMatrix(&dp, m_, n_);
         I_dp = (mccM(&dp) != 1 || mccN(&dp) != 1);
         p_dp = mccPR(&dp);
         I_dt = (mccM(&dt) != 1 || mccN(&dt) != 1);
         p_dt = mccPR(&dt);
         I_DP = (mccM(&DP) != 1 || mccN(&DP) != 1);
         p_DP = mccPR(&DP);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_dp+=I_dp, p_dt+=I_dt, p_DP+=I_DP)
               {
                  *p_dp = (*p_dt * (double) *p_DP);
               }
            }
         }
      }
      
      /* s = f .* d + (1-f) .* dp; */
      if(mccNOTSET(&f))
      {
         mexErrMsgTxt( "variable f undefined, line 208" );
      }
      if(mccNOTSET(&f))
      {
         mexErrMsgTxt( "variable f undefined, line 208" );
      }
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_s;
         int I_s=1;
         double *p_f;
         int I_f=1;
         double *p_d;
         int I_d=1;
         double *p_1f;
         int I_1f=1;
         double *p_dp;
         int I_dp=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&d), mccN(&d));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dp), mccN(&dp));
         mccAllocateMatrix(&s, m_, n_);
         I_s = (mccM(&s) != 1 || mccN(&s) != 1);
         p_s = mccPR(&s);
         I_f = (mccM(&f) != 1 || mccN(&f) != 1);
         p_f = mccPR(&f);
         I_d = (mccM(&d) != 1 || mccN(&d) != 1);
         p_d = mccPR(&d);
         I_1f = (mccM(&f) != 1 || mccN(&f) != 1);
         p_1f = mccPR(&f);
         I_dp = (mccM(&dp) != 1 || mccN(&dp) != 1);
         p_dp = mccPR(&dp);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_s+=I_s, p_f+=I_f, p_d+=I_d, p_1f+=I_1f, p_dp+=I_dp)
               {
                  *p_s = ((*p_f * (double) *p_d) + ((1 - *p_1f) * (double) *p_dp));
               }
            }
         }
      }
      
      /* % 2 - condition aux limites en x=x0 */
      /* % a - conditionnement des donnees */
      /* u = (T0(:,:,3) ./ (dx.^2) + T0(:,:,2) ./ (2.*dx)); */
      if(mccNOTSET(&T0))
      {
         mexErrMsgTxt( "variable T0 undefined, line 212" );
      }
      if(mccNOTSET(&T0))
      {
         mexErrMsgTxt( "variable T0 undefined, line 212" );
      }
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_u;
         int I_u=1;
         double *p_T0;
         int I_T0=1, J_T0;
         double *p_dx;
         int I_dx=1;
         double *p_1T0;
         int I_1T0=1, J_1T0;
         double *p_1dx;
         int I_1dx=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&T0), 1);
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dx), mccN(&dx));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&T0), 1);
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dx), mccN(&dx));
         mccAllocateMatrix(&u, m_, n_);
         mccCheckMatrixSize(&T0, mccM(&T0), 3);
         mccCheckMatrixSize(&T0, mccM(&T0), 2);
         I_u = (mccM(&u) != 1 || mccN(&u) != 1);
         p_u = mccPR(&u);
         if (mccM(&T0) == 1) { I_T0 = J_T0 = 0;}
         else { I_T0 = 1; J_T0=mccM(&T0)-m_; }
         p_T0 = mccPR(&T0) + 0 + mccM(&T0) * (3-1);
         I_dx = (mccM(&dx) != 1 || mccN(&dx) != 1);
         p_dx = mccPR(&dx);
         if (mccM(&T0) == 1) { I_1T0 = J_1T0 = 0;}
         else { I_1T0 = 1; J_1T0=mccM(&T0)-m_; }
         p_1T0 = mccPR(&T0) + 0 + mccM(&T0) * (2-1);
         I_1dx = (mccM(&dx) != 1 || mccN(&dx) != 1);
         p_1dx = mccPR(&dx);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_T0 += J_T0, p_1T0 += J_1T0)
            {
               for (i_=0; i_<m_; ++i_, p_u+=I_u, p_T0+=I_T0, p_dx+=I_dx, p_1T0+=I_1T0, p_1dx+=I_1dx)
               {
                  *p_u = ((*p_T0 / (double) mcmRealPowerInt(*p_dx, 2)) + (*p_1T0 / (double) (2 * (double) *p_1dx)));
               }
            }
         }
      }
      /* v = (T0(:,:,1) - 2 .* T0(:,:,3) ./ (dx.^2)); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_v;
         int I_v=1;
         double *p_T0;
         int I_T0=1, J_T0;
         double *p_1T0;
         int I_1T0=1, J_1T0;
         double *p_dx;
         int I_dx=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&T0), 1);
         m_ = mcmCalcResultSize(m_, &n_, mccM(&T0), 1);
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dx), mccN(&dx));
         mccAllocateMatrix(&v, m_, n_);
         mccCheckMatrixSize(&T0, mccM(&T0), 1);
         mccCheckMatrixSize(&T0, mccM(&T0), 3);
         I_v = (mccM(&v) != 1 || mccN(&v) != 1);
         p_v = mccPR(&v);
         if (mccM(&T0) == 1) { I_T0 = J_T0 = 0;}
         else { I_T0 = 1; J_T0=mccM(&T0)-m_; }
         p_T0 = mccPR(&T0) + 0 + mccM(&T0) * (1-1);
         if (mccM(&T0) == 1) { I_1T0 = J_1T0 = 0;}
         else { I_1T0 = 1; J_1T0=mccM(&T0)-m_; }
         p_1T0 = mccPR(&T0) + 0 + mccM(&T0) * (3-1);
         I_dx = (mccM(&dx) != 1 || mccN(&dx) != 1);
         p_dx = mccPR(&dx);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_T0 += J_T0, p_1T0 += J_1T0)
            {
               for (i_=0; i_<m_; ++i_, p_v+=I_v, p_T0+=I_T0, p_1T0+=I_1T0, p_dx+=I_dx)
               {
                  *p_v = (*p_T0 - ((2 * (double) *p_1T0) / (double) mcmRealPowerInt(*p_dx, 2)));
               }
            }
         }
      }
      /* w = (T0(:,:,3) ./ (dx.^2) - T0(:,:,2) ./ (2*dx)); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_w;
         int I_w=1;
         double *p_T0;
         int I_T0=1, J_T0;
         double *p_dx;
         int I_dx=1;
         double *p_1T0;
         int I_1T0=1, J_1T0;
         double *p_1dx;
         int I_1dx=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&T0), 1);
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dx), mccN(&dx));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&T0), 1);
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dx), mccN(&dx));
         mccAllocateMatrix(&w, m_, n_);
         mccCheckMatrixSize(&T0, mccM(&T0), 3);
         mccCheckMatrixSize(&T0, mccM(&T0), 2);
         I_w = (mccM(&w) != 1 || mccN(&w) != 1);
         p_w = mccPR(&w);
         if (mccM(&T0) == 1) { I_T0 = J_T0 = 0;}
         else { I_T0 = 1; J_T0=mccM(&T0)-m_; }
         p_T0 = mccPR(&T0) + 0 + mccM(&T0) * (3-1);
         I_dx = (mccM(&dx) != 1 || mccN(&dx) != 1);
         p_dx = mccPR(&dx);
         if (mccM(&T0) == 1) { I_1T0 = J_1T0 = 0;}
         else { I_1T0 = 1; J_1T0=mccM(&T0)-m_; }
         p_1T0 = mccPR(&T0) + 0 + mccM(&T0) * (2-1);
         I_1dx = (mccM(&dx) != 1 || mccN(&dx) != 1);
         p_1dx = mccPR(&dx);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_T0 += J_T0, p_1T0 += J_1T0)
            {
               for (i_=0; i_<m_; ++i_, p_w+=I_w, p_T0+=I_T0, p_dx+=I_dx, p_1T0+=I_1T0, p_1dx+=I_1dx)
               {
                  *p_w = ((*p_T0 / (double) mcmRealPowerInt(*p_dx, 2)) - (*p_1T0 / (double) (2 * (double) *p_1dx)));
               }
            }
         }
      }
      
      /* warning off */
      Mprhs_[0] = &S42_;
      Mplhs_[0] = 0;
      mccCallMATLAB(0, Mplhs_, 1, Mprhs_, "warning", 216);
      mccPrint(Mplhs_[0], 0);
      /* U = - u ./ w; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_U;
         int I_U=1;
         double *p_u;
         int I_u=1;
         double *p_w;
         int I_w=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&u), mccN(&u));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&w), mccN(&w));
         mccAllocateMatrix(&U, m_, n_);
         I_U = (mccM(&U) != 1 || mccN(&U) != 1);
         p_U = mccPR(&U);
         I_u = (mccM(&u) != 1 || mccN(&u) != 1);
         p_u = mccPR(&u);
         I_w = (mccM(&w) != 1 || mccN(&w) != 1);
         p_w = mccPR(&w);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_U+=I_U, p_u+=I_u, p_w+=I_w)
               {
                  *p_U = ((-*p_u) / (double) *p_w);
               }
            }
         }
      }
      /* V = - v ./ w; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_V;
         int I_V=1;
         double *p_v;
         int I_v=1;
         double *p_w;
         int I_w=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&v), mccN(&v));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&w), mccN(&w));
         mccAllocateMatrix(&V, m_, n_);
         I_V = (mccM(&V) != 1 || mccN(&V) != 1);
         p_V = mccPR(&V);
         I_v = (mccM(&v) != 1 || mccN(&v) != 1);
         p_v = mccPR(&v);
         I_w = (mccM(&w) != 1 || mccN(&w) != 1);
         p_w = mccPR(&w);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_V+=I_V, p_v+=I_v, p_w+=I_w)
               {
                  *p_V = ((-*p_v) / (double) *p_w);
               }
            }
         }
      }
      /* W = V0 ./ w; */
      if(mccNOTSET(&V0))
      {
         mexErrMsgTxt( "variable V0 undefined, line 219" );
      }
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_W;
         int I_W=1;
         double *p_V0;
         int I_V0=1;
         double *p_w;
         int I_w=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&V0), mccN(&V0));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&w), mccN(&w));
         mccAllocateMatrix(&W, m_, n_);
         I_W = (mccM(&W) != 1 || mccN(&W) != 1);
         p_W = mccPR(&W);
         I_V0 = (mccM(&V0) != 1 || mccN(&V0) != 1);
         p_V0 = mccPR(&V0);
         I_w = (mccM(&w) != 1 || mccN(&w) != 1);
         p_w = mccPR(&w);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_W+=I_W, p_V0+=I_V0, p_w+=I_w)
               {
                  *p_W = (*p_V0 / (double) *p_w);
               }
            }
         }
      }
      /* X = - u ./ v; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_X;
         int I_X=1;
         double *p_u;
         int I_u=1;
         double *p_v;
         int I_v=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&u), mccN(&u));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&v), mccN(&v));
         mccAllocateMatrix(&X, m_, n_);
         I_X = (mccM(&X) != 1 || mccN(&X) != 1);
         p_X = mccPR(&X);
         I_u = (mccM(&u) != 1 || mccN(&u) != 1);
         p_u = mccPR(&u);
         I_v = (mccM(&v) != 1 || mccN(&v) != 1);
         p_v = mccPR(&v);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_X+=I_X, p_u+=I_u, p_v+=I_v)
               {
                  *p_X = ((-*p_u) / (double) *p_v);
               }
            }
         }
      }
      /* Y = V0 ./ v; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_Y;
         int I_Y=1;
         double *p_V0;
         int I_V0=1;
         double *p_v;
         int I_v=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&V0), mccN(&V0));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&v), mccN(&v));
         mccAllocateMatrix(&Y, m_, n_);
         I_Y = (mccM(&Y) != 1 || mccN(&Y) != 1);
         p_Y = mccPR(&Y);
         I_V0 = (mccM(&V0) != 1 || mccN(&V0) != 1);
         p_V0 = mccPR(&V0);
         I_v = (mccM(&v) != 1 || mccN(&v) != 1);
         p_v = mccPR(&v);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_Y+=I_Y, p_V0+=I_V0, p_v+=I_v)
               {
                  *p_Y = (*p_V0 / (double) *p_v);
               }
            }
         }
      }
      /* Z = V0 ./ u; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_Z;
         int I_Z=1;
         double *p_V0;
         int I_V0=1;
         double *p_u;
         int I_u=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&V0), mccN(&V0));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&u), mccN(&u));
         mccAllocateMatrix(&Z, m_, n_);
         I_Z = (mccM(&Z) != 1 || mccN(&Z) != 1);
         p_Z = mccPR(&Z);
         I_V0 = (mccM(&V0) != 1 || mccN(&V0) != 1);
         p_V0 = mccPR(&V0);
         I_u = (mccM(&u) != 1 || mccN(&u) != 1);
         p_u = mccPR(&u);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_Z+=I_Z, p_V0+=I_V0, p_u+=I_u)
               {
                  *p_Z = (*p_V0 / (double) *p_u);
               }
            }
         }
      }
      /* warning on */
      Mprhs_[0] = &S43_;
      Mplhs_[0] = 0;
      mccCallMATLAB(0, Mplhs_, 1, Mprhs_, "warning", 223);
      mccPrint(Mplhs_[0], 0);
      
      /* comp = ones(1,M); */
      mccOnesMN(&M_comp, 1, M);
      
      /* % b - cas  w~=0 */
      /* i1 = find( w(1,:) ~= 0 ); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         unsigned short *p_BM0_;
         int I_BM0_=1;
         double *p_w;
         int I_w=1, J_w;
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&w));
         mccAllocateMatrix(&BM0_, m_, n_);
         mccCheckMatrixSize(&w, 1, mccN(&w));
         I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
         p_BM0_ = mccSPR(&BM0_);
         if (mccN(&w) == 1) { I_w = J_w = 0;}
         else { I_w = 1; J_w=mccM(&w)-m_; }
         p_w = mccPR(&w) + (1-1) + mccM(&w) * 0;
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_w += J_w)
            {
               for (i_=0; i_<m_; ++i_, p_BM0_+=I_BM0_, p_w+=I_w)
               {
                  *p_BM0_ = ( (*p_w != 0) || mccREL_NAN(*p_w) );
               }
            }
         }
      }
      mccFind(&i1, &BM0_);
      /* i2 = find( w(2,:) ~=0) ; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         unsigned short *p_BM0_;
         int I_BM0_=1;
         double *p_w;
         int I_w=1, J_w;
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&w));
         mccAllocateMatrix(&BM0_, m_, n_);
         mccCheckMatrixSize(&w, 2, mccN(&w));
         I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
         p_BM0_ = mccSPR(&BM0_);
         if (mccN(&w) == 1) { I_w = J_w = 0;}
         else { I_w = 1; J_w=mccM(&w)-m_; }
         p_w = mccPR(&w) + (2-1) + mccM(&w) * 0;
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_w += J_w)
            {
               for (i_=0; i_<m_; ++i_, p_BM0_+=I_BM0_, p_w+=I_w)
               {
                  *p_BM0_ = ( (*p_w != 0) || mccREL_NAN(*p_w) );
               }
            }
         }
      }
      mccFind(&i2, &BM0_);
      /* if ~ (isempty(i1) & isempty(i2)) */
      B0_ = mccIsEmpty(&i1);
      B1_ = mccIsEmpty(&i2);
      B2_ = (!(!!B0_ && !!B1_));
      if ((double)B2_)
      {
         /* b(1,i1,:) = b(1,i1,:) + a(1,i1,:) .* permute(V(1,i1)' * comp,[3 1 2]); */
         mccFindIndex(&IM0_, &i1);
         mccFindIndex(&IM2_, &i1);
         mccFindIndex(&IM3_, &i1);
         mccFindIndex(&IM4_, &IM3_);
         mccFindIndex(&IM5_, &IM4_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM0_;
            int I_RM0_=1;
            double *p_V;
            int I_V=1;
            int *p_IM5_;
            int I_IM5_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM5_) * mccN(&IM5_)));
            mccAllocateMatrix(&RM0_, m_, n_);
            mccCheckMatrixSize(&V, 1, mccGetMaxIndex(&IM5_ ,mccN(&V)));
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            I_IM5_ = (mccM(&IM5_) != 1 || mccN(&IM5_) != 1);
            p_IM5_ = mccIPR(&IM5_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM5_ += I_IM5_)
               {
                  p_V = mccPR(&V) + mccM(&V) * ((int)(*p_IM5_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_V+=I_V)
                  {
                     *p_RM0_ = *p_V;
                  }
               }
            }
         }
         mccConjTrans(&RM1_, &RM0_);
         mccRealMatrixMultiply(&RM2_, &RM1_, &M_comp);
         Mprhs_[0] = &RM2_;
         Mprhs_[1] = &S0_;
         Mplhs_[0] = &RM3_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 231);
         mccFindIndex(&IM1_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_b;
            int I_b=1, J_b;
            double *p_1b;
            int I_1b=1, J_1b;
            double *p_a;
            int I_a=1, J_a;
            double *p_RM3_;
            int I_RM3_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&a), mccN(&a));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            if (!mccNOTSET(&b)) m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            mccGrowMatrix(&b, m_, n_);
            mccCheckMatrixSize(&b, mccM(&b), mccN(&b));
            mccCheckMatrixSize(&a, mccM(&a), mccN(&a));
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_b = J_b = 0;}
            else { I_b = 1; J_b=mccM(&b)-m_; }
            p_b = mccPR(&b) + 0 + mccM(&b) * 0;
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_1b = J_1b = 0;}
            else { I_1b = 1; J_1b=mccM(&b)-m_; }
            p_1b = mccPR(&b) + 0 + mccM(&b) * 0;
            if (mccM(&a) == 1 && mccN(&a) == 1) { I_a = J_a = 0;}
            else { I_a = 1; J_a=mccM(&a)-m_; }
            p_a = mccPR(&a) + 0 + mccM(&a) * 0;
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_b += J_b, p_1b += J_1b, p_a += J_a)
               {
                  for (i_=0; i_<m_; ++i_, p_b+=I_b, p_1b+=I_1b, p_a+=I_a, p_RM3_+=I_RM3_)
                  {
                     *p_b = (*p_1b + (*p_a * (double) *p_RM3_));
                  }
               }
            }
         }
         /* bp(1,i2,:) = bp(1,i2,:) + ap(1,i2,:) .* permute(V(2,i2)' * comp,[3 1 2]); */
         mccFindIndex(&IM4_, &i2);
         mccFindIndex(&IM3_, &i2);
         mccFindIndex(&IM2_, &i2);
         mccFindIndex(&IM0_, &IM2_);
         mccFindIndex(&IM1_, &IM0_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM3_;
            int I_RM3_=1;
            double *p_V;
            int I_V=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM1_) * mccN(&IM1_)));
            mccAllocateMatrix(&RM3_, m_, n_);
            mccCheckMatrixSize(&V, 2, mccGetMaxIndex(&IM1_ ,mccN(&V)));
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_V = mccPR(&V) + mccM(&V) * ((int)(*p_IM1_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM3_+=I_RM3_, p_V+=I_V)
                  {
                     *p_RM3_ = *p_V;
                  }
               }
            }
         }
         mccConjTrans(&RM2_, &RM3_);
         mccRealMatrixMultiply(&RM1_, &RM2_, &M_comp);
         Mprhs_[0] = &RM1_;
         Mprhs_[1] = &S1_;
         Mplhs_[0] = &RM0_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 232);
         mccFindIndex(&IM5_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_bp;
            int I_bp=1, J_bp;
            double *p_1bp;
            int I_1bp=1, J_1bp;
            double *p_ap;
            int I_ap=1, J_ap;
            double *p_RM0_;
            int I_RM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), mccN(&ap));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM0_), mccN(&RM0_));
            if (!mccNOTSET(&bp)) m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            mccGrowMatrix(&bp, m_, n_);
            mccCheckMatrixSize(&bp, mccM(&bp), mccN(&bp));
            mccCheckMatrixSize(&ap, mccM(&ap), mccN(&ap));
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_bp = J_bp = 0;}
            else { I_bp = 1; J_bp=mccM(&bp)-m_; }
            p_bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_1bp = J_1bp = 0;}
            else { I_1bp = 1; J_1bp=mccM(&bp)-m_; }
            p_1bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            if (mccM(&ap) == 1 && mccN(&ap) == 1) { I_ap = J_ap = 0;}
            else { I_ap = 1; J_ap=mccM(&ap)-m_; }
            p_ap = mccPR(&ap) + 0 + mccM(&ap) * 0;
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_bp += J_bp, p_1bp += J_1bp, p_ap += J_ap)
               {
                  for (i_=0; i_<m_; ++i_, p_bp+=I_bp, p_1bp+=I_1bp, p_ap+=I_ap, p_RM0_+=I_RM0_)
                  {
                     *p_bp = (*p_1bp + (*p_ap * (double) *p_RM0_));
                  }
               }
            }
         }
         
         /* c(1,i1,:) = c(1,i1,:) + a(1,i1,:) .* permute(U(1,i1)' * comp,[3 1 2]); */
         mccFindIndex(&IM0_, &i1);
         mccFindIndex(&IM2_, &i1);
         mccFindIndex(&IM3_, &i1);
         mccFindIndex(&IM4_, &IM3_);
         mccFindIndex(&IM5_, &IM4_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM0_;
            int I_RM0_=1;
            double *p_U;
            int I_U=1;
            int *p_IM5_;
            int I_IM5_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM5_) * mccN(&IM5_)));
            mccAllocateMatrix(&RM0_, m_, n_);
            mccCheckMatrixSize(&U, 1, mccGetMaxIndex(&IM5_ ,mccN(&U)));
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            I_IM5_ = (mccM(&IM5_) != 1 || mccN(&IM5_) != 1);
            p_IM5_ = mccIPR(&IM5_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM5_ += I_IM5_)
               {
                  p_U = mccPR(&U) + mccM(&U) * ((int)(*p_IM5_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_U+=I_U)
                  {
                     *p_RM0_ = *p_U;
                  }
               }
            }
         }
         mccConjTrans(&RM1_, &RM0_);
         mccRealMatrixMultiply(&RM2_, &RM1_, &M_comp);
         Mprhs_[0] = &RM2_;
         Mprhs_[1] = &S2_;
         Mplhs_[0] = &RM3_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 234);
         mccFindIndex(&IM1_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_c;
            int I_c=1, J_c;
            double *p_1c;
            int I_1c=1, J_1c;
            double *p_a;
            int I_a=1, J_a;
            double *p_RM3_;
            int I_RM3_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&a), mccN(&a));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            if (!mccNOTSET(&c)) m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            mccGrowMatrix(&c, m_, n_);
            mccCheckMatrixSize(&c, mccM(&c), mccN(&c));
            mccCheckMatrixSize(&a, mccM(&a), mccN(&a));
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_c = J_c = 0;}
            else { I_c = 1; J_c=mccM(&c)-m_; }
            p_c = mccPR(&c) + 0 + mccM(&c) * 0;
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_1c = J_1c = 0;}
            else { I_1c = 1; J_1c=mccM(&c)-m_; }
            p_1c = mccPR(&c) + 0 + mccM(&c) * 0;
            if (mccM(&a) == 1 && mccN(&a) == 1) { I_a = J_a = 0;}
            else { I_a = 1; J_a=mccM(&a)-m_; }
            p_a = mccPR(&a) + 0 + mccM(&a) * 0;
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_c += J_c, p_1c += J_1c, p_a += J_a)
               {
                  for (i_=0; i_<m_; ++i_, p_c+=I_c, p_1c+=I_1c, p_a+=I_a, p_RM3_+=I_RM3_)
                  {
                     *p_c = (*p_1c + (*p_a * (double) *p_RM3_));
                  }
               }
            }
         }
         /* cp(1,i2,:) = cp(1,i2,:) + ap(1,i2,:) .* permute(U(2,i2)' * comp,[3 1 2]); */
         mccFindIndex(&IM4_, &i2);
         mccFindIndex(&IM3_, &i2);
         mccFindIndex(&IM2_, &i2);
         mccFindIndex(&IM0_, &IM2_);
         mccFindIndex(&IM1_, &IM0_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM3_;
            int I_RM3_=1;
            double *p_U;
            int I_U=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM1_) * mccN(&IM1_)));
            mccAllocateMatrix(&RM3_, m_, n_);
            mccCheckMatrixSize(&U, 2, mccGetMaxIndex(&IM1_ ,mccN(&U)));
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_U = mccPR(&U) + mccM(&U) * ((int)(*p_IM1_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM3_+=I_RM3_, p_U+=I_U)
                  {
                     *p_RM3_ = *p_U;
                  }
               }
            }
         }
         mccConjTrans(&RM2_, &RM3_);
         mccRealMatrixMultiply(&RM1_, &RM2_, &M_comp);
         Mprhs_[0] = &RM1_;
         Mprhs_[1] = &S3_;
         Mplhs_[0] = &RM0_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 235);
         mccFindIndex(&IM5_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_cp;
            int I_cp=1, J_cp;
            double *p_1cp;
            int I_1cp=1, J_1cp;
            double *p_ap;
            int I_ap=1, J_ap;
            double *p_RM0_;
            int I_RM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&cp), mccN(&cp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), mccN(&ap));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM0_), mccN(&RM0_));
            if (!mccNOTSET(&cp)) m_ = mcmCalcResultSize(m_, &n_, mccM(&cp), mccN(&cp));
            mccGrowMatrix(&cp, m_, n_);
            mccCheckMatrixSize(&cp, mccM(&cp), mccN(&cp));
            mccCheckMatrixSize(&ap, mccM(&ap), mccN(&ap));
            if (mccM(&cp) == 1 && mccN(&cp) == 1) { I_cp = J_cp = 0;}
            else { I_cp = 1; J_cp=mccM(&cp)-m_; }
            p_cp = mccPR(&cp) + 0 + mccM(&cp) * 0;
            if (mccM(&cp) == 1 && mccN(&cp) == 1) { I_1cp = J_1cp = 0;}
            else { I_1cp = 1; J_1cp=mccM(&cp)-m_; }
            p_1cp = mccPR(&cp) + 0 + mccM(&cp) * 0;
            if (mccM(&ap) == 1 && mccN(&ap) == 1) { I_ap = J_ap = 0;}
            else { I_ap = 1; J_ap=mccM(&ap)-m_; }
            p_ap = mccPR(&ap) + 0 + mccM(&ap) * 0;
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_cp += J_cp, p_1cp += J_1cp, p_ap += J_ap)
               {
                  for (i_=0; i_<m_; ++i_, p_cp+=I_cp, p_1cp+=I_1cp, p_ap+=I_ap, p_RM0_+=I_RM0_)
                  {
                     *p_cp = (*p_1cp + (*p_ap * (double) *p_RM0_));
                  }
               }
            }
         }
         
         /* s(1,:) = s(1,:) + squeeze( f .* sum( a(1,i1,:) .* permute(W(1,i1)' * comp,[3 1 2]),2) + ... */
         mccFindIndex(&IM1_, &i1);
         mccFindIndex(&IM0_, &IM1_);
         mccFindIndex(&IM2_, &i1);
         mccFindIndex(&IM3_, &IM2_);
         mccFindIndex(&IM4_, &IM3_);
         mccFindIndex(&IM5_, &IM4_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM0_;
            int I_RM0_=1;
            double *p_W;
            int I_W=1;
            int *p_IM5_;
            int I_IM5_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM5_) * mccN(&IM5_)));
            mccAllocateMatrix(&RM0_, m_, n_);
            mccCheckMatrixSize(&W, 1, mccGetMaxIndex(&IM5_ ,mccN(&W)));
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            I_IM5_ = (mccM(&IM5_) != 1 || mccN(&IM5_) != 1);
            p_IM5_ = mccIPR(&IM5_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM5_ += I_IM5_)
               {
                  p_W = mccPR(&W) + mccM(&W) * ((int)(*p_IM5_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_W+=I_W)
                  {
                     *p_RM0_ = *p_W;
                  }
               }
            }
         }
         mccConjTrans(&RM1_, &RM0_);
         mccRealMatrixMultiply(&RM2_, &RM1_, &M_comp);
         Mprhs_[0] = &RM2_;
         Mprhs_[1] = &S4_;
         Mplhs_[0] = &RM3_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 237);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM4_;
            int I_RM4_=1;
            double *p_a;
            int I_a=1, J_a;
            double *p_RM3_;
            int I_RM3_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&a), mccN(&a));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            mccAllocateMatrix(&RM4_, m_, n_);
            mccCheckMatrixSize(&a, mccM(&a), mccN(&a));
            I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
            p_RM4_ = mccPR(&RM4_);
            if (mccM(&a) == 1 && mccN(&a) == 1) { I_a = J_a = 0;}
            else { I_a = 1; J_a=mccM(&a)-m_; }
            p_a = mccPR(&a) + 0 + mccM(&a) * 0;
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_a += J_a)
               {
                  for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_a+=I_a, p_RM3_+=I_RM3_)
                  {
                     *p_RM4_ = (*p_a * (double) *p_RM3_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM4_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM5_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 237);
         mccFindIndex(&IM6_, &i2);
         mccFindIndex(&IM7_, &IM6_);
         mccFindIndex(&IM8_, &i2);
         mccFindIndex(&IM9_, &IM8_);
         mccFindIndex(&IM10_, &IM9_);
         mccFindIndex(&IM11_, &IM10_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM6_;
            int I_RM6_=1;
            double *p_W;
            int I_W=1;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM11_) * mccN(&IM11_)));
            mccAllocateMatrix(&RM6_, m_, n_);
            mccCheckMatrixSize(&W, 2, mccGetMaxIndex(&IM11_ ,mccN(&W)));
            I_RM6_ = (mccM(&RM6_) != 1 || mccN(&RM6_) != 1);
            p_RM6_ = mccPR(&RM6_);
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM11_ += I_IM11_)
               {
                  p_W = mccPR(&W) + mccM(&W) * ((int)(*p_IM11_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM6_+=I_RM6_, p_W+=I_W)
                  {
                     *p_RM6_ = *p_W;
                  }
               }
            }
         }
         mccConjTrans(&RM7_, &RM6_);
         mccRealMatrixMultiply(&RM8_, &RM7_, &M_comp);
         Mprhs_[0] = &RM8_;
         Mprhs_[1] = &S5_;
         Mplhs_[0] = &RM9_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 238);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM10_;
            int I_RM10_=1;
            double *p_ap;
            int I_ap=1, J_ap;
            double *p_RM9_;
            int I_RM9_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), mccN(&ap));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM9_), mccN(&RM9_));
            mccAllocateMatrix(&RM10_, m_, n_);
            mccCheckMatrixSize(&ap, mccM(&ap), mccN(&ap));
            I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
            p_RM10_ = mccPR(&RM10_);
            if (mccM(&ap) == 1 && mccN(&ap) == 1) { I_ap = J_ap = 0;}
            else { I_ap = 1; J_ap=mccM(&ap)-m_; }
            p_ap = mccPR(&ap) + 0 + mccM(&ap) * 0;
            I_RM9_ = (mccM(&RM9_) != 1 || mccN(&RM9_) != 1);
            p_RM9_ = mccPR(&RM9_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_ap += J_ap)
               {
                  for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_ap+=I_ap, p_RM9_+=I_RM9_)
                  {
                     *p_RM10_ = (*p_ap * (double) *p_RM9_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM10_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM11_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 238);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM12_;
            int I_RM12_=1;
            double *p_f;
            int I_f=1;
            double *p_RM5_;
            int I_RM5_=1;
            double *p_1f;
            int I_1f=1;
            double *p_RM11_;
            int I_RM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM5_), mccN(&RM5_));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM11_), mccN(&RM11_));
            mccAllocateMatrix(&RM12_, m_, n_);
            I_RM12_ = (mccM(&RM12_) != 1 || mccN(&RM12_) != 1);
            p_RM12_ = mccPR(&RM12_);
            I_f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_f = mccPR(&f);
            I_RM5_ = (mccM(&RM5_) != 1 || mccN(&RM5_) != 1);
            p_RM5_ = mccPR(&RM5_);
            I_1f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_1f = mccPR(&f);
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM12_+=I_RM12_, p_f+=I_f, p_RM5_+=I_RM5_, p_1f+=I_1f, p_RM11_+=I_RM11_)
                  {
                     *p_RM12_ = ((*p_f * (double) *p_RM5_) + ((1 - *p_1f) * (double) *p_RM11_));
                  }
               }
            }
         }
         mccLOG(&RM12_) = 0;
         mccSTRING(&RM12_) = 0;
         Mprhs_[0] = &RM12_;
         Mplhs_[0] = &RM13_;
         mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "squeeze", 238);
         mccConjTrans(&RM14_, &RM13_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_s;
            int I_s=1, J_s;
            double *p_1s;
            int I_1s=1, J_1s;
            double *p_RM14_;
            int I_RM14_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM14_), mccN(&RM14_));
            if (!mccNOTSET(&s)) m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            mccGrowMatrix(&s, 1, n_);
            mccCheckMatrixSize(&s, 1, mccN(&s));
            if (mccN(&s) == 1) { I_s = J_s = 0;}
            else { I_s = 1; J_s=mccM(&s)-m_; }
            p_s = mccPR(&s) + (1-1) + mccM(&s) * 0;
            if (mccN(&s) == 1) { I_1s = J_1s = 0;}
            else { I_1s = 1; J_1s=mccM(&s)-m_; }
            p_1s = mccPR(&s) + (1-1) + mccM(&s) * 0;
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_s += J_s, p_1s += J_1s)
               {
                  for (i_=0; i_<m_; ++i_, p_s+=I_s, p_1s+=I_1s, p_RM14_+=I_RM14_)
                  {
                     *p_s = (*p_1s + *p_RM14_);
                  }
               }
            }
         }
         
         /* a(1,i1,:)   = zeros(1,length(i1),M); */
         I0_ = mccGetLength(&i1);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM10_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 240);
         mccFindIndex(&IM11_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_a;
            int I_a=1, J_a;
            int *p_IM10_;
            int I_IM10_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM10_), mccN(&IM10_));
            if (!mccNOTSET(&a)) m_ = mcmCalcResultSize(m_, &n_, mccM(&a), mccN(&a));
            mccGrowMatrix(&a, m_, n_);
            if (mccM(&a) == 1 && mccN(&a) == 1) { I_a = J_a = 0;}
            else { I_a = 1; J_a=mccM(&a)-m_; }
            p_a = mccPR(&a) + 0 + mccM(&a) * 0;
            I_IM10_ = (mccM(&IM10_) != 1 || mccN(&IM10_) != 1);
            p_IM10_ = mccIPR(&IM10_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_a += J_a)
               {
                  for (i_=0; i_<m_; ++i_, p_a+=I_a, p_IM10_+=I_IM10_)
                  {
                     *p_a = ((int)*p_IM10_);
                  }
               }
            }
         }
         /* ap(1,i2,:)  = zeros(1,length(i2),M); */
         I0_ = mccGetLength(&i2);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM11_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 241);
         mccFindIndex(&IM10_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_ap;
            int I_ap=1, J_ap;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM11_), mccN(&IM11_));
            if (!mccNOTSET(&ap)) m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), mccN(&ap));
            mccGrowMatrix(&ap, m_, n_);
            if (mccM(&ap) == 1 && mccN(&ap) == 1) { I_ap = J_ap = 0;}
            else { I_ap = 1; J_ap=mccM(&ap)-m_; }
            p_ap = mccPR(&ap) + 0 + mccM(&ap) * 0;
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_ap += J_ap)
               {
                  for (i_=0; i_<m_; ++i_, p_ap+=I_ap, p_IM11_+=I_IM11_)
                  {
                     *p_ap = ((int)*p_IM11_);
                  }
               }
            }
         }
         /* end */
      }
      
      /* % c - cas w ==0 et v ~=0 */
      /* i1 = find(( w(1,:) == 0 ) & (v(1,:) ~= 0)); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         unsigned short *p_BM0_;
         int I_BM0_=1;
         double *p_w;
         int I_w=1, J_w;
         double *p_v;
         int I_v=1, J_v;
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&w));
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&v));
         mccAllocateMatrix(&BM0_, m_, n_);
         mccCheckMatrixSize(&w, 1, mccN(&w));
         mccCheckMatrixSize(&v, 1, mccN(&v));
         I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
         p_BM0_ = mccSPR(&BM0_);
         if (mccN(&w) == 1) { I_w = J_w = 0;}
         else { I_w = 1; J_w=mccM(&w)-m_; }
         p_w = mccPR(&w) + (1-1) + mccM(&w) * 0;
         if (mccN(&v) == 1) { I_v = J_v = 0;}
         else { I_v = 1; J_v=mccM(&v)-m_; }
         p_v = mccPR(&v) + (1-1) + mccM(&v) * 0;
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_w += J_w, p_v += J_v)
            {
               for (i_=0; i_<m_; ++i_, p_BM0_+=I_BM0_, p_w+=I_w, p_v+=I_v)
               {
                  *p_BM0_ = (!!( (*p_w == 0) && !mccREL_NAN(*p_w) ) && !!( (*p_v != 0) || mccREL_NAN(*p_v) ));
               }
            }
         }
      }
      mccFind(&i1, &BM0_);
      /* i2 = find(( w(2,:) == 0 ) & (v(2,:) ~= 0)); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         unsigned short *p_BM0_;
         int I_BM0_=1;
         double *p_w;
         int I_w=1, J_w;
         double *p_v;
         int I_v=1, J_v;
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&w));
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&v));
         mccAllocateMatrix(&BM0_, m_, n_);
         mccCheckMatrixSize(&w, 2, mccN(&w));
         mccCheckMatrixSize(&v, 2, mccN(&v));
         I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
         p_BM0_ = mccSPR(&BM0_);
         if (mccN(&w) == 1) { I_w = J_w = 0;}
         else { I_w = 1; J_w=mccM(&w)-m_; }
         p_w = mccPR(&w) + (2-1) + mccM(&w) * 0;
         if (mccN(&v) == 1) { I_v = J_v = 0;}
         else { I_v = 1; J_v=mccM(&v)-m_; }
         p_v = mccPR(&v) + (2-1) + mccM(&v) * 0;
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_w += J_w, p_v += J_v)
            {
               for (i_=0; i_<m_; ++i_, p_BM0_+=I_BM0_, p_w+=I_w, p_v+=I_v)
               {
                  *p_BM0_ = (!!( (*p_w == 0) && !mccREL_NAN(*p_w) ) && !!( (*p_v != 0) || mccREL_NAN(*p_v) ));
               }
            }
         }
      }
      mccFind(&i2, &BM0_);
      /* if ~ (isempty(i1) & isempty(i2)) */
      B2_ = mccIsEmpty(&i1);
      B1_ = mccIsEmpty(&i2);
      B0_ = (!(!!B2_ && !!B1_));
      if ((double)B0_)
      {
         /* c(1,i1,:) = c(1,i1,:) + b(1,i1,:) .* permute(X(1,i1)' * comp,[3 1 2]); */
         mccFindIndex(&IM10_, &i1);
         mccFindIndex(&IM9_, &i1);
         mccFindIndex(&IM8_, &i1);
         mccFindIndex(&IM7_, &IM8_);
         mccFindIndex(&IM6_, &IM7_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM14_;
            int I_RM14_=1;
            double *p_X;
            int I_X=1;
            int *p_IM6_;
            int I_IM6_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM6_) * mccN(&IM6_)));
            mccAllocateMatrix(&RM14_, m_, n_);
            mccCheckMatrixSize(&X, 1, mccGetMaxIndex(&IM6_ ,mccN(&X)));
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            I_IM6_ = (mccM(&IM6_) != 1 || mccN(&IM6_) != 1);
            p_IM6_ = mccIPR(&IM6_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM6_ += I_IM6_)
               {
                  p_X = mccPR(&X) + mccM(&X) * ((int)(*p_IM6_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM14_+=I_RM14_, p_X+=I_X)
                  {
                     *p_RM14_ = *p_X;
                  }
               }
            }
         }
         mccConjTrans(&RM13_, &RM14_);
         mccRealMatrixMultiply(&RM12_, &RM13_, &M_comp);
         Mprhs_[0] = &RM12_;
         Mprhs_[1] = &S6_;
         Mplhs_[0] = &RM11_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 248);
         mccFindIndex(&IM11_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_c;
            int I_c=1, J_c;
            double *p_1c;
            int I_1c=1, J_1c;
            double *p_b;
            int I_b=1, J_b;
            double *p_RM11_;
            int I_RM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM11_), mccN(&RM11_));
            if (!mccNOTSET(&c)) m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            mccGrowMatrix(&c, m_, n_);
            mccCheckMatrixSize(&c, mccM(&c), mccN(&c));
            mccCheckMatrixSize(&b, mccM(&b), mccN(&b));
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_c = J_c = 0;}
            else { I_c = 1; J_c=mccM(&c)-m_; }
            p_c = mccPR(&c) + 0 + mccM(&c) * 0;
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_1c = J_1c = 0;}
            else { I_1c = 1; J_1c=mccM(&c)-m_; }
            p_1c = mccPR(&c) + 0 + mccM(&c) * 0;
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_b = J_b = 0;}
            else { I_b = 1; J_b=mccM(&b)-m_; }
            p_b = mccPR(&b) + 0 + mccM(&b) * 0;
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_c += J_c, p_1c += J_1c, p_b += J_b)
               {
                  for (i_=0; i_<m_; ++i_, p_c+=I_c, p_1c+=I_1c, p_b+=I_b, p_RM11_+=I_RM11_)
                  {
                     *p_c = (*p_1c + (*p_b * (double) *p_RM11_));
                  }
               }
            }
         }
         /* cp(1,i2,:) = cp(1,i2,:) + bp(1,i2,:) .* permute(X(2,i2)' * comp,[3 1 2]); */
         mccFindIndex(&IM7_, &i2);
         mccFindIndex(&IM8_, &i2);
         mccFindIndex(&IM9_, &i2);
         mccFindIndex(&IM10_, &IM9_);
         mccFindIndex(&IM11_, &IM10_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM11_;
            int I_RM11_=1;
            double *p_X;
            int I_X=1;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM11_) * mccN(&IM11_)));
            mccAllocateMatrix(&RM11_, m_, n_);
            mccCheckMatrixSize(&X, 2, mccGetMaxIndex(&IM11_ ,mccN(&X)));
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM11_ += I_IM11_)
               {
                  p_X = mccPR(&X) + mccM(&X) * ((int)(*p_IM11_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM11_+=I_RM11_, p_X+=I_X)
                  {
                     *p_RM11_ = *p_X;
                  }
               }
            }
         }
         mccConjTrans(&RM12_, &RM11_);
         mccRealMatrixMultiply(&RM13_, &RM12_, &M_comp);
         Mprhs_[0] = &RM13_;
         Mprhs_[1] = &S7_;
         Mplhs_[0] = &RM14_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 249);
         mccFindIndex(&IM6_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_cp;
            int I_cp=1, J_cp;
            double *p_1cp;
            int I_1cp=1, J_1cp;
            double *p_bp;
            int I_bp=1, J_bp;
            double *p_RM14_;
            int I_RM14_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&cp), mccN(&cp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM14_), mccN(&RM14_));
            if (!mccNOTSET(&cp)) m_ = mcmCalcResultSize(m_, &n_, mccM(&cp), mccN(&cp));
            mccGrowMatrix(&cp, m_, n_);
            mccCheckMatrixSize(&cp, mccM(&cp), mccN(&cp));
            mccCheckMatrixSize(&bp, mccM(&bp), mccN(&bp));
            if (mccM(&cp) == 1 && mccN(&cp) == 1) { I_cp = J_cp = 0;}
            else { I_cp = 1; J_cp=mccM(&cp)-m_; }
            p_cp = mccPR(&cp) + 0 + mccM(&cp) * 0;
            if (mccM(&cp) == 1 && mccN(&cp) == 1) { I_1cp = J_1cp = 0;}
            else { I_1cp = 1; J_1cp=mccM(&cp)-m_; }
            p_1cp = mccPR(&cp) + 0 + mccM(&cp) * 0;
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_bp = J_bp = 0;}
            else { I_bp = 1; J_bp=mccM(&bp)-m_; }
            p_bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_cp += J_cp, p_1cp += J_1cp, p_bp += J_bp)
               {
                  for (i_=0; i_<m_; ++i_, p_cp+=I_cp, p_1cp+=I_1cp, p_bp+=I_bp, p_RM14_+=I_RM14_)
                  {
                     *p_cp = (*p_1cp + (*p_bp * (double) *p_RM14_));
                  }
               }
            }
         }
         
         /* b(2,i1,:) = b(2,i1,:) + a(2,i1,:) .* permute(X(1,i1)' * comp,[3 1 2]); */
         mccFindIndex(&IM10_, &i1);
         mccFindIndex(&IM9_, &i1);
         mccFindIndex(&IM8_, &i1);
         mccFindIndex(&IM7_, &IM8_);
         mccFindIndex(&IM6_, &IM7_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM14_;
            int I_RM14_=1;
            double *p_X;
            int I_X=1;
            int *p_IM6_;
            int I_IM6_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM6_) * mccN(&IM6_)));
            mccAllocateMatrix(&RM14_, m_, n_);
            mccCheckMatrixSize(&X, 1, mccGetMaxIndex(&IM6_ ,mccN(&X)));
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            I_IM6_ = (mccM(&IM6_) != 1 || mccN(&IM6_) != 1);
            p_IM6_ = mccIPR(&IM6_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM6_ += I_IM6_)
               {
                  p_X = mccPR(&X) + mccM(&X) * ((int)(*p_IM6_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM14_+=I_RM14_, p_X+=I_X)
                  {
                     *p_RM14_ = *p_X;
                  }
               }
            }
         }
         mccConjTrans(&RM13_, &RM14_);
         mccRealMatrixMultiply(&RM12_, &RM13_, &M_comp);
         Mprhs_[0] = &RM12_;
         Mprhs_[1] = &S8_;
         Mplhs_[0] = &RM11_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 251);
         mccFindIndex(&IM11_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_b;
            int I_b=1, J_b;
            double *p_1b;
            int I_1b=1, J_1b;
            double *p_a;
            int I_a=1, J_a;
            double *p_RM11_;
            int I_RM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&a), mccN(&a));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM11_), mccN(&RM11_));
            if (!mccNOTSET(&b)) m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            mccGrowMatrix(&b, m_, n_);
            mccCheckMatrixSize(&b, mccM(&b), mccN(&b));
            mccCheckMatrixSize(&a, mccM(&a), mccN(&a));
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_b = J_b = 0;}
            else { I_b = 1; J_b=mccM(&b)-m_; }
            p_b = mccPR(&b) + 0 + mccM(&b) * 0;
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_1b = J_1b = 0;}
            else { I_1b = 1; J_1b=mccM(&b)-m_; }
            p_1b = mccPR(&b) + 0 + mccM(&b) * 0;
            if (mccM(&a) == 1 && mccN(&a) == 1) { I_a = J_a = 0;}
            else { I_a = 1; J_a=mccM(&a)-m_; }
            p_a = mccPR(&a) + 0 + mccM(&a) * 0;
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_b += J_b, p_1b += J_1b, p_a += J_a)
               {
                  for (i_=0; i_<m_; ++i_, p_b+=I_b, p_1b+=I_1b, p_a+=I_a, p_RM11_+=I_RM11_)
                  {
                     *p_b = (*p_1b + (*p_a * (double) *p_RM11_));
                  }
               }
            }
         }
         /* bp(2,i2,:) = bp(2,i2,:) + ap(2,i2,:) .* permute(X(2,i2)' * comp,[3 1 2]); */
         mccFindIndex(&IM7_, &i2);
         mccFindIndex(&IM8_, &i2);
         mccFindIndex(&IM9_, &i2);
         mccFindIndex(&IM10_, &IM9_);
         mccFindIndex(&IM11_, &IM10_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM11_;
            int I_RM11_=1;
            double *p_X;
            int I_X=1;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM11_) * mccN(&IM11_)));
            mccAllocateMatrix(&RM11_, m_, n_);
            mccCheckMatrixSize(&X, 2, mccGetMaxIndex(&IM11_ ,mccN(&X)));
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM11_ += I_IM11_)
               {
                  p_X = mccPR(&X) + mccM(&X) * ((int)(*p_IM11_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM11_+=I_RM11_, p_X+=I_X)
                  {
                     *p_RM11_ = *p_X;
                  }
               }
            }
         }
         mccConjTrans(&RM12_, &RM11_);
         mccRealMatrixMultiply(&RM13_, &RM12_, &M_comp);
         Mprhs_[0] = &RM13_;
         Mprhs_[1] = &S9_;
         Mplhs_[0] = &RM14_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 252);
         mccFindIndex(&IM6_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_bp;
            int I_bp=1, J_bp;
            double *p_1bp;
            int I_1bp=1, J_1bp;
            double *p_ap;
            int I_ap=1, J_ap;
            double *p_RM14_;
            int I_RM14_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), mccN(&ap));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM14_), mccN(&RM14_));
            if (!mccNOTSET(&bp)) m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            mccGrowMatrix(&bp, m_, n_);
            mccCheckMatrixSize(&bp, mccM(&bp), mccN(&bp));
            mccCheckMatrixSize(&ap, mccM(&ap), mccN(&ap));
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_bp = J_bp = 0;}
            else { I_bp = 1; J_bp=mccM(&bp)-m_; }
            p_bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_1bp = J_1bp = 0;}
            else { I_1bp = 1; J_1bp=mccM(&bp)-m_; }
            p_1bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            if (mccM(&ap) == 1 && mccN(&ap) == 1) { I_ap = J_ap = 0;}
            else { I_ap = 1; J_ap=mccM(&ap)-m_; }
            p_ap = mccPR(&ap) + 0 + mccM(&ap) * 0;
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_bp += J_bp, p_1bp += J_1bp, p_ap += J_ap)
               {
                  for (i_=0; i_<m_; ++i_, p_bp+=I_bp, p_1bp+=I_1bp, p_ap+=I_ap, p_RM14_+=I_RM14_)
                  {
                     *p_bp = (*p_1bp + (*p_ap * (double) *p_RM14_));
                  }
               }
            }
         }
         
         /* s(1,:) = s(1,:) + squeeze( f .* sum( b(1,i1,:) .* permute(Y(1,i1)' * comp,[3 1 2]),2) + ... */
         mccFindIndex(&IM11_, &i1);
         mccFindIndex(&IM10_, &IM11_);
         mccFindIndex(&IM9_, &i1);
         mccFindIndex(&IM8_, &IM9_);
         mccFindIndex(&IM7_, &IM8_);
         mccFindIndex(&IM6_, &IM7_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM14_;
            int I_RM14_=1;
            double *p_Y;
            int I_Y=1;
            int *p_IM6_;
            int I_IM6_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM6_) * mccN(&IM6_)));
            mccAllocateMatrix(&RM14_, m_, n_);
            mccCheckMatrixSize(&Y, 1, mccGetMaxIndex(&IM6_ ,mccN(&Y)));
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            I_IM6_ = (mccM(&IM6_) != 1 || mccN(&IM6_) != 1);
            p_IM6_ = mccIPR(&IM6_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM6_ += I_IM6_)
               {
                  p_Y = mccPR(&Y) + mccM(&Y) * ((int)(*p_IM6_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM14_+=I_RM14_, p_Y+=I_Y)
                  {
                     *p_RM14_ = *p_Y;
                  }
               }
            }
         }
         mccConjTrans(&RM13_, &RM14_);
         mccRealMatrixMultiply(&RM12_, &RM13_, &M_comp);
         Mprhs_[0] = &RM12_;
         Mprhs_[1] = &S10_;
         Mplhs_[0] = &RM11_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 254);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM10_;
            int I_RM10_=1;
            double *p_b;
            int I_b=1, J_b;
            double *p_RM11_;
            int I_RM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM11_), mccN(&RM11_));
            mccAllocateMatrix(&RM10_, m_, n_);
            mccCheckMatrixSize(&b, mccM(&b), mccN(&b));
            I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
            p_RM10_ = mccPR(&RM10_);
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_b = J_b = 0;}
            else { I_b = 1; J_b=mccM(&b)-m_; }
            p_b = mccPR(&b) + 0 + mccM(&b) * 0;
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_b += J_b)
               {
                  for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_b+=I_b, p_RM11_+=I_RM11_)
                  {
                     *p_RM10_ = (*p_b * (double) *p_RM11_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM10_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM9_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 254);
         mccFindIndex(&IM5_, &i2);
         mccFindIndex(&IM4_, &IM5_);
         mccFindIndex(&IM3_, &i2);
         mccFindIndex(&IM2_, &IM3_);
         mccFindIndex(&IM0_, &IM2_);
         mccFindIndex(&IM1_, &IM0_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM8_;
            int I_RM8_=1;
            double *p_Y;
            int I_Y=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM1_) * mccN(&IM1_)));
            mccAllocateMatrix(&RM8_, m_, n_);
            mccCheckMatrixSize(&Y, 2, mccGetMaxIndex(&IM1_ ,mccN(&Y)));
            I_RM8_ = (mccM(&RM8_) != 1 || mccN(&RM8_) != 1);
            p_RM8_ = mccPR(&RM8_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_Y = mccPR(&Y) + mccM(&Y) * ((int)(*p_IM1_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM8_+=I_RM8_, p_Y+=I_Y)
                  {
                     *p_RM8_ = *p_Y;
                  }
               }
            }
         }
         mccConjTrans(&RM7_, &RM8_);
         mccRealMatrixMultiply(&RM6_, &RM7_, &M_comp);
         Mprhs_[0] = &RM6_;
         Mprhs_[1] = &S11_;
         Mplhs_[0] = &RM5_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 255);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM4_;
            int I_RM4_=1;
            double *p_bp;
            int I_bp=1, J_bp;
            double *p_RM5_;
            int I_RM5_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM5_), mccN(&RM5_));
            mccAllocateMatrix(&RM4_, m_, n_);
            mccCheckMatrixSize(&bp, mccM(&bp), mccN(&bp));
            I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
            p_RM4_ = mccPR(&RM4_);
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_bp = J_bp = 0;}
            else { I_bp = 1; J_bp=mccM(&bp)-m_; }
            p_bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            I_RM5_ = (mccM(&RM5_) != 1 || mccN(&RM5_) != 1);
            p_RM5_ = mccPR(&RM5_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_bp += J_bp)
               {
                  for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_bp+=I_bp, p_RM5_+=I_RM5_)
                  {
                     *p_RM4_ = (*p_bp * (double) *p_RM5_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM4_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM3_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 255);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM2_;
            int I_RM2_=1;
            double *p_f;
            int I_f=1;
            double *p_RM9_;
            int I_RM9_=1;
            double *p_1f;
            int I_1f=1;
            double *p_RM3_;
            int I_RM3_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM9_), mccN(&RM9_));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            mccAllocateMatrix(&RM2_, m_, n_);
            I_RM2_ = (mccM(&RM2_) != 1 || mccN(&RM2_) != 1);
            p_RM2_ = mccPR(&RM2_);
            I_f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_f = mccPR(&f);
            I_RM9_ = (mccM(&RM9_) != 1 || mccN(&RM9_) != 1);
            p_RM9_ = mccPR(&RM9_);
            I_1f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_1f = mccPR(&f);
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM2_+=I_RM2_, p_f+=I_f, p_RM9_+=I_RM9_, p_1f+=I_1f, p_RM3_+=I_RM3_)
                  {
                     *p_RM2_ = ((*p_f * (double) *p_RM9_) + ((1 - *p_1f) * (double) *p_RM3_));
                  }
               }
            }
         }
         mccLOG(&RM2_) = 0;
         mccSTRING(&RM2_) = 0;
         Mprhs_[0] = &RM2_;
         Mplhs_[0] = &RM1_;
         mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "squeeze", 255);
         mccConjTrans(&RM0_, &RM1_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_s;
            int I_s=1, J_s;
            double *p_1s;
            int I_1s=1, J_1s;
            double *p_RM0_;
            int I_RM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM0_), mccN(&RM0_));
            if (!mccNOTSET(&s)) m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            mccGrowMatrix(&s, 1, n_);
            mccCheckMatrixSize(&s, 1, mccN(&s));
            if (mccN(&s) == 1) { I_s = J_s = 0;}
            else { I_s = 1; J_s=mccM(&s)-m_; }
            p_s = mccPR(&s) + (1-1) + mccM(&s) * 0;
            if (mccN(&s) == 1) { I_1s = J_1s = 0;}
            else { I_1s = 1; J_1s=mccM(&s)-m_; }
            p_1s = mccPR(&s) + (1-1) + mccM(&s) * 0;
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_s += J_s, p_1s += J_1s)
               {
                  for (i_=0; i_<m_; ++i_, p_s+=I_s, p_1s+=I_1s, p_RM0_+=I_RM0_)
                  {
                     *p_s = (*p_1s + *p_RM0_);
                  }
               }
            }
         }
         
         /* s(2,:) = s(2,:) + squeeze( f .* sum( a(2,i1,:) .* permute(Y(1,i1)' * comp,[3 1 2]),2) + ... */
         mccFindIndex(&IM1_, &i1);
         mccFindIndex(&IM0_, &IM1_);
         mccFindIndex(&IM2_, &i1);
         mccFindIndex(&IM3_, &IM2_);
         mccFindIndex(&IM4_, &IM3_);
         mccFindIndex(&IM5_, &IM4_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM0_;
            int I_RM0_=1;
            double *p_Y;
            int I_Y=1;
            int *p_IM5_;
            int I_IM5_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM5_) * mccN(&IM5_)));
            mccAllocateMatrix(&RM0_, m_, n_);
            mccCheckMatrixSize(&Y, 1, mccGetMaxIndex(&IM5_ ,mccN(&Y)));
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            I_IM5_ = (mccM(&IM5_) != 1 || mccN(&IM5_) != 1);
            p_IM5_ = mccIPR(&IM5_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM5_ += I_IM5_)
               {
                  p_Y = mccPR(&Y) + mccM(&Y) * ((int)(*p_IM5_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_Y+=I_Y)
                  {
                     *p_RM0_ = *p_Y;
                  }
               }
            }
         }
         mccConjTrans(&RM1_, &RM0_);
         mccRealMatrixMultiply(&RM2_, &RM1_, &M_comp);
         Mprhs_[0] = &RM2_;
         Mprhs_[1] = &S12_;
         Mplhs_[0] = &RM3_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 257);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM4_;
            int I_RM4_=1;
            double *p_a;
            int I_a=1, J_a;
            double *p_RM3_;
            int I_RM3_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&a), mccN(&a));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            mccAllocateMatrix(&RM4_, m_, n_);
            mccCheckMatrixSize(&a, mccM(&a), mccN(&a));
            I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
            p_RM4_ = mccPR(&RM4_);
            if (mccM(&a) == 1 && mccN(&a) == 1) { I_a = J_a = 0;}
            else { I_a = 1; J_a=mccM(&a)-m_; }
            p_a = mccPR(&a) + 0 + mccM(&a) * 0;
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_a += J_a)
               {
                  for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_a+=I_a, p_RM3_+=I_RM3_)
                  {
                     *p_RM4_ = (*p_a * (double) *p_RM3_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM4_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM5_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 257);
         mccFindIndex(&IM6_, &i2);
         mccFindIndex(&IM7_, &IM6_);
         mccFindIndex(&IM8_, &i2);
         mccFindIndex(&IM9_, &IM8_);
         mccFindIndex(&IM10_, &IM9_);
         mccFindIndex(&IM11_, &IM10_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM6_;
            int I_RM6_=1;
            double *p_Y;
            int I_Y=1;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM11_) * mccN(&IM11_)));
            mccAllocateMatrix(&RM6_, m_, n_);
            mccCheckMatrixSize(&Y, 2, mccGetMaxIndex(&IM11_ ,mccN(&Y)));
            I_RM6_ = (mccM(&RM6_) != 1 || mccN(&RM6_) != 1);
            p_RM6_ = mccPR(&RM6_);
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM11_ += I_IM11_)
               {
                  p_Y = mccPR(&Y) + mccM(&Y) * ((int)(*p_IM11_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM6_+=I_RM6_, p_Y+=I_Y)
                  {
                     *p_RM6_ = *p_Y;
                  }
               }
            }
         }
         mccConjTrans(&RM7_, &RM6_);
         mccRealMatrixMultiply(&RM8_, &RM7_, &M_comp);
         Mprhs_[0] = &RM8_;
         Mprhs_[1] = &S13_;
         Mplhs_[0] = &RM9_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 258);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM10_;
            int I_RM10_=1;
            double *p_ap;
            int I_ap=1, J_ap;
            double *p_RM9_;
            int I_RM9_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), mccN(&ap));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM9_), mccN(&RM9_));
            mccAllocateMatrix(&RM10_, m_, n_);
            mccCheckMatrixSize(&ap, mccM(&ap), mccN(&ap));
            I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
            p_RM10_ = mccPR(&RM10_);
            if (mccM(&ap) == 1 && mccN(&ap) == 1) { I_ap = J_ap = 0;}
            else { I_ap = 1; J_ap=mccM(&ap)-m_; }
            p_ap = mccPR(&ap) + 0 + mccM(&ap) * 0;
            I_RM9_ = (mccM(&RM9_) != 1 || mccN(&RM9_) != 1);
            p_RM9_ = mccPR(&RM9_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_ap += J_ap)
               {
                  for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_ap+=I_ap, p_RM9_+=I_RM9_)
                  {
                     *p_RM10_ = (*p_ap * (double) *p_RM9_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM10_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM11_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 258);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM12_;
            int I_RM12_=1;
            double *p_f;
            int I_f=1;
            double *p_RM5_;
            int I_RM5_=1;
            double *p_1f;
            int I_1f=1;
            double *p_RM11_;
            int I_RM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM5_), mccN(&RM5_));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM11_), mccN(&RM11_));
            mccAllocateMatrix(&RM12_, m_, n_);
            I_RM12_ = (mccM(&RM12_) != 1 || mccN(&RM12_) != 1);
            p_RM12_ = mccPR(&RM12_);
            I_f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_f = mccPR(&f);
            I_RM5_ = (mccM(&RM5_) != 1 || mccN(&RM5_) != 1);
            p_RM5_ = mccPR(&RM5_);
            I_1f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_1f = mccPR(&f);
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM12_+=I_RM12_, p_f+=I_f, p_RM5_+=I_RM5_, p_1f+=I_1f, p_RM11_+=I_RM11_)
                  {
                     *p_RM12_ = ((*p_f * (double) *p_RM5_) + ((1 - *p_1f) * (double) *p_RM11_));
                  }
               }
            }
         }
         mccLOG(&RM12_) = 0;
         mccSTRING(&RM12_) = 0;
         Mprhs_[0] = &RM12_;
         Mplhs_[0] = &RM13_;
         mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "squeeze", 258);
         mccConjTrans(&RM14_, &RM13_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_s;
            int I_s=1, J_s;
            double *p_1s;
            int I_1s=1, J_1s;
            double *p_RM14_;
            int I_RM14_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM14_), mccN(&RM14_));
            if (!mccNOTSET(&s)) m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            mccGrowMatrix(&s, 2, n_);
            mccCheckMatrixSize(&s, 2, mccN(&s));
            if (mccN(&s) == 1) { I_s = J_s = 0;}
            else { I_s = 1; J_s=mccM(&s)-m_; }
            p_s = mccPR(&s) + (2-1) + mccM(&s) * 0;
            if (mccN(&s) == 1) { I_1s = J_1s = 0;}
            else { I_1s = 1; J_1s=mccM(&s)-m_; }
            p_1s = mccPR(&s) + (2-1) + mccM(&s) * 0;
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_s += J_s, p_1s += J_1s)
               {
                  for (i_=0; i_<m_; ++i_, p_s+=I_s, p_1s+=I_1s, p_RM14_+=I_RM14_)
                  {
                     *p_s = (*p_1s + *p_RM14_);
                  }
               }
            }
         }
         
         /* F(1,i1)  = X(1,i1) .* F(2,i1)  + Y(1,i1); */
         mccFindIndex(&IM10_, &i1);
         mccFindIndex(&IM9_, &i1);
         mccFindIndex(&IM8_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM14_;
            int I_RM14_=1;
            double *p_X;
            int I_X=1;
            int *p_IM10_;
            int I_IM10_=1;
            double *p_F;
            int I_F=1;
            int *p_IM9_;
            int I_IM9_=1;
            double *p_Y;
            int I_Y=1;
            int *p_IM8_;
            int I_IM8_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM10_) * mccN(&IM10_)));
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM9_) * mccN(&IM9_)));
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM8_) * mccN(&IM8_)));
            mccAllocateMatrix(&RM14_, m_, n_);
            mccCheckMatrixSize(&X, 1, mccGetMaxIndex(&IM10_ ,mccN(&X)));
            mccCheckMatrixSize(&F, 2, mccGetMaxIndex(&IM9_ ,mccN(&F)));
            mccCheckMatrixSize(&Y, 1, mccGetMaxIndex(&IM8_ ,mccN(&Y)));
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            I_IM10_ = (mccM(&IM10_) != 1 || mccN(&IM10_) != 1);
            p_IM10_ = mccIPR(&IM10_);
            I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
            p_IM9_ = mccIPR(&IM9_);
            I_IM8_ = (mccM(&IM8_) != 1 || mccN(&IM8_) != 1);
            p_IM8_ = mccIPR(&IM8_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM10_ += I_IM10_, p_IM9_ += I_IM9_, p_IM8_ += I_IM8_)
               {
                  p_X = mccPR(&X) + mccM(&X) * ((int)(*p_IM10_ - .5)) + (1-1);
                  p_F = mccPR(&F) + mccM(&F) * ((int)(*p_IM9_ - .5)) + (2-1);
                  p_Y = mccPR(&Y) + mccM(&Y) * ((int)(*p_IM8_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM14_+=I_RM14_, p_X+=I_X, p_F+=I_F, p_Y+=I_Y)
                  {
                     *p_RM14_ = ((*p_X * (double) *p_F) + *p_Y);
                  }
               }
            }
         }
         mccFindIndex(&IM11_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_F;
            int I_F=1;
            int *p_IM11_;
            int I_IM11_=1;
            double *p_RM14_;
            int I_RM14_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM11_) * mccN(&IM11_)));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM14_), mccN(&RM14_));
            mccGrowMatrix(&F, 1, mccGetMaxIndex(&IM11_ ,mccN(&F)));
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM11_ += I_IM11_)
               {
                  p_F = mccPR(&F) + mccM(&F) * ((int)(*p_IM11_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_F+=I_F, p_RM14_+=I_RM14_)
                  {
                     *p_F = *p_RM14_;
                  }
               }
            }
         }
         /* FP(1,i2) = X(2,i2) .* FP(2,i2) + Y(2,i2); */
         mccFindIndex(&IM9_, &i2);
         mccFindIndex(&IM10_, &i2);
         mccFindIndex(&IM11_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM14_;
            int I_RM14_=1;
            double *p_X;
            int I_X=1;
            int *p_IM9_;
            int I_IM9_=1;
            double *p_FP;
            int I_FP=1;
            int *p_IM10_;
            int I_IM10_=1;
            double *p_Y;
            int I_Y=1;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM9_) * mccN(&IM9_)));
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM10_) * mccN(&IM10_)));
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM11_) * mccN(&IM11_)));
            mccAllocateMatrix(&RM14_, m_, n_);
            mccCheckMatrixSize(&X, 2, mccGetMaxIndex(&IM9_ ,mccN(&X)));
            mccCheckMatrixSize(&FP, 2, mccGetMaxIndex(&IM10_ ,mccN(&FP)));
            mccCheckMatrixSize(&Y, 2, mccGetMaxIndex(&IM11_ ,mccN(&Y)));
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
            p_IM9_ = mccIPR(&IM9_);
            I_IM10_ = (mccM(&IM10_) != 1 || mccN(&IM10_) != 1);
            p_IM10_ = mccIPR(&IM10_);
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM9_ += I_IM9_, p_IM10_ += I_IM10_, p_IM11_ += I_IM11_)
               {
                  p_X = mccPR(&X) + mccM(&X) * ((int)(*p_IM9_ - .5)) + (2-1);
                  p_FP = mccPR(&FP) + mccM(&FP) * ((int)(*p_IM10_ - .5)) + (2-1);
                  p_Y = mccPR(&Y) + mccM(&Y) * ((int)(*p_IM11_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM14_+=I_RM14_, p_X+=I_X, p_FP+=I_FP, p_Y+=I_Y)
                  {
                     *p_RM14_ = ((*p_X * (double) *p_FP) + *p_Y);
                  }
               }
            }
         }
         mccFindIndex(&IM8_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_FP;
            int I_FP=1;
            int *p_IM8_;
            int I_IM8_=1;
            double *p_RM14_;
            int I_RM14_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM8_) * mccN(&IM8_)));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM14_), mccN(&RM14_));
            mccGrowMatrix(&FP, 1, mccGetMaxIndex(&IM8_ ,mccN(&FP)));
            I_IM8_ = (mccM(&IM8_) != 1 || mccN(&IM8_) != 1);
            p_IM8_ = mccIPR(&IM8_);
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM8_ += I_IM8_)
               {
                  p_FP = mccPR(&FP) + mccM(&FP) * ((int)(*p_IM8_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_FP+=I_FP, p_RM14_+=I_RM14_)
                  {
                     *p_FP = *p_RM14_;
                  }
               }
            }
         }
         
         /* a(2,i1,:)   = zeros(1,length(i1),M); */
         I0_ = mccGetLength(&i1);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM10_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 263);
         mccFindIndex(&IM11_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_a;
            int I_a=1, J_a;
            int *p_IM10_;
            int I_IM10_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM10_), mccN(&IM10_));
            if (!mccNOTSET(&a)) m_ = mcmCalcResultSize(m_, &n_, mccM(&a), mccN(&a));
            mccGrowMatrix(&a, m_, n_);
            if (mccM(&a) == 1 && mccN(&a) == 1) { I_a = J_a = 0;}
            else { I_a = 1; J_a=mccM(&a)-m_; }
            p_a = mccPR(&a) + 0 + mccM(&a) * 0;
            I_IM10_ = (mccM(&IM10_) != 1 || mccN(&IM10_) != 1);
            p_IM10_ = mccIPR(&IM10_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_a += J_a)
               {
                  for (i_=0; i_<m_; ++i_, p_a+=I_a, p_IM10_+=I_IM10_)
                  {
                     *p_a = ((int)*p_IM10_);
                  }
               }
            }
         }
         /* ap(2,i2,:)  = zeros(1,length(i2),M); */
         I0_ = mccGetLength(&i2);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM11_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 264);
         mccFindIndex(&IM10_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_ap;
            int I_ap=1, J_ap;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM11_), mccN(&IM11_));
            if (!mccNOTSET(&ap)) m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), mccN(&ap));
            mccGrowMatrix(&ap, m_, n_);
            if (mccM(&ap) == 1 && mccN(&ap) == 1) { I_ap = J_ap = 0;}
            else { I_ap = 1; J_ap=mccM(&ap)-m_; }
            p_ap = mccPR(&ap) + 0 + mccM(&ap) * 0;
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_ap += J_ap)
               {
                  for (i_=0; i_<m_; ++i_, p_ap+=I_ap, p_IM11_+=I_IM11_)
                  {
                     *p_ap = ((int)*p_IM11_);
                  }
               }
            }
         }
         /* b(1,i1,:)   = zeros(1,length(i1),M); */
         I0_ = mccGetLength(&i1);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM10_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 265);
         mccFindIndex(&IM11_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_b;
            int I_b=1, J_b;
            int *p_IM10_;
            int I_IM10_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM10_), mccN(&IM10_));
            if (!mccNOTSET(&b)) m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            mccGrowMatrix(&b, m_, n_);
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_b = J_b = 0;}
            else { I_b = 1; J_b=mccM(&b)-m_; }
            p_b = mccPR(&b) + 0 + mccM(&b) * 0;
            I_IM10_ = (mccM(&IM10_) != 1 || mccN(&IM10_) != 1);
            p_IM10_ = mccIPR(&IM10_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_b += J_b)
               {
                  for (i_=0; i_<m_; ++i_, p_b+=I_b, p_IM10_+=I_IM10_)
                  {
                     *p_b = ((int)*p_IM10_);
                  }
               }
            }
         }
         /* bp(1,i2,:)  = zeros(1,length(i2),M); */
         I0_ = mccGetLength(&i2);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM11_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 266);
         mccFindIndex(&IM10_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_bp;
            int I_bp=1, J_bp;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM11_), mccN(&IM11_));
            if (!mccNOTSET(&bp)) m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            mccGrowMatrix(&bp, m_, n_);
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_bp = J_bp = 0;}
            else { I_bp = 1; J_bp=mccM(&bp)-m_; }
            p_bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_bp += J_bp)
               {
                  for (i_=0; i_<m_; ++i_, p_bp+=I_bp, p_IM11_+=I_IM11_)
                  {
                     *p_bp = ((int)*p_IM11_);
                  }
               }
            }
         }
         /* end */
      }
      
      /* % d - cas w ==0, v ==0 et u~=0 */
      /* i1 = find(( w(1,:) == 0 ) & (v(1,:) == 0) & (u(1,:) ~= 0)); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         unsigned short *p_BM0_;
         int I_BM0_=1;
         double *p_w;
         int I_w=1, J_w;
         double *p_v;
         int I_v=1, J_v;
         double *p_u;
         int I_u=1, J_u;
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&w));
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&v));
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&u));
         mccAllocateMatrix(&BM0_, m_, n_);
         mccCheckMatrixSize(&w, 1, mccN(&w));
         mccCheckMatrixSize(&v, 1, mccN(&v));
         mccCheckMatrixSize(&u, 1, mccN(&u));
         I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
         p_BM0_ = mccSPR(&BM0_);
         if (mccN(&w) == 1) { I_w = J_w = 0;}
         else { I_w = 1; J_w=mccM(&w)-m_; }
         p_w = mccPR(&w) + (1-1) + mccM(&w) * 0;
         if (mccN(&v) == 1) { I_v = J_v = 0;}
         else { I_v = 1; J_v=mccM(&v)-m_; }
         p_v = mccPR(&v) + (1-1) + mccM(&v) * 0;
         if (mccN(&u) == 1) { I_u = J_u = 0;}
         else { I_u = 1; J_u=mccM(&u)-m_; }
         p_u = mccPR(&u) + (1-1) + mccM(&u) * 0;
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_w += J_w, p_v += J_v, p_u += J_u)
            {
               for (i_=0; i_<m_; ++i_, p_BM0_+=I_BM0_, p_w+=I_w, p_v+=I_v, p_u+=I_u)
               {
                  *p_BM0_ = (!!(!!( (*p_w == 0) && !mccREL_NAN(*p_w) ) && !!( (*p_v == 0) && !mccREL_NAN(*p_v) )) && !!( (*p_u != 0) || mccREL_NAN(*p_u) ));
               }
            }
         }
      }
      mccFind(&i1, &BM0_);
      /* i2 = find(( w(2,:) ==0) & (v(2,:) == 0) & (u(2,:) ~= 0)); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         unsigned short *p_BM0_;
         int I_BM0_=1;
         double *p_w;
         int I_w=1, J_w;
         double *p_v;
         int I_v=1, J_v;
         double *p_u;
         int I_u=1, J_u;
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&w));
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&v));
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&u));
         mccAllocateMatrix(&BM0_, m_, n_);
         mccCheckMatrixSize(&w, 2, mccN(&w));
         mccCheckMatrixSize(&v, 2, mccN(&v));
         mccCheckMatrixSize(&u, 2, mccN(&u));
         I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
         p_BM0_ = mccSPR(&BM0_);
         if (mccN(&w) == 1) { I_w = J_w = 0;}
         else { I_w = 1; J_w=mccM(&w)-m_; }
         p_w = mccPR(&w) + (2-1) + mccM(&w) * 0;
         if (mccN(&v) == 1) { I_v = J_v = 0;}
         else { I_v = 1; J_v=mccM(&v)-m_; }
         p_v = mccPR(&v) + (2-1) + mccM(&v) * 0;
         if (mccN(&u) == 1) { I_u = J_u = 0;}
         else { I_u = 1; J_u=mccM(&u)-m_; }
         p_u = mccPR(&u) + (2-1) + mccM(&u) * 0;
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_w += J_w, p_v += J_v, p_u += J_u)
            {
               for (i_=0; i_<m_; ++i_, p_BM0_+=I_BM0_, p_w+=I_w, p_v+=I_v, p_u+=I_u)
               {
                  *p_BM0_ = (!!(!!( (*p_w == 0) && !mccREL_NAN(*p_w) ) && !!( (*p_v == 0) && !mccREL_NAN(*p_v) )) && !!( (*p_u != 0) || mccREL_NAN(*p_u) ));
               }
            }
         }
      }
      mccFind(&i2, &BM0_);
      /* if ~ (isempty(i1) & isempty(i2)) */
      B0_ = mccIsEmpty(&i1);
      B1_ = mccIsEmpty(&i2);
      B2_ = (!(!!B0_ && !!B1_));
      if ((double)B2_)
      {
         /* s(1,:) = s(1,:) + squeeze( f .* sum( c(1,i1,:) .* permute(Z(1,i1)' * comp,[3 1 2]),2) + ... */
         mccFindIndex(&IM11_, &i1);
         mccFindIndex(&IM10_, &IM11_);
         mccFindIndex(&IM9_, &i1);
         mccFindIndex(&IM8_, &IM9_);
         mccFindIndex(&IM7_, &IM8_);
         mccFindIndex(&IM6_, &IM7_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM14_;
            int I_RM14_=1;
            double *p_Z;
            int I_Z=1;
            int *p_IM6_;
            int I_IM6_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM6_) * mccN(&IM6_)));
            mccAllocateMatrix(&RM14_, m_, n_);
            mccCheckMatrixSize(&Z, 1, mccGetMaxIndex(&IM6_ ,mccN(&Z)));
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            I_IM6_ = (mccM(&IM6_) != 1 || mccN(&IM6_) != 1);
            p_IM6_ = mccIPR(&IM6_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM6_ += I_IM6_)
               {
                  p_Z = mccPR(&Z) + mccM(&Z) * ((int)(*p_IM6_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM14_+=I_RM14_, p_Z+=I_Z)
                  {
                     *p_RM14_ = *p_Z;
                  }
               }
            }
         }
         mccConjTrans(&RM13_, &RM14_);
         mccRealMatrixMultiply(&RM12_, &RM13_, &M_comp);
         Mprhs_[0] = &RM12_;
         Mprhs_[1] = &S14_;
         Mplhs_[0] = &RM11_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 273);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM10_;
            int I_RM10_=1;
            double *p_c;
            int I_c=1, J_c;
            double *p_RM11_;
            int I_RM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM11_), mccN(&RM11_));
            mccAllocateMatrix(&RM10_, m_, n_);
            mccCheckMatrixSize(&c, mccM(&c), mccN(&c));
            I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
            p_RM10_ = mccPR(&RM10_);
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_c = J_c = 0;}
            else { I_c = 1; J_c=mccM(&c)-m_; }
            p_c = mccPR(&c) + 0 + mccM(&c) * 0;
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_c += J_c)
               {
                  for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_c+=I_c, p_RM11_+=I_RM11_)
                  {
                     *p_RM10_ = (*p_c * (double) *p_RM11_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM10_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM9_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 273);
         mccFindIndex(&IM5_, &i2);
         mccFindIndex(&IM4_, &IM5_);
         mccFindIndex(&IM3_, &i2);
         mccFindIndex(&IM2_, &IM3_);
         mccFindIndex(&IM0_, &IM2_);
         mccFindIndex(&IM1_, &IM0_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM8_;
            int I_RM8_=1;
            double *p_Z;
            int I_Z=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM1_) * mccN(&IM1_)));
            mccAllocateMatrix(&RM8_, m_, n_);
            mccCheckMatrixSize(&Z, 2, mccGetMaxIndex(&IM1_ ,mccN(&Z)));
            I_RM8_ = (mccM(&RM8_) != 1 || mccN(&RM8_) != 1);
            p_RM8_ = mccPR(&RM8_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_Z = mccPR(&Z) + mccM(&Z) * ((int)(*p_IM1_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM8_+=I_RM8_, p_Z+=I_Z)
                  {
                     *p_RM8_ = *p_Z;
                  }
               }
            }
         }
         mccConjTrans(&RM7_, &RM8_);
         mccRealMatrixMultiply(&RM6_, &RM7_, &M_comp);
         Mprhs_[0] = &RM6_;
         Mprhs_[1] = &S15_;
         Mplhs_[0] = &RM5_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 274);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM4_;
            int I_RM4_=1;
            double *p_cp;
            int I_cp=1, J_cp;
            double *p_RM5_;
            int I_RM5_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&cp), mccN(&cp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM5_), mccN(&RM5_));
            mccAllocateMatrix(&RM4_, m_, n_);
            mccCheckMatrixSize(&cp, mccM(&cp), mccN(&cp));
            I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
            p_RM4_ = mccPR(&RM4_);
            if (mccM(&cp) == 1 && mccN(&cp) == 1) { I_cp = J_cp = 0;}
            else { I_cp = 1; J_cp=mccM(&cp)-m_; }
            p_cp = mccPR(&cp) + 0 + mccM(&cp) * 0;
            I_RM5_ = (mccM(&RM5_) != 1 || mccN(&RM5_) != 1);
            p_RM5_ = mccPR(&RM5_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_cp += J_cp)
               {
                  for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_cp+=I_cp, p_RM5_+=I_RM5_)
                  {
                     *p_RM4_ = (*p_cp * (double) *p_RM5_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM4_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM3_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 274);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM2_;
            int I_RM2_=1;
            double *p_f;
            int I_f=1;
            double *p_RM9_;
            int I_RM9_=1;
            double *p_1f;
            int I_1f=1;
            double *p_RM3_;
            int I_RM3_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM9_), mccN(&RM9_));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            mccAllocateMatrix(&RM2_, m_, n_);
            I_RM2_ = (mccM(&RM2_) != 1 || mccN(&RM2_) != 1);
            p_RM2_ = mccPR(&RM2_);
            I_f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_f = mccPR(&f);
            I_RM9_ = (mccM(&RM9_) != 1 || mccN(&RM9_) != 1);
            p_RM9_ = mccPR(&RM9_);
            I_1f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_1f = mccPR(&f);
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM2_+=I_RM2_, p_f+=I_f, p_RM9_+=I_RM9_, p_1f+=I_1f, p_RM3_+=I_RM3_)
                  {
                     *p_RM2_ = ((*p_f * (double) *p_RM9_) + ((1 - *p_1f) * (double) *p_RM3_));
                  }
               }
            }
         }
         mccLOG(&RM2_) = 0;
         mccSTRING(&RM2_) = 0;
         Mprhs_[0] = &RM2_;
         Mplhs_[0] = &RM1_;
         mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "squeeze", 274);
         mccConjTrans(&RM0_, &RM1_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_s;
            int I_s=1, J_s;
            double *p_1s;
            int I_1s=1, J_1s;
            double *p_RM0_;
            int I_RM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM0_), mccN(&RM0_));
            if (!mccNOTSET(&s)) m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            mccGrowMatrix(&s, 1, n_);
            mccCheckMatrixSize(&s, 1, mccN(&s));
            if (mccN(&s) == 1) { I_s = J_s = 0;}
            else { I_s = 1; J_s=mccM(&s)-m_; }
            p_s = mccPR(&s) + (1-1) + mccM(&s) * 0;
            if (mccN(&s) == 1) { I_1s = J_1s = 0;}
            else { I_1s = 1; J_1s=mccM(&s)-m_; }
            p_1s = mccPR(&s) + (1-1) + mccM(&s) * 0;
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_s += J_s, p_1s += J_1s)
               {
                  for (i_=0; i_<m_; ++i_, p_s+=I_s, p_1s+=I_1s, p_RM0_+=I_RM0_)
                  {
                     *p_s = (*p_1s + *p_RM0_);
                  }
               }
            }
         }
         
         /* s(2,:) = s(2,:) + squeeze( f .* sum( b(2,i1,:) .* permute(Z(1,i1)' * comp,[3 1 2]),2) + ... */
         mccFindIndex(&IM1_, &i1);
         mccFindIndex(&IM0_, &IM1_);
         mccFindIndex(&IM2_, &i1);
         mccFindIndex(&IM3_, &IM2_);
         mccFindIndex(&IM4_, &IM3_);
         mccFindIndex(&IM5_, &IM4_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM0_;
            int I_RM0_=1;
            double *p_Z;
            int I_Z=1;
            int *p_IM5_;
            int I_IM5_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM5_) * mccN(&IM5_)));
            mccAllocateMatrix(&RM0_, m_, n_);
            mccCheckMatrixSize(&Z, 1, mccGetMaxIndex(&IM5_ ,mccN(&Z)));
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            I_IM5_ = (mccM(&IM5_) != 1 || mccN(&IM5_) != 1);
            p_IM5_ = mccIPR(&IM5_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM5_ += I_IM5_)
               {
                  p_Z = mccPR(&Z) + mccM(&Z) * ((int)(*p_IM5_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_Z+=I_Z)
                  {
                     *p_RM0_ = *p_Z;
                  }
               }
            }
         }
         mccConjTrans(&RM1_, &RM0_);
         mccRealMatrixMultiply(&RM2_, &RM1_, &M_comp);
         Mprhs_[0] = &RM2_;
         Mprhs_[1] = &S16_;
         Mplhs_[0] = &RM3_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 276);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM4_;
            int I_RM4_=1;
            double *p_b;
            int I_b=1, J_b;
            double *p_RM3_;
            int I_RM3_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            mccAllocateMatrix(&RM4_, m_, n_);
            mccCheckMatrixSize(&b, mccM(&b), mccN(&b));
            I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
            p_RM4_ = mccPR(&RM4_);
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_b = J_b = 0;}
            else { I_b = 1; J_b=mccM(&b)-m_; }
            p_b = mccPR(&b) + 0 + mccM(&b) * 0;
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_b += J_b)
               {
                  for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_b+=I_b, p_RM3_+=I_RM3_)
                  {
                     *p_RM4_ = (*p_b * (double) *p_RM3_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM4_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM5_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 276);
         mccFindIndex(&IM6_, &i2);
         mccFindIndex(&IM7_, &IM6_);
         mccFindIndex(&IM8_, &i2);
         mccFindIndex(&IM9_, &IM8_);
         mccFindIndex(&IM10_, &IM9_);
         mccFindIndex(&IM11_, &IM10_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM6_;
            int I_RM6_=1;
            double *p_Z;
            int I_Z=1;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM11_) * mccN(&IM11_)));
            mccAllocateMatrix(&RM6_, m_, n_);
            mccCheckMatrixSize(&Z, 2, mccGetMaxIndex(&IM11_ ,mccN(&Z)));
            I_RM6_ = (mccM(&RM6_) != 1 || mccN(&RM6_) != 1);
            p_RM6_ = mccPR(&RM6_);
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM11_ += I_IM11_)
               {
                  p_Z = mccPR(&Z) + mccM(&Z) * ((int)(*p_IM11_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM6_+=I_RM6_, p_Z+=I_Z)
                  {
                     *p_RM6_ = *p_Z;
                  }
               }
            }
         }
         mccConjTrans(&RM7_, &RM6_);
         mccRealMatrixMultiply(&RM8_, &RM7_, &M_comp);
         Mprhs_[0] = &RM8_;
         Mprhs_[1] = &S17_;
         Mplhs_[0] = &RM9_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 277);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM10_;
            int I_RM10_=1;
            double *p_bp;
            int I_bp=1, J_bp;
            double *p_RM9_;
            int I_RM9_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM9_), mccN(&RM9_));
            mccAllocateMatrix(&RM10_, m_, n_);
            mccCheckMatrixSize(&bp, mccM(&bp), mccN(&bp));
            I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
            p_RM10_ = mccPR(&RM10_);
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_bp = J_bp = 0;}
            else { I_bp = 1; J_bp=mccM(&bp)-m_; }
            p_bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            I_RM9_ = (mccM(&RM9_) != 1 || mccN(&RM9_) != 1);
            p_RM9_ = mccPR(&RM9_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_bp += J_bp)
               {
                  for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_bp+=I_bp, p_RM9_+=I_RM9_)
                  {
                     *p_RM10_ = (*p_bp * (double) *p_RM9_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM10_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM11_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 277);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM12_;
            int I_RM12_=1;
            double *p_f;
            int I_f=1;
            double *p_RM5_;
            int I_RM5_=1;
            double *p_1f;
            int I_1f=1;
            double *p_RM11_;
            int I_RM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM5_), mccN(&RM5_));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM11_), mccN(&RM11_));
            mccAllocateMatrix(&RM12_, m_, n_);
            I_RM12_ = (mccM(&RM12_) != 1 || mccN(&RM12_) != 1);
            p_RM12_ = mccPR(&RM12_);
            I_f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_f = mccPR(&f);
            I_RM5_ = (mccM(&RM5_) != 1 || mccN(&RM5_) != 1);
            p_RM5_ = mccPR(&RM5_);
            I_1f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_1f = mccPR(&f);
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM12_+=I_RM12_, p_f+=I_f, p_RM5_+=I_RM5_, p_1f+=I_1f, p_RM11_+=I_RM11_)
                  {
                     *p_RM12_ = ((*p_f * (double) *p_RM5_) + ((1 - *p_1f) * (double) *p_RM11_));
                  }
               }
            }
         }
         mccLOG(&RM12_) = 0;
         mccSTRING(&RM12_) = 0;
         Mprhs_[0] = &RM12_;
         Mplhs_[0] = &RM13_;
         mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "squeeze", 277);
         mccConjTrans(&RM14_, &RM13_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_s;
            int I_s=1, J_s;
            double *p_1s;
            int I_1s=1, J_1s;
            double *p_RM14_;
            int I_RM14_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM14_), mccN(&RM14_));
            if (!mccNOTSET(&s)) m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            mccGrowMatrix(&s, 2, n_);
            mccCheckMatrixSize(&s, 2, mccN(&s));
            if (mccN(&s) == 1) { I_s = J_s = 0;}
            else { I_s = 1; J_s=mccM(&s)-m_; }
            p_s = mccPR(&s) + (2-1) + mccM(&s) * 0;
            if (mccN(&s) == 1) { I_1s = J_1s = 0;}
            else { I_1s = 1; J_1s=mccM(&s)-m_; }
            p_1s = mccPR(&s) + (2-1) + mccM(&s) * 0;
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_s += J_s, p_1s += J_1s)
               {
                  for (i_=0; i_<m_; ++i_, p_s+=I_s, p_1s+=I_1s, p_RM14_+=I_RM14_)
                  {
                     *p_s = (*p_1s + *p_RM14_);
                  }
               }
            }
         }
         
         /* s(3,:) = s(3,:) + squeeze( f .* sum( a(3,i1,:) .* permute(Z(1,i1)' * comp,[3 1 2]),2) + ... */
         mccFindIndex(&IM11_, &i1);
         mccFindIndex(&IM10_, &IM11_);
         mccFindIndex(&IM9_, &i1);
         mccFindIndex(&IM8_, &IM9_);
         mccFindIndex(&IM7_, &IM8_);
         mccFindIndex(&IM6_, &IM7_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM14_;
            int I_RM14_=1;
            double *p_Z;
            int I_Z=1;
            int *p_IM6_;
            int I_IM6_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM6_) * mccN(&IM6_)));
            mccAllocateMatrix(&RM14_, m_, n_);
            mccCheckMatrixSize(&Z, 1, mccGetMaxIndex(&IM6_ ,mccN(&Z)));
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            I_IM6_ = (mccM(&IM6_) != 1 || mccN(&IM6_) != 1);
            p_IM6_ = mccIPR(&IM6_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM6_ += I_IM6_)
               {
                  p_Z = mccPR(&Z) + mccM(&Z) * ((int)(*p_IM6_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM14_+=I_RM14_, p_Z+=I_Z)
                  {
                     *p_RM14_ = *p_Z;
                  }
               }
            }
         }
         mccConjTrans(&RM13_, &RM14_);
         mccRealMatrixMultiply(&RM12_, &RM13_, &M_comp);
         Mprhs_[0] = &RM12_;
         Mprhs_[1] = &S18_;
         Mplhs_[0] = &RM11_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 279);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM10_;
            int I_RM10_=1;
            double *p_a;
            int I_a=1, J_a;
            double *p_RM11_;
            int I_RM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&a), mccN(&a));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM11_), mccN(&RM11_));
            mccAllocateMatrix(&RM10_, m_, n_);
            mccCheckMatrixSize(&a, mccM(&a), mccN(&a));
            I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
            p_RM10_ = mccPR(&RM10_);
            if (mccM(&a) == 1 && mccN(&a) == 1) { I_a = J_a = 0;}
            else { I_a = 1; J_a=mccM(&a)-m_; }
            p_a = mccPR(&a) + 0 + mccM(&a) * 0;
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_a += J_a)
               {
                  for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_a+=I_a, p_RM11_+=I_RM11_)
                  {
                     *p_RM10_ = (*p_a * (double) *p_RM11_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM10_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM9_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 279);
         mccFindIndex(&IM5_, &i2);
         mccFindIndex(&IM4_, &IM5_);
         mccFindIndex(&IM3_, &i2);
         mccFindIndex(&IM2_, &IM3_);
         mccFindIndex(&IM0_, &IM2_);
         mccFindIndex(&IM1_, &IM0_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM8_;
            int I_RM8_=1;
            double *p_Z;
            int I_Z=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM1_) * mccN(&IM1_)));
            mccAllocateMatrix(&RM8_, m_, n_);
            mccCheckMatrixSize(&Z, 2, mccGetMaxIndex(&IM1_ ,mccN(&Z)));
            I_RM8_ = (mccM(&RM8_) != 1 || mccN(&RM8_) != 1);
            p_RM8_ = mccPR(&RM8_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_Z = mccPR(&Z) + mccM(&Z) * ((int)(*p_IM1_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM8_+=I_RM8_, p_Z+=I_Z)
                  {
                     *p_RM8_ = *p_Z;
                  }
               }
            }
         }
         mccConjTrans(&RM7_, &RM8_);
         mccRealMatrixMultiply(&RM6_, &RM7_, &M_comp);
         Mprhs_[0] = &RM6_;
         Mprhs_[1] = &S19_;
         Mplhs_[0] = &RM5_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 280);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM4_;
            int I_RM4_=1;
            double *p_ap;
            int I_ap=1, J_ap;
            double *p_RM5_;
            int I_RM5_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), mccN(&ap));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM5_), mccN(&RM5_));
            mccAllocateMatrix(&RM4_, m_, n_);
            mccCheckMatrixSize(&ap, mccM(&ap), mccN(&ap));
            I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
            p_RM4_ = mccPR(&RM4_);
            if (mccM(&ap) == 1 && mccN(&ap) == 1) { I_ap = J_ap = 0;}
            else { I_ap = 1; J_ap=mccM(&ap)-m_; }
            p_ap = mccPR(&ap) + 0 + mccM(&ap) * 0;
            I_RM5_ = (mccM(&RM5_) != 1 || mccN(&RM5_) != 1);
            p_RM5_ = mccPR(&RM5_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_ap += J_ap)
               {
                  for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_ap+=I_ap, p_RM5_+=I_RM5_)
                  {
                     *p_RM4_ = (*p_ap * (double) *p_RM5_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM4_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM3_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 280);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM2_;
            int I_RM2_=1;
            double *p_f;
            int I_f=1;
            double *p_RM9_;
            int I_RM9_=1;
            double *p_1f;
            int I_1f=1;
            double *p_RM3_;
            int I_RM3_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM9_), mccN(&RM9_));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            mccAllocateMatrix(&RM2_, m_, n_);
            I_RM2_ = (mccM(&RM2_) != 1 || mccN(&RM2_) != 1);
            p_RM2_ = mccPR(&RM2_);
            I_f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_f = mccPR(&f);
            I_RM9_ = (mccM(&RM9_) != 1 || mccN(&RM9_) != 1);
            p_RM9_ = mccPR(&RM9_);
            I_1f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_1f = mccPR(&f);
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM2_+=I_RM2_, p_f+=I_f, p_RM9_+=I_RM9_, p_1f+=I_1f, p_RM3_+=I_RM3_)
                  {
                     *p_RM2_ = ((*p_f * (double) *p_RM9_) + ((1 - *p_1f) * (double) *p_RM3_));
                  }
               }
            }
         }
         mccLOG(&RM2_) = 0;
         mccSTRING(&RM2_) = 0;
         Mprhs_[0] = &RM2_;
         Mplhs_[0] = &RM1_;
         mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "squeeze", 280);
         mccConjTrans(&RM0_, &RM1_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_s;
            int I_s=1, J_s;
            double *p_1s;
            int I_1s=1, J_1s;
            double *p_RM0_;
            int I_RM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM0_), mccN(&RM0_));
            if (!mccNOTSET(&s)) m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            mccGrowMatrix(&s, 3, n_);
            mccCheckMatrixSize(&s, 3, mccN(&s));
            if (mccN(&s) == 1) { I_s = J_s = 0;}
            else { I_s = 1; J_s=mccM(&s)-m_; }
            p_s = mccPR(&s) + (3-1) + mccM(&s) * 0;
            if (mccN(&s) == 1) { I_1s = J_1s = 0;}
            else { I_1s = 1; J_1s=mccM(&s)-m_; }
            p_1s = mccPR(&s) + (3-1) + mccM(&s) * 0;
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_s += J_s, p_1s += J_1s)
               {
                  for (i_=0; i_<m_; ++i_, p_s+=I_s, p_1s+=I_1s, p_RM0_+=I_RM0_)
                  {
                     *p_s = (*p_1s + *p_RM0_);
                  }
               }
            }
         }
         
         /* F(2,i1)  = Z(1,i1);                      */
         mccFindIndex(&IM0_, &i1);
         mccFindIndex(&IM1_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_F;
            int I_F=1;
            int *p_IM1_;
            int I_IM1_=1;
            double *p_Z;
            int I_Z=1;
            int *p_IM0_;
            int I_IM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM1_) * mccN(&IM1_)));
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM0_) * mccN(&IM0_)));
            mccGrowMatrix(&F, 2, mccGetMaxIndex(&IM1_ ,mccN(&F)));
            mccCheckMatrixSize(&Z, 1, mccGetMaxIndex(&IM0_ ,mccN(&Z)));
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_, p_IM0_ += I_IM0_)
               {
                  p_F = mccPR(&F) + mccM(&F) * ((int)(*p_IM1_ - .5)) + (2-1);
                  p_Z = mccPR(&Z) + mccM(&Z) * ((int)(*p_IM0_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_F+=I_F, p_Z+=I_Z)
                  {
                     *p_F = *p_Z;
                  }
               }
            }
         }
         /* FP(2,i2) = Z(2,i2);                      */
         mccFindIndex(&IM1_, &i2);
         mccFindIndex(&IM0_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_FP;
            int I_FP=1;
            int *p_IM0_;
            int I_IM0_=1;
            double *p_Z;
            int I_Z=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM0_) * mccN(&IM0_)));
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM1_) * mccN(&IM1_)));
            mccGrowMatrix(&FP, 2, mccGetMaxIndex(&IM0_ ,mccN(&FP)));
            mccCheckMatrixSize(&Z, 2, mccGetMaxIndex(&IM1_ ,mccN(&Z)));
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM0_ += I_IM0_, p_IM1_ += I_IM1_)
               {
                  p_FP = mccPR(&FP) + mccM(&FP) * ((int)(*p_IM0_ - .5)) + (2-1);
                  p_Z = mccPR(&Z) + mccM(&Z) * ((int)(*p_IM1_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_FP+=I_FP, p_Z+=I_Z)
                  {
                     *p_FP = *p_Z;
                  }
               }
            }
         }
         
         /* a(3,i1,:)   = zeros(1,length(i1),M); */
         I0_ = mccGetLength(&i1);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM0_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 285);
         mccFindIndex(&IM1_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_a;
            int I_a=1, J_a;
            int *p_IM0_;
            int I_IM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM0_), mccN(&IM0_));
            if (!mccNOTSET(&a)) m_ = mcmCalcResultSize(m_, &n_, mccM(&a), mccN(&a));
            mccGrowMatrix(&a, m_, n_);
            if (mccM(&a) == 1 && mccN(&a) == 1) { I_a = J_a = 0;}
            else { I_a = 1; J_a=mccM(&a)-m_; }
            p_a = mccPR(&a) + 0 + mccM(&a) * 0;
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_a += J_a)
               {
                  for (i_=0; i_<m_; ++i_, p_a+=I_a, p_IM0_+=I_IM0_)
                  {
                     *p_a = ((int)*p_IM0_);
                  }
               }
            }
         }
         /* ap(3,i2,:)  = zeros(1,length(i2),M); */
         I0_ = mccGetLength(&i2);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM1_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 286);
         mccFindIndex(&IM0_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_ap;
            int I_ap=1, J_ap;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM1_), mccN(&IM1_));
            if (!mccNOTSET(&ap)) m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), mccN(&ap));
            mccGrowMatrix(&ap, m_, n_);
            if (mccM(&ap) == 1 && mccN(&ap) == 1) { I_ap = J_ap = 0;}
            else { I_ap = 1; J_ap=mccM(&ap)-m_; }
            p_ap = mccPR(&ap) + 0 + mccM(&ap) * 0;
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_ap += J_ap)
               {
                  for (i_=0; i_<m_; ++i_, p_ap+=I_ap, p_IM1_+=I_IM1_)
                  {
                     *p_ap = ((int)*p_IM1_);
                  }
               }
            }
         }
         /* b(2,i1,:)   = zeros(1,length(i1),M); */
         I0_ = mccGetLength(&i1);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM0_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 287);
         mccFindIndex(&IM1_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_b;
            int I_b=1, J_b;
            int *p_IM0_;
            int I_IM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM0_), mccN(&IM0_));
            if (!mccNOTSET(&b)) m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            mccGrowMatrix(&b, m_, n_);
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_b = J_b = 0;}
            else { I_b = 1; J_b=mccM(&b)-m_; }
            p_b = mccPR(&b) + 0 + mccM(&b) * 0;
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_b += J_b)
               {
                  for (i_=0; i_<m_; ++i_, p_b+=I_b, p_IM0_+=I_IM0_)
                  {
                     *p_b = ((int)*p_IM0_);
                  }
               }
            }
         }
         /* bp(2,i2,:)  = zeros(1,length(i2),M); */
         I0_ = mccGetLength(&i2);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM1_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 288);
         mccFindIndex(&IM0_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_bp;
            int I_bp=1, J_bp;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM1_), mccN(&IM1_));
            if (!mccNOTSET(&bp)) m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            mccGrowMatrix(&bp, m_, n_);
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_bp = J_bp = 0;}
            else { I_bp = 1; J_bp=mccM(&bp)-m_; }
            p_bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_bp += J_bp)
               {
                  for (i_=0; i_<m_; ++i_, p_bp+=I_bp, p_IM1_+=I_IM1_)
                  {
                     *p_bp = ((int)*p_IM1_);
                  }
               }
            }
         }
         /* c(1,i1,:)   = zeros(1,length(i1),M); */
         I0_ = mccGetLength(&i1);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM0_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 289);
         mccFindIndex(&IM1_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_c;
            int I_c=1, J_c;
            int *p_IM0_;
            int I_IM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM0_), mccN(&IM0_));
            if (!mccNOTSET(&c)) m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            mccGrowMatrix(&c, m_, n_);
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_c = J_c = 0;}
            else { I_c = 1; J_c=mccM(&c)-m_; }
            p_c = mccPR(&c) + 0 + mccM(&c) * 0;
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_c += J_c)
               {
                  for (i_=0; i_<m_; ++i_, p_c+=I_c, p_IM0_+=I_IM0_)
                  {
                     *p_c = ((int)*p_IM0_);
                  }
               }
            }
         }
         /* cp(1,i2,:)  = zeros(1,length(i2),M); */
         I0_ = mccGetLength(&i2);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM1_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 290);
         mccFindIndex(&IM0_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_cp;
            int I_cp=1, J_cp;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM1_), mccN(&IM1_));
            if (!mccNOTSET(&cp)) m_ = mcmCalcResultSize(m_, &n_, mccM(&cp), mccN(&cp));
            mccGrowMatrix(&cp, m_, n_);
            if (mccM(&cp) == 1 && mccN(&cp) == 1) { I_cp = J_cp = 0;}
            else { I_cp = 1; J_cp=mccM(&cp)-m_; }
            p_cp = mccPR(&cp) + 0 + mccM(&cp) * 0;
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_cp += J_cp)
               {
                  for (i_=0; i_<m_; ++i_, p_cp+=I_cp, p_IM1_+=I_IM1_)
                  {
                     *p_cp = ((int)*p_IM1_);
                  }
               }
            }
         }
         /* end */
      }
      
      /* % 3 - condition aux limites en x=x1 */
      /* % a - conditionnement des donnees */
      /* u = (T1(:,:,3) ./ (dx.^2) + T1(:,:,2) ./ (2.*dx)); */
      if(mccNOTSET(&T1))
      {
         mexErrMsgTxt( "variable T1 undefined, line 295" );
      }
      if(mccNOTSET(&T1))
      {
         mexErrMsgTxt( "variable T1 undefined, line 295" );
      }
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_u;
         int I_u=1;
         double *p_T1;
         int I_T1=1, J_T1;
         double *p_dx;
         int I_dx=1;
         double *p_1T1;
         int I_1T1=1, J_1T1;
         double *p_1dx;
         int I_1dx=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&T1), 1);
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dx), mccN(&dx));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&T1), 1);
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dx), mccN(&dx));
         mccAllocateMatrix(&u, m_, n_);
         mccCheckMatrixSize(&T1, mccM(&T1), 3);
         mccCheckMatrixSize(&T1, mccM(&T1), 2);
         I_u = (mccM(&u) != 1 || mccN(&u) != 1);
         p_u = mccPR(&u);
         if (mccM(&T1) == 1) { I_T1 = J_T1 = 0;}
         else { I_T1 = 1; J_T1=mccM(&T1)-m_; }
         p_T1 = mccPR(&T1) + 0 + mccM(&T1) * (3-1);
         I_dx = (mccM(&dx) != 1 || mccN(&dx) != 1);
         p_dx = mccPR(&dx);
         if (mccM(&T1) == 1) { I_1T1 = J_1T1 = 0;}
         else { I_1T1 = 1; J_1T1=mccM(&T1)-m_; }
         p_1T1 = mccPR(&T1) + 0 + mccM(&T1) * (2-1);
         I_1dx = (mccM(&dx) != 1 || mccN(&dx) != 1);
         p_1dx = mccPR(&dx);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_T1 += J_T1, p_1T1 += J_1T1)
            {
               for (i_=0; i_<m_; ++i_, p_u+=I_u, p_T1+=I_T1, p_dx+=I_dx, p_1T1+=I_1T1, p_1dx+=I_1dx)
               {
                  *p_u = ((*p_T1 / (double) mcmRealPowerInt(*p_dx, 2)) + (*p_1T1 / (double) (2 * (double) *p_1dx)));
               }
            }
         }
      }
      /* v = (T1(:,:,1) - 2 .* T1(:,:,3) ./ (dx.^2)); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_v;
         int I_v=1;
         double *p_T1;
         int I_T1=1, J_T1;
         double *p_1T1;
         int I_1T1=1, J_1T1;
         double *p_dx;
         int I_dx=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&T1), 1);
         m_ = mcmCalcResultSize(m_, &n_, mccM(&T1), 1);
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dx), mccN(&dx));
         mccAllocateMatrix(&v, m_, n_);
         mccCheckMatrixSize(&T1, mccM(&T1), 1);
         mccCheckMatrixSize(&T1, mccM(&T1), 3);
         I_v = (mccM(&v) != 1 || mccN(&v) != 1);
         p_v = mccPR(&v);
         if (mccM(&T1) == 1) { I_T1 = J_T1 = 0;}
         else { I_T1 = 1; J_T1=mccM(&T1)-m_; }
         p_T1 = mccPR(&T1) + 0 + mccM(&T1) * (1-1);
         if (mccM(&T1) == 1) { I_1T1 = J_1T1 = 0;}
         else { I_1T1 = 1; J_1T1=mccM(&T1)-m_; }
         p_1T1 = mccPR(&T1) + 0 + mccM(&T1) * (3-1);
         I_dx = (mccM(&dx) != 1 || mccN(&dx) != 1);
         p_dx = mccPR(&dx);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_T1 += J_T1, p_1T1 += J_1T1)
            {
               for (i_=0; i_<m_; ++i_, p_v+=I_v, p_T1+=I_T1, p_1T1+=I_1T1, p_dx+=I_dx)
               {
                  *p_v = (*p_T1 - ((2 * (double) *p_1T1) / (double) mcmRealPowerInt(*p_dx, 2)));
               }
            }
         }
      }
      /* w = (T1(:,:,3) ./ (dx.^2) - T1(:,:,2) ./ (2*dx)); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_w;
         int I_w=1;
         double *p_T1;
         int I_T1=1, J_T1;
         double *p_dx;
         int I_dx=1;
         double *p_1T1;
         int I_1T1=1, J_1T1;
         double *p_1dx;
         int I_1dx=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&T1), 1);
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dx), mccN(&dx));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&T1), 1);
         m_ = mcmCalcResultSize(m_, &n_, mccM(&dx), mccN(&dx));
         mccAllocateMatrix(&w, m_, n_);
         mccCheckMatrixSize(&T1, mccM(&T1), 3);
         mccCheckMatrixSize(&T1, mccM(&T1), 2);
         I_w = (mccM(&w) != 1 || mccN(&w) != 1);
         p_w = mccPR(&w);
         if (mccM(&T1) == 1) { I_T1 = J_T1 = 0;}
         else { I_T1 = 1; J_T1=mccM(&T1)-m_; }
         p_T1 = mccPR(&T1) + 0 + mccM(&T1) * (3-1);
         I_dx = (mccM(&dx) != 1 || mccN(&dx) != 1);
         p_dx = mccPR(&dx);
         if (mccM(&T1) == 1) { I_1T1 = J_1T1 = 0;}
         else { I_1T1 = 1; J_1T1=mccM(&T1)-m_; }
         p_1T1 = mccPR(&T1) + 0 + mccM(&T1) * (2-1);
         I_1dx = (mccM(&dx) != 1 || mccN(&dx) != 1);
         p_1dx = mccPR(&dx);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_T1 += J_T1, p_1T1 += J_1T1)
            {
               for (i_=0; i_<m_; ++i_, p_w+=I_w, p_T1+=I_T1, p_dx+=I_dx, p_1T1+=I_1T1, p_1dx+=I_1dx)
               {
                  *p_w = ((*p_T1 / (double) mcmRealPowerInt(*p_dx, 2)) - (*p_1T1 / (double) (2 * (double) *p_1dx)));
               }
            }
         }
      }
      
      /* warning off */
      Mprhs_[0] = &S44_;
      Mplhs_[0] = 0;
      mccCallMATLAB(0, Mplhs_, 1, Mprhs_, "warning", 299);
      mccPrint(Mplhs_[0], 0);
      /* U = - w ./ u; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_U;
         int I_U=1;
         double *p_w;
         int I_w=1;
         double *p_u;
         int I_u=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&w), mccN(&w));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&u), mccN(&u));
         mccAllocateMatrix(&U, m_, n_);
         I_U = (mccM(&U) != 1 || mccN(&U) != 1);
         p_U = mccPR(&U);
         I_w = (mccM(&w) != 1 || mccN(&w) != 1);
         p_w = mccPR(&w);
         I_u = (mccM(&u) != 1 || mccN(&u) != 1);
         p_u = mccPR(&u);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_U+=I_U, p_w+=I_w, p_u+=I_u)
               {
                  *p_U = ((-*p_w) / (double) *p_u);
               }
            }
         }
      }
      /* V = - v ./ u; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_V;
         int I_V=1;
         double *p_v;
         int I_v=1;
         double *p_u;
         int I_u=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&v), mccN(&v));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&u), mccN(&u));
         mccAllocateMatrix(&V, m_, n_);
         I_V = (mccM(&V) != 1 || mccN(&V) != 1);
         p_V = mccPR(&V);
         I_v = (mccM(&v) != 1 || mccN(&v) != 1);
         p_v = mccPR(&v);
         I_u = (mccM(&u) != 1 || mccN(&u) != 1);
         p_u = mccPR(&u);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_V+=I_V, p_v+=I_v, p_u+=I_u)
               {
                  *p_V = ((-*p_v) / (double) *p_u);
               }
            }
         }
      }
      /* W = V1 ./ u; */
      if(mccNOTSET(&V1))
      {
         mexErrMsgTxt( "variable V1 undefined, line 302" );
      }
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_W;
         int I_W=1;
         double *p_V1;
         int I_V1=1;
         double *p_u;
         int I_u=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&V1), mccN(&V1));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&u), mccN(&u));
         mccAllocateMatrix(&W, m_, n_);
         I_W = (mccM(&W) != 1 || mccN(&W) != 1);
         p_W = mccPR(&W);
         I_V1 = (mccM(&V1) != 1 || mccN(&V1) != 1);
         p_V1 = mccPR(&V1);
         I_u = (mccM(&u) != 1 || mccN(&u) != 1);
         p_u = mccPR(&u);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_W+=I_W, p_V1+=I_V1, p_u+=I_u)
               {
                  *p_W = (*p_V1 / (double) *p_u);
               }
            }
         }
      }
      /* X = - w ./ v; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_X;
         int I_X=1;
         double *p_w;
         int I_w=1;
         double *p_v;
         int I_v=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&w), mccN(&w));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&v), mccN(&v));
         mccAllocateMatrix(&X, m_, n_);
         I_X = (mccM(&X) != 1 || mccN(&X) != 1);
         p_X = mccPR(&X);
         I_w = (mccM(&w) != 1 || mccN(&w) != 1);
         p_w = mccPR(&w);
         I_v = (mccM(&v) != 1 || mccN(&v) != 1);
         p_v = mccPR(&v);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_X+=I_X, p_w+=I_w, p_v+=I_v)
               {
                  *p_X = ((-*p_w) / (double) *p_v);
               }
            }
         }
      }
      /* Y = V1 ./ v; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_Y;
         int I_Y=1;
         double *p_V1;
         int I_V1=1;
         double *p_v;
         int I_v=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&V1), mccN(&V1));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&v), mccN(&v));
         mccAllocateMatrix(&Y, m_, n_);
         I_Y = (mccM(&Y) != 1 || mccN(&Y) != 1);
         p_Y = mccPR(&Y);
         I_V1 = (mccM(&V1) != 1 || mccN(&V1) != 1);
         p_V1 = mccPR(&V1);
         I_v = (mccM(&v) != 1 || mccN(&v) != 1);
         p_v = mccPR(&v);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_Y+=I_Y, p_V1+=I_V1, p_v+=I_v)
               {
                  *p_Y = (*p_V1 / (double) *p_v);
               }
            }
         }
      }
      /* Z = V1 ./ w; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_Z;
         int I_Z=1;
         double *p_V1;
         int I_V1=1;
         double *p_w;
         int I_w=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&V1), mccN(&V1));
         m_ = mcmCalcResultSize(m_, &n_, mccM(&w), mccN(&w));
         mccAllocateMatrix(&Z, m_, n_);
         I_Z = (mccM(&Z) != 1 || mccN(&Z) != 1);
         p_Z = mccPR(&Z);
         I_V1 = (mccM(&V1) != 1 || mccN(&V1) != 1);
         p_V1 = mccPR(&V1);
         I_w = (mccM(&w) != 1 || mccN(&w) != 1);
         p_w = mccPR(&w);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_Z+=I_Z, p_V1+=I_V1, p_w+=I_w)
               {
                  *p_Z = (*p_V1 / (double) *p_w);
               }
            }
         }
      }
      /* warning on  */
      Mprhs_[0] = &S45_;
      Mplhs_[0] = 0;
      mccCallMATLAB(0, Mplhs_, 1, Mprhs_, "warning", 306);
      mccPrint(Mplhs_[0], 0);
      
      /* comp = ones(1,M); */
      mccOnesMN(&M_comp, 1, M);
      
      /* % b - cas  u~=0 */
      /* i1 = find( u(1,:) ~= 0 ); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         unsigned short *p_BM0_;
         int I_BM0_=1;
         double *p_u;
         int I_u=1, J_u;
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&u));
         mccAllocateMatrix(&BM0_, m_, n_);
         mccCheckMatrixSize(&u, 1, mccN(&u));
         I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
         p_BM0_ = mccSPR(&BM0_);
         if (mccN(&u) == 1) { I_u = J_u = 0;}
         else { I_u = 1; J_u=mccM(&u)-m_; }
         p_u = mccPR(&u) + (1-1) + mccM(&u) * 0;
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_u += J_u)
            {
               for (i_=0; i_<m_; ++i_, p_BM0_+=I_BM0_, p_u+=I_u)
               {
                  *p_BM0_ = ( (*p_u != 0) || mccREL_NAN(*p_u) );
               }
            }
         }
      }
      mccFind(&i1, &BM0_);
      /* i2 = find( u(2,:) ~=0) ; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         unsigned short *p_BM0_;
         int I_BM0_=1;
         double *p_u;
         int I_u=1, J_u;
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&u));
         mccAllocateMatrix(&BM0_, m_, n_);
         mccCheckMatrixSize(&u, 2, mccN(&u));
         I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
         p_BM0_ = mccSPR(&BM0_);
         if (mccN(&u) == 1) { I_u = J_u = 0;}
         else { I_u = 1; J_u=mccM(&u)-m_; }
         p_u = mccPR(&u) + (2-1) + mccM(&u) * 0;
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_u += J_u)
            {
               for (i_=0; i_<m_; ++i_, p_BM0_+=I_BM0_, p_u+=I_u)
               {
                  *p_BM0_ = ( (*p_u != 0) || mccREL_NAN(*p_u) );
               }
            }
         }
      }
      mccFind(&i2, &BM0_);
      /* if ~ (isempty(i1) & isempty(i2)) */
      B2_ = mccIsEmpty(&i1);
      B1_ = mccIsEmpty(&i2);
      B0_ = (!(!!B2_ && !!B1_));
      if ((double)B0_)
      {
         /* a(K,i1,:) = a(K,i1,:) + c(K,i1,:) .* permute(U(1,i1)' * comp,[3 1 2]); */
         mccFindIndex(&IM0_, &i1);
         mccFindIndex(&IM2_, &i1);
         mccFindIndex(&IM3_, &i1);
         mccFindIndex(&IM4_, &IM3_);
         mccFindIndex(&IM5_, &IM4_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM0_;
            int I_RM0_=1;
            double *p_U;
            int I_U=1;
            int *p_IM5_;
            int I_IM5_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM5_) * mccN(&IM5_)));
            mccAllocateMatrix(&RM0_, m_, n_);
            mccCheckMatrixSize(&U, 1, mccGetMaxIndex(&IM5_ ,mccN(&U)));
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            I_IM5_ = (mccM(&IM5_) != 1 || mccN(&IM5_) != 1);
            p_IM5_ = mccIPR(&IM5_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM5_ += I_IM5_)
               {
                  p_U = mccPR(&U) + mccM(&U) * ((int)(*p_IM5_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_U+=I_U)
                  {
                     *p_RM0_ = *p_U;
                  }
               }
            }
         }
         mccConjTrans(&RM1_, &RM0_);
         mccRealMatrixMultiply(&RM2_, &RM1_, &M_comp);
         Mprhs_[0] = &RM2_;
         Mprhs_[1] = &S20_;
         Mplhs_[0] = &RM3_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 314);
         mccFindIndex(&IM1_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_a;
            int I_a=1, J_a;
            double *p_1a;
            int I_1a=1, J_1a;
            double *p_c;
            int I_c=1, J_c;
            double *p_RM3_;
            int I_RM3_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&a), mccN(&a));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            if (!mccNOTSET(&a)) m_ = mcmCalcResultSize(m_, &n_, mccM(&a), mccN(&a));
            mccGrowMatrix(&a, m_, n_);
            mccCheckMatrixSize(&a, mccM(&a), mccN(&a));
            mccCheckMatrixSize(&c, mccM(&c), mccN(&c));
            if (mccM(&a) == 1 && mccN(&a) == 1) { I_a = J_a = 0;}
            else { I_a = 1; J_a=mccM(&a)-m_; }
            p_a = mccPR(&a) + 0 + mccM(&a) * 0;
            if (mccM(&a) == 1 && mccN(&a) == 1) { I_1a = J_1a = 0;}
            else { I_1a = 1; J_1a=mccM(&a)-m_; }
            p_1a = mccPR(&a) + 0 + mccM(&a) * 0;
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_c = J_c = 0;}
            else { I_c = 1; J_c=mccM(&c)-m_; }
            p_c = mccPR(&c) + 0 + mccM(&c) * 0;
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_a += J_a, p_1a += J_1a, p_c += J_c)
               {
                  for (i_=0; i_<m_; ++i_, p_a+=I_a, p_1a+=I_1a, p_c+=I_c, p_RM3_+=I_RM3_)
                  {
                     *p_a = (*p_1a + (*p_c * (double) *p_RM3_));
                  }
               }
            }
         }
         /* ap(K,i2,:) = ap(K,i2,:) + cp(K,i2,:) .* permute(U(2,i2)' * comp,[3 1 2]); */
         mccFindIndex(&IM4_, &i2);
         mccFindIndex(&IM3_, &i2);
         mccFindIndex(&IM2_, &i2);
         mccFindIndex(&IM0_, &IM2_);
         mccFindIndex(&IM1_, &IM0_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM3_;
            int I_RM3_=1;
            double *p_U;
            int I_U=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM1_) * mccN(&IM1_)));
            mccAllocateMatrix(&RM3_, m_, n_);
            mccCheckMatrixSize(&U, 2, mccGetMaxIndex(&IM1_ ,mccN(&U)));
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_U = mccPR(&U) + mccM(&U) * ((int)(*p_IM1_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM3_+=I_RM3_, p_U+=I_U)
                  {
                     *p_RM3_ = *p_U;
                  }
               }
            }
         }
         mccConjTrans(&RM2_, &RM3_);
         mccRealMatrixMultiply(&RM1_, &RM2_, &M_comp);
         Mprhs_[0] = &RM1_;
         Mprhs_[1] = &S21_;
         Mplhs_[0] = &RM0_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 315);
         mccFindIndex(&IM5_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_ap;
            int I_ap=1, J_ap;
            double *p_1ap;
            int I_1ap=1, J_1ap;
            double *p_cp;
            int I_cp=1, J_cp;
            double *p_RM0_;
            int I_RM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), mccN(&ap));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&cp), mccN(&cp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM0_), mccN(&RM0_));
            if (!mccNOTSET(&ap)) m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), mccN(&ap));
            mccGrowMatrix(&ap, m_, n_);
            mccCheckMatrixSize(&ap, mccM(&ap), mccN(&ap));
            mccCheckMatrixSize(&cp, mccM(&cp), mccN(&cp));
            if (mccM(&ap) == 1 && mccN(&ap) == 1) { I_ap = J_ap = 0;}
            else { I_ap = 1; J_ap=mccM(&ap)-m_; }
            p_ap = mccPR(&ap) + 0 + mccM(&ap) * 0;
            if (mccM(&ap) == 1 && mccN(&ap) == 1) { I_1ap = J_1ap = 0;}
            else { I_1ap = 1; J_1ap=mccM(&ap)-m_; }
            p_1ap = mccPR(&ap) + 0 + mccM(&ap) * 0;
            if (mccM(&cp) == 1 && mccN(&cp) == 1) { I_cp = J_cp = 0;}
            else { I_cp = 1; J_cp=mccM(&cp)-m_; }
            p_cp = mccPR(&cp) + 0 + mccM(&cp) * 0;
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_ap += J_ap, p_1ap += J_1ap, p_cp += J_cp)
               {
                  for (i_=0; i_<m_; ++i_, p_ap+=I_ap, p_1ap+=I_1ap, p_cp+=I_cp, p_RM0_+=I_RM0_)
                  {
                     *p_ap = (*p_1ap + (*p_cp * (double) *p_RM0_));
                  }
               }
            }
         }
         
         /* b(K,i1,:) = b(K,i1,:) + c(K,i1,:) .* permute(V(1,i1)' * comp,[3 1 2]); */
         mccFindIndex(&IM0_, &i1);
         mccFindIndex(&IM2_, &i1);
         mccFindIndex(&IM3_, &i1);
         mccFindIndex(&IM4_, &IM3_);
         mccFindIndex(&IM5_, &IM4_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM0_;
            int I_RM0_=1;
            double *p_V;
            int I_V=1;
            int *p_IM5_;
            int I_IM5_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM5_) * mccN(&IM5_)));
            mccAllocateMatrix(&RM0_, m_, n_);
            mccCheckMatrixSize(&V, 1, mccGetMaxIndex(&IM5_ ,mccN(&V)));
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            I_IM5_ = (mccM(&IM5_) != 1 || mccN(&IM5_) != 1);
            p_IM5_ = mccIPR(&IM5_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM5_ += I_IM5_)
               {
                  p_V = mccPR(&V) + mccM(&V) * ((int)(*p_IM5_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_V+=I_V)
                  {
                     *p_RM0_ = *p_V;
                  }
               }
            }
         }
         mccConjTrans(&RM1_, &RM0_);
         mccRealMatrixMultiply(&RM2_, &RM1_, &M_comp);
         Mprhs_[0] = &RM2_;
         Mprhs_[1] = &S22_;
         Mplhs_[0] = &RM3_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 317);
         mccFindIndex(&IM1_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_b;
            int I_b=1, J_b;
            double *p_1b;
            int I_1b=1, J_1b;
            double *p_c;
            int I_c=1, J_c;
            double *p_RM3_;
            int I_RM3_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            if (!mccNOTSET(&b)) m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            mccGrowMatrix(&b, m_, n_);
            mccCheckMatrixSize(&b, mccM(&b), mccN(&b));
            mccCheckMatrixSize(&c, mccM(&c), mccN(&c));
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_b = J_b = 0;}
            else { I_b = 1; J_b=mccM(&b)-m_; }
            p_b = mccPR(&b) + 0 + mccM(&b) * 0;
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_1b = J_1b = 0;}
            else { I_1b = 1; J_1b=mccM(&b)-m_; }
            p_1b = mccPR(&b) + 0 + mccM(&b) * 0;
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_c = J_c = 0;}
            else { I_c = 1; J_c=mccM(&c)-m_; }
            p_c = mccPR(&c) + 0 + mccM(&c) * 0;
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_b += J_b, p_1b += J_1b, p_c += J_c)
               {
                  for (i_=0; i_<m_; ++i_, p_b+=I_b, p_1b+=I_1b, p_c+=I_c, p_RM3_+=I_RM3_)
                  {
                     *p_b = (*p_1b + (*p_c * (double) *p_RM3_));
                  }
               }
            }
         }
         /* bp(K,i2,:) = bp(K,i2,:) + cp(K,i2,:) .* permute(V(2,i2)' * comp,[3 1 2]); */
         mccFindIndex(&IM4_, &i2);
         mccFindIndex(&IM3_, &i2);
         mccFindIndex(&IM2_, &i2);
         mccFindIndex(&IM0_, &IM2_);
         mccFindIndex(&IM1_, &IM0_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM3_;
            int I_RM3_=1;
            double *p_V;
            int I_V=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM1_) * mccN(&IM1_)));
            mccAllocateMatrix(&RM3_, m_, n_);
            mccCheckMatrixSize(&V, 2, mccGetMaxIndex(&IM1_ ,mccN(&V)));
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_V = mccPR(&V) + mccM(&V) * ((int)(*p_IM1_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM3_+=I_RM3_, p_V+=I_V)
                  {
                     *p_RM3_ = *p_V;
                  }
               }
            }
         }
         mccConjTrans(&RM2_, &RM3_);
         mccRealMatrixMultiply(&RM1_, &RM2_, &M_comp);
         Mprhs_[0] = &RM1_;
         Mprhs_[1] = &S23_;
         Mplhs_[0] = &RM0_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 318);
         mccFindIndex(&IM5_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_bp;
            int I_bp=1, J_bp;
            double *p_1bp;
            int I_1bp=1, J_1bp;
            double *p_cp;
            int I_cp=1, J_cp;
            double *p_RM0_;
            int I_RM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&cp), mccN(&cp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM0_), mccN(&RM0_));
            if (!mccNOTSET(&bp)) m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            mccGrowMatrix(&bp, m_, n_);
            mccCheckMatrixSize(&bp, mccM(&bp), mccN(&bp));
            mccCheckMatrixSize(&cp, mccM(&cp), mccN(&cp));
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_bp = J_bp = 0;}
            else { I_bp = 1; J_bp=mccM(&bp)-m_; }
            p_bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_1bp = J_1bp = 0;}
            else { I_1bp = 1; J_1bp=mccM(&bp)-m_; }
            p_1bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            if (mccM(&cp) == 1 && mccN(&cp) == 1) { I_cp = J_cp = 0;}
            else { I_cp = 1; J_cp=mccM(&cp)-m_; }
            p_cp = mccPR(&cp) + 0 + mccM(&cp) * 0;
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_bp += J_bp, p_1bp += J_1bp, p_cp += J_cp)
               {
                  for (i_=0; i_<m_; ++i_, p_bp+=I_bp, p_1bp+=I_1bp, p_cp+=I_cp, p_RM0_+=I_RM0_)
                  {
                     *p_bp = (*p_1bp + (*p_cp * (double) *p_RM0_));
                  }
               }
            }
         }
         
         /* s(K,:) = s(K,:) + squeeze( f .* sum( c(K,i1,:) .* permute(W(1,i1)' * comp,[3 1 2]),2) + ... */
         mccFindIndex(&IM1_, &i1);
         mccFindIndex(&IM0_, &IM1_);
         mccFindIndex(&IM2_, &i1);
         mccFindIndex(&IM3_, &IM2_);
         mccFindIndex(&IM4_, &IM3_);
         mccFindIndex(&IM5_, &IM4_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM0_;
            int I_RM0_=1;
            double *p_W;
            int I_W=1;
            int *p_IM5_;
            int I_IM5_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM5_) * mccN(&IM5_)));
            mccAllocateMatrix(&RM0_, m_, n_);
            mccCheckMatrixSize(&W, 1, mccGetMaxIndex(&IM5_ ,mccN(&W)));
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            I_IM5_ = (mccM(&IM5_) != 1 || mccN(&IM5_) != 1);
            p_IM5_ = mccIPR(&IM5_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM5_ += I_IM5_)
               {
                  p_W = mccPR(&W) + mccM(&W) * ((int)(*p_IM5_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_W+=I_W)
                  {
                     *p_RM0_ = *p_W;
                  }
               }
            }
         }
         mccConjTrans(&RM1_, &RM0_);
         mccRealMatrixMultiply(&RM2_, &RM1_, &M_comp);
         Mprhs_[0] = &RM2_;
         Mprhs_[1] = &S24_;
         Mplhs_[0] = &RM3_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 320);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM4_;
            int I_RM4_=1;
            double *p_c;
            int I_c=1, J_c;
            double *p_RM3_;
            int I_RM3_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            mccAllocateMatrix(&RM4_, m_, n_);
            mccCheckMatrixSize(&c, mccM(&c), mccN(&c));
            I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
            p_RM4_ = mccPR(&RM4_);
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_c = J_c = 0;}
            else { I_c = 1; J_c=mccM(&c)-m_; }
            p_c = mccPR(&c) + 0 + mccM(&c) * 0;
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_c += J_c)
               {
                  for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_c+=I_c, p_RM3_+=I_RM3_)
                  {
                     *p_RM4_ = (*p_c * (double) *p_RM3_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM4_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM5_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 320);
         mccFindIndex(&IM6_, &i2);
         mccFindIndex(&IM7_, &IM6_);
         mccFindIndex(&IM8_, &i2);
         mccFindIndex(&IM9_, &IM8_);
         mccFindIndex(&IM10_, &IM9_);
         mccFindIndex(&IM11_, &IM10_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM6_;
            int I_RM6_=1;
            double *p_W;
            int I_W=1;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM11_) * mccN(&IM11_)));
            mccAllocateMatrix(&RM6_, m_, n_);
            mccCheckMatrixSize(&W, 2, mccGetMaxIndex(&IM11_ ,mccN(&W)));
            I_RM6_ = (mccM(&RM6_) != 1 || mccN(&RM6_) != 1);
            p_RM6_ = mccPR(&RM6_);
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM11_ += I_IM11_)
               {
                  p_W = mccPR(&W) + mccM(&W) * ((int)(*p_IM11_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM6_+=I_RM6_, p_W+=I_W)
                  {
                     *p_RM6_ = *p_W;
                  }
               }
            }
         }
         mccConjTrans(&RM7_, &RM6_);
         mccRealMatrixMultiply(&RM8_, &RM7_, &M_comp);
         Mprhs_[0] = &RM8_;
         Mprhs_[1] = &S25_;
         Mplhs_[0] = &RM9_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 321);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM10_;
            int I_RM10_=1;
            double *p_cp;
            int I_cp=1, J_cp;
            double *p_RM9_;
            int I_RM9_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&cp), mccN(&cp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM9_), mccN(&RM9_));
            mccAllocateMatrix(&RM10_, m_, n_);
            mccCheckMatrixSize(&cp, mccM(&cp), mccN(&cp));
            I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
            p_RM10_ = mccPR(&RM10_);
            if (mccM(&cp) == 1 && mccN(&cp) == 1) { I_cp = J_cp = 0;}
            else { I_cp = 1; J_cp=mccM(&cp)-m_; }
            p_cp = mccPR(&cp) + 0 + mccM(&cp) * 0;
            I_RM9_ = (mccM(&RM9_) != 1 || mccN(&RM9_) != 1);
            p_RM9_ = mccPR(&RM9_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_cp += J_cp)
               {
                  for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_cp+=I_cp, p_RM9_+=I_RM9_)
                  {
                     *p_RM10_ = (*p_cp * (double) *p_RM9_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM10_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM11_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 321);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM12_;
            int I_RM12_=1;
            double *p_f;
            int I_f=1;
            double *p_RM5_;
            int I_RM5_=1;
            double *p_1f;
            int I_1f=1;
            double *p_RM11_;
            int I_RM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM5_), mccN(&RM5_));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM11_), mccN(&RM11_));
            mccAllocateMatrix(&RM12_, m_, n_);
            I_RM12_ = (mccM(&RM12_) != 1 || mccN(&RM12_) != 1);
            p_RM12_ = mccPR(&RM12_);
            I_f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_f = mccPR(&f);
            I_RM5_ = (mccM(&RM5_) != 1 || mccN(&RM5_) != 1);
            p_RM5_ = mccPR(&RM5_);
            I_1f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_1f = mccPR(&f);
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM12_+=I_RM12_, p_f+=I_f, p_RM5_+=I_RM5_, p_1f+=I_1f, p_RM11_+=I_RM11_)
                  {
                     *p_RM12_ = ((*p_f * (double) *p_RM5_) + ((1 - *p_1f) * (double) *p_RM11_));
                  }
               }
            }
         }
         mccLOG(&RM12_) = 0;
         mccSTRING(&RM12_) = 0;
         Mprhs_[0] = &RM12_;
         Mplhs_[0] = &RM13_;
         mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "squeeze", 321);
         mccConjTrans(&RM14_, &RM13_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_s;
            int I_s=1, J_s;
            double *p_1s;
            int I_1s=1, J_1s;
            double *p_RM14_;
            int I_RM14_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM14_), mccN(&RM14_));
            if (!mccNOTSET(&s)) m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            mccGrowMatrix(&s, K, n_);
            mccCheckMatrixSize(&s, K, mccN(&s));
            if (mccN(&s) == 1) { I_s = J_s = 0;}
            else { I_s = 1; J_s=mccM(&s)-m_; }
            p_s = mccPR(&s) + (K-1) + mccM(&s) * 0;
            if (mccN(&s) == 1) { I_1s = J_1s = 0;}
            else { I_1s = 1; J_1s=mccM(&s)-m_; }
            p_1s = mccPR(&s) + (K-1) + mccM(&s) * 0;
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_s += J_s, p_1s += J_1s)
               {
                  for (i_=0; i_<m_; ++i_, p_s+=I_s, p_1s+=I_1s, p_RM14_+=I_RM14_)
                  {
                     *p_s = (*p_1s + *p_RM14_);
                  }
               }
            }
         }
         
         /* c(K,i1,:)   = zeros(1,length(i1),M); */
         I0_ = mccGetLength(&i1);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM10_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 323);
         mccFindIndex(&IM11_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_c;
            int I_c=1, J_c;
            int *p_IM10_;
            int I_IM10_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM10_), mccN(&IM10_));
            if (!mccNOTSET(&c)) m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            mccGrowMatrix(&c, m_, n_);
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_c = J_c = 0;}
            else { I_c = 1; J_c=mccM(&c)-m_; }
            p_c = mccPR(&c) + 0 + mccM(&c) * 0;
            I_IM10_ = (mccM(&IM10_) != 1 || mccN(&IM10_) != 1);
            p_IM10_ = mccIPR(&IM10_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_c += J_c)
               {
                  for (i_=0; i_<m_; ++i_, p_c+=I_c, p_IM10_+=I_IM10_)
                  {
                     *p_c = ((int)*p_IM10_);
                  }
               }
            }
         }
         /* cp(K,i2,:)  = zeros(1,length(i2),M); */
         I0_ = mccGetLength(&i2);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM11_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 324);
         mccFindIndex(&IM10_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_cp;
            int I_cp=1, J_cp;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM11_), mccN(&IM11_));
            if (!mccNOTSET(&cp)) m_ = mcmCalcResultSize(m_, &n_, mccM(&cp), mccN(&cp));
            mccGrowMatrix(&cp, m_, n_);
            if (mccM(&cp) == 1 && mccN(&cp) == 1) { I_cp = J_cp = 0;}
            else { I_cp = 1; J_cp=mccM(&cp)-m_; }
            p_cp = mccPR(&cp) + 0 + mccM(&cp) * 0;
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_cp += J_cp)
               {
                  for (i_=0; i_<m_; ++i_, p_cp+=I_cp, p_IM11_+=I_IM11_)
                  {
                     *p_cp = ((int)*p_IM11_);
                  }
               }
            }
         }
         /* end                      */
      }
      /* % c - cas u ==0 et v ~=0 */
      /* i1 = find(( u(1,:) == 0 ) & (v(1,:) ~= 0)); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         unsigned short *p_BM0_;
         int I_BM0_=1;
         double *p_u;
         int I_u=1, J_u;
         double *p_v;
         int I_v=1, J_v;
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&u));
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&v));
         mccAllocateMatrix(&BM0_, m_, n_);
         mccCheckMatrixSize(&u, 1, mccN(&u));
         mccCheckMatrixSize(&v, 1, mccN(&v));
         I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
         p_BM0_ = mccSPR(&BM0_);
         if (mccN(&u) == 1) { I_u = J_u = 0;}
         else { I_u = 1; J_u=mccM(&u)-m_; }
         p_u = mccPR(&u) + (1-1) + mccM(&u) * 0;
         if (mccN(&v) == 1) { I_v = J_v = 0;}
         else { I_v = 1; J_v=mccM(&v)-m_; }
         p_v = mccPR(&v) + (1-1) + mccM(&v) * 0;
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_u += J_u, p_v += J_v)
            {
               for (i_=0; i_<m_; ++i_, p_BM0_+=I_BM0_, p_u+=I_u, p_v+=I_v)
               {
                  *p_BM0_ = (!!( (*p_u == 0) && !mccREL_NAN(*p_u) ) && !!( (*p_v != 0) || mccREL_NAN(*p_v) ));
               }
            }
         }
      }
      mccFind(&i1, &BM0_);
      /* i2 = find(( u(2,:) ==0) & (v(2,:) ~= 0)); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         unsigned short *p_BM0_;
         int I_BM0_=1;
         double *p_u;
         int I_u=1, J_u;
         double *p_v;
         int I_v=1, J_v;
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&u));
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&v));
         mccAllocateMatrix(&BM0_, m_, n_);
         mccCheckMatrixSize(&u, 2, mccN(&u));
         mccCheckMatrixSize(&v, 2, mccN(&v));
         I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
         p_BM0_ = mccSPR(&BM0_);
         if (mccN(&u) == 1) { I_u = J_u = 0;}
         else { I_u = 1; J_u=mccM(&u)-m_; }
         p_u = mccPR(&u) + (2-1) + mccM(&u) * 0;
         if (mccN(&v) == 1) { I_v = J_v = 0;}
         else { I_v = 1; J_v=mccM(&v)-m_; }
         p_v = mccPR(&v) + (2-1) + mccM(&v) * 0;
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_u += J_u, p_v += J_v)
            {
               for (i_=0; i_<m_; ++i_, p_BM0_+=I_BM0_, p_u+=I_u, p_v+=I_v)
               {
                  *p_BM0_ = (!!( (*p_u == 0) && !mccREL_NAN(*p_u) ) && !!( (*p_v != 0) || mccREL_NAN(*p_v) ));
               }
            }
         }
      }
      mccFind(&i2, &BM0_);
      /* if ~ (isempty(i1) & isempty(i2)) */
      B0_ = mccIsEmpty(&i1);
      B1_ = mccIsEmpty(&i2);
      B2_ = (!(!!B0_ && !!B1_));
      if ((double)B2_)
      {
         /* a(K,i1,:) = a(K,i1,:) + b(K,i1,:) .* permute(X(1,i1)' * comp,[3 1 2]); */
         mccFindIndex(&IM10_, &i1);
         mccFindIndex(&IM9_, &i1);
         mccFindIndex(&IM8_, &i1);
         mccFindIndex(&IM7_, &IM8_);
         mccFindIndex(&IM6_, &IM7_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM14_;
            int I_RM14_=1;
            double *p_X;
            int I_X=1;
            int *p_IM6_;
            int I_IM6_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM6_) * mccN(&IM6_)));
            mccAllocateMatrix(&RM14_, m_, n_);
            mccCheckMatrixSize(&X, 1, mccGetMaxIndex(&IM6_ ,mccN(&X)));
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            I_IM6_ = (mccM(&IM6_) != 1 || mccN(&IM6_) != 1);
            p_IM6_ = mccIPR(&IM6_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM6_ += I_IM6_)
               {
                  p_X = mccPR(&X) + mccM(&X) * ((int)(*p_IM6_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM14_+=I_RM14_, p_X+=I_X)
                  {
                     *p_RM14_ = *p_X;
                  }
               }
            }
         }
         mccConjTrans(&RM13_, &RM14_);
         mccRealMatrixMultiply(&RM12_, &RM13_, &M_comp);
         Mprhs_[0] = &RM12_;
         Mprhs_[1] = &S26_;
         Mplhs_[0] = &RM11_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 330);
         mccFindIndex(&IM11_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_a;
            int I_a=1, J_a;
            double *p_1a;
            int I_1a=1, J_1a;
            double *p_b;
            int I_b=1, J_b;
            double *p_RM11_;
            int I_RM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&a), mccN(&a));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM11_), mccN(&RM11_));
            if (!mccNOTSET(&a)) m_ = mcmCalcResultSize(m_, &n_, mccM(&a), mccN(&a));
            mccGrowMatrix(&a, m_, n_);
            mccCheckMatrixSize(&a, mccM(&a), mccN(&a));
            mccCheckMatrixSize(&b, mccM(&b), mccN(&b));
            if (mccM(&a) == 1 && mccN(&a) == 1) { I_a = J_a = 0;}
            else { I_a = 1; J_a=mccM(&a)-m_; }
            p_a = mccPR(&a) + 0 + mccM(&a) * 0;
            if (mccM(&a) == 1 && mccN(&a) == 1) { I_1a = J_1a = 0;}
            else { I_1a = 1; J_1a=mccM(&a)-m_; }
            p_1a = mccPR(&a) + 0 + mccM(&a) * 0;
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_b = J_b = 0;}
            else { I_b = 1; J_b=mccM(&b)-m_; }
            p_b = mccPR(&b) + 0 + mccM(&b) * 0;
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_a += J_a, p_1a += J_1a, p_b += J_b)
               {
                  for (i_=0; i_<m_; ++i_, p_a+=I_a, p_1a+=I_1a, p_b+=I_b, p_RM11_+=I_RM11_)
                  {
                     *p_a = (*p_1a + (*p_b * (double) *p_RM11_));
                  }
               }
            }
         }
         /* ap(K,i2,:) = ap(K,i2,:) + bp(K,i2,:) .* permute(X(2,i2)' * comp,[3 1 2]); */
         mccFindIndex(&IM7_, &i2);
         mccFindIndex(&IM8_, &i2);
         mccFindIndex(&IM9_, &i2);
         mccFindIndex(&IM10_, &IM9_);
         mccFindIndex(&IM11_, &IM10_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM11_;
            int I_RM11_=1;
            double *p_X;
            int I_X=1;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM11_) * mccN(&IM11_)));
            mccAllocateMatrix(&RM11_, m_, n_);
            mccCheckMatrixSize(&X, 2, mccGetMaxIndex(&IM11_ ,mccN(&X)));
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM11_ += I_IM11_)
               {
                  p_X = mccPR(&X) + mccM(&X) * ((int)(*p_IM11_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM11_+=I_RM11_, p_X+=I_X)
                  {
                     *p_RM11_ = *p_X;
                  }
               }
            }
         }
         mccConjTrans(&RM12_, &RM11_);
         mccRealMatrixMultiply(&RM13_, &RM12_, &M_comp);
         Mprhs_[0] = &RM13_;
         Mprhs_[1] = &S27_;
         Mplhs_[0] = &RM14_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 331);
         mccFindIndex(&IM6_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_ap;
            int I_ap=1, J_ap;
            double *p_1ap;
            int I_1ap=1, J_1ap;
            double *p_bp;
            int I_bp=1, J_bp;
            double *p_RM14_;
            int I_RM14_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), mccN(&ap));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM14_), mccN(&RM14_));
            if (!mccNOTSET(&ap)) m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), mccN(&ap));
            mccGrowMatrix(&ap, m_, n_);
            mccCheckMatrixSize(&ap, mccM(&ap), mccN(&ap));
            mccCheckMatrixSize(&bp, mccM(&bp), mccN(&bp));
            if (mccM(&ap) == 1 && mccN(&ap) == 1) { I_ap = J_ap = 0;}
            else { I_ap = 1; J_ap=mccM(&ap)-m_; }
            p_ap = mccPR(&ap) + 0 + mccM(&ap) * 0;
            if (mccM(&ap) == 1 && mccN(&ap) == 1) { I_1ap = J_1ap = 0;}
            else { I_1ap = 1; J_1ap=mccM(&ap)-m_; }
            p_1ap = mccPR(&ap) + 0 + mccM(&ap) * 0;
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_bp = J_bp = 0;}
            else { I_bp = 1; J_bp=mccM(&bp)-m_; }
            p_bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_ap += J_ap, p_1ap += J_1ap, p_bp += J_bp)
               {
                  for (i_=0; i_<m_; ++i_, p_ap+=I_ap, p_1ap+=I_1ap, p_bp+=I_bp, p_RM14_+=I_RM14_)
                  {
                     *p_ap = (*p_1ap + (*p_bp * (double) *p_RM14_));
                  }
               }
            }
         }
         
         /* % bug corrige le 07/08/2000 */
         /* b(K-1,i1,:) = b(K-1,i1,:) + c(K-1,i1,:) .* permute(X(1,i1)' * comp,[3 1 2]); */
         mccFindIndex(&IM10_, &i1);
         mccFindIndex(&IM9_, &i1);
         mccFindIndex(&IM8_, &i1);
         mccFindIndex(&IM7_, &IM8_);
         mccFindIndex(&IM6_, &IM7_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM14_;
            int I_RM14_=1;
            double *p_X;
            int I_X=1;
            int *p_IM6_;
            int I_IM6_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM6_) * mccN(&IM6_)));
            mccAllocateMatrix(&RM14_, m_, n_);
            mccCheckMatrixSize(&X, 1, mccGetMaxIndex(&IM6_ ,mccN(&X)));
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            I_IM6_ = (mccM(&IM6_) != 1 || mccN(&IM6_) != 1);
            p_IM6_ = mccIPR(&IM6_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM6_ += I_IM6_)
               {
                  p_X = mccPR(&X) + mccM(&X) * ((int)(*p_IM6_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM14_+=I_RM14_, p_X+=I_X)
                  {
                     *p_RM14_ = *p_X;
                  }
               }
            }
         }
         mccConjTrans(&RM13_, &RM14_);
         mccRealMatrixMultiply(&RM12_, &RM13_, &M_comp);
         Mprhs_[0] = &RM12_;
         Mprhs_[1] = &S28_;
         Mplhs_[0] = &RM11_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 334);
         mccFindIndex(&IM11_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_b;
            int I_b=1, J_b;
            double *p_1b;
            int I_1b=1, J_1b;
            double *p_c;
            int I_c=1, J_c;
            double *p_RM11_;
            int I_RM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM11_), mccN(&RM11_));
            if (!mccNOTSET(&b)) m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            mccGrowMatrix(&b, m_, n_);
            mccCheckMatrixSize(&b, mccM(&b), mccN(&b));
            mccCheckMatrixSize(&c, mccM(&c), mccN(&c));
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_b = J_b = 0;}
            else { I_b = 1; J_b=mccM(&b)-m_; }
            p_b = mccPR(&b) + 0 + mccM(&b) * 0;
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_1b = J_1b = 0;}
            else { I_1b = 1; J_1b=mccM(&b)-m_; }
            p_1b = mccPR(&b) + 0 + mccM(&b) * 0;
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_c = J_c = 0;}
            else { I_c = 1; J_c=mccM(&c)-m_; }
            p_c = mccPR(&c) + 0 + mccM(&c) * 0;
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_b += J_b, p_1b += J_1b, p_c += J_c)
               {
                  for (i_=0; i_<m_; ++i_, p_b+=I_b, p_1b+=I_1b, p_c+=I_c, p_RM11_+=I_RM11_)
                  {
                     *p_b = (*p_1b + (*p_c * (double) *p_RM11_));
                  }
               }
            }
         }
         /* bp(K-1,i2,:) = bp(K-1,i2,:) + c(K-1,i2,:) .* permute(X(2,i2)' * comp,[3 1 2]); */
         mccFindIndex(&IM7_, &i2);
         mccFindIndex(&IM8_, &i2);
         mccFindIndex(&IM9_, &i2);
         mccFindIndex(&IM10_, &IM9_);
         mccFindIndex(&IM11_, &IM10_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM11_;
            int I_RM11_=1;
            double *p_X;
            int I_X=1;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM11_) * mccN(&IM11_)));
            mccAllocateMatrix(&RM11_, m_, n_);
            mccCheckMatrixSize(&X, 2, mccGetMaxIndex(&IM11_ ,mccN(&X)));
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM11_ += I_IM11_)
               {
                  p_X = mccPR(&X) + mccM(&X) * ((int)(*p_IM11_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM11_+=I_RM11_, p_X+=I_X)
                  {
                     *p_RM11_ = *p_X;
                  }
               }
            }
         }
         mccConjTrans(&RM12_, &RM11_);
         mccRealMatrixMultiply(&RM13_, &RM12_, &M_comp);
         Mprhs_[0] = &RM13_;
         Mprhs_[1] = &S29_;
         Mplhs_[0] = &RM14_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 335);
         mccFindIndex(&IM6_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_bp;
            int I_bp=1, J_bp;
            double *p_1bp;
            int I_1bp=1, J_1bp;
            double *p_c;
            int I_c=1, J_c;
            double *p_RM14_;
            int I_RM14_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM14_), mccN(&RM14_));
            if (!mccNOTSET(&bp)) m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            mccGrowMatrix(&bp, m_, n_);
            mccCheckMatrixSize(&bp, mccM(&bp), mccN(&bp));
            mccCheckMatrixSize(&c, mccM(&c), mccN(&c));
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_bp = J_bp = 0;}
            else { I_bp = 1; J_bp=mccM(&bp)-m_; }
            p_bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_1bp = J_1bp = 0;}
            else { I_1bp = 1; J_1bp=mccM(&bp)-m_; }
            p_1bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_c = J_c = 0;}
            else { I_c = 1; J_c=mccM(&c)-m_; }
            p_c = mccPR(&c) + 0 + mccM(&c) * 0;
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_bp += J_bp, p_1bp += J_1bp, p_c += J_c)
               {
                  for (i_=0; i_<m_; ++i_, p_bp+=I_bp, p_1bp+=I_1bp, p_c+=I_c, p_RM14_+=I_RM14_)
                  {
                     *p_bp = (*p_1bp + (*p_c * (double) *p_RM14_));
                  }
               }
            }
         }
         
         /* s(K,:) = s(K,:) + squeeze( f .* sum( b(K,i1,:) .* permute(Y(1,i1)' * comp,[3 1 2]),2) + ... */
         mccFindIndex(&IM11_, &i1);
         mccFindIndex(&IM10_, &IM11_);
         mccFindIndex(&IM9_, &i1);
         mccFindIndex(&IM8_, &IM9_);
         mccFindIndex(&IM7_, &IM8_);
         mccFindIndex(&IM6_, &IM7_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM14_;
            int I_RM14_=1;
            double *p_Y;
            int I_Y=1;
            int *p_IM6_;
            int I_IM6_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM6_) * mccN(&IM6_)));
            mccAllocateMatrix(&RM14_, m_, n_);
            mccCheckMatrixSize(&Y, 1, mccGetMaxIndex(&IM6_ ,mccN(&Y)));
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            I_IM6_ = (mccM(&IM6_) != 1 || mccN(&IM6_) != 1);
            p_IM6_ = mccIPR(&IM6_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM6_ += I_IM6_)
               {
                  p_Y = mccPR(&Y) + mccM(&Y) * ((int)(*p_IM6_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM14_+=I_RM14_, p_Y+=I_Y)
                  {
                     *p_RM14_ = *p_Y;
                  }
               }
            }
         }
         mccConjTrans(&RM13_, &RM14_);
         mccRealMatrixMultiply(&RM12_, &RM13_, &M_comp);
         Mprhs_[0] = &RM12_;
         Mprhs_[1] = &S30_;
         Mplhs_[0] = &RM11_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 337);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM10_;
            int I_RM10_=1;
            double *p_b;
            int I_b=1, J_b;
            double *p_RM11_;
            int I_RM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM11_), mccN(&RM11_));
            mccAllocateMatrix(&RM10_, m_, n_);
            mccCheckMatrixSize(&b, mccM(&b), mccN(&b));
            I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
            p_RM10_ = mccPR(&RM10_);
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_b = J_b = 0;}
            else { I_b = 1; J_b=mccM(&b)-m_; }
            p_b = mccPR(&b) + 0 + mccM(&b) * 0;
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_b += J_b)
               {
                  for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_b+=I_b, p_RM11_+=I_RM11_)
                  {
                     *p_RM10_ = (*p_b * (double) *p_RM11_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM10_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM9_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 337);
         mccFindIndex(&IM5_, &i2);
         mccFindIndex(&IM4_, &IM5_);
         mccFindIndex(&IM3_, &i2);
         mccFindIndex(&IM2_, &IM3_);
         mccFindIndex(&IM0_, &IM2_);
         mccFindIndex(&IM1_, &IM0_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM8_;
            int I_RM8_=1;
            double *p_Y;
            int I_Y=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM1_) * mccN(&IM1_)));
            mccAllocateMatrix(&RM8_, m_, n_);
            mccCheckMatrixSize(&Y, 2, mccGetMaxIndex(&IM1_ ,mccN(&Y)));
            I_RM8_ = (mccM(&RM8_) != 1 || mccN(&RM8_) != 1);
            p_RM8_ = mccPR(&RM8_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_Y = mccPR(&Y) + mccM(&Y) * ((int)(*p_IM1_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM8_+=I_RM8_, p_Y+=I_Y)
                  {
                     *p_RM8_ = *p_Y;
                  }
               }
            }
         }
         mccConjTrans(&RM7_, &RM8_);
         mccRealMatrixMultiply(&RM6_, &RM7_, &M_comp);
         Mprhs_[0] = &RM6_;
         Mprhs_[1] = &S31_;
         Mplhs_[0] = &RM5_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 338);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM4_;
            int I_RM4_=1;
            double *p_bp;
            int I_bp=1, J_bp;
            double *p_RM5_;
            int I_RM5_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM5_), mccN(&RM5_));
            mccAllocateMatrix(&RM4_, m_, n_);
            mccCheckMatrixSize(&bp, mccM(&bp), mccN(&bp));
            I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
            p_RM4_ = mccPR(&RM4_);
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_bp = J_bp = 0;}
            else { I_bp = 1; J_bp=mccM(&bp)-m_; }
            p_bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            I_RM5_ = (mccM(&RM5_) != 1 || mccN(&RM5_) != 1);
            p_RM5_ = mccPR(&RM5_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_bp += J_bp)
               {
                  for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_bp+=I_bp, p_RM5_+=I_RM5_)
                  {
                     *p_RM4_ = (*p_bp * (double) *p_RM5_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM4_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM3_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 338);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM2_;
            int I_RM2_=1;
            double *p_f;
            int I_f=1;
            double *p_RM9_;
            int I_RM9_=1;
            double *p_1f;
            int I_1f=1;
            double *p_RM3_;
            int I_RM3_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM9_), mccN(&RM9_));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            mccAllocateMatrix(&RM2_, m_, n_);
            I_RM2_ = (mccM(&RM2_) != 1 || mccN(&RM2_) != 1);
            p_RM2_ = mccPR(&RM2_);
            I_f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_f = mccPR(&f);
            I_RM9_ = (mccM(&RM9_) != 1 || mccN(&RM9_) != 1);
            p_RM9_ = mccPR(&RM9_);
            I_1f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_1f = mccPR(&f);
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM2_+=I_RM2_, p_f+=I_f, p_RM9_+=I_RM9_, p_1f+=I_1f, p_RM3_+=I_RM3_)
                  {
                     *p_RM2_ = ((*p_f * (double) *p_RM9_) + ((1 - *p_1f) * (double) *p_RM3_));
                  }
               }
            }
         }
         mccLOG(&RM2_) = 0;
         mccSTRING(&RM2_) = 0;
         Mprhs_[0] = &RM2_;
         Mplhs_[0] = &RM1_;
         mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "squeeze", 338);
         mccConjTrans(&RM0_, &RM1_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_s;
            int I_s=1, J_s;
            double *p_1s;
            int I_1s=1, J_1s;
            double *p_RM0_;
            int I_RM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM0_), mccN(&RM0_));
            if (!mccNOTSET(&s)) m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            mccGrowMatrix(&s, K, n_);
            mccCheckMatrixSize(&s, K, mccN(&s));
            if (mccN(&s) == 1) { I_s = J_s = 0;}
            else { I_s = 1; J_s=mccM(&s)-m_; }
            p_s = mccPR(&s) + (K-1) + mccM(&s) * 0;
            if (mccN(&s) == 1) { I_1s = J_1s = 0;}
            else { I_1s = 1; J_1s=mccM(&s)-m_; }
            p_1s = mccPR(&s) + (K-1) + mccM(&s) * 0;
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_s += J_s, p_1s += J_1s)
               {
                  for (i_=0; i_<m_; ++i_, p_s+=I_s, p_1s+=I_1s, p_RM0_+=I_RM0_)
                  {
                     *p_s = (*p_1s + *p_RM0_);
                  }
               }
            }
         }
         
         /* s(K-1,:) = s(K-1,:) + squeeze( f .* sum( c(K-1,i1,:) .* permute(Y(1,i1)' * comp,[3 1 2]),2) + ... */
         mccFindIndex(&IM1_, &i1);
         mccFindIndex(&IM0_, &IM1_);
         mccFindIndex(&IM2_, &i1);
         mccFindIndex(&IM3_, &IM2_);
         mccFindIndex(&IM4_, &IM3_);
         mccFindIndex(&IM5_, &IM4_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM0_;
            int I_RM0_=1;
            double *p_Y;
            int I_Y=1;
            int *p_IM5_;
            int I_IM5_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM5_) * mccN(&IM5_)));
            mccAllocateMatrix(&RM0_, m_, n_);
            mccCheckMatrixSize(&Y, 1, mccGetMaxIndex(&IM5_ ,mccN(&Y)));
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            I_IM5_ = (mccM(&IM5_) != 1 || mccN(&IM5_) != 1);
            p_IM5_ = mccIPR(&IM5_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM5_ += I_IM5_)
               {
                  p_Y = mccPR(&Y) + mccM(&Y) * ((int)(*p_IM5_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_Y+=I_Y)
                  {
                     *p_RM0_ = *p_Y;
                  }
               }
            }
         }
         mccConjTrans(&RM1_, &RM0_);
         mccRealMatrixMultiply(&RM2_, &RM1_, &M_comp);
         Mprhs_[0] = &RM2_;
         Mprhs_[1] = &S32_;
         Mplhs_[0] = &RM3_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 340);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM4_;
            int I_RM4_=1;
            double *p_c;
            int I_c=1, J_c;
            double *p_RM3_;
            int I_RM3_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            mccAllocateMatrix(&RM4_, m_, n_);
            mccCheckMatrixSize(&c, mccM(&c), mccN(&c));
            I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
            p_RM4_ = mccPR(&RM4_);
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_c = J_c = 0;}
            else { I_c = 1; J_c=mccM(&c)-m_; }
            p_c = mccPR(&c) + 0 + mccM(&c) * 0;
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_c += J_c)
               {
                  for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_c+=I_c, p_RM3_+=I_RM3_)
                  {
                     *p_RM4_ = (*p_c * (double) *p_RM3_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM4_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM5_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 340);
         mccFindIndex(&IM6_, &i2);
         mccFindIndex(&IM7_, &IM6_);
         mccFindIndex(&IM8_, &i2);
         mccFindIndex(&IM9_, &IM8_);
         mccFindIndex(&IM10_, &IM9_);
         mccFindIndex(&IM11_, &IM10_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM6_;
            int I_RM6_=1;
            double *p_Y;
            int I_Y=1;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM11_) * mccN(&IM11_)));
            mccAllocateMatrix(&RM6_, m_, n_);
            mccCheckMatrixSize(&Y, 2, mccGetMaxIndex(&IM11_ ,mccN(&Y)));
            I_RM6_ = (mccM(&RM6_) != 1 || mccN(&RM6_) != 1);
            p_RM6_ = mccPR(&RM6_);
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM11_ += I_IM11_)
               {
                  p_Y = mccPR(&Y) + mccM(&Y) * ((int)(*p_IM11_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM6_+=I_RM6_, p_Y+=I_Y)
                  {
                     *p_RM6_ = *p_Y;
                  }
               }
            }
         }
         mccConjTrans(&RM7_, &RM6_);
         mccRealMatrixMultiply(&RM8_, &RM7_, &M_comp);
         Mprhs_[0] = &RM8_;
         Mprhs_[1] = &S33_;
         Mplhs_[0] = &RM9_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 341);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM10_;
            int I_RM10_=1;
            double *p_cp;
            int I_cp=1, J_cp;
            double *p_RM9_;
            int I_RM9_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&cp), mccN(&cp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM9_), mccN(&RM9_));
            mccAllocateMatrix(&RM10_, m_, n_);
            mccCheckMatrixSize(&cp, mccM(&cp), mccN(&cp));
            I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
            p_RM10_ = mccPR(&RM10_);
            if (mccM(&cp) == 1 && mccN(&cp) == 1) { I_cp = J_cp = 0;}
            else { I_cp = 1; J_cp=mccM(&cp)-m_; }
            p_cp = mccPR(&cp) + 0 + mccM(&cp) * 0;
            I_RM9_ = (mccM(&RM9_) != 1 || mccN(&RM9_) != 1);
            p_RM9_ = mccPR(&RM9_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_cp += J_cp)
               {
                  for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_cp+=I_cp, p_RM9_+=I_RM9_)
                  {
                     *p_RM10_ = (*p_cp * (double) *p_RM9_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM10_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM11_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 341);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM12_;
            int I_RM12_=1;
            double *p_f;
            int I_f=1;
            double *p_RM5_;
            int I_RM5_=1;
            double *p_1f;
            int I_1f=1;
            double *p_RM11_;
            int I_RM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM5_), mccN(&RM5_));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM11_), mccN(&RM11_));
            mccAllocateMatrix(&RM12_, m_, n_);
            I_RM12_ = (mccM(&RM12_) != 1 || mccN(&RM12_) != 1);
            p_RM12_ = mccPR(&RM12_);
            I_f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_f = mccPR(&f);
            I_RM5_ = (mccM(&RM5_) != 1 || mccN(&RM5_) != 1);
            p_RM5_ = mccPR(&RM5_);
            I_1f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_1f = mccPR(&f);
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM12_+=I_RM12_, p_f+=I_f, p_RM5_+=I_RM5_, p_1f+=I_1f, p_RM11_+=I_RM11_)
                  {
                     *p_RM12_ = ((*p_f * (double) *p_RM5_) + ((1 - *p_1f) * (double) *p_RM11_));
                  }
               }
            }
         }
         mccLOG(&RM12_) = 0;
         mccSTRING(&RM12_) = 0;
         Mprhs_[0] = &RM12_;
         Mplhs_[0] = &RM13_;
         mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "squeeze", 341);
         mccConjTrans(&RM14_, &RM13_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_s;
            int I_s=1, J_s;
            double *p_1s;
            int I_1s=1, J_1s;
            double *p_RM14_;
            int I_RM14_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM14_), mccN(&RM14_));
            if (!mccNOTSET(&s)) m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            mccGrowMatrix(&s, (K-1), n_);
            mccCheckMatrixSize(&s, (K-1), mccN(&s));
            if (mccN(&s) == 1) { I_s = J_s = 0;}
            else { I_s = 1; J_s=mccM(&s)-m_; }
            p_s = mccPR(&s) + ((K-1)-1) + mccM(&s) * 0;
            if (mccN(&s) == 1) { I_1s = J_1s = 0;}
            else { I_1s = 1; J_1s=mccM(&s)-m_; }
            p_1s = mccPR(&s) + ((K-1)-1) + mccM(&s) * 0;
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_s += J_s, p_1s += J_1s)
               {
                  for (i_=0; i_<m_; ++i_, p_s+=I_s, p_1s+=I_1s, p_RM14_+=I_RM14_)
                  {
                     *p_s = (*p_1s + *p_RM14_);
                  }
               }
            }
         }
         
         /* F(K,i1)  = X(1,i1) .* F(K-1,i1)  + Y(1,i1); */
         mccFindIndex(&IM10_, &i1);
         mccFindIndex(&IM9_, &i1);
         mccFindIndex(&IM8_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM14_;
            int I_RM14_=1;
            double *p_X;
            int I_X=1;
            int *p_IM10_;
            int I_IM10_=1;
            double *p_F;
            int I_F=1;
            int *p_IM9_;
            int I_IM9_=1;
            double *p_Y;
            int I_Y=1;
            int *p_IM8_;
            int I_IM8_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM10_) * mccN(&IM10_)));
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM9_) * mccN(&IM9_)));
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM8_) * mccN(&IM8_)));
            mccAllocateMatrix(&RM14_, m_, n_);
            mccCheckMatrixSize(&X, 1, mccGetMaxIndex(&IM10_ ,mccN(&X)));
            mccCheckMatrixSize(&F, (K-1), mccGetMaxIndex(&IM9_ ,mccN(&F)));
            mccCheckMatrixSize(&Y, 1, mccGetMaxIndex(&IM8_ ,mccN(&Y)));
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            I_IM10_ = (mccM(&IM10_) != 1 || mccN(&IM10_) != 1);
            p_IM10_ = mccIPR(&IM10_);
            I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
            p_IM9_ = mccIPR(&IM9_);
            I_IM8_ = (mccM(&IM8_) != 1 || mccN(&IM8_) != 1);
            p_IM8_ = mccIPR(&IM8_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM10_ += I_IM10_, p_IM9_ += I_IM9_, p_IM8_ += I_IM8_)
               {
                  p_X = mccPR(&X) + mccM(&X) * ((int)(*p_IM10_ - .5)) + (1-1);
                  p_F = mccPR(&F) + mccM(&F) * ((int)(*p_IM9_ - .5)) + ((K-1)-1);
                  p_Y = mccPR(&Y) + mccM(&Y) * ((int)(*p_IM8_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM14_+=I_RM14_, p_X+=I_X, p_F+=I_F, p_Y+=I_Y)
                  {
                     *p_RM14_ = ((*p_X * (double) *p_F) + *p_Y);
                  }
               }
            }
         }
         mccFindIndex(&IM11_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_F;
            int I_F=1;
            int *p_IM11_;
            int I_IM11_=1;
            double *p_RM14_;
            int I_RM14_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM11_) * mccN(&IM11_)));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM14_), mccN(&RM14_));
            mccGrowMatrix(&F, K, mccGetMaxIndex(&IM11_ ,mccN(&F)));
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM11_ += I_IM11_)
               {
                  p_F = mccPR(&F) + mccM(&F) * ((int)(*p_IM11_ - .5)) + (K-1);
                  for (i_=0; i_<m_; ++i_, p_F+=I_F, p_RM14_+=I_RM14_)
                  {
                     *p_F = *p_RM14_;
                  }
               }
            }
         }
         /* FP(K,i2) = X(2,i2) .* FP(K-1,i2) + Y(2,i2); */
         mccFindIndex(&IM9_, &i2);
         mccFindIndex(&IM10_, &i2);
         mccFindIndex(&IM11_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM14_;
            int I_RM14_=1;
            double *p_X;
            int I_X=1;
            int *p_IM9_;
            int I_IM9_=1;
            double *p_FP;
            int I_FP=1;
            int *p_IM10_;
            int I_IM10_=1;
            double *p_Y;
            int I_Y=1;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM9_) * mccN(&IM9_)));
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM10_) * mccN(&IM10_)));
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM11_) * mccN(&IM11_)));
            mccAllocateMatrix(&RM14_, m_, n_);
            mccCheckMatrixSize(&X, 2, mccGetMaxIndex(&IM9_ ,mccN(&X)));
            mccCheckMatrixSize(&FP, (K-1), mccGetMaxIndex(&IM10_ ,mccN(&FP)));
            mccCheckMatrixSize(&Y, 2, mccGetMaxIndex(&IM11_ ,mccN(&Y)));
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
            p_IM9_ = mccIPR(&IM9_);
            I_IM10_ = (mccM(&IM10_) != 1 || mccN(&IM10_) != 1);
            p_IM10_ = mccIPR(&IM10_);
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM9_ += I_IM9_, p_IM10_ += I_IM10_, p_IM11_ += I_IM11_)
               {
                  p_X = mccPR(&X) + mccM(&X) * ((int)(*p_IM9_ - .5)) + (2-1);
                  p_FP = mccPR(&FP) + mccM(&FP) * ((int)(*p_IM10_ - .5)) + ((K-1)-1);
                  p_Y = mccPR(&Y) + mccM(&Y) * ((int)(*p_IM11_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM14_+=I_RM14_, p_X+=I_X, p_FP+=I_FP, p_Y+=I_Y)
                  {
                     *p_RM14_ = ((*p_X * (double) *p_FP) + *p_Y);
                  }
               }
            }
         }
         mccFindIndex(&IM8_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_FP;
            int I_FP=1;
            int *p_IM8_;
            int I_IM8_=1;
            double *p_RM14_;
            int I_RM14_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM8_) * mccN(&IM8_)));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM14_), mccN(&RM14_));
            mccGrowMatrix(&FP, K, mccGetMaxIndex(&IM8_ ,mccN(&FP)));
            I_IM8_ = (mccM(&IM8_) != 1 || mccN(&IM8_) != 1);
            p_IM8_ = mccIPR(&IM8_);
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM8_ += I_IM8_)
               {
                  p_FP = mccPR(&FP) + mccM(&FP) * ((int)(*p_IM8_ - .5)) + (K-1);
                  for (i_=0; i_<m_; ++i_, p_FP+=I_FP, p_RM14_+=I_RM14_)
                  {
                     *p_FP = *p_RM14_;
                  }
               }
            }
         }
         
         /* c(K-1,i1,:)   = zeros(1,length(i1),M); */
         I0_ = mccGetLength(&i1);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM10_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 346);
         mccFindIndex(&IM11_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_c;
            int I_c=1, J_c;
            int *p_IM10_;
            int I_IM10_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM10_), mccN(&IM10_));
            if (!mccNOTSET(&c)) m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            mccGrowMatrix(&c, m_, n_);
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_c = J_c = 0;}
            else { I_c = 1; J_c=mccM(&c)-m_; }
            p_c = mccPR(&c) + 0 + mccM(&c) * 0;
            I_IM10_ = (mccM(&IM10_) != 1 || mccN(&IM10_) != 1);
            p_IM10_ = mccIPR(&IM10_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_c += J_c)
               {
                  for (i_=0; i_<m_; ++i_, p_c+=I_c, p_IM10_+=I_IM10_)
                  {
                     *p_c = ((int)*p_IM10_);
                  }
               }
            }
         }
         /* cp(K-1,i2,:)  = zeros(1,length(i2),M); */
         I0_ = mccGetLength(&i2);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM11_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 347);
         mccFindIndex(&IM10_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_cp;
            int I_cp=1, J_cp;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM11_), mccN(&IM11_));
            if (!mccNOTSET(&cp)) m_ = mcmCalcResultSize(m_, &n_, mccM(&cp), mccN(&cp));
            mccGrowMatrix(&cp, m_, n_);
            if (mccM(&cp) == 1 && mccN(&cp) == 1) { I_cp = J_cp = 0;}
            else { I_cp = 1; J_cp=mccM(&cp)-m_; }
            p_cp = mccPR(&cp) + 0 + mccM(&cp) * 0;
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_cp += J_cp)
               {
                  for (i_=0; i_<m_; ++i_, p_cp+=I_cp, p_IM11_+=I_IM11_)
                  {
                     *p_cp = ((int)*p_IM11_);
                  }
               }
            }
         }
         /* b(K,i1,:)   = zeros(1,length(i1),M); */
         I0_ = mccGetLength(&i1);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM10_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 348);
         mccFindIndex(&IM11_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_b;
            int I_b=1, J_b;
            int *p_IM10_;
            int I_IM10_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM10_), mccN(&IM10_));
            if (!mccNOTSET(&b)) m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            mccGrowMatrix(&b, m_, n_);
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_b = J_b = 0;}
            else { I_b = 1; J_b=mccM(&b)-m_; }
            p_b = mccPR(&b) + 0 + mccM(&b) * 0;
            I_IM10_ = (mccM(&IM10_) != 1 || mccN(&IM10_) != 1);
            p_IM10_ = mccIPR(&IM10_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_b += J_b)
               {
                  for (i_=0; i_<m_; ++i_, p_b+=I_b, p_IM10_+=I_IM10_)
                  {
                     *p_b = ((int)*p_IM10_);
                  }
               }
            }
         }
         /* bp(K,i2,:)  = zeros(1,length(i2),M); */
         I0_ = mccGetLength(&i2);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM11_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 349);
         mccFindIndex(&IM10_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_bp;
            int I_bp=1, J_bp;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM11_), mccN(&IM11_));
            if (!mccNOTSET(&bp)) m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            mccGrowMatrix(&bp, m_, n_);
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_bp = J_bp = 0;}
            else { I_bp = 1; J_bp=mccM(&bp)-m_; }
            p_bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_bp += J_bp)
               {
                  for (i_=0; i_<m_; ++i_, p_bp+=I_bp, p_IM11_+=I_IM11_)
                  {
                     *p_bp = ((int)*p_IM11_);
                  }
               }
            }
         }
         /* end */
      }
      /* % d - cas u ==0, v ==0 et w ~=0 */
      /* i1 = find(( u(1,:) == 0 ) & (v(1,:) == 0) & (w(1,:) ~= 0)); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         unsigned short *p_BM0_;
         int I_BM0_=1;
         double *p_u;
         int I_u=1, J_u;
         double *p_v;
         int I_v=1, J_v;
         double *p_w;
         int I_w=1, J_w;
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&u));
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&v));
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&w));
         mccAllocateMatrix(&BM0_, m_, n_);
         mccCheckMatrixSize(&u, 1, mccN(&u));
         mccCheckMatrixSize(&v, 1, mccN(&v));
         mccCheckMatrixSize(&w, 1, mccN(&w));
         I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
         p_BM0_ = mccSPR(&BM0_);
         if (mccN(&u) == 1) { I_u = J_u = 0;}
         else { I_u = 1; J_u=mccM(&u)-m_; }
         p_u = mccPR(&u) + (1-1) + mccM(&u) * 0;
         if (mccN(&v) == 1) { I_v = J_v = 0;}
         else { I_v = 1; J_v=mccM(&v)-m_; }
         p_v = mccPR(&v) + (1-1) + mccM(&v) * 0;
         if (mccN(&w) == 1) { I_w = J_w = 0;}
         else { I_w = 1; J_w=mccM(&w)-m_; }
         p_w = mccPR(&w) + (1-1) + mccM(&w) * 0;
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_u += J_u, p_v += J_v, p_w += J_w)
            {
               for (i_=0; i_<m_; ++i_, p_BM0_+=I_BM0_, p_u+=I_u, p_v+=I_v, p_w+=I_w)
               {
                  *p_BM0_ = (!!(!!( (*p_u == 0) && !mccREL_NAN(*p_u) ) && !!( (*p_v == 0) && !mccREL_NAN(*p_v) )) && !!( (*p_w != 0) || mccREL_NAN(*p_w) ));
               }
            }
         }
      }
      mccFind(&i1, &BM0_);
      /* i2 = find(( u(2,:) ==0) & (v(2,:) == 0) & (w(2,:) ~= 0)); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         unsigned short *p_BM0_;
         int I_BM0_=1;
         double *p_u;
         int I_u=1, J_u;
         double *p_v;
         int I_v=1, J_v;
         double *p_w;
         int I_w=1, J_w;
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&u));
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&v));
         m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&w));
         mccAllocateMatrix(&BM0_, m_, n_);
         mccCheckMatrixSize(&u, 2, mccN(&u));
         mccCheckMatrixSize(&v, 2, mccN(&v));
         mccCheckMatrixSize(&w, 2, mccN(&w));
         I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
         p_BM0_ = mccSPR(&BM0_);
         if (mccN(&u) == 1) { I_u = J_u = 0;}
         else { I_u = 1; J_u=mccM(&u)-m_; }
         p_u = mccPR(&u) + (2-1) + mccM(&u) * 0;
         if (mccN(&v) == 1) { I_v = J_v = 0;}
         else { I_v = 1; J_v=mccM(&v)-m_; }
         p_v = mccPR(&v) + (2-1) + mccM(&v) * 0;
         if (mccN(&w) == 1) { I_w = J_w = 0;}
         else { I_w = 1; J_w=mccM(&w)-m_; }
         p_w = mccPR(&w) + (2-1) + mccM(&w) * 0;
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_, p_u += J_u, p_v += J_v, p_w += J_w)
            {
               for (i_=0; i_<m_; ++i_, p_BM0_+=I_BM0_, p_u+=I_u, p_v+=I_v, p_w+=I_w)
               {
                  *p_BM0_ = (!!(!!( (*p_u == 0) && !mccREL_NAN(*p_u) ) && !!( (*p_v == 0) && !mccREL_NAN(*p_v) )) && !!( (*p_w != 0) || mccREL_NAN(*p_w) ));
               }
            }
         }
      }
      mccFind(&i2, &BM0_);
      /* if ~ (isempty(i1) & isempty(i2)) */
      B2_ = mccIsEmpty(&i1);
      B1_ = mccIsEmpty(&i2);
      B0_ = (!(!!B2_ && !!B1_));
      if ((double)B0_)
      {
         /* s(K,:) = s(K,:) + squeeze( f .* sum( a(K,i1,:) .* permute(Z(1,i1)' * comp,[3 1 2]),2) + ... */
         mccFindIndex(&IM11_, &i1);
         mccFindIndex(&IM10_, &IM11_);
         mccFindIndex(&IM9_, &i1);
         mccFindIndex(&IM8_, &IM9_);
         mccFindIndex(&IM7_, &IM8_);
         mccFindIndex(&IM6_, &IM7_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM14_;
            int I_RM14_=1;
            double *p_Z;
            int I_Z=1;
            int *p_IM6_;
            int I_IM6_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM6_) * mccN(&IM6_)));
            mccAllocateMatrix(&RM14_, m_, n_);
            mccCheckMatrixSize(&Z, 1, mccGetMaxIndex(&IM6_ ,mccN(&Z)));
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            I_IM6_ = (mccM(&IM6_) != 1 || mccN(&IM6_) != 1);
            p_IM6_ = mccIPR(&IM6_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM6_ += I_IM6_)
               {
                  p_Z = mccPR(&Z) + mccM(&Z) * ((int)(*p_IM6_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM14_+=I_RM14_, p_Z+=I_Z)
                  {
                     *p_RM14_ = *p_Z;
                  }
               }
            }
         }
         mccConjTrans(&RM13_, &RM14_);
         mccRealMatrixMultiply(&RM12_, &RM13_, &M_comp);
         Mprhs_[0] = &RM12_;
         Mprhs_[1] = &S34_;
         Mplhs_[0] = &RM11_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 355);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM10_;
            int I_RM10_=1;
            double *p_a;
            int I_a=1, J_a;
            double *p_RM11_;
            int I_RM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&a), mccN(&a));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM11_), mccN(&RM11_));
            mccAllocateMatrix(&RM10_, m_, n_);
            mccCheckMatrixSize(&a, mccM(&a), mccN(&a));
            I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
            p_RM10_ = mccPR(&RM10_);
            if (mccM(&a) == 1 && mccN(&a) == 1) { I_a = J_a = 0;}
            else { I_a = 1; J_a=mccM(&a)-m_; }
            p_a = mccPR(&a) + 0 + mccM(&a) * 0;
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_a += J_a)
               {
                  for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_a+=I_a, p_RM11_+=I_RM11_)
                  {
                     *p_RM10_ = (*p_a * (double) *p_RM11_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM10_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM9_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 355);
         mccFindIndex(&IM5_, &i2);
         mccFindIndex(&IM4_, &IM5_);
         mccFindIndex(&IM3_, &i2);
         mccFindIndex(&IM2_, &IM3_);
         mccFindIndex(&IM0_, &IM2_);
         mccFindIndex(&IM1_, &IM0_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM8_;
            int I_RM8_=1;
            double *p_Z;
            int I_Z=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM1_) * mccN(&IM1_)));
            mccAllocateMatrix(&RM8_, m_, n_);
            mccCheckMatrixSize(&Z, 2, mccGetMaxIndex(&IM1_ ,mccN(&Z)));
            I_RM8_ = (mccM(&RM8_) != 1 || mccN(&RM8_) != 1);
            p_RM8_ = mccPR(&RM8_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_Z = mccPR(&Z) + mccM(&Z) * ((int)(*p_IM1_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM8_+=I_RM8_, p_Z+=I_Z)
                  {
                     *p_RM8_ = *p_Z;
                  }
               }
            }
         }
         mccConjTrans(&RM7_, &RM8_);
         mccRealMatrixMultiply(&RM6_, &RM7_, &M_comp);
         Mprhs_[0] = &RM6_;
         Mprhs_[1] = &S35_;
         Mplhs_[0] = &RM5_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 356);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM4_;
            int I_RM4_=1;
            double *p_ap;
            int I_ap=1, J_ap;
            double *p_RM5_;
            int I_RM5_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), mccN(&ap));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM5_), mccN(&RM5_));
            mccAllocateMatrix(&RM4_, m_, n_);
            mccCheckMatrixSize(&ap, mccM(&ap), mccN(&ap));
            I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
            p_RM4_ = mccPR(&RM4_);
            if (mccM(&ap) == 1 && mccN(&ap) == 1) { I_ap = J_ap = 0;}
            else { I_ap = 1; J_ap=mccM(&ap)-m_; }
            p_ap = mccPR(&ap) + 0 + mccM(&ap) * 0;
            I_RM5_ = (mccM(&RM5_) != 1 || mccN(&RM5_) != 1);
            p_RM5_ = mccPR(&RM5_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_ap += J_ap)
               {
                  for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_ap+=I_ap, p_RM5_+=I_RM5_)
                  {
                     *p_RM4_ = (*p_ap * (double) *p_RM5_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM4_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM3_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 356);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM2_;
            int I_RM2_=1;
            double *p_f;
            int I_f=1;
            double *p_RM9_;
            int I_RM9_=1;
            double *p_1f;
            int I_1f=1;
            double *p_RM3_;
            int I_RM3_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM9_), mccN(&RM9_));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            mccAllocateMatrix(&RM2_, m_, n_);
            I_RM2_ = (mccM(&RM2_) != 1 || mccN(&RM2_) != 1);
            p_RM2_ = mccPR(&RM2_);
            I_f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_f = mccPR(&f);
            I_RM9_ = (mccM(&RM9_) != 1 || mccN(&RM9_) != 1);
            p_RM9_ = mccPR(&RM9_);
            I_1f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_1f = mccPR(&f);
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM2_+=I_RM2_, p_f+=I_f, p_RM9_+=I_RM9_, p_1f+=I_1f, p_RM3_+=I_RM3_)
                  {
                     *p_RM2_ = ((*p_f * (double) *p_RM9_) + ((1 - *p_1f) * (double) *p_RM3_));
                  }
               }
            }
         }
         mccLOG(&RM2_) = 0;
         mccSTRING(&RM2_) = 0;
         Mprhs_[0] = &RM2_;
         Mplhs_[0] = &RM1_;
         mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "squeeze", 356);
         mccConjTrans(&RM0_, &RM1_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_s;
            int I_s=1, J_s;
            double *p_1s;
            int I_1s=1, J_1s;
            double *p_RM0_;
            int I_RM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM0_), mccN(&RM0_));
            if (!mccNOTSET(&s)) m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            mccGrowMatrix(&s, K, n_);
            mccCheckMatrixSize(&s, K, mccN(&s));
            if (mccN(&s) == 1) { I_s = J_s = 0;}
            else { I_s = 1; J_s=mccM(&s)-m_; }
            p_s = mccPR(&s) + (K-1) + mccM(&s) * 0;
            if (mccN(&s) == 1) { I_1s = J_1s = 0;}
            else { I_1s = 1; J_1s=mccM(&s)-m_; }
            p_1s = mccPR(&s) + (K-1) + mccM(&s) * 0;
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_s += J_s, p_1s += J_1s)
               {
                  for (i_=0; i_<m_; ++i_, p_s+=I_s, p_1s+=I_1s, p_RM0_+=I_RM0_)
                  {
                     *p_s = (*p_1s + *p_RM0_);
                  }
               }
            }
         }
         
         /* s(K-1,:) = s(K-1,:) + squeeze( f .* sum( b(K-1,i1,:) .* permute(Z(1,i1)' * comp,[3 1 2]),2) + ... */
         mccFindIndex(&IM1_, &i1);
         mccFindIndex(&IM0_, &IM1_);
         mccFindIndex(&IM2_, &i1);
         mccFindIndex(&IM3_, &IM2_);
         mccFindIndex(&IM4_, &IM3_);
         mccFindIndex(&IM5_, &IM4_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM0_;
            int I_RM0_=1;
            double *p_Z;
            int I_Z=1;
            int *p_IM5_;
            int I_IM5_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM5_) * mccN(&IM5_)));
            mccAllocateMatrix(&RM0_, m_, n_);
            mccCheckMatrixSize(&Z, 1, mccGetMaxIndex(&IM5_ ,mccN(&Z)));
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            I_IM5_ = (mccM(&IM5_) != 1 || mccN(&IM5_) != 1);
            p_IM5_ = mccIPR(&IM5_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM5_ += I_IM5_)
               {
                  p_Z = mccPR(&Z) + mccM(&Z) * ((int)(*p_IM5_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_Z+=I_Z)
                  {
                     *p_RM0_ = *p_Z;
                  }
               }
            }
         }
         mccConjTrans(&RM1_, &RM0_);
         mccRealMatrixMultiply(&RM2_, &RM1_, &M_comp);
         Mprhs_[0] = &RM2_;
         Mprhs_[1] = &S36_;
         Mplhs_[0] = &RM3_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 358);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM4_;
            int I_RM4_=1;
            double *p_b;
            int I_b=1, J_b;
            double *p_RM3_;
            int I_RM3_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            mccAllocateMatrix(&RM4_, m_, n_);
            mccCheckMatrixSize(&b, mccM(&b), mccN(&b));
            I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
            p_RM4_ = mccPR(&RM4_);
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_b = J_b = 0;}
            else { I_b = 1; J_b=mccM(&b)-m_; }
            p_b = mccPR(&b) + 0 + mccM(&b) * 0;
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_b += J_b)
               {
                  for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_b+=I_b, p_RM3_+=I_RM3_)
                  {
                     *p_RM4_ = (*p_b * (double) *p_RM3_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM4_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM5_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 358);
         mccFindIndex(&IM6_, &i2);
         mccFindIndex(&IM7_, &IM6_);
         mccFindIndex(&IM8_, &i2);
         mccFindIndex(&IM9_, &IM8_);
         mccFindIndex(&IM10_, &IM9_);
         mccFindIndex(&IM11_, &IM10_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM6_;
            int I_RM6_=1;
            double *p_Z;
            int I_Z=1;
            int *p_IM11_;
            int I_IM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM11_) * mccN(&IM11_)));
            mccAllocateMatrix(&RM6_, m_, n_);
            mccCheckMatrixSize(&Z, 2, mccGetMaxIndex(&IM11_ ,mccN(&Z)));
            I_RM6_ = (mccM(&RM6_) != 1 || mccN(&RM6_) != 1);
            p_RM6_ = mccPR(&RM6_);
            I_IM11_ = (mccM(&IM11_) != 1 || mccN(&IM11_) != 1);
            p_IM11_ = mccIPR(&IM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM11_ += I_IM11_)
               {
                  p_Z = mccPR(&Z) + mccM(&Z) * ((int)(*p_IM11_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM6_+=I_RM6_, p_Z+=I_Z)
                  {
                     *p_RM6_ = *p_Z;
                  }
               }
            }
         }
         mccConjTrans(&RM7_, &RM6_);
         mccRealMatrixMultiply(&RM8_, &RM7_, &M_comp);
         Mprhs_[0] = &RM8_;
         Mprhs_[1] = &S37_;
         Mplhs_[0] = &RM9_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 359);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM10_;
            int I_RM10_=1;
            double *p_bp;
            int I_bp=1, J_bp;
            double *p_RM9_;
            int I_RM9_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM9_), mccN(&RM9_));
            mccAllocateMatrix(&RM10_, m_, n_);
            mccCheckMatrixSize(&bp, mccM(&bp), mccN(&bp));
            I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
            p_RM10_ = mccPR(&RM10_);
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_bp = J_bp = 0;}
            else { I_bp = 1; J_bp=mccM(&bp)-m_; }
            p_bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            I_RM9_ = (mccM(&RM9_) != 1 || mccN(&RM9_) != 1);
            p_RM9_ = mccPR(&RM9_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_bp += J_bp)
               {
                  for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_bp+=I_bp, p_RM9_+=I_RM9_)
                  {
                     *p_RM10_ = (*p_bp * (double) *p_RM9_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM10_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM11_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 359);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM12_;
            int I_RM12_=1;
            double *p_f;
            int I_f=1;
            double *p_RM5_;
            int I_RM5_=1;
            double *p_1f;
            int I_1f=1;
            double *p_RM11_;
            int I_RM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM5_), mccN(&RM5_));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM11_), mccN(&RM11_));
            mccAllocateMatrix(&RM12_, m_, n_);
            I_RM12_ = (mccM(&RM12_) != 1 || mccN(&RM12_) != 1);
            p_RM12_ = mccPR(&RM12_);
            I_f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_f = mccPR(&f);
            I_RM5_ = (mccM(&RM5_) != 1 || mccN(&RM5_) != 1);
            p_RM5_ = mccPR(&RM5_);
            I_1f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_1f = mccPR(&f);
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM12_+=I_RM12_, p_f+=I_f, p_RM5_+=I_RM5_, p_1f+=I_1f, p_RM11_+=I_RM11_)
                  {
                     *p_RM12_ = ((*p_f * (double) *p_RM5_) + ((1 - *p_1f) * (double) *p_RM11_));
                  }
               }
            }
         }
         mccLOG(&RM12_) = 0;
         mccSTRING(&RM12_) = 0;
         Mprhs_[0] = &RM12_;
         Mplhs_[0] = &RM13_;
         mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "squeeze", 359);
         mccConjTrans(&RM14_, &RM13_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_s;
            int I_s=1, J_s;
            double *p_1s;
            int I_1s=1, J_1s;
            double *p_RM14_;
            int I_RM14_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM14_), mccN(&RM14_));
            if (!mccNOTSET(&s)) m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            mccGrowMatrix(&s, (K-1), n_);
            mccCheckMatrixSize(&s, (K-1), mccN(&s));
            if (mccN(&s) == 1) { I_s = J_s = 0;}
            else { I_s = 1; J_s=mccM(&s)-m_; }
            p_s = mccPR(&s) + ((K-1)-1) + mccM(&s) * 0;
            if (mccN(&s) == 1) { I_1s = J_1s = 0;}
            else { I_1s = 1; J_1s=mccM(&s)-m_; }
            p_1s = mccPR(&s) + ((K-1)-1) + mccM(&s) * 0;
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_s += J_s, p_1s += J_1s)
               {
                  for (i_=0; i_<m_; ++i_, p_s+=I_s, p_1s+=I_1s, p_RM14_+=I_RM14_)
                  {
                     *p_s = (*p_1s + *p_RM14_);
                  }
               }
            }
         }
         
         /* s(K-2,:) = s(K-2,:) + squeeze( f .* sum( c(K-2,i1,:) .* permute(Z(1,i1)' * comp,[3 1 2]),2) + ... */
         mccFindIndex(&IM11_, &i1);
         mccFindIndex(&IM10_, &IM11_);
         mccFindIndex(&IM9_, &i1);
         mccFindIndex(&IM8_, &IM9_);
         mccFindIndex(&IM7_, &IM8_);
         mccFindIndex(&IM6_, &IM7_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM14_;
            int I_RM14_=1;
            double *p_Z;
            int I_Z=1;
            int *p_IM6_;
            int I_IM6_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM6_) * mccN(&IM6_)));
            mccAllocateMatrix(&RM14_, m_, n_);
            mccCheckMatrixSize(&Z, 1, mccGetMaxIndex(&IM6_ ,mccN(&Z)));
            I_RM14_ = (mccM(&RM14_) != 1 || mccN(&RM14_) != 1);
            p_RM14_ = mccPR(&RM14_);
            I_IM6_ = (mccM(&IM6_) != 1 || mccN(&IM6_) != 1);
            p_IM6_ = mccIPR(&IM6_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM6_ += I_IM6_)
               {
                  p_Z = mccPR(&Z) + mccM(&Z) * ((int)(*p_IM6_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_RM14_+=I_RM14_, p_Z+=I_Z)
                  {
                     *p_RM14_ = *p_Z;
                  }
               }
            }
         }
         mccConjTrans(&RM13_, &RM14_);
         mccRealMatrixMultiply(&RM12_, &RM13_, &M_comp);
         Mprhs_[0] = &RM12_;
         Mprhs_[1] = &S38_;
         Mplhs_[0] = &RM11_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 361);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM10_;
            int I_RM10_=1;
            double *p_c;
            int I_c=1, J_c;
            double *p_RM11_;
            int I_RM11_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM11_), mccN(&RM11_));
            mccAllocateMatrix(&RM10_, m_, n_);
            mccCheckMatrixSize(&c, mccM(&c), mccN(&c));
            I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
            p_RM10_ = mccPR(&RM10_);
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_c = J_c = 0;}
            else { I_c = 1; J_c=mccM(&c)-m_; }
            p_c = mccPR(&c) + 0 + mccM(&c) * 0;
            I_RM11_ = (mccM(&RM11_) != 1 || mccN(&RM11_) != 1);
            p_RM11_ = mccPR(&RM11_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_c += J_c)
               {
                  for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_c+=I_c, p_RM11_+=I_RM11_)
                  {
                     *p_RM10_ = (*p_c * (double) *p_RM11_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM10_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM9_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 361);
         mccFindIndex(&IM5_, &i2);
         mccFindIndex(&IM4_, &IM5_);
         mccFindIndex(&IM3_, &i2);
         mccFindIndex(&IM2_, &IM3_);
         mccFindIndex(&IM0_, &IM2_);
         mccFindIndex(&IM1_, &IM0_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM8_;
            int I_RM8_=1;
            double *p_Z;
            int I_Z=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM1_) * mccN(&IM1_)));
            mccAllocateMatrix(&RM8_, m_, n_);
            mccCheckMatrixSize(&Z, 2, mccGetMaxIndex(&IM1_ ,mccN(&Z)));
            I_RM8_ = (mccM(&RM8_) != 1 || mccN(&RM8_) != 1);
            p_RM8_ = mccPR(&RM8_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_Z = mccPR(&Z) + mccM(&Z) * ((int)(*p_IM1_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_RM8_+=I_RM8_, p_Z+=I_Z)
                  {
                     *p_RM8_ = *p_Z;
                  }
               }
            }
         }
         mccConjTrans(&RM7_, &RM8_);
         mccRealMatrixMultiply(&RM6_, &RM7_, &M_comp);
         Mprhs_[0] = &RM6_;
         Mprhs_[1] = &S39_;
         Mplhs_[0] = &RM5_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "permute", 362);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM4_;
            int I_RM4_=1;
            double *p_cp;
            int I_cp=1, J_cp;
            double *p_RM5_;
            int I_RM5_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&cp), mccN(&cp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM5_), mccN(&RM5_));
            mccAllocateMatrix(&RM4_, m_, n_);
            mccCheckMatrixSize(&cp, mccM(&cp), mccN(&cp));
            I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
            p_RM4_ = mccPR(&RM4_);
            if (mccM(&cp) == 1 && mccN(&cp) == 1) { I_cp = J_cp = 0;}
            else { I_cp = 1; J_cp=mccM(&cp)-m_; }
            p_cp = mccPR(&cp) + 0 + mccM(&cp) * 0;
            I_RM5_ = (mccM(&RM5_) != 1 || mccN(&RM5_) != 1);
            p_RM5_ = mccPR(&RM5_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_cp += J_cp)
               {
                  for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_cp+=I_cp, p_RM5_+=I_RM5_)
                  {
                     *p_RM4_ = (*p_cp * (double) *p_RM5_);
                  }
               }
            }
         }
         Mprhs_[0] = &RM4_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM3_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 362);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM2_;
            int I_RM2_=1;
            double *p_f;
            int I_f=1;
            double *p_RM9_;
            int I_RM9_=1;
            double *p_1f;
            int I_1f=1;
            double *p_RM3_;
            int I_RM3_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM9_), mccN(&RM9_));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            mccAllocateMatrix(&RM2_, m_, n_);
            I_RM2_ = (mccM(&RM2_) != 1 || mccN(&RM2_) != 1);
            p_RM2_ = mccPR(&RM2_);
            I_f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_f = mccPR(&f);
            I_RM9_ = (mccM(&RM9_) != 1 || mccN(&RM9_) != 1);
            p_RM9_ = mccPR(&RM9_);
            I_1f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_1f = mccPR(&f);
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM2_+=I_RM2_, p_f+=I_f, p_RM9_+=I_RM9_, p_1f+=I_1f, p_RM3_+=I_RM3_)
                  {
                     *p_RM2_ = ((*p_f * (double) *p_RM9_) + ((1 - *p_1f) * (double) *p_RM3_));
                  }
               }
            }
         }
         mccLOG(&RM2_) = 0;
         mccSTRING(&RM2_) = 0;
         Mprhs_[0] = &RM2_;
         Mplhs_[0] = &RM1_;
         mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "squeeze", 362);
         mccConjTrans(&RM0_, &RM1_);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_s;
            int I_s=1, J_s;
            double *p_1s;
            int I_1s=1, J_1s;
            double *p_RM0_;
            int I_RM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM0_), mccN(&RM0_));
            if (!mccNOTSET(&s)) m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&s));
            mccGrowMatrix(&s, (K-2), n_);
            mccCheckMatrixSize(&s, (K-2), mccN(&s));
            if (mccN(&s) == 1) { I_s = J_s = 0;}
            else { I_s = 1; J_s=mccM(&s)-m_; }
            p_s = mccPR(&s) + ((K-2)-1) + mccM(&s) * 0;
            if (mccN(&s) == 1) { I_1s = J_1s = 0;}
            else { I_1s = 1; J_1s=mccM(&s)-m_; }
            p_1s = mccPR(&s) + ((K-2)-1) + mccM(&s) * 0;
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_s += J_s, p_1s += J_1s)
               {
                  for (i_=0; i_<m_; ++i_, p_s+=I_s, p_1s+=I_1s, p_RM0_+=I_RM0_)
                  {
                     *p_s = (*p_1s + *p_RM0_);
                  }
               }
            }
         }
         
         /* F(K-1,i1)  = Z(1,i1);                      */
         mccFindIndex(&IM0_, &i1);
         mccFindIndex(&IM1_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_F;
            int I_F=1;
            int *p_IM1_;
            int I_IM1_=1;
            double *p_Z;
            int I_Z=1;
            int *p_IM0_;
            int I_IM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM1_) * mccN(&IM1_)));
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM0_) * mccN(&IM0_)));
            mccGrowMatrix(&F, (K-1), mccGetMaxIndex(&IM1_ ,mccN(&F)));
            mccCheckMatrixSize(&Z, 1, mccGetMaxIndex(&IM0_ ,mccN(&Z)));
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_, p_IM0_ += I_IM0_)
               {
                  p_F = mccPR(&F) + mccM(&F) * ((int)(*p_IM1_ - .5)) + ((K-1)-1);
                  p_Z = mccPR(&Z) + mccM(&Z) * ((int)(*p_IM0_ - .5)) + (1-1);
                  for (i_=0; i_<m_; ++i_, p_F+=I_F, p_Z+=I_Z)
                  {
                     *p_F = *p_Z;
                  }
               }
            }
         }
         /* FP(K-1,i2) = Z(2,i2);    */
         mccFindIndex(&IM1_, &i2);
         mccFindIndex(&IM0_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_FP;
            int I_FP=1;
            int *p_IM0_;
            int I_IM0_=1;
            double *p_Z;
            int I_Z=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM0_) * mccN(&IM0_)));
            m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM1_) * mccN(&IM1_)));
            mccGrowMatrix(&FP, (K-1), mccGetMaxIndex(&IM0_ ,mccN(&FP)));
            mccCheckMatrixSize(&Z, 2, mccGetMaxIndex(&IM1_ ,mccN(&Z)));
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM0_ += I_IM0_, p_IM1_ += I_IM1_)
               {
                  p_FP = mccPR(&FP) + mccM(&FP) * ((int)(*p_IM0_ - .5)) + ((K-1)-1);
                  p_Z = mccPR(&Z) + mccM(&Z) * ((int)(*p_IM1_ - .5)) + (2-1);
                  for (i_=0; i_<m_; ++i_, p_FP+=I_FP, p_Z+=I_Z)
                  {
                     *p_FP = *p_Z;
                  }
               }
            }
         }
         
         /* a(K,i1,:)     = zeros(1,length(i1),M); */
         I0_ = mccGetLength(&i1);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM0_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 367);
         mccFindIndex(&IM1_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_a;
            int I_a=1, J_a;
            int *p_IM0_;
            int I_IM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM0_), mccN(&IM0_));
            if (!mccNOTSET(&a)) m_ = mcmCalcResultSize(m_, &n_, mccM(&a), mccN(&a));
            mccGrowMatrix(&a, m_, n_);
            if (mccM(&a) == 1 && mccN(&a) == 1) { I_a = J_a = 0;}
            else { I_a = 1; J_a=mccM(&a)-m_; }
            p_a = mccPR(&a) + 0 + mccM(&a) * 0;
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_a += J_a)
               {
                  for (i_=0; i_<m_; ++i_, p_a+=I_a, p_IM0_+=I_IM0_)
                  {
                     *p_a = ((int)*p_IM0_);
                  }
               }
            }
         }
         /* ap(K,i2,:)    = zeros(1,length(i2),M); */
         I0_ = mccGetLength(&i2);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM1_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 368);
         mccFindIndex(&IM0_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_ap;
            int I_ap=1, J_ap;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM1_), mccN(&IM1_));
            if (!mccNOTSET(&ap)) m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), mccN(&ap));
            mccGrowMatrix(&ap, m_, n_);
            if (mccM(&ap) == 1 && mccN(&ap) == 1) { I_ap = J_ap = 0;}
            else { I_ap = 1; J_ap=mccM(&ap)-m_; }
            p_ap = mccPR(&ap) + 0 + mccM(&ap) * 0;
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_ap += J_ap)
               {
                  for (i_=0; i_<m_; ++i_, p_ap+=I_ap, p_IM1_+=I_IM1_)
                  {
                     *p_ap = ((int)*p_IM1_);
                  }
               }
            }
         }
         /* b(K-1,i1,:)   = zeros(1,length(i1),M); */
         I0_ = mccGetLength(&i1);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM0_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 369);
         mccFindIndex(&IM1_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_b;
            int I_b=1, J_b;
            int *p_IM0_;
            int I_IM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM0_), mccN(&IM0_));
            if (!mccNOTSET(&b)) m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            mccGrowMatrix(&b, m_, n_);
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_b = J_b = 0;}
            else { I_b = 1; J_b=mccM(&b)-m_; }
            p_b = mccPR(&b) + 0 + mccM(&b) * 0;
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_b += J_b)
               {
                  for (i_=0; i_<m_; ++i_, p_b+=I_b, p_IM0_+=I_IM0_)
                  {
                     *p_b = ((int)*p_IM0_);
                  }
               }
            }
         }
         /* bp(K-1,i2,:)  = zeros(1,length(i2),M); */
         I0_ = mccGetLength(&i2);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM1_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 370);
         mccFindIndex(&IM0_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_bp;
            int I_bp=1, J_bp;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM1_), mccN(&IM1_));
            if (!mccNOTSET(&bp)) m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            mccGrowMatrix(&bp, m_, n_);
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_bp = J_bp = 0;}
            else { I_bp = 1; J_bp=mccM(&bp)-m_; }
            p_bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_bp += J_bp)
               {
                  for (i_=0; i_<m_; ++i_, p_bp+=I_bp, p_IM1_+=I_IM1_)
                  {
                     *p_bp = ((int)*p_IM1_);
                  }
               }
            }
         }
         /* c(K-2,i1,:)   = zeros(1,length(i1),M); */
         I0_ = mccGetLength(&i1);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM0_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 371);
         mccFindIndex(&IM1_, &i1);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_c;
            int I_c=1, J_c;
            int *p_IM0_;
            int I_IM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM0_), mccN(&IM0_));
            if (!mccNOTSET(&c)) m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            mccGrowMatrix(&c, m_, n_);
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_c = J_c = 0;}
            else { I_c = 1; J_c=mccM(&c)-m_; }
            p_c = mccPR(&c) + 0 + mccM(&c) * 0;
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_c += J_c)
               {
                  for (i_=0; i_<m_; ++i_, p_c+=I_c, p_IM0_+=I_IM0_)
                  {
                     *p_c = ((int)*p_IM0_);
                  }
               }
            }
         }
         /* cp(K-2,i2,:)  = zeros(1,length(i2),M); */
         I0_ = mccGetLength(&i2);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM1_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 372);
         mccFindIndex(&IM0_, &i2);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_cp;
            int I_cp=1, J_cp;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&IM1_), mccN(&IM1_));
            if (!mccNOTSET(&cp)) m_ = mcmCalcResultSize(m_, &n_, mccM(&cp), mccN(&cp));
            mccGrowMatrix(&cp, m_, n_);
            if (mccM(&cp) == 1 && mccN(&cp) == 1) { I_cp = J_cp = 0;}
            else { I_cp = 1; J_cp=mccM(&cp)-m_; }
            p_cp = mccPR(&cp) + 0 + mccM(&cp) * 0;
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_cp += J_cp)
               {
                  for (i_=0; i_<m_; ++i_, p_cp+=I_cp, p_IM1_+=I_IM1_)
                  {
                     *p_cp = ((int)*p_IM1_);
                  }
               }
            }
         }
         /* end */
      }
      /* % 4 - equation predictive et interpretative */
      /* n = find(mode);            % c'est du matlab 5 ! */
      mccFind(&n, &mode);
      /* if ~isempty(n) */
      B0_ = mccIsEmpty(&n);
      B1_ = (!B0_);
      if ((double)B1_)
      {
         
         /* % decalage d'indice pour ce terme -1 */
         /* Fa=F(1:(K-1),n);   */
         mccFindIndex(&IM1_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_Fa;
            int I_Fa=1;
            double *p_F;
            int I_F=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, ((int)((K - 1) - 1) + 1), (mccM(&IM1_) * mccN(&IM1_)));
            mccAllocateMatrix(&Fa, m_, n_);
            mccCheckMatrixSize(&F, (K - 1), mccGetMaxIndex(&IM1_ ,mccN(&F)));
            I_Fa = (mccM(&Fa) != 1 || mccN(&Fa) != 1);
            p_Fa = mccPR(&Fa);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_F = mccPR(&F) + mccM(&F) * ((int)(*p_IM1_ - .5)) + ((int)(1 - .5));
                  for (i_=0; i_<m_; ++i_, p_Fa+=I_Fa, p_F+=I_F)
                  {
                     *p_Fa = *p_F;
                  }
               }
            }
         }
         /* FPa=FP(1:(K-1),n);   */
         mccFindIndex(&IM1_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_FPa;
            int I_FPa=1;
            double *p_FP;
            int I_FP=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, ((int)((K - 1) - 1) + 1), (mccM(&IM1_) * mccN(&IM1_)));
            mccAllocateMatrix(&FPa, m_, n_);
            mccCheckMatrixSize(&FP, (K - 1), mccGetMaxIndex(&IM1_ ,mccN(&FP)));
            I_FPa = (mccM(&FPa) != 1 || mccN(&FPa) != 1);
            p_FPa = mccPR(&FPa);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_FP = mccPR(&FP) + mccM(&FP) * ((int)(*p_IM1_ - .5)) + ((int)(1 - .5));
                  for (i_=0; i_<m_; ++i_, p_FPa+=I_FPa, p_FP+=I_FP)
                  {
                     *p_FPa = *p_FP;
                  }
               }
            }
         }
         /* % on cree un matrice de la forme [F,F,F...F] (M fois)  */
         /* % sur la troisieme dimension */
         /* Fa=reshape(Fa(:)*ones(1,M),K-1,length(n),M);      */
         mccCopy(&RM0_, &Fa);
         mxSetM( &RM0_, mccM(&RM0_) * mccN(&RM0_));
         mxSetN( &RM0_, 1);
         mccOnesMN(&IM1_, 1, M);
         mccRealMatrixMultiply(&RM1_, &RM0_, &IM1_);
         I0_ = mccGetLength(&n);
         Mprhs_[0] = &RM1_;
         Mprhs_[1] = mccTempMatrix((K - 1), 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[3] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &Fa;
         mccCallMATLAB(1, Mplhs_, 4, Mprhs_, "reshape", 383);
         /* FPa=reshape(FPa(:)*ones(1,M),K-1,length(n),M);      */
         mccCopy(&RM1_, &FPa);
         mxSetM( &RM1_, mccM(&RM1_) * mccN(&RM1_));
         mxSetN( &RM1_, 1);
         mccOnesMN(&IM1_, 1, M);
         mccRealMatrixMultiply(&RM0_, &RM1_, &IM1_);
         I0_ = mccGetLength(&n);
         Mprhs_[0] = &RM0_;
         Mprhs_[1] = mccTempMatrix((K - 1), 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[3] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &FPa;
         mccCallMATLAB(1, Mplhs_, 4, Mprhs_, "reshape", 384);
         
         /* % decalage d'indice pour ce terme 0 */
         /* Fb=F(:,n);   */
         mccFindIndex(&IM1_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_Fb;
            int I_Fb=1;
            double *p_F;
            int I_F=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&F), (mccM(&IM1_) * mccN(&IM1_)));
            mccAllocateMatrix(&Fb, m_, n_);
            mccCheckMatrixSize(&F, mccM(&F), mccGetMaxIndex(&IM1_ ,mccN(&F)));
            I_Fb = (mccM(&Fb) != 1 || mccN(&Fb) != 1);
            p_Fb = mccPR(&Fb);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_F = mccPR(&F) + mccM(&F) * ((int)(*p_IM1_ - .5)) + 0;
                  for (i_=0; i_<m_; ++i_, p_Fb+=I_Fb, p_F+=I_F)
                  {
                     *p_Fb = *p_F;
                  }
               }
            }
         }
         /* FPb=FP(:,n);   */
         mccFindIndex(&IM1_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_FPb;
            int I_FPb=1;
            double *p_FP;
            int I_FP=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&FP), (mccM(&IM1_) * mccN(&IM1_)));
            mccAllocateMatrix(&FPb, m_, n_);
            mccCheckMatrixSize(&FP, mccM(&FP), mccGetMaxIndex(&IM1_ ,mccN(&FP)));
            I_FPb = (mccM(&FPb) != 1 || mccN(&FPb) != 1);
            p_FPb = mccPR(&FPb);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_FP = mccPR(&FP) + mccM(&FP) * ((int)(*p_IM1_ - .5)) + 0;
                  for (i_=0; i_<m_; ++i_, p_FPb+=I_FPb, p_FP+=I_FP)
                  {
                     *p_FPb = *p_FP;
                  }
               }
            }
         }
         /* % on cree un matrice de la forme [F,F,F...F] (M fois)  */
         /* % sur la troisieme dimension */
         /* Fb=reshape(Fb(:)*ones(1,M),K,length(n),M);      */
         mccCopy(&RM0_, &Fb);
         mxSetM( &RM0_, mccM(&RM0_) * mccN(&RM0_));
         mxSetN( &RM0_, 1);
         mccOnesMN(&IM1_, 1, M);
         mccRealMatrixMultiply(&RM1_, &RM0_, &IM1_);
         I0_ = mccGetLength(&n);
         Mprhs_[0] = &RM1_;
         Mprhs_[1] = mccTempMatrix(K, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[3] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &Fb;
         mccCallMATLAB(1, Mplhs_, 4, Mprhs_, "reshape", 391);
         /* FPb=reshape(FPb(:)*ones(1,M),K,length(n),M);      */
         mccCopy(&RM1_, &FPb);
         mxSetM( &RM1_, mccM(&RM1_) * mccN(&RM1_));
         mxSetN( &RM1_, 1);
         mccOnesMN(&IM1_, 1, M);
         mccRealMatrixMultiply(&RM0_, &RM1_, &IM1_);
         I0_ = mccGetLength(&n);
         Mprhs_[0] = &RM0_;
         Mprhs_[1] = mccTempMatrix(K, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[3] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &FPb;
         mccCallMATLAB(1, Mplhs_, 4, Mprhs_, "reshape", 392);
         
         /* % decalage d'indice pour ce terme +1 */
         /* Fc=F(2:K,n);   */
         mccFindIndex(&IM1_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_Fc;
            int I_Fc=1;
            double *p_F;
            int I_F=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, ((int)(K - 2) + 1), (mccM(&IM1_) * mccN(&IM1_)));
            mccAllocateMatrix(&Fc, m_, n_);
            mccCheckMatrixSize(&F, K, mccGetMaxIndex(&IM1_ ,mccN(&F)));
            I_Fc = (mccM(&Fc) != 1 || mccN(&Fc) != 1);
            p_Fc = mccPR(&Fc);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_F = mccPR(&F) + mccM(&F) * ((int)(*p_IM1_ - .5)) + ((int)(2 - .5));
                  for (i_=0; i_<m_; ++i_, p_Fc+=I_Fc, p_F+=I_F)
                  {
                     *p_Fc = *p_F;
                  }
               }
            }
         }
         /* FPc=FP(2:K,n);   */
         mccFindIndex(&IM1_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_FPc;
            int I_FPc=1;
            double *p_FP;
            int I_FP=1;
            int *p_IM1_;
            int I_IM1_=1;
            m_ = mcmCalcResultSize(m_, &n_, ((int)(K - 2) + 1), (mccM(&IM1_) * mccN(&IM1_)));
            mccAllocateMatrix(&FPc, m_, n_);
            mccCheckMatrixSize(&FP, K, mccGetMaxIndex(&IM1_ ,mccN(&FP)));
            I_FPc = (mccM(&FPc) != 1 || mccN(&FPc) != 1);
            p_FPc = mccPR(&FPc);
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_)
               {
                  p_FP = mccPR(&FP) + mccM(&FP) * ((int)(*p_IM1_ - .5)) + ((int)(2 - .5));
                  for (i_=0; i_<m_; ++i_, p_FPc+=I_FPc, p_FP+=I_FP)
                  {
                     *p_FPc = *p_FP;
                  }
               }
            }
         }
         /* % on cree un matrice de la forme [F,F,F...F] (M fois)  */
         /* % sur la troisieme dimension */
         /* Fc=reshape(Fc(:)*ones(1,M),K-1,length(n),M);      */
         mccCopy(&RM0_, &Fc);
         mxSetM( &RM0_, mccM(&RM0_) * mccN(&RM0_));
         mxSetN( &RM0_, 1);
         mccOnesMN(&IM1_, 1, M);
         mccRealMatrixMultiply(&RM1_, &RM0_, &IM1_);
         I0_ = mccGetLength(&n);
         Mprhs_[0] = &RM1_;
         Mprhs_[1] = mccTempMatrix((K - 1), 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[3] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &Fc;
         mccCallMATLAB(1, Mplhs_, 4, Mprhs_, "reshape", 399);
         /* FPc=reshape(FPc(:)*ones(1,M),K-1,length(n),M);      */
         mccCopy(&RM1_, &FPc);
         mxSetM( &RM1_, mccM(&RM1_) * mccN(&RM1_));
         mxSetN( &RM1_, 1);
         mccOnesMN(&IM1_, 1, M);
         mccRealMatrixMultiply(&RM0_, &RM1_, &IM1_);
         I0_ = mccGetLength(&n);
         Mprhs_[0] = &RM0_;
         Mprhs_[1] = mccTempMatrix((K - 1), 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[3] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &FPc;
         mccCallMATLAB(1, Mplhs_, 4, Mprhs_, "reshape", 400);
         
         /* % modification du terme source     */
         /* s     = s + squeeze( sum( ... */
         I0_ = mccGetLength(&n);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I0_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM1_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 404);
         mccFindIndex(&IM0_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM0_;
            int I_RM0_=1;
            double *p_a;
            int I_a=1, J_a;
            double *p_Fa;
            int I_Fa=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&a), mccN(&a));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&Fa), mccN(&Fa));
            mccAllocateMatrix(&RM0_, m_, n_);
            mccCheckMatrixSize(&a, mccM(&a), mccN(&a));
            I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
            p_RM0_ = mccPR(&RM0_);
            if (mccM(&a) == 1 && mccN(&a) == 1) { I_a = J_a = 0;}
            else { I_a = 1; J_a=mccM(&a)-m_; }
            p_a = mccPR(&a) + 0 + mccM(&a) * 0;
            I_Fa = (mccM(&Fa) != 1 || mccN(&Fa) != 1);
            p_Fa = mccPR(&Fa);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_a += J_a)
               {
                  for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_a+=I_a, p_Fa+=I_Fa)
                  {
                     *p_RM0_ = (*p_a * (double) *p_Fa);
                  }
               }
            }
         }
         mccSTRING(&RM0_) = 0;
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = &IM1_;
         Mprhs_[2] = &RM0_;
         Mplhs_[0] = &RM1_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "cat", 404);
         mccFindIndex(&IM2_, &n);
         mccFindIndex(&IM3_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM2_;
            int I_RM2_=1;
            double *p_c;
            int I_c=1, J_c;
            double *p_Fc;
            int I_Fc=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&Fc), mccN(&Fc));
            mccAllocateMatrix(&RM2_, m_, n_);
            mccCheckMatrixSize(&c, mccM(&c), mccN(&c));
            I_RM2_ = (mccM(&RM2_) != 1 || mccN(&RM2_) != 1);
            p_RM2_ = mccPR(&RM2_);
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_c = J_c = 0;}
            else { I_c = 1; J_c=mccM(&c)-m_; }
            p_c = mccPR(&c) + 0 + mccM(&c) * 0;
            I_Fc = (mccM(&Fc) != 1 || mccN(&Fc) != 1);
            p_Fc = mccPR(&Fc);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_c += J_c)
               {
                  for (i_=0; i_<m_; ++i_, p_RM2_+=I_RM2_, p_c+=I_c, p_Fc+=I_Fc)
                  {
                     *p_RM2_ = (*p_c * (double) *p_Fc);
                  }
               }
            }
         }
         mccSTRING(&RM2_) = 0;
         I1_ = mccGetLength(&n);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I1_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM4_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 406);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = &RM2_;
         Mprhs_[2] = &IM4_;
         Mplhs_[0] = &RM3_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "cat", 406);
         I2_ = mccGetLength(&n);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I2_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM5_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 407);
         mccFindIndex(&IM6_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM4_;
            int I_RM4_=1;
            double *p_ap;
            int I_ap=1, J_ap;
            double *p_FPa;
            int I_FPa=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), mccN(&ap));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&FPa), mccN(&FPa));
            mccAllocateMatrix(&RM4_, m_, n_);
            mccCheckMatrixSize(&ap, mccM(&ap), mccN(&ap));
            I_RM4_ = (mccM(&RM4_) != 1 || mccN(&RM4_) != 1);
            p_RM4_ = mccPR(&RM4_);
            if (mccM(&ap) == 1 && mccN(&ap) == 1) { I_ap = J_ap = 0;}
            else { I_ap = 1; J_ap=mccM(&ap)-m_; }
            p_ap = mccPR(&ap) + 0 + mccM(&ap) * 0;
            I_FPa = (mccM(&FPa) != 1 || mccN(&FPa) != 1);
            p_FPa = mccPR(&FPa);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_ap += J_ap)
               {
                  for (i_=0; i_<m_; ++i_, p_RM4_+=I_RM4_, p_ap+=I_ap, p_FPa+=I_FPa)
                  {
                     *p_RM4_ = (*p_ap * (double) *p_FPa);
                  }
               }
            }
         }
         mccSTRING(&RM4_) = 0;
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = &IM5_;
         Mprhs_[2] = &RM4_;
         Mplhs_[0] = &RM5_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "cat", 407);
         mccFindIndex(&IM7_, &n);
         mccFindIndex(&IM8_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM6_;
            int I_RM6_=1;
            double *p_cp;
            int I_cp=1, J_cp;
            double *p_FPc;
            int I_FPc=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&cp), mccN(&cp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&FPc), mccN(&FPc));
            mccAllocateMatrix(&RM6_, m_, n_);
            mccCheckMatrixSize(&cp, mccM(&cp), mccN(&cp));
            I_RM6_ = (mccM(&RM6_) != 1 || mccN(&RM6_) != 1);
            p_RM6_ = mccPR(&RM6_);
            if (mccM(&cp) == 1 && mccN(&cp) == 1) { I_cp = J_cp = 0;}
            else { I_cp = 1; J_cp=mccM(&cp)-m_; }
            p_cp = mccPR(&cp) + 0 + mccM(&cp) * 0;
            I_FPc = (mccM(&FPc) != 1 || mccN(&FPc) != 1);
            p_FPc = mccPR(&FPc);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_cp += J_cp)
               {
                  for (i_=0; i_<m_; ++i_, p_RM6_+=I_RM6_, p_cp+=I_cp, p_FPc+=I_FPc)
                  {
                     *p_RM6_ = (*p_cp * (double) *p_FPc);
                  }
               }
            }
         }
         mccSTRING(&RM6_) = 0;
         I3_ = mccGetLength(&n);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I3_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &IM9_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 409);
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = &RM6_;
         Mprhs_[2] = &IM9_;
         Mplhs_[0] = &RM7_;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "cat", 409);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM8_;
            int I_RM8_=1;
            double *p_f;
            int I_f=1;
            double *p_RM1_;
            int I_RM1_=1;
            double *p_b;
            int I_b=1, J_b;
            double *p_Fb;
            int I_Fb=1;
            double *p_RM3_;
            int I_RM3_=1;
            double *p_1f;
            int I_1f=1;
            double *p_RM5_;
            int I_RM5_=1;
            double *p_bp;
            int I_bp=1, J_bp;
            double *p_FPb;
            int I_FPb=1;
            double *p_RM7_;
            int I_RM7_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM1_), mccN(&RM1_));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&Fb), mccN(&Fb));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM3_), mccN(&RM3_));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&f), mccN(&f));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM5_), mccN(&RM5_));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&FPb), mccN(&FPb));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM7_), mccN(&RM7_));
            mccAllocateMatrix(&RM8_, m_, n_);
            mccCheckMatrixSize(&b, mccM(&b), mccN(&b));
            mccCheckMatrixSize(&bp, mccM(&bp), mccN(&bp));
            I_RM8_ = (mccM(&RM8_) != 1 || mccN(&RM8_) != 1);
            p_RM8_ = mccPR(&RM8_);
            I_f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_f = mccPR(&f);
            I_RM1_ = (mccM(&RM1_) != 1 || mccN(&RM1_) != 1);
            p_RM1_ = mccPR(&RM1_);
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_b = J_b = 0;}
            else { I_b = 1; J_b=mccM(&b)-m_; }
            p_b = mccPR(&b) + 0 + mccM(&b) * 0;
            I_Fb = (mccM(&Fb) != 1 || mccN(&Fb) != 1);
            p_Fb = mccPR(&Fb);
            I_RM3_ = (mccM(&RM3_) != 1 || mccN(&RM3_) != 1);
            p_RM3_ = mccPR(&RM3_);
            I_1f = (mccM(&f) != 1 || mccN(&f) != 1);
            p_1f = mccPR(&f);
            I_RM5_ = (mccM(&RM5_) != 1 || mccN(&RM5_) != 1);
            p_RM5_ = mccPR(&RM5_);
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_bp = J_bp = 0;}
            else { I_bp = 1; J_bp=mccM(&bp)-m_; }
            p_bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            I_FPb = (mccM(&FPb) != 1 || mccN(&FPb) != 1);
            p_FPb = mccPR(&FPb);
            I_RM7_ = (mccM(&RM7_) != 1 || mccN(&RM7_) != 1);
            p_RM7_ = mccPR(&RM7_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_b += J_b, p_bp += J_bp)
               {
                  for (i_=0; i_<m_; ++i_, p_RM8_+=I_RM8_, p_f+=I_f, p_RM1_+=I_RM1_, p_b+=I_b, p_Fb+=I_Fb, p_RM3_+=I_RM3_, p_1f+=I_1f, p_RM5_+=I_RM5_, p_bp+=I_bp, p_FPb+=I_FPb, p_RM7_+=I_RM7_)
                  {
                     *p_RM8_ = ((*p_f * (double) ((*p_RM1_ + (*p_b * (double) *p_Fb)) + *p_RM3_)) + ((1 - *p_1f) * (double) ((*p_RM5_ + (*p_bp * (double) *p_FPb)) + *p_RM7_)));
                  }
               }
            }
         }
         Mprhs_[0] = &RM8_;
         Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
         Mplhs_[0] = &RM9_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sum", 409);
         Mprhs_[0] = &RM9_;
         Mplhs_[0] = &RM10_;
         mccCallMATLAB(1, Mplhs_, 1, Mprhs_, "squeeze", 409);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_s;
            int I_s=1;
            double *p_1s;
            int I_1s=1;
            double *p_RM10_;
            int I_RM10_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&s), mccN(&s));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&RM10_), mccN(&RM10_));
            mccAllocateMatrix(&s, m_, n_);
            I_s = (mccM(&s) != 1 || mccN(&s) != 1);
            p_s = mccPR(&s);
            I_1s = (mccM(&s) != 1 || mccN(&s) != 1);
            p_1s = mccPR(&s);
            I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
            p_RM10_ = mccPR(&RM10_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_s+=I_s, p_1s+=I_1s, p_RM10_+=I_RM10_)
                  {
                     *p_s = (*p_1s + *p_RM10_);
                  }
               }
            }
         }
         
         /* % annulation des coefficients deja utilises */
         /* comp=zeros(K,length(n),M); */
         I3_ = mccGetLength(&n);
         Mprhs_[0] = mccTempMatrix(K, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix(I3_, 0., mccINT, 0 );
         Mprhs_[2] = mccTempMatrix(M, 0., mccINT, 0 );
         Mplhs_[0] = &M_comp;
         mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "zeros", 412);
         /* a(:,n,:)=comp; */
         mccFindIndex(&IM9_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_a;
            int I_a=1, J_a;
            int *p_M_comp;
            int I_M_comp=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&M_comp), mccN(&M_comp));
            if (!mccNOTSET(&a)) m_ = mcmCalcResultSize(m_, &n_, mccM(&a), mccN(&a));
            mccGrowMatrix(&a, m_, n_);
            if (mccM(&a) == 1 && mccN(&a) == 1) { I_a = J_a = 0;}
            else { I_a = 1; J_a=mccM(&a)-m_; }
            p_a = mccPR(&a) + 0 + mccM(&a) * 0;
            I_M_comp = (mccM(&M_comp) != 1 || mccN(&M_comp) != 1);
            p_M_comp = mccIPR(&M_comp);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_a += J_a)
               {
                  for (i_=0; i_<m_; ++i_, p_a+=I_a, p_M_comp+=I_M_comp)
                  {
                     *p_a = ((int)*p_M_comp);
                  }
               }
            }
         }
         /* ap(:,n,:)=comp; */
         mccFindIndex(&IM9_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_ap;
            int I_ap=1, J_ap;
            int *p_M_comp;
            int I_M_comp=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&M_comp), mccN(&M_comp));
            if (!mccNOTSET(&ap)) m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), mccN(&ap));
            mccGrowMatrix(&ap, m_, n_);
            if (mccM(&ap) == 1 && mccN(&ap) == 1) { I_ap = J_ap = 0;}
            else { I_ap = 1; J_ap=mccM(&ap)-m_; }
            p_ap = mccPR(&ap) + 0 + mccM(&ap) * 0;
            I_M_comp = (mccM(&M_comp) != 1 || mccN(&M_comp) != 1);
            p_M_comp = mccIPR(&M_comp);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_ap += J_ap)
               {
                  for (i_=0; i_<m_; ++i_, p_ap+=I_ap, p_M_comp+=I_M_comp)
                  {
                     *p_ap = ((int)*p_M_comp);
                  }
               }
            }
         }
         /* b(:,n,:)=comp; */
         mccFindIndex(&IM9_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_b;
            int I_b=1, J_b;
            int *p_M_comp;
            int I_M_comp=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&M_comp), mccN(&M_comp));
            if (!mccNOTSET(&b)) m_ = mcmCalcResultSize(m_, &n_, mccM(&b), mccN(&b));
            mccGrowMatrix(&b, m_, n_);
            if (mccM(&b) == 1 && mccN(&b) == 1) { I_b = J_b = 0;}
            else { I_b = 1; J_b=mccM(&b)-m_; }
            p_b = mccPR(&b) + 0 + mccM(&b) * 0;
            I_M_comp = (mccM(&M_comp) != 1 || mccN(&M_comp) != 1);
            p_M_comp = mccIPR(&M_comp);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_b += J_b)
               {
                  for (i_=0; i_<m_; ++i_, p_b+=I_b, p_M_comp+=I_M_comp)
                  {
                     *p_b = ((int)*p_M_comp);
                  }
               }
            }
         }
         /* bp(:,n,:)=comp; */
         mccFindIndex(&IM9_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_bp;
            int I_bp=1, J_bp;
            int *p_M_comp;
            int I_M_comp=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&M_comp), mccN(&M_comp));
            if (!mccNOTSET(&bp)) m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), mccN(&bp));
            mccGrowMatrix(&bp, m_, n_);
            if (mccM(&bp) == 1 && mccN(&bp) == 1) { I_bp = J_bp = 0;}
            else { I_bp = 1; J_bp=mccM(&bp)-m_; }
            p_bp = mccPR(&bp) + 0 + mccM(&bp) * 0;
            I_M_comp = (mccM(&M_comp) != 1 || mccN(&M_comp) != 1);
            p_M_comp = mccIPR(&M_comp);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_bp += J_bp)
               {
                  for (i_=0; i_<m_; ++i_, p_bp+=I_bp, p_M_comp+=I_M_comp)
                  {
                     *p_bp = ((int)*p_M_comp);
                  }
               }
            }
         }
         /* c(:,n,:)=comp; */
         mccFindIndex(&IM9_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_c;
            int I_c=1, J_c;
            int *p_M_comp;
            int I_M_comp=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&M_comp), mccN(&M_comp));
            if (!mccNOTSET(&c)) m_ = mcmCalcResultSize(m_, &n_, mccM(&c), mccN(&c));
            mccGrowMatrix(&c, m_, n_);
            if (mccM(&c) == 1 && mccN(&c) == 1) { I_c = J_c = 0;}
            else { I_c = 1; J_c=mccM(&c)-m_; }
            p_c = mccPR(&c) + 0 + mccM(&c) * 0;
            I_M_comp = (mccM(&M_comp) != 1 || mccN(&M_comp) != 1);
            p_M_comp = mccIPR(&M_comp);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_c += J_c)
               {
                  for (i_=0; i_<m_; ++i_, p_c+=I_c, p_M_comp+=I_M_comp)
                  {
                     *p_c = ((int)*p_M_comp);
                  }
               }
            }
         }
         /* cp(:,n,:)=comp; */
         mccFindIndex(&IM9_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_cp;
            int I_cp=1, J_cp;
            int *p_M_comp;
            int I_M_comp=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&M_comp), mccN(&M_comp));
            if (!mccNOTSET(&cp)) m_ = mcmCalcResultSize(m_, &n_, mccM(&cp), mccN(&cp));
            mccGrowMatrix(&cp, m_, n_);
            if (mccM(&cp) == 1 && mccN(&cp) == 1) { I_cp = J_cp = 0;}
            else { I_cp = 1; J_cp=mccM(&cp)-m_; }
            p_cp = mccPR(&cp) + 0 + mccM(&cp) * 0;
            I_M_comp = (mccM(&M_comp) != 1 || mccN(&M_comp) != 1);
            p_M_comp = mccIPR(&M_comp);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_cp += J_cp)
               {
                  for (i_=0; i_<m_; ++i_, p_cp+=I_cp, p_M_comp+=I_M_comp)
                  {
                     *p_cp = ((int)*p_M_comp);
                  }
               }
            }
         }
         
         /* end */
      }
      
      /* % 5 - creation des matrices creuses et des vecteurs pour */
      /* % la resolution du systeme lineaire (au sens des moindres carres) */
      
      /* % a - Les vecteurs de donnees au temps t */
      /* FS = F(:); */
      mccCopy(&FS, &F);
      mxSetM( &FS, mccM(&FS) * mccN(&FS));
      mxSetN( &FS, 1);
      /* S  = s(:) + FS; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_S;
         int I_S=1;
         double *p_s;
         int I_s=1;
         double *p_FS;
         int I_FS=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&s)*mccN(&s), 1);
         m_ = mcmCalcResultSize(m_, &n_, mccM(&FS), mccN(&FS));
         mccAllocateMatrix(&S, m_, n_);
         mccCheckVectorSize(&s, mccM(&s)*mccN(&s));
         I_S = (mccM(&S) != 1 || mccN(&S) != 1);
         p_S = mccPR(&S);
         I_s = (mccM(&s) != 1 || mccN(&s) != 1);
         p_s = mccPR(&s) + 0;
         I_FS = (mccM(&FS) != 1 || mccN(&FS) != 1);
         p_FS = mccPR(&FS);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_S+=I_S, p_s+=I_s, p_FS+=I_FS)
               {
                  *p_S = (*p_s + *p_FS);
               }
            }
         }
      }
      
      /* % vecteur d'indice */
      /* ia  = 2:K; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         int *p_ia;
         int I_ia=1;
         m_ = mcmCalcResultSize(m_, &n_, 1, ((int)(K - 2) + 1));
         m_ = mcmCalcResultSize(m_, &n_, 1, ((int)(K - 2) + 1));
         mccAllocateMatrix(&ia, m_, n_);
         I_ia = (mccM(&ia) != 1 || mccN(&ia) != 1);
         p_ia = mccIPR(&ia);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_ia+=I_ia)
               {
                  *p_ia = ((int)(2 + i_+j_*m_));
               }
            }
         }
      }
      /* ib  = 1:K; */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         int *p_ib;
         int I_ib=1;
         m_ = mcmCalcResultSize(m_, &n_, 1, ((int)(K - 1) + 1));
         m_ = mcmCalcResultSize(m_, &n_, 1, ((int)(K - 1) + 1));
         mccAllocateMatrix(&ib, m_, n_);
         I_ib = (mccM(&ib) != 1 || mccN(&ib) != 1);
         p_ib = mccIPR(&ib);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_ib+=I_ib)
               {
                  *p_ib = ((int)(1 + i_+j_*m_));
               }
            }
         }
      }
      /* ic  = 1:(K - 1); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         int *p_ic;
         int I_ic=1;
         m_ = mcmCalcResultSize(m_, &n_, 1, ((int)((K - 1) - 1) + 1));
         m_ = mcmCalcResultSize(m_, &n_, 1, ((int)((K - 1) - 1) + 1));
         mccAllocateMatrix(&ic, m_, n_);
         I_ic = (mccM(&ic) != 1 || mccN(&ic) != 1);
         p_ic = mccIPR(&ic);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_ic+=I_ic)
               {
                  *p_ic = ((int)(1 + i_+j_*m_));
               }
            }
         }
      }
      
      /* % initialisation des buffer a vide */
      /* alpha  = []; */
      mccCreateEmpty(&alpha);
      /* alphap = []; */
      mccCreateEmpty(&alphap);
      /* p      = []; */
      mccCreateEmpty(&p);
      /* q      = []; */
      mccCreateEmpty(&q);
      
      /* % boucle de remplissage des bufffers */
      /* for m = 1:M */
      for (I3_ = 1; I3_ <= M; I3_ = I3_ + 1)
      {
         m = I3_;
         /* for n = 1:M */
         for (I2_ = 1; I2_ <= M; I2_ = I2_ + 1)
         {
            n_1 = I2_;
            
            /* % pour a */
            /* alpha    = cat(1,alpha,a(ia,n,m)); */
            mccFindIndex(&IM9_, &ia);
            {
               int i_, j_;
               int m_=1, n_=1, cx_ = 0;
               double *p_RM10_;
               int I_RM10_=1;
               double *p_a;
               int I_a=1, J_a;
               m_ = mcmCalcResultSize(m_, &n_, mccM(&a), 1);
               mccAllocateMatrix(&RM10_, m_, n_);
               mccCheckMatrixSize(&a, mccM(&a), m);
               I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
               p_RM10_ = mccPR(&RM10_);
               if (mccM(&a) == 1) { I_a = J_a = 0;}
               else { I_a = 1; J_a=mccM(&a)-m_; }
               p_a = mccPR(&a) + 0 + mccM(&a) * (m-1);
               if (m_ != 0)
               {
                  for (j_=0; j_<n_; ++j_, p_a += J_a)
                  {
                     for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_a+=I_a)
                     {
                        *p_RM10_ = *p_a;
                     }
                  }
               }
            }
            mccSTRING(&RM10_) = 0;
            Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
            Mprhs_[1] = &alpha;
            Mprhs_[2] = &RM10_;
            Mplhs_[0] = &alpha;
            mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "cat", 445);
            /* alphap   = cat(1,alphap,ap(ia,n,m)); */
            mccFindIndex(&IM9_, &ia);
            {
               int i_, j_;
               int m_=1, n_=1, cx_ = 0;
               double *p_RM10_;
               int I_RM10_=1;
               double *p_ap;
               int I_ap=1, J_ap;
               m_ = mcmCalcResultSize(m_, &n_, mccM(&ap), 1);
               mccAllocateMatrix(&RM10_, m_, n_);
               mccCheckMatrixSize(&ap, mccM(&ap), m);
               I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
               p_RM10_ = mccPR(&RM10_);
               if (mccM(&ap) == 1) { I_ap = J_ap = 0;}
               else { I_ap = 1; J_ap=mccM(&ap)-m_; }
               p_ap = mccPR(&ap) + 0 + mccM(&ap) * (m-1);
               if (m_ != 0)
               {
                  for (j_=0; j_<n_; ++j_, p_ap += J_ap)
                  {
                     for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_ap+=I_ap)
                     {
                        *p_RM10_ = *p_ap;
                     }
                  }
               }
            }
            mccSTRING(&RM10_) = 0;
            Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
            Mprhs_[1] = &alphap;
            Mprhs_[2] = &RM10_;
            Mplhs_[0] = &alphap;
            mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "cat", 446);
            /* p        = cat(1,p,ia' + K .* (m - 1));  */
            mccConjTrans(&IM9_, &ia);
            R0_ = (K * (double) (m - 1));
            {
               int i_, j_;
               int m_=1, n_=1, cx_ = 0;
               double *p_RM10_;
               int I_RM10_=1;
               int *p_IM9_;
               int I_IM9_=1;
               m_ = mcmCalcResultSize(m_, &n_, mccM(&IM9_), mccN(&IM9_));
               mccAllocateMatrix(&RM10_, m_, n_);
               I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
               p_RM10_ = mccPR(&RM10_);
               I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
               p_IM9_ = mccIPR(&IM9_);
               if (m_ != 0)
               {
                  for (j_=0; j_<n_; ++j_)
                  {
                     for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_IM9_+=I_IM9_)
                     {
                        *p_RM10_ = (((int)*p_IM9_) + R0_);
                     }
                  }
               }
            }
            mccSTRING(&RM10_) = 0;
            Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
            Mprhs_[1] = &p;
            Mprhs_[2] = &RM10_;
            Mplhs_[0] = &p;
            mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "cat", 447);
            /* q        = cat(1,q,ia' - 1 + K .* (n - 1));  */
            mccConjTrans(&IM9_, &ia);
            R0_ = (K * (double) (n_1 - 1));
            {
               int i_, j_;
               int m_=1, n_=1, cx_ = 0;
               double *p_RM10_;
               int I_RM10_=1;
               int *p_IM9_;
               int I_IM9_=1;
               m_ = mcmCalcResultSize(m_, &n_, mccM(&IM9_), mccN(&IM9_));
               mccAllocateMatrix(&RM10_, m_, n_);
               I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
               p_RM10_ = mccPR(&RM10_);
               I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
               p_IM9_ = mccIPR(&IM9_);
               if (m_ != 0)
               {
                  for (j_=0; j_<n_; ++j_)
                  {
                     for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_IM9_+=I_IM9_)
                     {
                        *p_RM10_ = ((((int)*p_IM9_) - 1) + R0_);
                     }
                  }
               }
            }
            mccSTRING(&RM10_) = 0;
            Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
            Mprhs_[1] = &q;
            Mprhs_[2] = &RM10_;
            Mplhs_[0] = &q;
            mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "cat", 448);
            
            /* % pour b */
            /* alpha    = cat(1,alpha,b(ib,n,m)); */
            mccFindIndex(&IM9_, &ib);
            {
               int i_, j_;
               int m_=1, n_=1, cx_ = 0;
               double *p_RM10_;
               int I_RM10_=1;
               double *p_b;
               int I_b=1, J_b;
               m_ = mcmCalcResultSize(m_, &n_, mccM(&b), 1);
               mccAllocateMatrix(&RM10_, m_, n_);
               mccCheckMatrixSize(&b, mccM(&b), m);
               I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
               p_RM10_ = mccPR(&RM10_);
               if (mccM(&b) == 1) { I_b = J_b = 0;}
               else { I_b = 1; J_b=mccM(&b)-m_; }
               p_b = mccPR(&b) + 0 + mccM(&b) * (m-1);
               if (m_ != 0)
               {
                  for (j_=0; j_<n_; ++j_, p_b += J_b)
                  {
                     for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_b+=I_b)
                     {
                        *p_RM10_ = *p_b;
                     }
                  }
               }
            }
            mccSTRING(&RM10_) = 0;
            Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
            Mprhs_[1] = &alpha;
            Mprhs_[2] = &RM10_;
            Mplhs_[0] = &alpha;
            mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "cat", 451);
            /* alphap   = cat(1,alphap,bp(ib,n,m)); */
            mccFindIndex(&IM9_, &ib);
            {
               int i_, j_;
               int m_=1, n_=1, cx_ = 0;
               double *p_RM10_;
               int I_RM10_=1;
               double *p_bp;
               int I_bp=1, J_bp;
               m_ = mcmCalcResultSize(m_, &n_, mccM(&bp), 1);
               mccAllocateMatrix(&RM10_, m_, n_);
               mccCheckMatrixSize(&bp, mccM(&bp), m);
               I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
               p_RM10_ = mccPR(&RM10_);
               if (mccM(&bp) == 1) { I_bp = J_bp = 0;}
               else { I_bp = 1; J_bp=mccM(&bp)-m_; }
               p_bp = mccPR(&bp) + 0 + mccM(&bp) * (m-1);
               if (m_ != 0)
               {
                  for (j_=0; j_<n_; ++j_, p_bp += J_bp)
                  {
                     for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_bp+=I_bp)
                     {
                        *p_RM10_ = *p_bp;
                     }
                  }
               }
            }
            mccSTRING(&RM10_) = 0;
            Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
            Mprhs_[1] = &alphap;
            Mprhs_[2] = &RM10_;
            Mplhs_[0] = &alphap;
            mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "cat", 452);
            /* p        = cat(1,p,ib' + K .* (m - 1));  */
            mccConjTrans(&IM9_, &ib);
            R0_ = (K * (double) (m - 1));
            {
               int i_, j_;
               int m_=1, n_=1, cx_ = 0;
               double *p_RM10_;
               int I_RM10_=1;
               int *p_IM9_;
               int I_IM9_=1;
               m_ = mcmCalcResultSize(m_, &n_, mccM(&IM9_), mccN(&IM9_));
               mccAllocateMatrix(&RM10_, m_, n_);
               I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
               p_RM10_ = mccPR(&RM10_);
               I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
               p_IM9_ = mccIPR(&IM9_);
               if (m_ != 0)
               {
                  for (j_=0; j_<n_; ++j_)
                  {
                     for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_IM9_+=I_IM9_)
                     {
                        *p_RM10_ = (((int)*p_IM9_) + R0_);
                     }
                  }
               }
            }
            mccSTRING(&RM10_) = 0;
            Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
            Mprhs_[1] = &p;
            Mprhs_[2] = &RM10_;
            Mplhs_[0] = &p;
            mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "cat", 453);
            /* q        = cat(1,q,ib' + K .* (n - 1));  */
            mccConjTrans(&IM9_, &ib);
            R0_ = (K * (double) (n_1 - 1));
            {
               int i_, j_;
               int m_=1, n_=1, cx_ = 0;
               double *p_RM10_;
               int I_RM10_=1;
               int *p_IM9_;
               int I_IM9_=1;
               m_ = mcmCalcResultSize(m_, &n_, mccM(&IM9_), mccN(&IM9_));
               mccAllocateMatrix(&RM10_, m_, n_);
               I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
               p_RM10_ = mccPR(&RM10_);
               I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
               p_IM9_ = mccIPR(&IM9_);
               if (m_ != 0)
               {
                  for (j_=0; j_<n_; ++j_)
                  {
                     for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_IM9_+=I_IM9_)
                     {
                        *p_RM10_ = (((int)*p_IM9_) + R0_);
                     }
                  }
               }
            }
            mccSTRING(&RM10_) = 0;
            Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
            Mprhs_[1] = &q;
            Mprhs_[2] = &RM10_;
            Mplhs_[0] = &q;
            mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "cat", 454);
            
            /* % pour c */
            /* alpha    = cat(1,alpha,c(ic,n,m)); */
            mccFindIndex(&IM9_, &ic);
            {
               int i_, j_;
               int m_=1, n_=1, cx_ = 0;
               double *p_RM10_;
               int I_RM10_=1;
               double *p_c;
               int I_c=1, J_c;
               m_ = mcmCalcResultSize(m_, &n_, mccM(&c), 1);
               mccAllocateMatrix(&RM10_, m_, n_);
               mccCheckMatrixSize(&c, mccM(&c), m);
               I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
               p_RM10_ = mccPR(&RM10_);
               if (mccM(&c) == 1) { I_c = J_c = 0;}
               else { I_c = 1; J_c=mccM(&c)-m_; }
               p_c = mccPR(&c) + 0 + mccM(&c) * (m-1);
               if (m_ != 0)
               {
                  for (j_=0; j_<n_; ++j_, p_c += J_c)
                  {
                     for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_c+=I_c)
                     {
                        *p_RM10_ = *p_c;
                     }
                  }
               }
            }
            mccSTRING(&RM10_) = 0;
            Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
            Mprhs_[1] = &alpha;
            Mprhs_[2] = &RM10_;
            Mplhs_[0] = &alpha;
            mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "cat", 457);
            mccLOG(&alpha) = 0;
            /* alphap   = cat(1,alphap,cp(ic,n,m)); */
            mccFindIndex(&IM9_, &ic);
            {
               int i_, j_;
               int m_=1, n_=1, cx_ = 0;
               double *p_RM10_;
               int I_RM10_=1;
               double *p_cp;
               int I_cp=1, J_cp;
               m_ = mcmCalcResultSize(m_, &n_, mccM(&cp), 1);
               mccAllocateMatrix(&RM10_, m_, n_);
               mccCheckMatrixSize(&cp, mccM(&cp), m);
               I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
               p_RM10_ = mccPR(&RM10_);
               if (mccM(&cp) == 1) { I_cp = J_cp = 0;}
               else { I_cp = 1; J_cp=mccM(&cp)-m_; }
               p_cp = mccPR(&cp) + 0 + mccM(&cp) * (m-1);
               if (m_ != 0)
               {
                  for (j_=0; j_<n_; ++j_, p_cp += J_cp)
                  {
                     for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_cp+=I_cp)
                     {
                        *p_RM10_ = *p_cp;
                     }
                  }
               }
            }
            mccSTRING(&RM10_) = 0;
            Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
            Mprhs_[1] = &alphap;
            Mprhs_[2] = &RM10_;
            Mplhs_[0] = &alphap;
            mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "cat", 458);
            mccLOG(&alphap) = 0;
            /* p        = cat(1,p,ic' + K .* (m - 1));  */
            mccConjTrans(&IM9_, &ic);
            R0_ = (K * (double) (m - 1));
            {
               int i_, j_;
               int m_=1, n_=1, cx_ = 0;
               double *p_RM10_;
               int I_RM10_=1;
               int *p_IM9_;
               int I_IM9_=1;
               m_ = mcmCalcResultSize(m_, &n_, mccM(&IM9_), mccN(&IM9_));
               mccAllocateMatrix(&RM10_, m_, n_);
               I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
               p_RM10_ = mccPR(&RM10_);
               I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
               p_IM9_ = mccIPR(&IM9_);
               if (m_ != 0)
               {
                  for (j_=0; j_<n_; ++j_)
                  {
                     for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_IM9_+=I_IM9_)
                     {
                        *p_RM10_ = (((int)*p_IM9_) + R0_);
                     }
                  }
               }
            }
            mccSTRING(&RM10_) = 0;
            Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
            Mprhs_[1] = &p;
            Mprhs_[2] = &RM10_;
            Mplhs_[0] = &p;
            mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "cat", 459);
            mccLOG(&p) = 0;
            /* q        = cat(1,q,ic' + 1 + K .* (n - 1));  */
            mccConjTrans(&IM9_, &ic);
            R0_ = (K * (double) (n_1 - 1));
            {
               int i_, j_;
               int m_=1, n_=1, cx_ = 0;
               double *p_RM10_;
               int I_RM10_=1;
               int *p_IM9_;
               int I_IM9_=1;
               m_ = mcmCalcResultSize(m_, &n_, mccM(&IM9_), mccN(&IM9_));
               mccAllocateMatrix(&RM10_, m_, n_);
               I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
               p_RM10_ = mccPR(&RM10_);
               I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
               p_IM9_ = mccIPR(&IM9_);
               if (m_ != 0)
               {
                  for (j_=0; j_<n_; ++j_)
                  {
                     for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_IM9_+=I_IM9_)
                     {
                        *p_RM10_ = ((((int)*p_IM9_) + 1) + R0_);
                     }
                  }
               }
            }
            mccSTRING(&RM10_) = 0;
            Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
            Mprhs_[1] = &q;
            Mprhs_[2] = &RM10_;
            Mplhs_[0] = &q;
            mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "cat", 460);
            mccLOG(&q) = 0;
            
            /* end */
         }
         /* end */
      }
      
      /* % duplication des indice avnat retrait des elements non nul	 */
      /* pp = p; */
      mccCopy(&pp, &p);
      mccLOG(&pp) = mccLOG(&p);
      mccSTRING(&pp) = mccSTRING(&p);
      /* qp = q; */
      mccCopy(&qp, &q);
      mccLOG(&qp) = mccLOG(&q);
      mccSTRING(&qp) = mccSTRING(&q);
      
      /* % g - retrait des elements nuls */
      /* iz = find(alpha  ~= 0); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         unsigned short *p_BM0_;
         int I_BM0_=1;
         double *p_alpha;
         int I_alpha=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&alpha), mccN(&alpha));
         mccAllocateMatrix(&BM0_, m_, n_);
         I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
         p_BM0_ = mccSPR(&BM0_);
         I_alpha = (mccM(&alpha) != 1 || mccN(&alpha) != 1);
         p_alpha = mccPR(&alpha);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_BM0_+=I_BM0_, p_alpha+=I_alpha)
               {
                  *p_BM0_ = ( (*p_alpha != 0) || mccREL_NAN(*p_alpha) );
               }
            }
         }
      }
      mccFind(&iz, &BM0_);
      /* alpha = alpha(iz); */
      mccFindIndex(&IM9_, &iz);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM10_;
         int I_RM10_=1;
         double *p_alpha;
         int *p_IM9_;
         int I_IM9_=1;
         m_ = mccCalcSubscriptDimensions(m_, &n_, mccM(&IM9_), mccN(&IM9_), &alpha);
         mccAllocateMatrix(&RM10_, m_, n_);
         mccCheckVectorSize(&alpha, mccGetMaxIndex(&IM9_ ,mccM(&alpha)*mccN(&alpha)));
         I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
         p_RM10_ = mccPR(&RM10_);
         I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
         p_IM9_ = mccIPR(&IM9_);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_IM9_+=I_IM9_)
               {
                  *p_RM10_ = mccPR(&alpha)[((int)(*p_IM9_ - .5))];
               }
            }
         }
      }
      mccLOG(&RM10_) = mccLOG(&alpha);
      mccSTRING(&RM10_) = mccSTRING(&alpha);
      mccCopy(&alpha, &RM10_);
      mccLOG(&alpha) = mccLOG(&RM10_);
      mccSTRING(&alpha) = mccSTRING(&RM10_);
      /* p  =p(iz); */
      mccFindIndex(&IM9_, &iz);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM10_;
         int I_RM10_=1;
         double *p_p;
         int *p_IM9_;
         int I_IM9_=1;
         m_ = mccCalcSubscriptDimensions(m_, &n_, mccM(&IM9_), mccN(&IM9_), &p);
         mccAllocateMatrix(&RM10_, m_, n_);
         mccCheckVectorSize(&p, mccGetMaxIndex(&IM9_ ,mccM(&p)*mccN(&p)));
         I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
         p_RM10_ = mccPR(&RM10_);
         I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
         p_IM9_ = mccIPR(&IM9_);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_IM9_+=I_IM9_)
               {
                  *p_RM10_ = mccPR(&p)[((int)(*p_IM9_ - .5))];
               }
            }
         }
      }
      mccLOG(&RM10_) = mccLOG(&p);
      mccSTRING(&RM10_) = mccSTRING(&p);
      mccCopy(&p, &RM10_);
      mccLOG(&p) = mccLOG(&RM10_);
      mccSTRING(&p) = mccSTRING(&RM10_);
      /* q = q(iz); */
      mccFindIndex(&IM9_, &iz);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM10_;
         int I_RM10_=1;
         double *p_q;
         int *p_IM9_;
         int I_IM9_=1;
         m_ = mccCalcSubscriptDimensions(m_, &n_, mccM(&IM9_), mccN(&IM9_), &q);
         mccAllocateMatrix(&RM10_, m_, n_);
         mccCheckVectorSize(&q, mccGetMaxIndex(&IM9_ ,mccM(&q)*mccN(&q)));
         I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
         p_RM10_ = mccPR(&RM10_);
         I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
         p_IM9_ = mccIPR(&IM9_);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_IM9_+=I_IM9_)
               {
                  *p_RM10_ = mccPR(&q)[((int)(*p_IM9_ - .5))];
               }
            }
         }
      }
      mccLOG(&RM10_) = mccLOG(&q);
      mccSTRING(&RM10_) = mccSTRING(&q);
      mccCopy(&q, &RM10_);
      mccLOG(&q) = mccLOG(&RM10_);
      mccSTRING(&q) = mccSTRING(&RM10_);
      
      /* izp = find(alphap ~= 0); */
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         unsigned short *p_BM0_;
         int I_BM0_=1;
         double *p_alphap;
         int I_alphap=1;
         m_ = mcmCalcResultSize(m_, &n_, mccM(&alphap), mccN(&alphap));
         mccAllocateMatrix(&BM0_, m_, n_);
         I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
         p_BM0_ = mccSPR(&BM0_);
         I_alphap = (mccM(&alphap) != 1 || mccN(&alphap) != 1);
         p_alphap = mccPR(&alphap);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_BM0_+=I_BM0_, p_alphap+=I_alphap)
               {
                  *p_BM0_ = ( (*p_alphap != 0) || mccREL_NAN(*p_alphap) );
               }
            }
         }
      }
      mccFind(&izp, &BM0_);
      /* alphap = alphap(izp); */
      mccFindIndex(&IM9_, &izp);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM10_;
         int I_RM10_=1;
         double *p_alphap;
         int *p_IM9_;
         int I_IM9_=1;
         m_ = mccCalcSubscriptDimensions(m_, &n_, mccM(&IM9_), mccN(&IM9_), &alphap);
         mccAllocateMatrix(&RM10_, m_, n_);
         mccCheckVectorSize(&alphap, mccGetMaxIndex(&IM9_ ,mccM(&alphap)*mccN(&alphap)));
         I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
         p_RM10_ = mccPR(&RM10_);
         I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
         p_IM9_ = mccIPR(&IM9_);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_IM9_+=I_IM9_)
               {
                  *p_RM10_ = mccPR(&alphap)[((int)(*p_IM9_ - .5))];
               }
            }
         }
      }
      mccLOG(&RM10_) = mccLOG(&alphap);
      mccSTRING(&RM10_) = mccSTRING(&alphap);
      mccCopy(&alphap, &RM10_);
      mccLOG(&alphap) = mccLOG(&RM10_);
      mccSTRING(&alphap) = mccSTRING(&RM10_);
      /* pp = pp(izp); */
      mccFindIndex(&IM9_, &izp);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM10_;
         int I_RM10_=1;
         double *p_pp;
         int *p_IM9_;
         int I_IM9_=1;
         m_ = mccCalcSubscriptDimensions(m_, &n_, mccM(&IM9_), mccN(&IM9_), &pp);
         mccAllocateMatrix(&RM10_, m_, n_);
         mccCheckVectorSize(&pp, mccGetMaxIndex(&IM9_ ,mccM(&pp)*mccN(&pp)));
         I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
         p_RM10_ = mccPR(&RM10_);
         I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
         p_IM9_ = mccIPR(&IM9_);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_IM9_+=I_IM9_)
               {
                  *p_RM10_ = mccPR(&pp)[((int)(*p_IM9_ - .5))];
               }
            }
         }
      }
      mccLOG(&RM10_) = mccLOG(&pp);
      mccSTRING(&RM10_) = mccSTRING(&pp);
      mccCopy(&pp, &RM10_);
      mccLOG(&pp) = mccLOG(&RM10_);
      mccSTRING(&pp) = mccSTRING(&RM10_);
      /* qp = qp(izp); */
      mccFindIndex(&IM9_, &izp);
      {
         int i_, j_;
         int m_=1, n_=1, cx_ = 0;
         double *p_RM10_;
         int I_RM10_=1;
         double *p_qp;
         int *p_IM9_;
         int I_IM9_=1;
         m_ = mccCalcSubscriptDimensions(m_, &n_, mccM(&IM9_), mccN(&IM9_), &qp);
         mccAllocateMatrix(&RM10_, m_, n_);
         mccCheckVectorSize(&qp, mccGetMaxIndex(&IM9_ ,mccM(&qp)*mccN(&qp)));
         I_RM10_ = (mccM(&RM10_) != 1 || mccN(&RM10_) != 1);
         p_RM10_ = mccPR(&RM10_);
         I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
         p_IM9_ = mccIPR(&IM9_);
         if (m_ != 0)
         {
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_RM10_+=I_RM10_, p_IM9_+=I_IM9_)
               {
                  *p_RM10_ = mccPR(&qp)[((int)(*p_IM9_ - .5))];
               }
            }
         }
      }
      mccLOG(&RM10_) = mccLOG(&qp);
      mccSTRING(&RM10_) = mccSTRING(&qp);
      mccCopy(&qp, &RM10_);
      mccLOG(&qp) = mccLOG(&RM10_);
      mccSTRING(&qp) = mccSTRING(&RM10_);
      
      /* % h - creation des matrices creuses  */
      /* km = K .* M ; */
      km = (K * (double) M);
      
      /* % 6 - resolution du systeme  */
      /* [fpv,ALPHA,ALPHAP,IDENTITE] = rpde1dsolver_creux(km,f,p,q,alpha,pp,qp,alphap,FS,S,nargout > 5); */
      I3_ = mccNargout();
      Mprhs_[0] = mccTempMatrix(km, 0., mccREAL, 0 );
      Mprhs_[1] = &f;
      Mprhs_[2] = &p;
      Mprhs_[3] = &q;
      Mprhs_[4] = &alpha;
      Mprhs_[5] = &pp;
      Mprhs_[6] = &qp;
      Mprhs_[7] = &alphap;
      Mprhs_[8] = &FS;
      Mprhs_[9] = &S;
      Mprhs_[10] = mccTempMatrix((I3_ > 5), 0., mccBOOL, MCC_LOGICAL );
      Mplhs_[0] = &fpv;
      Mplhs_[1] = &ALPHA;
      Mplhs_[2] = &ALPHAP;
      Mplhs_[3] = &IDENTITE;
      mccCallMATLAB(4, Mplhs_, 11, Mprhs_, "rpde1dsolver_creux", 484);
      
      /* % 7 - mise en forme de la sortie */
      /* fp = reshape(fpv,[K,M]); */
      mccCatenateColumns(&IM9_, mccTempMatrix(K, 0., mccINT, 0 ), mccTempMatrix(M, 0., mccINT, 0 ));
      mccReshape2(&fp, &fpv, &IM9_);
      
      /* % 8  - on recupere les donnees interpretees (modifier par les conditions aux limites) */
      /* n = find(FP);             */
      mccFind(&n, &FP);
      /* if ~isempty(n)	 */
      B1_ = mccIsEmpty(&n);
      B0_ = (!B1_);
      if ((double)B0_)
      {
         /* fp(n)=FP(n); */
         mccFindIndex(&IM8_, &n);
         mccFindIndex(&IM9_, &n);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_fp;
            int *p_IM9_;
            int I_IM9_=1;
            double *p_FP;
            int *p_IM8_;
            int I_IM8_=1;
            m_ = mccCalcSubscriptDimensions(m_, &n_, mccM(&IM8_), mccN(&IM8_), &FP);
            if (m_ == 1 && n_ == 1) m_ = mccCalcSubscriptDimensions(m_, &n_, mccM(&IM9_), mccN(&IM9_), &fp);
            mccGrowVector(&fp, mccGetMaxIndex(&IM9_ ,mccM(&fp)*mccN(&fp)));
            mccCheckVectorSize(&FP, mccGetMaxIndex(&IM8_ ,mccM(&FP)*mccN(&FP)));
            I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
            p_IM9_ = mccIPR(&IM9_);
            I_IM8_ = (mccM(&IM8_) != 1 || mccN(&IM8_) != 1);
            p_IM8_ = mccIPR(&IM8_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_IM9_+=I_IM9_, p_IM8_+=I_IM8_)
                  {
                     mccPR(&fp)[((int)(*p_IM9_ - .5))] = mccPR(&FP)[((int)(*p_IM8_ - .5))];
                  }
               }
            }
         }
         /* end */
      }
      
      /* % 9 - calcul des deriveees */
      /* if nargout >1  */
      I3_ = mccNargout();
      B0_ = (I3_ > 1);
      if ((double)B0_)
      {
         
         /* % a - prolongation de la derivee 2  */
         /* %     par continiute de la derive 3. */
         /* fp0   = 4 .* fp(1,:)  - 6 .* fp(2,:)   + 4 .* fp(3,:)    - fp(4,:); */
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_fp0;
            int I_fp0=1;
            double *p_fp;
            int I_fp=1, J_fp;
            double *p_1fp;
            int I_1fp=1, J_1fp;
            double *p_2fp;
            int I_2fp=1, J_2fp;
            double *p_3fp;
            int I_3fp=1, J_3fp;
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&fp));
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&fp));
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&fp));
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&fp));
            mccAllocateMatrix(&fp0, m_, n_);
            mccCheckMatrixSize(&fp, 1, mccN(&fp));
            mccCheckMatrixSize(&fp, 2, mccN(&fp));
            mccCheckMatrixSize(&fp, 3, mccN(&fp));
            mccCheckMatrixSize(&fp, 4, mccN(&fp));
            I_fp0 = (mccM(&fp0) != 1 || mccN(&fp0) != 1);
            p_fp0 = mccPR(&fp0);
            if (mccN(&fp) == 1) { I_fp = J_fp = 0;}
            else { I_fp = 1; J_fp=mccM(&fp)-m_; }
            p_fp = mccPR(&fp) + (1-1) + mccM(&fp) * 0;
            if (mccN(&fp) == 1) { I_1fp = J_1fp = 0;}
            else { I_1fp = 1; J_1fp=mccM(&fp)-m_; }
            p_1fp = mccPR(&fp) + (2-1) + mccM(&fp) * 0;
            if (mccN(&fp) == 1) { I_2fp = J_2fp = 0;}
            else { I_2fp = 1; J_2fp=mccM(&fp)-m_; }
            p_2fp = mccPR(&fp) + (3-1) + mccM(&fp) * 0;
            if (mccN(&fp) == 1) { I_3fp = J_3fp = 0;}
            else { I_3fp = 1; J_3fp=mccM(&fp)-m_; }
            p_3fp = mccPR(&fp) + (4-1) + mccM(&fp) * 0;
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_fp += J_fp, p_1fp += J_1fp, p_2fp += J_2fp, p_3fp += J_3fp)
               {
                  for (i_=0; i_<m_; ++i_, p_fp0+=I_fp0, p_fp+=I_fp, p_1fp+=I_1fp, p_2fp+=I_2fp, p_3fp+=I_3fp)
                  {
                     *p_fp0 = ((((4 * (double) *p_fp) - (6 * (double) *p_1fp)) + (4 * (double) *p_2fp)) - *p_3fp);
                  }
               }
            }
         }
         /* fpkp1 = 4 .* fp(K,:)  - 6 .* fp(K-1,:) + 4 .* fp(K-2,:)  - fp(K-3,:); */
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_fpkp1;
            int I_fpkp1=1;
            double *p_fp;
            int I_fp=1, J_fp;
            double *p_1fp;
            int I_1fp=1, J_1fp;
            double *p_2fp;
            int I_2fp=1, J_2fp;
            double *p_3fp;
            int I_3fp=1, J_3fp;
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&fp));
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&fp));
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&fp));
            m_ = mcmCalcResultSize(m_, &n_, 1, mccN(&fp));
            mccAllocateMatrix(&fpkp1, m_, n_);
            mccCheckMatrixSize(&fp, K, mccN(&fp));
            mccCheckMatrixSize(&fp, (K-1), mccN(&fp));
            mccCheckMatrixSize(&fp, (K-2), mccN(&fp));
            mccCheckMatrixSize(&fp, (K-3), mccN(&fp));
            I_fpkp1 = (mccM(&fpkp1) != 1 || mccN(&fpkp1) != 1);
            p_fpkp1 = mccPR(&fpkp1);
            if (mccN(&fp) == 1) { I_fp = J_fp = 0;}
            else { I_fp = 1; J_fp=mccM(&fp)-m_; }
            p_fp = mccPR(&fp) + (K-1) + mccM(&fp) * 0;
            if (mccN(&fp) == 1) { I_1fp = J_1fp = 0;}
            else { I_1fp = 1; J_1fp=mccM(&fp)-m_; }
            p_1fp = mccPR(&fp) + ((K-1)-1) + mccM(&fp) * 0;
            if (mccN(&fp) == 1) { I_2fp = J_2fp = 0;}
            else { I_2fp = 1; J_2fp=mccM(&fp)-m_; }
            p_2fp = mccPR(&fp) + ((K-2)-1) + mccM(&fp) * 0;
            if (mccN(&fp) == 1) { I_3fp = J_3fp = 0;}
            else { I_3fp = 1; J_3fp=mccM(&fp)-m_; }
            p_3fp = mccPR(&fp) + ((K-3)-1) + mccM(&fp) * 0;
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_fp += J_fp, p_1fp += J_1fp, p_2fp += J_2fp, p_3fp += J_3fp)
               {
                  for (i_=0; i_<m_; ++i_, p_fpkp1+=I_fpkp1, p_fp+=I_fp, p_1fp+=I_1fp, p_2fp+=I_2fp, p_3fp+=I_3fp)
                  {
                     *p_fpkp1 = ((((4 * (double) *p_fp) - (6 * (double) *p_1fp)) + (4 * (double) *p_2fp)) - *p_3fp);
                  }
               }
            }
         }
         
         /* % b - matrice pour le calcul des derivees */
         /* fpp = cat(1,fp0,fp,fpkp1); */
         Mprhs_[0] = mccTempMatrix(1, 0., mccINT, 0 );
         Mprhs_[1] = &fp0;
         Mprhs_[2] = &fp;
         Mprhs_[3] = &fpkp1;
         Mplhs_[0] = &fpp;
         mccCallMATLAB(1, Mplhs_, 4, Mprhs_, "cat", 504);
         
         /* % c - calcule des derivees */
         /* im = 1:K; */
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            int *p_im;
            int I_im=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, ((int)(K - 1) + 1));
            m_ = mcmCalcResultSize(m_, &n_, 1, ((int)(K - 1) + 1));
            mccAllocateMatrix(&im, m_, n_);
            I_im = (mccM(&im) != 1 || mccN(&im) != 1);
            p_im = mccIPR(&im);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_im+=I_im)
                  {
                     *p_im = ((int)(1 + i_+j_*m_));
                  }
               }
            }
         }
         /* i0 = 2:(K+1); */
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            int *p_i0;
            int I_i0=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, ((int)((K + 1) - 2) + 1));
            m_ = mcmCalcResultSize(m_, &n_, 1, ((int)((K + 1) - 2) + 1));
            mccAllocateMatrix(&i0, m_, n_);
            I_i0 = (mccM(&i0) != 1 || mccN(&i0) != 1);
            p_i0 = mccIPR(&i0);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_i0+=I_i0)
                  {
                     *p_i0 = ((int)(2 + i_+j_*m_));
                  }
               }
            }
         }
         /* ip = 3:(K +2); */
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            int *p_ip;
            int I_ip=1;
            m_ = mcmCalcResultSize(m_, &n_, 1, ((int)((K + 2) - 3) + 1));
            m_ = mcmCalcResultSize(m_, &n_, 1, ((int)((K + 2) - 3) + 1));
            mccAllocateMatrix(&ip, m_, n_);
            I_ip = (mccM(&ip) != 1 || mccN(&ip) != 1);
            p_ip = mccIPR(&ip);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_ip+=I_ip)
                  {
                     *p_ip = ((int)(3 + i_+j_*m_));
                  }
               }
            }
         }
         
         /* dfpdx   = (fpp(ip,:) - fpp(im,:)) ./ (2 .* dx); */
         mccFindIndex(&IM8_, &ip);
         mccFindIndex(&IM9_, &im);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_dfpdx;
            int I_dfpdx=1;
            double *p_fpp;
            int I_fpp=1, J_fpp;
            int *p_IM8_;
            int I_IM8_=1;
            double *p_1fpp;
            int I_1fpp=1, J_1fpp;
            int *p_IM9_;
            int I_IM9_=1;
            double *p_dx;
            int I_dx=1;
            m_ = mcmCalcResultSize(m_, &n_, (mccM(&IM8_) * mccN(&IM8_)), mccN(&fpp));
            m_ = mcmCalcResultSize(m_, &n_, (mccM(&IM9_) * mccN(&IM9_)), mccN(&fpp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&dx), mccN(&dx));
            mccAllocateMatrix(&dfpdx, m_, n_);
            mccCheckMatrixSize(&fpp, mccGetMaxIndex(&IM8_ ,mccM(&fpp)), mccN(&fpp));
            mccCheckMatrixSize(&fpp, mccGetMaxIndex(&IM9_ ,mccM(&fpp)), mccN(&fpp));
            I_dfpdx = (mccM(&dfpdx) != 1 || mccN(&dfpdx) != 1);
            p_dfpdx = mccPR(&dfpdx);
            J_fpp = ((mccM(&fpp) != 1 || mccN(&fpp) != 1) ? mccM(&fpp) : 0);
            p_fpp = mccPR(&fpp) + mccM(&fpp) * 0;
            I_IM8_ = (mccM(&IM8_) != 1 || mccN(&IM8_) != 1);
            J_1fpp = ((mccM(&fpp) != 1 || mccN(&fpp) != 1) ? mccM(&fpp) : 0);
            p_1fpp = mccPR(&fpp) + mccM(&fpp) * 0;
            I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
            I_dx = (mccM(&dx) != 1 || mccN(&dx) != 1);
            p_dx = mccPR(&dx);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_fpp += J_fpp, p_1fpp += J_1fpp)
               {
                  p_IM8_ = mccIPR(&IM8_);
                  p_IM9_ = mccIPR(&IM9_);
                  for (i_=0; i_<m_; ++i_, p_dfpdx+=I_dfpdx, p_IM8_+=I_IM8_, p_IM9_+=I_IM9_, p_dx+=I_dx)
                  {
                     *p_dfpdx = ((p_fpp[((int)(*p_IM8_ - .5))] - p_1fpp[((int)(*p_IM9_ - .5))]) / (double) (2 * (double) *p_dx));
                  }
               }
            }
         }
         /* d2fpdx2 = (fpp(ip,:) - 2 .* fpp(i0,:) + fpp(im,:)) ./ (dx .^ 2);   */
         mccFindIndex(&IM9_, &ip);
         mccFindIndex(&IM8_, &i0);
         mccFindIndex(&IM7_, &im);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_d2fpdx2;
            int I_d2fpdx2=1;
            double *p_fpp;
            int I_fpp=1, J_fpp;
            int *p_IM9_;
            int I_IM9_=1;
            double *p_1fpp;
            int I_1fpp=1, J_1fpp;
            int *p_IM8_;
            int I_IM8_=1;
            double *p_2fpp;
            int I_2fpp=1, J_2fpp;
            int *p_IM7_;
            int I_IM7_=1;
            double *p_dx;
            int I_dx=1;
            m_ = mcmCalcResultSize(m_, &n_, (mccM(&IM9_) * mccN(&IM9_)), mccN(&fpp));
            m_ = mcmCalcResultSize(m_, &n_, (mccM(&IM8_) * mccN(&IM8_)), mccN(&fpp));
            m_ = mcmCalcResultSize(m_, &n_, (mccM(&IM7_) * mccN(&IM7_)), mccN(&fpp));
            m_ = mcmCalcResultSize(m_, &n_, mccM(&dx), mccN(&dx));
            mccAllocateMatrix(&d2fpdx2, m_, n_);
            mccCheckMatrixSize(&fpp, mccGetMaxIndex(&IM9_ ,mccM(&fpp)), mccN(&fpp));
            mccCheckMatrixSize(&fpp, mccGetMaxIndex(&IM8_ ,mccM(&fpp)), mccN(&fpp));
            mccCheckMatrixSize(&fpp, mccGetMaxIndex(&IM7_ ,mccM(&fpp)), mccN(&fpp));
            I_d2fpdx2 = (mccM(&d2fpdx2) != 1 || mccN(&d2fpdx2) != 1);
            p_d2fpdx2 = mccPR(&d2fpdx2);
            J_fpp = ((mccM(&fpp) != 1 || mccN(&fpp) != 1) ? mccM(&fpp) : 0);
            p_fpp = mccPR(&fpp) + mccM(&fpp) * 0;
            I_IM9_ = (mccM(&IM9_) != 1 || mccN(&IM9_) != 1);
            J_1fpp = ((mccM(&fpp) != 1 || mccN(&fpp) != 1) ? mccM(&fpp) : 0);
            p_1fpp = mccPR(&fpp) + mccM(&fpp) * 0;
            I_IM8_ = (mccM(&IM8_) != 1 || mccN(&IM8_) != 1);
            J_2fpp = ((mccM(&fpp) != 1 || mccN(&fpp) != 1) ? mccM(&fpp) : 0);
            p_2fpp = mccPR(&fpp) + mccM(&fpp) * 0;
            I_IM7_ = (mccM(&IM7_) != 1 || mccN(&IM7_) != 1);
            I_dx = (mccM(&dx) != 1 || mccN(&dx) != 1);
            p_dx = mccPR(&dx);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_fpp += J_fpp, p_1fpp += J_1fpp, p_2fpp += J_2fpp)
               {
                  p_IM9_ = mccIPR(&IM9_);
                  p_IM8_ = mccIPR(&IM8_);
                  p_IM7_ = mccIPR(&IM7_);
                  for (i_=0; i_<m_; ++i_, p_d2fpdx2+=I_d2fpdx2, p_IM9_+=I_IM9_, p_IM8_+=I_IM8_, p_IM7_+=I_IM7_, p_dx+=I_dx)
                  {
                     *p_d2fpdx2 = (((p_fpp[((int)(*p_IM9_ - .5))] - (2 * (double) p_1fpp[((int)(*p_IM8_ - .5))])) + p_2fpp[((int)(*p_IM7_ - .5))]) / (double) mcmRealPowerInt(*p_dx, 2));
                  }
               }
            }
         }
         
         /* end */
      }
      
      mccReturnFirstValue(&plhs_[0], &fp);
      mccReturnValue(&plhs_[1], &dfpdx);
      mccReturnValue(&plhs_[2], &d2fpdx2);
      mccReturnValue(&plhs_[3], &FS);
      mccReturnValue(&plhs_[4], &S);
      mccReturnValue(&plhs_[5], &ALPHA);
      mccReturnValue(&plhs_[6], &ALPHAP);
      mccReturnValue(&plhs_[7], &IDENTITE);
   }
   return;
}
