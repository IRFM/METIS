/*
 * MATLAB Compiler: 2.0.1
 * Date: Thu May 19 15:32:55 2005
 * Arguments: "-d" "/home/sccp/gsem/cgc/matlab5/zineb/v3.0/op" "-x" "zregular" 
 */
#include "zregular.h"

/*
 * The function "Mzregular" is the implementation version of the "zregular"
 * M-function from file "/home/sccp/gsem/cgc/matlab5/zineb/v3.0/op/zregular.m"
 * (lines 1-71). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * % ZREGULAR regularise un profil au centre
 * %---------------------------------------------------------------------
 * % fichier zregular.m ->  zregular
 * %
 * %
 * % fonction Matlab 5 : Attention fonction compile avec mcc
 * %
 * % Cette fonction regularise les profils au centre. Elle modifie les 4
 * % points les plus centraux. Elle utilise un polynome d'odre 3 de derivee
 * % nulle au centre et passant par les 3 points adjacants aux points modifies.
 * % Elle est utiliser pour les densite de courants calculer a partir des formules
 * % analytique ou au centre le numerateur et le denominateur d'un quotient tendent 
 * % vers 0. 
 * %
 * % Cette fonction ne marche que sur les donnees reels 
 * %  
 * % syntaxe  :
 * %  
 * %     yy=zregular(x,y);
 * %    
 * % entree :
 * %
 * %     x  = variable
 * %     y  = y(x) [vecteur]
 * %
 * % sorties :
 * % 
 * %     yy =  profil regularise au centre
 * % 
 * % ordre de compilation : mcc -V1.2 -riv zregular
 * %                        ( chemin vers le compilateur ->addpath /usr/local/matlab5/toolbox/compiler)
 * %     
 * % fonction ecrite par J-F Artaud , poste 46-78
 * % version 1.0, du 12/06/2001.
 * %
 * % liste des modifications : 
 * %
 * %    * 12/06/2001 -> optimisation et compilation
 * %
 * %--------------------------------------------------------------
 * %
 * function [y,a,b,c,d] = zregular(x,y)
 */
static mxArray * Mzregular(mxArray * * a,
                           mxArray * * b,
                           mxArray * * c,
                           mxArray * * d,
                           int nargout_,
                           mxArray * x,
                           mxArray * y_) {
    mxArray * y = mclGetUninitializedArray();
    mxArray * n = mclGetUninitializedArray();
    mxArray * x1 = mclGetUninitializedArray();
    mxArray * x2 = mclGetUninitializedArray();
    mxArray * x3 = mclGetUninitializedArray();
    mxArray * y1 = mclGetUninitializedArray();
    mxArray * y2 = mclGetUninitializedArray();
    mxArray * y3 = mclGetUninitializedArray();
    mclValidateInputs("zregular", 2, &x, &y_);
    mclCopyInputArg(&y, y_);
    /*
     * 
     * 
     * n =5;
     */
    mlfAssign(&n, mlfScalar(5.0));
    /*
     * 
     * x1 = x(n);
     */
    mlfAssign(&x1, mlfIndexRef(x, "(?)", n));
    /*
     * x2 = x(n+1);
     */
    mlfAssign(&x2, mlfIndexRef(x, "(?)", mlfPlus(n, mlfScalar(1.0))));
    /*
     * x3 = x(n+2);
     */
    mlfAssign(&x3, mlfIndexRef(x, "(?)", mlfPlus(n, mlfScalar(2.0))));
    /*
     * 
     * y1 = y(n);
     */
    mlfAssign(&y1, mlfIndexRef(y, "(?)", n));
    /*
     * y2 = y(n+1);
     */
    mlfAssign(&y2, mlfIndexRef(y, "(?)", mlfPlus(n, mlfScalar(1.0))));
    /*
     * y3 = y(n+2);
     */
    mlfAssign(&y3, mlfIndexRef(y, "(?)", mlfPlus(n, mlfScalar(2.0))));
    /*
     * 
     * 
     * a = -(x2.^2.*y3-x2.^2.*y1+x3.^2.*y1+y2.*x1.^2-y3.*x1.^2-y2.*x3.^2)./ ...
     */
    mlfAssign(
      a,
      mlfRdivide(
        mlfRdivide(
          mlfUminus(
            mlfMinus(
              mlfMinus(
                mlfPlus(
                  mlfPlus(
                    mlfMinus(
                      mlfTimes(mlfPower(x2, mlfScalar(2.0)), y3),
                      mlfTimes(mlfPower(x2, mlfScalar(2.0)), y1)),
                    mlfTimes(mlfPower(x3, mlfScalar(2.0)), y1)),
                  mlfTimes(y2, mlfPower(x1, mlfScalar(2.0)))),
                mlfTimes(y3, mlfPower(x1, mlfScalar(2.0)))),
              mlfTimes(y2, mlfPower(x3, mlfScalar(2.0))))),
          mlfMinus(x2, x3)),
        mlfMinus(
          mlfPlus(
            mlfMinus(
              mlfPlus(
                mlfMinus(
                  mlfTimes(
                    mlfPower(x2, mlfScalar(2.0)), mlfPower(x3, mlfScalar(2.0))),
                  mlfTimes(
                    mlfPower(x1, mlfScalar(2.0)),
                    mlfPower(x2, mlfScalar(2.0)))),
                mlfTimes(x2, mlfPower(x1, mlfScalar(3.0)))),
              mlfTimes(mlfTimes(x2, x3), mlfPower(x1, mlfScalar(2.0)))),
            mlfTimes(x3, mlfPower(x1, mlfScalar(3.0)))),
          mlfTimes(
            mlfPower(x3, mlfScalar(2.0)), mlfPower(x1, mlfScalar(2.0))))));
    /*
     * (x2-x3)./(x2.^2.*x3.^2-x1.^2.*x2.^2+x2.*x1.^3-x2.*x3.*x1.^2+x3.*x1.^3-x3.^2.*x1.^2);
     * b = (-x1.^3.*y3+x1.^3.*y2-x3.^3.*y2+y3.*x2.^3-y1.*x2.^3+y1.*x3.^3)./ ...
     */
    mlfAssign(
      b,
      mlfRdivide(
        mlfPlus(
          mlfMinus(
            mlfPlus(
              mlfMinus(
                mlfPlus(
                  mlfTimes(mlfUminus(mlfPower(x1, mlfScalar(3.0))), y3),
                  mlfTimes(mlfPower(x1, mlfScalar(3.0)), y2)),
                mlfTimes(mlfPower(x3, mlfScalar(3.0)), y2)),
              mlfTimes(y3, mlfPower(x2, mlfScalar(3.0)))),
            mlfTimes(y1, mlfPower(x2, mlfScalar(3.0)))),
          mlfTimes(y1, mlfPower(x3, mlfScalar(3.0)))),
        mlfPlus(
          mlfMinus(
            mlfPlus(
              mlfMinus(
                mlfMinus(
                  mlfTimes(
                    mlfPower(x1, mlfScalar(3.0)), mlfPower(x2, mlfScalar(2.0))),
                  mlfTimes(
                    mlfPower(x1, mlfScalar(3.0)),
                    mlfPower(x3, mlfScalar(2.0)))),
                mlfTimes(
                  mlfPower(x1, mlfScalar(2.0)), mlfPower(x2, mlfScalar(3.0)))),
              mlfTimes(
                mlfPower(x1, mlfScalar(2.0)), mlfPower(x3, mlfScalar(3.0)))),
            mlfTimes(
              mlfPower(x3, mlfScalar(3.0)), mlfPower(x2, mlfScalar(2.0)))),
          mlfTimes(
            mlfPower(x3, mlfScalar(2.0)), mlfPower(x2, mlfScalar(3.0))))));
    /*
     * (x1.^3.*x2.^2-x1.^3.*x3.^2-x1.^2.*x2.^3+x1.^2.*x3.^3-x3.^3.*x2.^2+x3.^2.*x2.^3);
     * d = (x2.^2.*x1.^3.*y3-x3.^2.*x1.^3.*y2-y3.*x1.^2.*x2.^3+y2.*x1.^2.*x3.^3+ ...
     */
    mlfAssign(
      d,
      mlfRdivide(
        mlfMinus(
          mlfPlus(
            mlfPlus(
              mlfMinus(
                mlfMinus(
                  mlfTimes(
                    mlfTimes(
                      mlfPower(x2, mlfScalar(2.0)),
                      mlfPower(x1, mlfScalar(3.0))),
                    y3),
                  mlfTimes(
                    mlfTimes(
                      mlfPower(x3, mlfScalar(2.0)),
                      mlfPower(x1, mlfScalar(3.0))),
                    y2)),
                mlfTimes(
                  mlfTimes(y3, mlfPower(x1, mlfScalar(2.0))),
                  mlfPower(x2, mlfScalar(3.0)))),
              mlfTimes(
                mlfTimes(y2, mlfPower(x1, mlfScalar(2.0))),
                mlfPower(x3, mlfScalar(3.0)))),
            mlfTimes(
              mlfTimes(mlfPower(x3, mlfScalar(2.0)), y1),
              mlfPower(x2, mlfScalar(3.0)))),
          mlfTimes(
            mlfTimes(mlfPower(x2, mlfScalar(2.0)), y1),
            mlfPower(x3, mlfScalar(3.0)))),
        mlfPlus(
          mlfMinus(
            mlfPlus(
              mlfMinus(
                mlfMinus(
                  mlfTimes(
                    mlfPower(x1, mlfScalar(3.0)), mlfPower(x2, mlfScalar(2.0))),
                  mlfTimes(
                    mlfPower(x1, mlfScalar(3.0)),
                    mlfPower(x3, mlfScalar(2.0)))),
                mlfTimes(
                  mlfPower(x1, mlfScalar(2.0)), mlfPower(x2, mlfScalar(3.0)))),
              mlfTimes(
                mlfPower(x1, mlfScalar(2.0)), mlfPower(x3, mlfScalar(3.0)))),
            mlfTimes(
              mlfPower(x3, mlfScalar(3.0)), mlfPower(x2, mlfScalar(2.0)))),
          mlfTimes(
            mlfPower(x3, mlfScalar(2.0)), mlfPower(x2, mlfScalar(3.0))))));
    /*
     * x3.^2.*y1.*x2.^3-x2.^2.*y1.*x3.^3)./ ...
     * (x1.^3.*x2.^2-x1.^3.*x3.^2-x1.^2.*x2.^3+x1.^2.*x3.^3-x3.^3.*x2.^2+x3.^2.*x2.^3);
     * %c = 0;
     * 
     * %p = [a,b,c,d];
     * 
     * %yy = y;
     * y(1:(n-1)) = a .* x(1:(n-1)) .^ 3 + b .* x(1:(n-1)) .^ 2 + d;
     */
    mlfIndexAssign(
      &y,
      "(?)",
      mlfColon(mlfScalar(1.0), mlfMinus(n, mlfScalar(1.0)), NULL),
      mlfPlus(
        mlfPlus(
          mlfTimes(
            *a,
            mlfPower(
              mlfIndexRef(
                x,
                "(?)",
                mlfColon(mlfScalar(1.0), mlfMinus(n, mlfScalar(1.0)), NULL)),
              mlfScalar(3.0))),
          mlfTimes(
            *b,
            mlfPower(
              mlfIndexRef(
                x,
                "(?)",
                mlfColon(mlfScalar(1.0), mlfMinus(n, mlfScalar(1.0)), NULL)),
              mlfScalar(2.0)))),
        *d));
    mclValidateOutputs("zregular", 5, nargout_, &y, a, b, c, d);
    mxDestroyArray(n);
    mxDestroyArray(x1);
    mxDestroyArray(x2);
    mxDestroyArray(x3);
    mxDestroyArray(y1);
    mxDestroyArray(y2);
    mxDestroyArray(y3);
    /*
     * 
     * %yy(1:(n-1)) = polyval(p,x(1:(n-1)));
     */
    return y;
}

/*
 * The function "mlfZregular" contains the normal interface for the "zregular"
 * M-function from file "/home/sccp/gsem/cgc/matlab5/zineb/v3.0/op/zregular.m"
 * (lines 1-71). This function processes any input arguments and passes them to
 * the implementation version of the function, appearing above.
 */
mxArray * mlfZregular(mxArray * * a,
                      mxArray * * b,
                      mxArray * * c,
                      mxArray * * d,
                      mxArray * x,
                      mxArray * y_) {
    int nargout = 1;
    mxArray * y = mclGetUninitializedArray();
    mxArray * a__ = mclGetUninitializedArray();
    mxArray * b__ = mclGetUninitializedArray();
    mxArray * c__ = mclGetUninitializedArray();
    mxArray * d__ = mclGetUninitializedArray();
    mlfEnterNewContext(4, 2, a, b, c, d, x, y_);
    if (a != NULL) {
        ++nargout;
    }
    if (b != NULL) {
        ++nargout;
    }
    if (c != NULL) {
        ++nargout;
    }
    if (d != NULL) {
        ++nargout;
    }
    y = Mzregular(&a__, &b__, &c__, &d__, nargout, x, y_);
    mlfRestorePreviousContext(4, 2, a, b, c, d, x, y_);
    if (a != NULL) {
        mclCopyOutputArg(a, a__);
    } else {
        mxDestroyArray(a__);
    }
    if (b != NULL) {
        mclCopyOutputArg(b, b__);
    } else {
        mxDestroyArray(b__);
    }
    if (c != NULL) {
        mclCopyOutputArg(c, c__);
    } else {
        mxDestroyArray(c__);
    }
    if (d != NULL) {
        mclCopyOutputArg(d, d__);
    } else {
        mxDestroyArray(d__);
    }
    return mlfReturnValue(y);
}

/*
 * The function "mlxZregular" contains the feval interface for the "zregular"
 * M-function from file "/home/sccp/gsem/cgc/matlab5/zineb/v3.0/op/zregular.m"
 * (lines 1-71). The feval function calls the implementation version of
 * zregular through this function. This function processes any input arguments
 * and passes them to the implementation version of the function, appearing
 * above.
 */
void mlxZregular(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[5];
    int i;
    if (nlhs > 5) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: zregular Line: 42 Column:"
            " 0 The function \"zregular\" was called with mo"
            "re than the declared number of outputs (5)"));
    }
    if (nrhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: zregular Line: 42 Column"
            ": 0 The function \"zregular\" was called with "
            "more than the declared number of inputs (2)"));
    }
    for (i = 0; i < 5; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 2 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 2; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 2, mprhs[0], mprhs[1]);
    mplhs[0]
      = Mzregular(
          &mplhs[1],
          &mplhs[2],
          &mplhs[3],
          &mplhs[4],
          nlhs,
          mprhs[0],
          mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 5 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 5; ++i) {
        mxDestroyArray(mplhs[i]);
    }
}
