/*
 * MATLAB Compiler: 2.0.1
 * Date: Thu May 19 15:32:59 2005
 * Arguments: "-d" "/home/sccp/gsem/cgc/matlab5/zineb/v3.0/solver" "-x"
 * "rpdederive" 
 */
#include "rpdederive.h"

/*
 * The function "Mrpdederive" is the implementation version of the "rpdederive"
 * M-function from file
 * "/home/sccp/gsem/cgc/matlab5/zineb/v3.0/solver/rpdederive.m" (lines 1-349).
 * It contains the actual compiled code for that M-function. It is a static
 * function and must only be called from one of the interface functions,
 * appearing below.
 */
/*
 * % RPDEDERIVE calcule la derivee d'une grandeur y (reel) par rapport a x (reel)
 * %----------------------------------------------------------------
 * % fichier rpdederive.m -> rpdederive
 * %
 * % Fonction Matlab 5 : attention compilee avec mcc -> ne marche que sur des reel
 * %
 * % Cette fonction calcule la derivee de y (reel) par rapport a x (reel) sur 3 points
 * % en fixant des conditions de prolongation aux bornes de l'intervalle. 
 * % La derivee est effectuee selon la dimension d des matrices.
 * % 
 * % voir aussi pdederive.m
 * %
 * % syntaxe :
 * % 
 * %    dydx = rpdederive(x,y,c0,c1,d,o); 
 * %
 * % entree :
 * %        
 * %   x  = matrice reelle  des X  de memes dimensions que y ou verifiant
 * %        length(x)=K. (les X, le long de la dimension de derivation
 * %        doivent etre equidistants pour obtenir une prolongation au 
 * %        bord de pra continuite precise)   
 * % 
 * %   y  = matrice  reelle des Y de dimension dim = size (y), 
 * %        dim etant un vecteur aynt au moins d elements.
 * %        et dim(d)=K
 * %   
 * %   c0 = matrice donnant les conditions de prolongation k=1,
 * %        k etant l'indice de x et y sur la dimension d. c0 
 * %        peut etre un scalaire (condition uniforme pour tous
 * %        les points) ou une matrice de memes dimension que x, 
 * %        sauf pour la dimsion d :
 * %            dim0 =size(c0), dim0(d)=1.
 * %        Cette matrice code pour des conditions de prolongation
 * %        variables selon les points.Les valeur de c0 possible sont :
 * %          * pour la derive premiere :
 * %                 0 -> derivee nulle en k=1
 * %                 1 -> derivee est calculee sur deux points  en k=1 
 * %                 2 -> prolongation par continuite de la derivee 
 * %                      d'ordre 3 en k=1 
 * %          * pour la derive seconde :
 * %                 0 -> derivee nulle en k=1
 * %                 1 -> derivee 1ere nulle en k=K 
 * %                 2 -> prolongation par continuite de la derivee 
 * %                      d'ordre 3  en k=1      
 * %        
 * %   c1 = matrice donnant les conditions de prolongation k=K,
 * %        k etant l'indice de x et y sur la dimension d. c1
 * %        peut etre un scalaire (condition uniforme pour tous
 * %        les points) ou une matrice de memes dimension que x, 
 * %        sauf pour la dimsion d :
 * %            dim1 =size(c1), dim1(d)=1.
 * %        Cette matrice code pour des conditions de prolongation
 * %        variables selon les points. Les valeur de c1 possible sont :
 * %          * pour la derive premiere :
 * %                 0 -> derivee nulle en k=K
 * %                 1 -> derivee est calculee sur deux points  en k=K 
 * %                 2 -> prolongation par continuite de la derivee 
 * %                      d'ordre 3 en k=K 
 * %          * pour la derive seconde :
 * %                 0 -> derivee 2ieme nulle en k=K
 * %                 1 -> derivee 1ere nulle en k=K 
 * %                 2 -> prolongation par continuite de la derivee 
 * %                      d'ordre 3  en k=K      
 * %
 * %   d = dimension d'espace sur laquelle porte la derivee
 * %
 * %   o = ordre de la derivee (1 ou 2)
 * %   
 * % sortie : 
 * %
 * %   dydx = derivee de y par rapport a x (memes dimensions que x)
 * %
 * % exemple :
 * % 
 * %   x=linspace(0,pi/2);
 * %   y=cos(x);
 * %   dydx=rpdederive(x,y,0,1,2);
 * %   plot(x,-sin(x),'o',x,dydx,'+');
 * %
 * % ordre de compilation : mcc -V1.2 -r rpdederive
 * %                        ( chemin vers le compilateur ->addpath /usr/local/matlab5/toolbox/compiler)
 * %
 * % fonction ecrite par J-F Artaud , poste 46-78
 * % version 1.5, du 14/06/2001.
 * % 
 * % 
 * % liste des modifications : 
 * %
 * %--------------------------------------------------------------
 * %
 * 
 * function dydx=derive(x,y,c0,c1,d,o)
 */
static mxArray * Mrpdederive(int nargout_,
                             mxArray * x,
                             mxArray * y,
                             mxArray * c0,
                             mxArray * c1,
                             mxArray * d,
                             mxArray * o) {
    mxArray * dydx = mclGetUninitializedArray();
    mxArray * K = mclGetUninitializedArray();
    mxArray * comp = mclGetUninitializedArray();
    mxArray * dx = mclGetUninitializedArray();
    mxArray * dy = mclGetUninitializedArray();
    mxArray * k0 = mclGetUninitializedArray();
    mxArray * k1 = mclGetUninitializedArray();
    mxArray * mode = mclGetUninitializedArray();
    mxArray * n = mclGetUninitializedArray();
    mxArray * nargin_ = mclGetUninitializedArray();
    mxArray * ts = mclGetUninitializedArray();
    mxArray * x0 = mclGetUninitializedArray();
    mxArray * x1 = mclGetUninitializedArray();
    mxArray * y0 = mclGetUninitializedArray();
    mxArray * y1 = mclGetUninitializedArray();
    mlfAssign(&nargin_, mlfNargin(0, x, y, c0, c1, d, o, NULL));
    mclValidateInputs("rpdederive", 6, &x, &y, &c0, &c1, &d, &o);
    mclCopyArray(&x);
    mclCopyArray(&y);
    mclCopyArray(&c0);
    mclCopyArray(&c1);
    mclCopyArray(&o);
    /*
     * 
     * % pas de verification des dimensions
     * if nargin == 5
     */
    if (mlfTobool(mlfEq(nargin_, mlfScalar(5.0)))) {
        /*
         * o=1;
         */
        mlfAssign(&o, mlfScalar(1.0));
    /*
     * end
     */
    }
    /*
     * % nb dimension des matrices
     * ts = size(y);
     */
    mlfAssign(&ts, mlfSize(mclValueVarargout(), y, NULL));
    /*
     * n = length(ts);
     */
    mlfAssign(&n, mlfLength(ts));
    /*
     * 
     * % mode 1 (x est un vecteur)
     * if length(x) == prod(size(x))
     */
    if (mlfTobool(
          mlfEq(
            mlfLength(x),
            mlfProd(mlfSize(mclValueVarargout(), x, NULL), NULL)))) {
        /*
         * mode = 1;
         */
        mlfAssign(&mode, mlfScalar(1.0));
        /*
         * x    = x(:);
         */
        mlfAssign(&x, mlfIndexRef(x, "(?)", mlfCreateColonIndex()));
    /*
     * else
     */
    } else {
        /*
         * mode = 0;
         */
        mlfAssign(&mode, mlfScalar(0.0));
    /*
     * end
     */
    }
    /*
     * 
     * % permutation si on ne derive pas selon la dimension 1 
     * if n == 2 & d == 2
     */
    {
        mxArray * a_ = mclInitialize(mlfEq(n, mlfScalar(2.0)));
        if (mlfTobool(a_) && mlfTobool(mlfAnd(a_, mlfEq(d, mlfScalar(2.0))))) {
            mxDestroyArray(a_);
            /*
             * if mode == 0
             */
            if (mlfTobool(mlfEq(mode, mlfScalar(0.0)))) {
                /*
                 * x = x.';
                 */
                mlfAssign(&x, mlfTranspose(x));
            /*
             * end
             */
            }
            /*
             * y  = y.';
             */
            mlfAssign(&y, mlfTranspose(y));
            /*
             * c0 = c0.';
             */
            mlfAssign(&c0, mlfTranspose(c0));
            /*
             * c1 = c1.';
             */
            mlfAssign(&c1, mlfTranspose(c1));
        /*
         * 
         * elseif d>1
         */
        } else {
            mxDestroyArray(a_);
            if (mlfTobool(mlfGt(d, mlfScalar(1.0)))) {
                /*
                 * if mode == 0
                 */
                if (mlfTobool(mlfEq(mode, mlfScalar(0.0)))) {
                    /*
                     * x = shiftdim(x,d-1);
                     */
                    mlfAssign(
                      &x, mlfShiftdim(NULL, x, mlfMinus(d, mlfScalar(1.0))));
                /*
                 * end
                 */
                }
                /*
                 * y = shiftdim(y,d-1);
                 */
                mlfAssign(
                  &y, mlfShiftdim(NULL, y, mlfMinus(d, mlfScalar(1.0))));
                /*
                 * c0 = shiftdim(c0,d-1);
                 */
                mlfAssign(
                  &c0, mlfShiftdim(NULL, c0, mlfMinus(d, mlfScalar(1.0))));
                /*
                 * c1 = shiftdim(c1,d-1);
                 */
                mlfAssign(
                  &c1, mlfShiftdim(NULL, c1, mlfMinus(d, mlfScalar(1.0))));
            }
        }
    /*
     * 
     * end
     */
    }
    /*
     * 
     * % dimension d'espace 
     * K  = size(x,1);
     */
    mlfAssign(&K, mlfSize(mclValueVarargout(), x, mlfScalar(1.0)));
    /*
     * k0 = 1:(K-1);
     */
    mlfAssign(&k0, mlfColon(mlfScalar(1.0), mlfMinus(K, mlfScalar(1.0)), NULL));
    /*
     * k1 = 2:K;
     */
    mlfAssign(&k1, mlfColon(mlfScalar(2.0), K, NULL));
    /*
     * 
     * % calcule de dx
     * % si x est un vecteur
     * if mode == 1 
     */
    if (mlfTobool(mlfEq(mode, mlfScalar(1.0)))) {
        /*
         * % prolongement au bord
         * x1 = 2 .* x(K) - x(K-1);
         */
        mlfAssign(
          &x1,
          mlfMinus(
            mlfTimes(mlfScalar(2.0), mlfIndexRef(x, "(?)", K)),
            mlfIndexRef(x, "(?)", mlfMinus(K, mlfScalar(1.0)))));
        /*
         * % prolongement au centre
         * x0 = 2 .* x(1) - x(2);
         */
        mlfAssign(
          &x0,
          mlfMinus(
            mlfTimes(mlfScalar(2.0), mlfIndexRef(x, "(?)", mlfScalar(1.0))),
            mlfIndexRef(x, "(?)", mlfScalar(2.0))));
        /*
         * % concatenation des matrices
         * dx= (cat(1,x(k1),x1) - cat(1,x0,x(k0))) ./ 2;
         */
        mlfAssign(
          &dx,
          mlfRdivide(
            mlfMinus(
              mlfCat(mlfScalar(1.0), mlfIndexRef(x, "(?)", k1), x1, NULL),
              mlfCat(mlfScalar(1.0), x0, mlfIndexRef(x, "(?)", k0), NULL)),
            mlfScalar(2.0)));
        /*
         * 
         * % complement pour la division
         * if n == 1
         */
        if (mlfTobool(mlfEq(n, mlfScalar(1.0)))) {
            /*
             * comp = 1;
             */
            mlfAssign(&comp, mlfScalar(1.0));
        /*
         * elseif n == 2
         */
        } else if (mlfTobool(mlfEq(n, mlfScalar(2.0)))) {
            /*
             * comp = ones(1,size(y,2));
             */
            mlfAssign(
              &comp,
              mlfOnes(
                mlfScalar(1.0),
                mlfSize(mclValueVarargout(), y, mlfScalar(2.0)),
                NULL));
        /*
         * elseif n == 3 
         */
        } else if (mlfTobool(mlfEq(n, mlfScalar(3.0)))) {
            /*
             * comp = ones(1,size(y,2),size(y,3));
             */
            mlfAssign(
              &comp,
              mlfOnes(
                mlfScalar(1.0),
                mlfSize(mclValueVarargout(), y, mlfScalar(2.0)),
                mlfSize(mclValueVarargout(), y, mlfScalar(3.0)),
                NULL));
        /*
         * 
         * elseif n == 4 
         */
        } else if (mlfTobool(mlfEq(n, mlfScalar(4.0)))) {
            /*
             * comp = ones(1,size(y,2),size(y,3),size(y,4));
             */
            mlfAssign(
              &comp,
              mlfOnes(
                mlfScalar(1.0),
                mlfSize(mclValueVarargout(), y, mlfScalar(2.0)),
                mlfSize(mclValueVarargout(), y, mlfScalar(3.0)),
                mlfSize(mclValueVarargout(), y, mlfScalar(4.0)),
                NULL));
        /*
         * 
         * else 
         */
        } else {
            /*
             * error('Pas encore implante - a vouds de l''ecrire')
             */
            mlfError(
              mxCreateString("Pas encore implante - a vouds de l'ecrire"));
        /*
         * 
         * end
         */
        }
    /*
     * 
     * else
     */
    } else {
        /*
         * % selon le nombre de dimension (1 a 4)
         * if n == 1
         */
        if (mlfTobool(mlfEq(n, mlfScalar(1.0)))) {
            /*
             * % prolongement au bord
             * x1 = 2 .* x(K) - x(K-1);
             */
            mlfAssign(
              &x1,
              mlfMinus(
                mlfTimes(mlfScalar(2.0), mlfIndexRef(x, "(?)", K)),
                mlfIndexRef(x, "(?)", mlfMinus(K, mlfScalar(1.0)))));
            /*
             * % prolongement au centre
             * x0 = 2 .* x(1) - x(2);
             */
            mlfAssign(
              &x0,
              mlfMinus(
                mlfTimes(mlfScalar(2.0), mlfIndexRef(x, "(?)", mlfScalar(1.0))),
                mlfIndexRef(x, "(?)", mlfScalar(2.0))));
            /*
             * % concatenation des matrices
             * dx= (cat(1,x(k1),x1) - cat(1,x0,x(k0))) ./ 2;
             */
            mlfAssign(
              &dx,
              mlfRdivide(
                mlfMinus(
                  mlfCat(mlfScalar(1.0), mlfIndexRef(x, "(?)", k1), x1, NULL),
                  mlfCat(mlfScalar(1.0), x0, mlfIndexRef(x, "(?)", k0), NULL)),
                mlfScalar(2.0)));
        /*
         * 
         * elseif n == 2
         */
        } else if (mlfTobool(mlfEq(n, mlfScalar(2.0)))) {
            /*
             * % prolongement au bord
             * x1 = 2 .* x(K,:) - x(K-1,:);
             */
            mlfAssign(
              &x1,
              mlfMinus(
                mlfTimes(
                  mlfScalar(2.0),
                  mlfIndexRef(x, "(?,?)", K, mlfCreateColonIndex())),
                mlfIndexRef(
                  x,
                  "(?,?)",
                  mlfMinus(K, mlfScalar(1.0)),
                  mlfCreateColonIndex())));
            /*
             * % prolongement au centre
             * x0 = 2 .* x(1,:) - x(2,:);
             */
            mlfAssign(
              &x0,
              mlfMinus(
                mlfTimes(
                  mlfScalar(2.0),
                  mlfIndexRef(
                    x, "(?,?)", mlfScalar(1.0), mlfCreateColonIndex())),
                mlfIndexRef(
                  x, "(?,?)", mlfScalar(2.0), mlfCreateColonIndex())));
            /*
             * % concatenation des matrices
             * dx= (cat(1,x(k1,:),x1) - cat(1,x0,x(k0,:))) ./ 2;
             */
            mlfAssign(
              &dx,
              mlfRdivide(
                mlfMinus(
                  mlfCat(
                    mlfScalar(1.0),
                    mlfIndexRef(x, "(?,?)", k1, mlfCreateColonIndex()),
                    x1,
                    NULL),
                  mlfCat(
                    mlfScalar(1.0),
                    x0,
                    mlfIndexRef(x, "(?,?)", k0, mlfCreateColonIndex()),
                    NULL)),
                mlfScalar(2.0)));
        /*
         * 
         * elseif n == 3
         */
        } else if (mlfTobool(mlfEq(n, mlfScalar(3.0)))) {
            /*
             * % prolongement au bord
             * x1 = 2 .* x(K,:,:) - x(K-1,:,:);
             */
            mlfAssign(
              &x1,
              mlfMinus(
                mlfTimes(
                  mlfScalar(2.0),
                  mlfIndexRef(
                    x,
                    "(?,?,?)",
                    K,
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex())),
                mlfIndexRef(
                  x,
                  "(?,?,?)",
                  mlfMinus(K, mlfScalar(1.0)),
                  mlfCreateColonIndex(),
                  mlfCreateColonIndex())));
            /*
             * % prolongement au centre
             * x0 = 2 .* x(1,:,:) - x(2,:,:);
             */
            mlfAssign(
              &x0,
              mlfMinus(
                mlfTimes(
                  mlfScalar(2.0),
                  mlfIndexRef(
                    x,
                    "(?,?,?)",
                    mlfScalar(1.0),
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex())),
                mlfIndexRef(
                  x,
                  "(?,?,?)",
                  mlfScalar(2.0),
                  mlfCreateColonIndex(),
                  mlfCreateColonIndex())));
            /*
             * % concatenation des matrices
             * dx= (cat(1,x(k1,:,:),x1) - cat(1,x0,x(k0,:,:))) ./ 2;
             */
            mlfAssign(
              &dx,
              mlfRdivide(
                mlfMinus(
                  mlfCat(
                    mlfScalar(1.0),
                    mlfIndexRef(
                      x,
                      "(?,?,?)",
                      k1,
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex()),
                    x1,
                    NULL),
                  mlfCat(
                    mlfScalar(1.0),
                    x0,
                    mlfIndexRef(
                      x,
                      "(?,?,?)",
                      k0,
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex()),
                    NULL)),
                mlfScalar(2.0)));
        /*
         * 
         * elseif n == 4
         */
        } else if (mlfTobool(mlfEq(n, mlfScalar(4.0)))) {
            /*
             * % prolongement au bord
             * x1 = 2 .* x(K,:,:,:) - x(K-1,:,:,:);
             */
            mlfAssign(
              &x1,
              mlfMinus(
                mlfTimes(
                  mlfScalar(2.0),
                  mlfIndexRef(
                    x,
                    "(?,?,?,?)",
                    K,
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex())),
                mlfIndexRef(
                  x,
                  "(?,?,?,?)",
                  mlfMinus(K, mlfScalar(1.0)),
                  mlfCreateColonIndex(),
                  mlfCreateColonIndex(),
                  mlfCreateColonIndex())));
            /*
             * % prolongement au centre
             * x0 = 2 .* x(1,:,:,:) - x(2,:,:,:);
             */
            mlfAssign(
              &x0,
              mlfMinus(
                mlfTimes(
                  mlfScalar(2.0),
                  mlfIndexRef(
                    x,
                    "(?,?,?,?)",
                    mlfScalar(1.0),
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex())),
                mlfIndexRef(
                  x,
                  "(?,?,?,?)",
                  mlfScalar(2.0),
                  mlfCreateColonIndex(),
                  mlfCreateColonIndex(),
                  mlfCreateColonIndex())));
            /*
             * % concatenation des matrices
             * dx= (cat(1,x(k1,:,:,:),x1) - cat(1,x0,x(k0,:,:,:))) ./ 2;
             */
            mlfAssign(
              &dx,
              mlfRdivide(
                mlfMinus(
                  mlfCat(
                    mlfScalar(1.0),
                    mlfIndexRef(
                      x,
                      "(?,?,?,?)",
                      k1,
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex()),
                    x1,
                    NULL),
                  mlfCat(
                    mlfScalar(1.0),
                    x0,
                    mlfIndexRef(
                      x,
                      "(?,?,?,?)",
                      k0,
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex()),
                    NULL)),
                mlfScalar(2.0)));
        /*
         * 
         * else 
         */
        } else {
            /*
             * error('Pas encore implante - a vouds de l''ecrire')
             */
            mlfError(
              mxCreateString("Pas encore implante - a vouds de l'ecrire"));
        /*
         * 
         * end	
         */
        }
    /*
     * 
     * end
     */
    }
    /*
     * 
     * % calcule de dy 
     * % en fonction de l'ordre de la derivee
     * if o ==1
     */
    if (mlfTobool(mlfEq(o, mlfScalar(1.0)))) {
        /*
         * % selon le nombre de dimension (1 a 4)
         * if n == 1
         */
        if (mlfTobool(mlfEq(n, mlfScalar(1.0)))) {
            /*
             * % prolongement au bord
             * y1 = (c1~=2) .* (2 .* y(K) - y(K-1)) + ...
             */
            mlfAssign(
              &y1,
              mlfPlus(
                mlfTimes(
                  mlfNe(c1, mlfScalar(2.0)),
                  mlfMinus(
                    mlfTimes(mlfScalar(2.0), mlfIndexRef(y, "(?)", K)),
                    mlfIndexRef(y, "(?)", mlfMinus(K, mlfScalar(1.0))))),
                mlfTimes(
                  mlfEq(c1, mlfScalar(2.0)),
                  mlfMinus(
                    mlfPlus(
                      mlfMinus(
                        mlfTimes(mlfScalar(4.0), mlfIndexRef(y, "(?)", K)),
                        mlfTimes(
                          mlfScalar(6.0),
                          mlfIndexRef(y, "(?)", mlfMinus(K, mlfScalar(1.0))))),
                      mlfTimes(
                        mlfScalar(4.0),
                        mlfIndexRef(y, "(?)", mlfMinus(K, mlfScalar(2.0))))),
                    mlfIndexRef(y, "(?)", mlfMinus(K, mlfScalar(3.0)))))));
            /*
             * (c1==2) .* (4 .* y(K) - 6 .* y(K-1) + 4 .* y(K-2) - y(K-3));
             * % prolongement au centre
             * y0 = (c0~=2)  .* (2 .* y(1) - y(2)) +  ...
             */
            mlfAssign(
              &y0,
              mlfPlus(
                mlfTimes(
                  mlfNe(c0, mlfScalar(2.0)),
                  mlfMinus(
                    mlfTimes(
                      mlfScalar(2.0), mlfIndexRef(y, "(?)", mlfScalar(1.0))),
                    mlfIndexRef(y, "(?)", mlfScalar(2.0)))),
                mlfTimes(
                  mlfEq(c0, mlfScalar(2.0)),
                  mlfMinus(
                    mlfPlus(
                      mlfMinus(
                        mlfTimes(
                          mlfScalar(4.0),
                          mlfIndexRef(y, "(?)", mlfScalar(1.0))),
                        mlfTimes(
                          mlfScalar(6.0),
                          mlfIndexRef(y, "(?)", mlfScalar(2.0)))),
                      mlfTimes(
                        mlfScalar(4.0), mlfIndexRef(y, "(?)", mlfScalar(3.0)))),
                    mlfIndexRef(y, "(?)", mlfScalar(4.0))))));
            /*
             * (c0==2) .* (4 .* y(1) - 6 .* y(2) + 4 .* y(3) - y(4));
             * % concatenation des matrices
             * dy=cat(1,y(k1),y1)-cat(1,y0,y(k0));
             */
            mlfAssign(
              &dy,
              mlfMinus(
                mlfCat(mlfScalar(1.0), mlfIndexRef(y, "(?)", k1), y1, NULL),
                mlfCat(mlfScalar(1.0), y0, mlfIndexRef(y, "(?)", k0), NULL)));
        /*
         * 
         * elseif n == 2
         */
        } else if (mlfTobool(mlfEq(n, mlfScalar(2.0)))) {
            /*
             * % prolongement au bord
             * y1 = (c1~=2) .* (2 .* y(K,:) - y(K-1,:)) + ...
             */
            mlfAssign(
              &y1,
              mlfPlus(
                mlfTimes(
                  mlfNe(c1, mlfScalar(2.0)),
                  mlfMinus(
                    mlfTimes(
                      mlfScalar(2.0),
                      mlfIndexRef(y, "(?,?)", K, mlfCreateColonIndex())),
                    mlfIndexRef(
                      y,
                      "(?,?)",
                      mlfMinus(K, mlfScalar(1.0)),
                      mlfCreateColonIndex()))),
                mlfTimes(
                  mlfEq(c1, mlfScalar(2.0)),
                  mlfMinus(
                    mlfPlus(
                      mlfMinus(
                        mlfTimes(
                          mlfScalar(4.0),
                          mlfIndexRef(y, "(?,?)", K, mlfCreateColonIndex())),
                        mlfTimes(
                          mlfScalar(6.0),
                          mlfIndexRef(
                            y,
                            "(?,?)",
                            mlfMinus(K, mlfScalar(1.0)),
                            mlfCreateColonIndex()))),
                      mlfTimes(
                        mlfScalar(4.0),
                        mlfIndexRef(
                          y,
                          "(?,?)",
                          mlfMinus(K, mlfScalar(2.0)),
                          mlfCreateColonIndex()))),
                    mlfIndexRef(
                      y,
                      "(?,?)",
                      mlfMinus(K, mlfScalar(3.0)),
                      mlfCreateColonIndex())))));
            /*
             * (c1==2) .* (4 .* y(K,:) - 6 .* y(K-1,:) + 4 .* y(K-2,:) - y(K-3,:));
             * % prolongement au centre
             * y0 = (c0~=2)  .* (2 .* y(1,:) - y(2,:)) +  ...
             */
            mlfAssign(
              &y0,
              mlfPlus(
                mlfTimes(
                  mlfNe(c0, mlfScalar(2.0)),
                  mlfMinus(
                    mlfTimes(
                      mlfScalar(2.0),
                      mlfIndexRef(
                        y, "(?,?)", mlfScalar(1.0), mlfCreateColonIndex())),
                    mlfIndexRef(
                      y, "(?,?)", mlfScalar(2.0), mlfCreateColonIndex()))),
                mlfTimes(
                  mlfEq(c0, mlfScalar(2.0)),
                  mlfMinus(
                    mlfPlus(
                      mlfMinus(
                        mlfTimes(
                          mlfScalar(4.0),
                          mlfIndexRef(
                            y, "(?,?)", mlfScalar(1.0), mlfCreateColonIndex())),
                        mlfTimes(
                          mlfScalar(6.0),
                          mlfIndexRef(
                            y,
                            "(?,?)",
                            mlfScalar(2.0),
                            mlfCreateColonIndex()))),
                      mlfTimes(
                        mlfScalar(4.0),
                        mlfIndexRef(
                          y, "(?,?)", mlfScalar(3.0), mlfCreateColonIndex()))),
                    mlfIndexRef(
                      y, "(?,?)", mlfScalar(4.0), mlfCreateColonIndex())))));
            /*
             * (c0==2) .* (4 .* y(1,:) - 6 .* y(2,:) + 4 .* y(3,:) - y(4,:));
             * % concatenation des matrices
             * dy=cat(1,y(k1,:),y1)-cat(1,y0,y(k0,:));
             */
            mlfAssign(
              &dy,
              mlfMinus(
                mlfCat(
                  mlfScalar(1.0),
                  mlfIndexRef(y, "(?,?)", k1, mlfCreateColonIndex()), y1, NULL),
                mlfCat(
                  mlfScalar(1.0),
                  y0,
                  mlfIndexRef(y, "(?,?)", k0, mlfCreateColonIndex()),
                  NULL)));
        /*
         * 
         * elseif n == 3
         */
        } else if (mlfTobool(mlfEq(n, mlfScalar(3.0)))) {
            /*
             * % prolongement au bord
             * y1 = (c1~=2) .* (2 .* y(K,:,:) - y(K-1,:,:)) + ...
             */
            mlfAssign(
              &y1,
              mlfPlus(
                mlfTimes(
                  mlfNe(c1, mlfScalar(2.0)),
                  mlfMinus(
                    mlfTimes(
                      mlfScalar(2.0),
                      mlfIndexRef(
                        y,
                        "(?,?,?)",
                        K,
                        mlfCreateColonIndex(),
                        mlfCreateColonIndex())),
                    mlfIndexRef(
                      y,
                      "(?,?,?)",
                      mlfMinus(K, mlfScalar(1.0)),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex()))),
                mlfTimes(
                  mlfEq(c1, mlfScalar(2.0)),
                  mlfMinus(
                    mlfPlus(
                      mlfMinus(
                        mlfTimes(
                          mlfScalar(4.0),
                          mlfIndexRef(
                            y,
                            "(?,?,?)",
                            K,
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex())),
                        mlfTimes(
                          mlfScalar(6.0),
                          mlfIndexRef(
                            y,
                            "(?,?,?)",
                            mlfMinus(K, mlfScalar(1.0)),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex()))),
                      mlfTimes(
                        mlfScalar(4.0),
                        mlfIndexRef(
                          y,
                          "(?,?,?)",
                          mlfMinus(K, mlfScalar(2.0)),
                          mlfCreateColonIndex(),
                          mlfCreateColonIndex()))),
                    mlfIndexRef(
                      y,
                      "(?,?,?)",
                      mlfMinus(K, mlfScalar(3.0)),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex())))));
            /*
             * (c1==2) .* (4 .* y(K,:,:) - 6 .* y(K-1,:,:) + 4 .* y(K-2,:,:) - y(K-3,:,:));
             * % prolongement au centre
             * y0 = (c0~=2)  .* (2 .* y(1,:,:) - y(2,:,:)) +  ...
             */
            mlfAssign(
              &y0,
              mlfPlus(
                mlfTimes(
                  mlfNe(c0, mlfScalar(2.0)),
                  mlfMinus(
                    mlfTimes(
                      mlfScalar(2.0),
                      mlfIndexRef(
                        y,
                        "(?,?,?)",
                        mlfScalar(1.0),
                        mlfCreateColonIndex(),
                        mlfCreateColonIndex())),
                    mlfIndexRef(
                      y,
                      "(?,?,?)",
                      mlfScalar(2.0),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex()))),
                mlfTimes(
                  mlfEq(c0, mlfScalar(2.0)),
                  mlfMinus(
                    mlfPlus(
                      mlfMinus(
                        mlfTimes(
                          mlfScalar(4.0),
                          mlfIndexRef(
                            y,
                            "(?,?,?)",
                            mlfScalar(1.0),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex())),
                        mlfTimes(
                          mlfScalar(6.0),
                          mlfIndexRef(
                            y,
                            "(?,?,?)",
                            mlfScalar(2.0),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex()))),
                      mlfTimes(
                        mlfScalar(4.0),
                        mlfIndexRef(
                          y,
                          "(?,?,?)",
                          mlfScalar(3.0),
                          mlfCreateColonIndex(),
                          mlfCreateColonIndex()))),
                    mlfIndexRef(
                      y,
                      "(?,?,?)",
                      mlfScalar(4.0),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex())))));
            /*
             * (c0==2) .* (4 .* y(1,:,:) - 6 .* y(2,:,:) + 4 .* y(3,:,:) - y(4,:,:));
             * % concatenation des matrices
             * dy=cat(1,y(k1,:,:),y1)-cat(1,y0,y(k0,:,:));
             */
            mlfAssign(
              &dy,
              mlfMinus(
                mlfCat(
                  mlfScalar(1.0),
                  mlfIndexRef(
                    y,
                    "(?,?,?)",
                    k1,
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex()),
                  y1,
                  NULL),
                mlfCat(
                  mlfScalar(1.0),
                  y0,
                  mlfIndexRef(
                    y,
                    "(?,?,?)",
                    k0,
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex()),
                  NULL)));
        /*
         * 
         * elseif n == 4
         */
        } else if (mlfTobool(mlfEq(n, mlfScalar(4.0)))) {
            /*
             * % prolongement au bord
             * y1 = (c1~=2) .* (2 .* y(K,:,:,:) - y(K-1,:,:,:)) + ...
             */
            mlfAssign(
              &y1,
              mlfPlus(
                mlfTimes(
                  mlfNe(c1, mlfScalar(2.0)),
                  mlfMinus(
                    mlfTimes(
                      mlfScalar(2.0),
                      mlfIndexRef(
                        y,
                        "(?,?,?,?)",
                        K,
                        mlfCreateColonIndex(),
                        mlfCreateColonIndex(),
                        mlfCreateColonIndex())),
                    mlfIndexRef(
                      y,
                      "(?,?,?,?)",
                      mlfMinus(K, mlfScalar(1.0)),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex()))),
                mlfTimes(
                  mlfEq(c1, mlfScalar(2.0)),
                  mlfMinus(
                    mlfPlus(
                      mlfMinus(
                        mlfTimes(
                          mlfScalar(4.0),
                          mlfIndexRef(
                            y,
                            "(?,?,?,?)",
                            K,
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex())),
                        mlfTimes(
                          mlfScalar(6.0),
                          mlfIndexRef(
                            y,
                            "(?,?,?,?)",
                            mlfMinus(K, mlfScalar(1.0)),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex()))),
                      mlfTimes(
                        mlfScalar(4.0),
                        mlfIndexRef(
                          y,
                          "(?,?,?,?)",
                          mlfMinus(K, mlfScalar(2.0)),
                          mlfCreateColonIndex(),
                          mlfCreateColonIndex(),
                          mlfCreateColonIndex()))),
                    mlfIndexRef(
                      y,
                      "(?,?,?,?)",
                      mlfMinus(K, mlfScalar(3.0)),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex())))));
            /*
             * (c1==2) .* (4 .* y(K,:,:,:) - 6 .* y(K-1,:,:,:) + 4 .* y(K-2,:,:,:) - y(K-3,:,:,:));
             * % prolongement au centre
             * y0 = (c0~=2)  .* (2 .* y(1,:,:,:) - y(2,:,:,:)) +  ...
             */
            mlfAssign(
              &y0,
              mlfPlus(
                mlfTimes(
                  mlfNe(c0, mlfScalar(2.0)),
                  mlfMinus(
                    mlfTimes(
                      mlfScalar(2.0),
                      mlfIndexRef(
                        y,
                        "(?,?,?,?)",
                        mlfScalar(1.0),
                        mlfCreateColonIndex(),
                        mlfCreateColonIndex(),
                        mlfCreateColonIndex())),
                    mlfIndexRef(
                      y,
                      "(?,?,?,?)",
                      mlfScalar(2.0),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex()))),
                mlfTimes(
                  mlfEq(c0, mlfScalar(2.0)),
                  mlfMinus(
                    mlfPlus(
                      mlfMinus(
                        mlfTimes(
                          mlfScalar(4.0),
                          mlfIndexRef(
                            y,
                            "(?,?,?,?)",
                            mlfScalar(1.0),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex())),
                        mlfTimes(
                          mlfScalar(6.0),
                          mlfIndexRef(
                            y,
                            "(?,?,?,?)",
                            mlfScalar(2.0),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex()))),
                      mlfTimes(
                        mlfScalar(4.0),
                        mlfIndexRef(
                          y,
                          "(?,?,?,?)",
                          mlfScalar(3.0),
                          mlfCreateColonIndex(),
                          mlfCreateColonIndex(),
                          mlfCreateColonIndex()))),
                    mlfIndexRef(
                      y,
                      "(?,?,?,?)",
                      mlfScalar(4.0),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex())))));
            /*
             * (c0==2) .* (4 .* y(1,:,:,:) - 6 .* y(2,:,:,:) + 4 .* y(3,:,:,:) - y(4,:,:,:));
             * % concatenation des matrices
             * dy=cat(1,y(k1,:,:,:),y1)-cat(1,y0,y(k0,:,:,:));
             */
            mlfAssign(
              &dy,
              mlfMinus(
                mlfCat(
                  mlfScalar(1.0),
                  mlfIndexRef(
                    y,
                    "(?,?,?,?)",
                    k1,
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex()),
                  y1,
                  NULL),
                mlfCat(
                  mlfScalar(1.0),
                  y0,
                  mlfIndexRef(
                    y,
                    "(?,?,?,?)",
                    k0,
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex()),
                  NULL)));
        /*
         * 
         * else 
         */
        } else {
            /*
             * error('Pas encore implante - a vous de l''ecrire')
             */
            mlfError(
              mxCreateString("Pas encore implante - a vous de l'ecrire"));
        /*
         * 
         * end	
         */
        }
    /*
     * else
     */
    } else {
        /*
         * % selon le nombre de dimension (1 a 4)
         * if n == 1
         */
        if (mlfTobool(mlfEq(n, mlfScalar(1.0)))) {
            /*
             * % prolongement au bord
             * y1 = (c1~=2) .* y(K-1) + ...
             */
            mlfAssign(
              &y1,
              mlfPlus(
                mlfTimes(
                  mlfNe(c1, mlfScalar(2.0)),
                  mlfIndexRef(y, "(?)", mlfMinus(K, mlfScalar(1.0)))),
                mlfTimes(
                  mlfEq(c1, mlfScalar(2.0)),
                  mlfMinus(
                    mlfPlus(
                      mlfMinus(
                        mlfTimes(mlfScalar(4.0), mlfIndexRef(y, "(?)", K)),
                        mlfTimes(
                          mlfScalar(6.0),
                          mlfIndexRef(y, "(?)", mlfMinus(K, mlfScalar(1.0))))),
                      mlfTimes(
                        mlfScalar(4.0),
                        mlfIndexRef(y, "(?)", mlfMinus(K, mlfScalar(2.0))))),
                    mlfIndexRef(y, "(?)", mlfMinus(K, mlfScalar(3.0)))))));
            /*
             * (c1==2) .* (4 .* y(K) - 6 .* y(K-1) + 4 .* y(K-2) - y(K-3));
             * % prolongement au centre
             * y0 = (c0~=2)  .* y(2) +  ...
             */
            mlfAssign(
              &y0,
              mlfPlus(
                mlfTimes(
                  mlfNe(c0, mlfScalar(2.0)),
                  mlfIndexRef(y, "(?)", mlfScalar(2.0))),
                mlfTimes(
                  mlfEq(c0, mlfScalar(2.0)),
                  mlfMinus(
                    mlfPlus(
                      mlfMinus(
                        mlfTimes(
                          mlfScalar(4.0),
                          mlfIndexRef(y, "(?)", mlfScalar(1.0))),
                        mlfTimes(
                          mlfScalar(6.0),
                          mlfIndexRef(y, "(?)", mlfScalar(2.0)))),
                      mlfTimes(
                        mlfScalar(4.0), mlfIndexRef(y, "(?)", mlfScalar(3.0)))),
                    mlfIndexRef(y, "(?)", mlfScalar(4.0))))));
            /*
             * (c0==2) .* (4 .* y(1) - 6 .* y(2) + 4 .* y(3) - y(4));
             * % concatenation des matrices
             * dy=cat(1,y(k1),y1) - 2 .* y + cat(1,y0,y(k0));
             */
            mlfAssign(
              &dy,
              mlfPlus(
                mlfMinus(
                  mlfCat(mlfScalar(1.0), mlfIndexRef(y, "(?)", k1), y1, NULL),
                  mlfTimes(mlfScalar(2.0), y)),
                mlfCat(mlfScalar(1.0), y0, mlfIndexRef(y, "(?)", k0), NULL)));
        /*
         * 
         * elseif n == 2
         */
        } else if (mlfTobool(mlfEq(n, mlfScalar(2.0)))) {
            /*
             * % prolongement au bord
             * y1 = (c1~=2) .* y(K-1,:) + ...
             */
            mlfAssign(
              &y1,
              mlfPlus(
                mlfTimes(
                  mlfNe(c1, mlfScalar(2.0)),
                  mlfIndexRef(
                    y,
                    "(?,?)",
                    mlfMinus(K, mlfScalar(1.0)),
                    mlfCreateColonIndex())),
                mlfTimes(
                  mlfEq(c1, mlfScalar(2.0)),
                  mlfMinus(
                    mlfPlus(
                      mlfMinus(
                        mlfTimes(
                          mlfScalar(4.0),
                          mlfIndexRef(y, "(?,?)", K, mlfCreateColonIndex())),
                        mlfTimes(
                          mlfScalar(6.0),
                          mlfIndexRef(
                            y,
                            "(?,?)",
                            mlfMinus(K, mlfScalar(1.0)),
                            mlfCreateColonIndex()))),
                      mlfTimes(
                        mlfScalar(4.0),
                        mlfIndexRef(
                          y,
                          "(?,?)",
                          mlfMinus(K, mlfScalar(2.0)),
                          mlfCreateColonIndex()))),
                    mlfIndexRef(
                      y,
                      "(?,?)",
                      mlfMinus(K, mlfScalar(3.0)),
                      mlfCreateColonIndex())))));
            /*
             * (c1==2) .* (4 .* y(K,:) - 6 .* y(K-1,:) + 4 .* y(K-2,:) - y(K-3,:));
             * % prolongement au centre
             * y0 = (c0~=2)  .* y(2,:) +  ...
             */
            mlfAssign(
              &y0,
              mlfPlus(
                mlfTimes(
                  mlfNe(c0, mlfScalar(2.0)),
                  mlfIndexRef(
                    y, "(?,?)", mlfScalar(2.0), mlfCreateColonIndex())),
                mlfTimes(
                  mlfEq(c0, mlfScalar(2.0)),
                  mlfMinus(
                    mlfPlus(
                      mlfMinus(
                        mlfTimes(
                          mlfScalar(4.0),
                          mlfIndexRef(
                            y, "(?,?)", mlfScalar(1.0), mlfCreateColonIndex())),
                        mlfTimes(
                          mlfScalar(6.0),
                          mlfIndexRef(
                            y,
                            "(?,?)",
                            mlfScalar(2.0),
                            mlfCreateColonIndex()))),
                      mlfTimes(
                        mlfScalar(4.0),
                        mlfIndexRef(
                          y, "(?,?)", mlfScalar(3.0), mlfCreateColonIndex()))),
                    mlfIndexRef(
                      y, "(?,?)", mlfScalar(4.0), mlfCreateColonIndex())))));
            /*
             * (c0==2) .* (4 .* y(1,:) - 6 .* y(2,:) + 4 .* y(3,:) - y(4,:));
             * % concatenation des matrices
             * dy=cat(1,y(k1,:),y1) - 2 .* y + cat(1,y0,y(k0,:));
             */
            mlfAssign(
              &dy,
              mlfPlus(
                mlfMinus(
                  mlfCat(
                    mlfScalar(1.0),
                    mlfIndexRef(y, "(?,?)", k1, mlfCreateColonIndex()),
                    y1,
                    NULL),
                  mlfTimes(mlfScalar(2.0), y)),
                mlfCat(
                  mlfScalar(1.0),
                  y0,
                  mlfIndexRef(y, "(?,?)", k0, mlfCreateColonIndex()),
                  NULL)));
        /*
         * 
         * elseif n == 3
         */
        } else if (mlfTobool(mlfEq(n, mlfScalar(3.0)))) {
            /*
             * % prolongement au bord
             * y1 = (c1~=2) .* y(K-1,:,:) + ...
             */
            mlfAssign(
              &y1,
              mlfPlus(
                mlfTimes(
                  mlfNe(c1, mlfScalar(2.0)),
                  mlfIndexRef(
                    y,
                    "(?,?,?)",
                    mlfMinus(K, mlfScalar(1.0)),
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex())),
                mlfTimes(
                  mlfEq(c1, mlfScalar(2.0)),
                  mlfMinus(
                    mlfPlus(
                      mlfMinus(
                        mlfTimes(
                          mlfScalar(4.0),
                          mlfIndexRef(
                            y,
                            "(?,?,?)",
                            K,
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex())),
                        mlfTimes(
                          mlfScalar(6.0),
                          mlfIndexRef(
                            y,
                            "(?,?,?)",
                            mlfMinus(K, mlfScalar(1.0)),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex()))),
                      mlfTimes(
                        mlfScalar(4.0),
                        mlfIndexRef(
                          y,
                          "(?,?,?)",
                          mlfMinus(K, mlfScalar(2.0)),
                          mlfCreateColonIndex(),
                          mlfCreateColonIndex()))),
                    mlfIndexRef(
                      y,
                      "(?,?,?)",
                      mlfMinus(K, mlfScalar(3.0)),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex())))));
            /*
             * (c1==2) .* (4 .* y(K,:,:) - 6 .* y(K-1,:,:) + 4 .* y(K-2,:,:) - y(K-3,:,:));
             * % prolongement au centre
             * y0 = (c0~=2)  .* y(2,:,:) +  ...
             */
            mlfAssign(
              &y0,
              mlfPlus(
                mlfTimes(
                  mlfNe(c0, mlfScalar(2.0)),
                  mlfIndexRef(
                    y,
                    "(?,?,?)",
                    mlfScalar(2.0),
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex())),
                mlfTimes(
                  mlfEq(c0, mlfScalar(2.0)),
                  mlfMinus(
                    mlfPlus(
                      mlfMinus(
                        mlfTimes(
                          mlfScalar(4.0),
                          mlfIndexRef(
                            y,
                            "(?,?,?)",
                            mlfScalar(1.0),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex())),
                        mlfTimes(
                          mlfScalar(6.0),
                          mlfIndexRef(
                            y,
                            "(?,?,?)",
                            mlfScalar(2.0),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex()))),
                      mlfTimes(
                        mlfScalar(4.0),
                        mlfIndexRef(
                          y,
                          "(?,?,?)",
                          mlfScalar(3.0),
                          mlfCreateColonIndex(),
                          mlfCreateColonIndex()))),
                    mlfIndexRef(
                      y,
                      "(?,?,?)",
                      mlfScalar(4.0),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex())))));
            /*
             * (c0==2) .* (4 .* y(1,:,:) - 6 .* y(2,:,:) + 4 .* y(3,:,:) - y(4,:,:));
             * % concatenation des matrices
             * dy=cat(1,y(k1,:,:),y1) - 2 .* y + cat(1,y0,y(k0,:,:));
             */
            mlfAssign(
              &dy,
              mlfPlus(
                mlfMinus(
                  mlfCat(
                    mlfScalar(1.0),
                    mlfIndexRef(
                      y,
                      "(?,?,?)",
                      k1,
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex()),
                    y1,
                    NULL),
                  mlfTimes(mlfScalar(2.0), y)),
                mlfCat(
                  mlfScalar(1.0),
                  y0,
                  mlfIndexRef(
                    y,
                    "(?,?,?)",
                    k0,
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex()),
                  NULL)));
        /*
         * 
         * elseif n == 4
         */
        } else if (mlfTobool(mlfEq(n, mlfScalar(4.0)))) {
            /*
             * % prolongement au bord
             * y1 = (c1~=2) .* y(K-1,:,:,:) + ...
             */
            mlfAssign(
              &y1,
              mlfPlus(
                mlfTimes(
                  mlfNe(c1, mlfScalar(2.0)),
                  mlfIndexRef(
                    y,
                    "(?,?,?,?)",
                    mlfMinus(K, mlfScalar(1.0)),
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex())),
                mlfTimes(
                  mlfEq(c1, mlfScalar(2.0)),
                  mlfMinus(
                    mlfPlus(
                      mlfMinus(
                        mlfTimes(
                          mlfScalar(4.0),
                          mlfIndexRef(
                            y,
                            "(?,?,?,?)",
                            K,
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex())),
                        mlfTimes(
                          mlfScalar(6.0),
                          mlfIndexRef(
                            y,
                            "(?,?,?,?)",
                            mlfMinus(K, mlfScalar(1.0)),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex()))),
                      mlfTimes(
                        mlfScalar(4.0),
                        mlfIndexRef(
                          y,
                          "(?,?,?,?)",
                          mlfMinus(K, mlfScalar(2.0)),
                          mlfCreateColonIndex(),
                          mlfCreateColonIndex(),
                          mlfCreateColonIndex()))),
                    mlfIndexRef(
                      y,
                      "(?,?,?,?)",
                      mlfMinus(K, mlfScalar(3.0)),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex())))));
            /*
             * (c1==2) .* (4 .* y(K,:,:,:) - 6 .* y(K-1,:,:,:) + 4 .* y(K-2,:,:,:) - y(K-3,:,:,:));
             * % prolongement au centre
             * y0 = (c0~=2)  .* y(2,:,:,:) +  ...
             */
            mlfAssign(
              &y0,
              mlfPlus(
                mlfTimes(
                  mlfNe(c0, mlfScalar(2.0)),
                  mlfIndexRef(
                    y,
                    "(?,?,?,?)",
                    mlfScalar(2.0),
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex())),
                mlfTimes(
                  mlfEq(c0, mlfScalar(2.0)),
                  mlfMinus(
                    mlfPlus(
                      mlfMinus(
                        mlfTimes(
                          mlfScalar(4.0),
                          mlfIndexRef(
                            y,
                            "(?,?,?,?)",
                            mlfScalar(1.0),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex())),
                        mlfTimes(
                          mlfScalar(6.0),
                          mlfIndexRef(
                            y,
                            "(?,?,?,?)",
                            mlfScalar(2.0),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex(),
                            mlfCreateColonIndex()))),
                      mlfTimes(
                        mlfScalar(4.0),
                        mlfIndexRef(
                          y,
                          "(?,?,?,?)",
                          mlfScalar(3.0),
                          mlfCreateColonIndex(),
                          mlfCreateColonIndex(),
                          mlfCreateColonIndex()))),
                    mlfIndexRef(
                      y,
                      "(?,?,?,?)",
                      mlfScalar(4.0),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex())))));
            /*
             * (c0==2) .* (4 .* y(1,:,:,:) - 6 .* y(2,:,:,:) + 4 .* y(3,:,:,:) - y(4,:,:,:));
             * % concatenation des matrices
             * dy=cat(1,y(k1,:,:,:),y1) - 2 .* y + cat(1,y0,y(k0,:,:,:));
             */
            mlfAssign(
              &dy,
              mlfPlus(
                mlfMinus(
                  mlfCat(
                    mlfScalar(1.0),
                    mlfIndexRef(
                      y,
                      "(?,?,?,?)",
                      k1,
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex(),
                      mlfCreateColonIndex()),
                    y1,
                    NULL),
                  mlfTimes(mlfScalar(2.0), y)),
                mlfCat(
                  mlfScalar(1.0),
                  y0,
                  mlfIndexRef(
                    y,
                    "(?,?,?,?)",
                    k0,
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex(),
                    mlfCreateColonIndex()),
                  NULL)));
        /*
         * 
         * else 
         */
        } else {
            /*
             * error('Pas encore implante - a vous de l''ecrire')
             */
            mlfError(
              mxCreateString("Pas encore implante - a vous de l'ecrire"));
        /*
         * 
         * end	
         */
        }
    /*
     * end
     */
    }
    /*
     * 
     * 
     * % calcul de la derivee
     * if mode ==1
     */
    if (mlfTobool(mlfEq(mode, mlfScalar(1.0)))) {
        /*
         * dx = dx * comp;
         */
        mlfAssign(&dx, mlfMtimes(dx, comp));
    /*
     * end
     */
    }
    /*
     * 
     * if o ==1
     */
    if (mlfTobool(mlfEq(o, mlfScalar(1.0)))) {
        /*
         * dydx = dy ./ (2 .* dx) ;
         */
        mlfAssign(&dydx, mlfRdivide(dy, mlfTimes(mlfScalar(2.0), dx)));
    /*
     * else
     */
    } else {
        /*
         * dydx = dy ./ (dx .^ 2);
         */
        mlfAssign(&dydx, mlfRdivide(dy, mlfPower(dx, mlfScalar(2.0))));
    /*
     * end
     */
    }
    /*
     * 
     * % mise a zeros selon les demandes des bornes
     * % selon le nombre de dimension (1 a 4)
     * if n == 1
     */
    if (mlfTobool(mlfEq(n, mlfScalar(1.0)))) {
        /*
         * % 0 au bord
         * dydx(K) =(c1~=0) .* dydx(K);
         */
        mlfIndexAssign(
          &dydx,
          "(?)",
          K,
          mlfTimes(mlfNe(c1, mlfScalar(0.0)), mlfIndexRef(dydx, "(?)", K)));
        /*
         * % 0 au centre
         * dydx(1) = (c0~=0) .* dydx(1);
         */
        mlfIndexAssign(
          &dydx,
          "(?)",
          mlfScalar(1.0),
          mlfTimes(
            mlfNe(c0, mlfScalar(0.0)),
            mlfIndexRef(dydx, "(?)", mlfScalar(1.0))));
    /*
     * 
     * elseif n == 2
     */
    } else if (mlfTobool(mlfEq(n, mlfScalar(2.0)))) {
        /*
         * % 0 au bord
         * dydx(K,:) =(c1~=0) .* dydx(K,:);
         */
        mlfIndexAssign(
          &dydx,
          "(?,?)",
          K,
          mlfCreateColonIndex(),
          mlfTimes(
            mlfNe(c1, mlfScalar(0.0)),
            mlfIndexRef(dydx, "(?,?)", K, mlfCreateColonIndex())));
        /*
         * % 0 au centre
         * dydx(1,:) = (c0~=0) .* dydx(1,:);
         */
        mlfIndexAssign(
          &dydx,
          "(?,?)",
          mlfScalar(1.0),
          mlfCreateColonIndex(),
          mlfTimes(
            mlfNe(c0, mlfScalar(0.0)),
            mlfIndexRef(dydx, "(?,?)", mlfScalar(1.0), mlfCreateColonIndex())));
    /*
     * 
     * elseif n == 3
     */
    } else if (mlfTobool(mlfEq(n, mlfScalar(3.0)))) {
        /*
         * % 0 au bord
         * dydx(K,:,:) =(c1~=0) .* dydx(K,:,:);
         */
        mlfIndexAssign(
          &dydx,
          "(?,?,?)",
          K,
          mlfCreateColonIndex(),
          mlfCreateColonIndex(),
          mlfTimes(
            mlfNe(c1, mlfScalar(0.0)),
            mlfIndexRef(
              dydx,
              "(?,?,?)",
              K,
              mlfCreateColonIndex(),
              mlfCreateColonIndex())));
        /*
         * % 0 au centre
         * dydx(1,:,:) = (c0~=0) .* dydx(1,:,:);
         */
        mlfIndexAssign(
          &dydx,
          "(?,?,?)",
          mlfScalar(1.0),
          mlfCreateColonIndex(),
          mlfCreateColonIndex(),
          mlfTimes(
            mlfNe(c0, mlfScalar(0.0)),
            mlfIndexRef(
              dydx,
              "(?,?,?)",
              mlfScalar(1.0),
              mlfCreateColonIndex(),
              mlfCreateColonIndex())));
    /*
     * 
     * elseif n == 4
     */
    } else if (mlfTobool(mlfEq(n, mlfScalar(4.0)))) {
        /*
         * % 0 au bord
         * dydx(K,:,:,:) =(c1~=0) .* dydx(K,:,:,:);
         */
        mlfIndexAssign(
          &dydx,
          "(?,?,?,?)",
          K,
          mlfCreateColonIndex(),
          mlfCreateColonIndex(),
          mlfCreateColonIndex(),
          mlfTimes(
            mlfNe(c1, mlfScalar(0.0)),
            mlfIndexRef(
              dydx,
              "(?,?,?,?)",
              K,
              mlfCreateColonIndex(),
              mlfCreateColonIndex(),
              mlfCreateColonIndex())));
        /*
         * % 0 au centre
         * dydx(1,:,:,:) = (c0~=0) .* dydx(1,:,:,:);
         */
        mlfIndexAssign(
          &dydx,
          "(?,?,?,?)",
          mlfScalar(1.0),
          mlfCreateColonIndex(),
          mlfCreateColonIndex(),
          mlfCreateColonIndex(),
          mlfTimes(
            mlfNe(c0, mlfScalar(0.0)),
            mlfIndexRef(
              dydx,
              "(?,?,?,?)",
              mlfScalar(1.0),
              mlfCreateColonIndex(),
              mlfCreateColonIndex(),
              mlfCreateColonIndex())));
    /*
     * 
     * else 
     */
    } else {
        /*
         * error('Pas encore implante - a vous de l''ecrire')
         */
        mlfError(mxCreateString("Pas encore implante - a vous de l'ecrire"));
    /*
     * 
     * end	
     */
    }
    /*
     * 
     * 
     * % permutation inverse si on ne derive pas selon la dimension 1 
     * if n == 2 & d ==2
     */
    {
        mxArray * a_ = mclInitialize(mlfEq(n, mlfScalar(2.0)));
        if (mlfTobool(a_) && mlfTobool(mlfAnd(a_, mlfEq(d, mlfScalar(2.0))))) {
            mxDestroyArray(a_);
            /*
             * dydx  = dydx.';
             */
            mlfAssign(&dydx, mlfTranspose(dydx));
        /*
         * 
         * elseif d>1
         */
        } else {
            mxDestroyArray(a_);
            if (mlfTobool(mlfGt(d, mlfScalar(1.0)))) {
                /*
                 * dydx=shiftdim(dydx,n-(d-1));
                 */
                mlfAssign(
                  &dydx,
                  mlfShiftdim(
                    NULL, dydx, mlfMinus(n, mlfMinus(d, mlfScalar(1.0)))));
            }
        }
    /*
     * end
     */
    }
    mclValidateOutputs("rpdederive", 1, nargout_, &dydx);
    mxDestroyArray(K);
    mxDestroyArray(c0);
    mxDestroyArray(c1);
    mxDestroyArray(comp);
    mxDestroyArray(dx);
    mxDestroyArray(dy);
    mxDestroyArray(k0);
    mxDestroyArray(k1);
    mxDestroyArray(mode);
    mxDestroyArray(n);
    mxDestroyArray(nargin_);
    mxDestroyArray(o);
    mxDestroyArray(ts);
    mxDestroyArray(x);
    mxDestroyArray(x0);
    mxDestroyArray(x1);
    mxDestroyArray(y);
    mxDestroyArray(y0);
    mxDestroyArray(y1);
    return dydx;
}

/*
 * The function "mlfRpdederive" contains the normal interface for the
 * "rpdederive" M-function from file
 * "/home/sccp/gsem/cgc/matlab5/zineb/v3.0/solver/rpdederive.m" (lines 1-349).
 * This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
mxArray * mlfRpdederive(mxArray * x,
                        mxArray * y,
                        mxArray * c0,
                        mxArray * c1,
                        mxArray * d,
                        mxArray * o) {
    int nargout = 1;
    mxArray * dydx = mclGetUninitializedArray();
    mlfEnterNewContext(0, 6, x, y, c0, c1, d, o);
    dydx = Mrpdederive(nargout, x, y, c0, c1, d, o);
    mlfRestorePreviousContext(0, 6, x, y, c0, c1, d, o);
    return mlfReturnValue(dydx);
}

/*
 * The function "mlxRpdederive" contains the feval interface for the
 * "rpdederive" M-function from file
 * "/home/sccp/gsem/cgc/matlab5/zineb/v3.0/solver/rpdederive.m" (lines 1-349).
 * The feval function calls the implementation version of rpdederive through
 * this function. This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
void mlxRpdederive(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[6];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: rpdederive Line: 93 Column"
            ": 0 The function \"rpdederive\" was called with "
            "more than the declared number of outputs (1)"));
    }
    if (nrhs > 6) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: rpdederive Line: 93 Column"
            ": 0 The function \"rpdederive\" was called with "
            "more than the declared number of inputs (6)"));
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 6 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 6; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(
      0, 6, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4], mprhs[5]);
    mplhs[0]
      = Mrpdederive(
          nlhs, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4], mprhs[5]);
    mlfRestorePreviousContext(
      0, 6, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4], mprhs[5]);
    plhs[0] = mplhs[0];
}
