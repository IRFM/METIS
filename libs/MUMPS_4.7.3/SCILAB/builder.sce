// $Id: builder_source.sce,v 1.11 2007/05/02 11:17:44 jylexcel Exp $


//******************* VARIABLE PART TO COSTUMIZE ***************************//

// -- MUMPS:

MUMPS_DIR = home + "/MUMPS_4.7.3/";
MUMPS_INC_DIR = MUMPS_DIR+"include/"; //path until dmumps_c.h and zmumps_c.h
MUMPS_LIB_DIR = MUMPS_DIR+"lib/";     //path until libdmumps.a, libzmumps.a and libpord.a
DMUMPS_LIB = MUMPS_LIB_DIR+"libdmumps.a";
ZMUMPS_LIB = MUMPS_LIB_DIR+"libzmumps.a";
LIB_MPISEQ = MUMPS_LIB_DIR+"../libseq/libmpiseq.a";

// -- SCILAB: Path to scilab routines

SCI_DIR_INC = SCI + "/routines/";

// -- BLAS library, if not already included in Scilab:

BLAS_LIB = "";

// -- ORDERINGS (should correspond to the ones defined MUMPS's Makefile.inc):

PORD_LIB =  MUMPS_LIB_DIR+"libpord.a";
METIS_LIB = "/home/afevre/metis-4.0/libmetis.a";  
ORDERINGS_LIB = PORD_LIB+" "+METIS_LIB;

// -- PTHREAD lib required by MUMPS versions > 4.6
PTHREAD_LIB="-lpthread";

// -- COMPILER FOR THE INTERFACE
COMPILER_= "gcc -c -O";

// -- FORTRAN RUNTIME LIBRARIES

// -- g95 
FORT_LIB_DIR = "/home/afevre/compil/g95-install/lib/gcc-lib/i686-pc-linux-gnu/4.0.1/";
FORT_LIB = FORT_LIB_DIR+"libf95.a"+" "+FORT_LIB_DIR+"libgcc.a";

// -- gfortran compiler
//FORT_LIB_DIR= "/usr/lib/";
//FORT_LIB=FORT_LIB_DIR+"libgfortran.a";

// --  ifort compiler 
//FORT_LIB_DIR = "/opt/intel/fc/9.0/lib/";
//FORT_LIB = FORT_LIB_DIR+"libifcore.a"+" "+FORT_LIB_DIR+"libimf.a"+" "+FORT_LIB_DIR+"libguide.a"+" "+FORT_LIB_DIR+"libirc.a";

//**************************************************************************//
//******************* DON't EDIT BELOW (Normally) **************************//

//---- Build the Makefile 
fd=mopen("Makefile","w");
mfprintf(fd,"SCIDIRINC = %s\n",SCI_DIR_INC);
mfprintf(fd,"MUMPSINCDIR = %s\n",MUMPS_INC_DIR);
mfprintf(fd,"CC = %s\n", COMPILER_);
mfprintf(fd,"all : \n");
mfprintf(fd,"\t$(CC) intdmumpsc.c -I${MUMPSINCDIR} -I${SCIDIRINC}\n");
mfprintf(fd,"\t$(CC) intzmumpsc.c -I${MUMPSINCDIR} -I${SCIDIRINC}\n");
mfprintf(fd,"clean:\n");
mfprintf(fd,"\trm *.o loader_inc.sce\n");
mclose(fd);

//---- Compile
unix("make");

//---- Build the Loader_inc.sce
fd=mopen("loader_inc.sce","w");
mfprintf(fd,"objects = [ path+\""intzmumpsc.o \"" ; \n")
mfprintf(fd,"		path+\""intdmumpsc.o\"" ; \n")
mfprintf(fd,"		\""%s\"" ; \n",DMUMPS_LIB)
mfprintf(fd,"		\""%s\"" ; \n",ZMUMPS_LIB)
mfprintf(fd,"		\""%s\"" ; \n",ORDERINGS_LIB)
mfprintf(fd,"		\""%s\"" ; \n",LIB_MPISEQ)
mfprintf(fd,"		\""%s\"" ; \n",PORD_LIB)
mfprintf(fd,"		\""%s\"" ; \n",METIS_LIB)
mfprintf(fd,"		\""%s\"" ; \n",BLAS_LIB)
mfprintf(fd,"		\""%s\"" ; \n",FORT_LIB)
mfprintf(fd,"		\""%s\"" ]; \n",PTHREAD_LIB)
mclose(fd);
