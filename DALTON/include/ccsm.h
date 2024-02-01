      INTEGER MXSMOP, MXSMSEL

      PARAMETER ( MXSMOP = 125 , MXSMSEL = 50 )

      INTEGER IPRSM, NSMSEL, NSMOPER

      INTEGER IASMOP, IBSMOP, ICSMOP, IDSMOP

      INTEGER NSMSELX, ISMSELX, ISMSL, ISMSEL

#if defined (SYS_CRAY)
      REAL EXSMFR, BSMFR
#else
      DOUBLE PRECISION EXSMFR, BSMFR
#endif
      COMMON /INFSMCC/  EXSMFR(MXSMSEL),  BSMFR(MXSMSEL),
     *                   ISMSEL(MXSMSEL,2),
     *                   IASMOP(MXSMOP), IBSMOP(MXSMOP), ICSMOP(MXSMOP),
     *                   IDSMOP(MXSMOP), NSMSEL, IPRSM, NSMOPER,
     *                   NSMSELX(8), ISMSELX(8)
