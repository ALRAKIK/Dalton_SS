      INTEGER MAXX1LBL
      PARAMETER ( MAXX1LBL = 120 )

      LOGICAL LX1OPN
      INTEGER NX1LBL,ISYOFX1

      LOGICAL LORXX1(MAXX1LBL)
      INTEGER ISYX1(MAXX1LBL)
      CHARACTER*8 LBLX1(MAXX1LBL)

#if defined (SYS_CRAY)
      REAL FRQX1(MAXX1LBL), AVEX1(MAXX1LBL)
#else
      DOUBLE PRECISION FRQX1(MAXX1LBL), AVEX1(MAXX1LBL)
#endif


      COMMON/IX1RSP/ ISYX1, NX1LBL, ISYOFX1(8), LX1OPN, LORXX1
      COMMON/CX1RSP/ LBLX1
      COMMON/RX1RSP/ FRQX1,AVEX1

