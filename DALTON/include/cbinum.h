#if defined (SYS_CRAY)
      REAL THRMID
#else
      DOUBLE PRECISION THRMID
#endif
      LOGICAL NRMCRD, PGMTST, HTEST, PREHES, REUHES, NUMVIB,
     &        ANALZ1, NPRPDR, HARMON, SPECTR, MIDAS, MINOUT
      COMMON /CBINUM/ NRMCRD, PGMTST, HTEST, PREHES, REUHES,
     &                NUMVIB, ANALZ1, NPRPDR, HARMON, SPECTR,
     &                MIDAS,  MINOUT
      COMMON /MIDINF / THRMID
