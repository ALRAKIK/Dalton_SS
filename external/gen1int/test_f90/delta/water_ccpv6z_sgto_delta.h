  ! number of test cases
  integer, parameter :: NUM_TESTS = 3
  ! ID of block of sub-shells on bra center
  integer, parameter :: BRA_BLOCK(NUM_TESTS) = (/8, 8, 8/)
  ! ID of block of sub-shells on ket center
  integer, parameter :: KET_BLOCK(NUM_TESTS) = (/6, 6, 6/)
  ! orders of electronic derivatives
  integer, parameter :: ORDER_ELEC(NUM_TESTS) = (/0, 0, 0/)
  ! indices of Dirac delta function origin
  integer, parameter :: IDX_DELORG(NUM_TESTS) = (/1, 2, 3/)
  ! coordinates of Dirac delta function origin
  real(REALK), parameter :: DEL_ORIGIN(3*NUM_TESTS) = &
    (/0.0_REALK, -0.224905893_REALK, 0.0_REALK,       &
      1.45235_REALK, 0.899623_REALK, 0.0_REALK,       &
     -1.45235_REALK, 0.899623_REALK, 0.0_REALK/)
  ! indices of dipole origin
  integer, parameter :: IDX_DIPORG(NUM_TESTS) = (/0, 0, 0/)
  ! coordinates of dipole origin
  real(REALK), parameter :: DIP_ORIGIN(3*NUM_TESTS) = &
    (/0.0_REALK, 0.0_REALK, 0.0_REALK,                &
      0.0_REALK, 0.0_REALK, 0.0_REALK,                &
      0.0_REALK, 0.0_REALK, 0.0_REALK/)
  ! scale constants
  real(REALK), parameter :: SCAL_CONST(NUM_TESTS) =       &
    (/8.3872954891254174_REALK, 8.3872954891254174_REALK, &
      8.3872954891254174_REALK/)
  ! orders of Cartesian multipole moments
  integer, parameter :: ORDER_MOM(NUM_TESTS) = (/0, 0, 0/)
  ! orders of partial geometric derivatives on bra center
  integer, parameter :: ORDER_GEO_BRA(NUM_TESTS) = (/0, 0, 0/)
  ! orders of partial geometric derivatives on ket center
  integer, parameter :: ORDER_GEO_KET(NUM_TESTS) = (/0, 0, 0/)
  ! orders of geometric derivatives on Dirac delta function origin
  integer, parameter :: ORDER_GEO_POT(NUM_TESTS) = (/0, 0, 0/)
  ! orders of geometric derivatives on dipole origin
  integer, parameter :: ORDER_GEO_MOM(NUM_TESTS) = (/0, 0, 0/)
  ! numbers of differentiated centers of total geometric derivatives
  integer, parameter :: NUM_CENTS(NUM_TESTS) = (/0, 0, 0/)
  ! indices of differentiated centers
  integer, parameter :: IDX_CENT(3*NUM_TESTS) = &
    (/0, 0, 0,    0, 0, 0,    0, 0, 0/)
  ! orders of geometric derivatives of differentiated centers
  integer, parameter :: ORDER_CENT(3*NUM_TESTS) = &
    (/0, 0, 0,    0, 0, 0,    0, 0, 0/)
  ! referenced results from Dalton
  real(REALK) REF_CONTR_INTS(396)
  ! results of test 1, ket-major order, FC O1 01
  data REF_CONTR_INTS(1:132) /                                             &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK/
  ! results of test 1, ket-major order, FC H1 02
  data REF_CONTR_INTS(133:264) /                                           &
     -0.0013437403505582395068_REALK,     -0.0026503147049221388587_REALK, &
     -0.0013229242189971954398_REALK,     -0.0007010940107028099284_REALK, &
     -0.0003870622065568072234_REALK,     -0.0002127553215110628004_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
     -0.0060492831258405904693_REALK,     -0.0119312514623768142763_REALK, &
     -0.0059555725564250294388_REALK,     -0.0031562021389106150945_REALK, &
     -0.0017424860939853404652_REALK,     -0.0009577871021101596119_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0037318052037916263006_REALK,      0.0073603938464786543283_REALK, &
      0.0036739951156670642464_REALK,      0.0019470623743649314707_REALK, &
      0.0010749403752143694805_REALK,      0.0005908592501665269276_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0048196958935111762143_REALK,      0.0095060856768339821454_REALK, &
      0.0047450330884820247218_REALK,      0.0025146673038030226133_REALK, &
      0.0013883055061196987196_REALK,      0.0007631057212678974376_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0025989713373866783243_REALK,      0.0051260587287458641778_REALK, &
      0.0025587101892713156985_REALK,      0.0013560084266823545056_REALK, &
      0.0007486294358942700107_REALK,      0.0004114968954039577580_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
     -0.0087304995789890850466_REALK,     -0.0172195256367041193979_REALK, &
     -0.0085952537870811890341_REALK,     -0.0045551218006736251404_REALK, &
     -0.0025148061007343614631_REALK,     -0.0013823059225008529557_REALK, &
     -0.1831426265675174525072_REALK,     -0.3612197818487098421869_REALK, &
     -0.1803055300946158623354_REALK,     -0.0955543223343152392513_REALK, &
     -0.0527539335440684598044_REALK,     -0.0289970963375219166480_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
     -0.8244759488369548527942_REALK,     -1.6261480353324506697987_REALK, &
     -0.8117038386500607938245_REALK,     -0.4301687818319735456818_REALK, &
     -0.2374889463409160239582_REALK,     -0.1305398364349724116273_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.5086195458644937605897_REALK,      1.0031713798394767955813_REALK, &
      0.5007404259312160377249_REALK,      0.2653712952683455750247_REALK, &
      0.1465070269255860369029_REALK,      0.0805300778251216153869_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.6568916121537995422841_REALK,      1.2956145125120093197779_REALK, &
      0.6467155820790473352844_REALK,      0.3427319680998020534624_REALK, &
      0.1892165527092195820291_REALK,      0.1040060946920599721510_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.3542220317377152150051_REALK,      0.6986467728916923025651_REALK, &
      0.3487347123970339080046_REALK,      0.1848146815023621325480_REALK, &
      0.1020330759275604737280_REALK,      0.0560842146455980961051_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
     -1.1899074277842527358473_REALK,     -2.3469036648654304499928_REALK, &
     -1.1714743506262070837920_REALK,     -0.6208319714175101200482_REALK, &
     -0.3427508851729956784737_REALK,     -0.1883988504635355221151_REALK/
  ! results of test 1, ket-major order, FC H2 03
  data REF_CONTR_INTS(265:396) /                                           &
     -0.0000741109745488961248_REALK,     -0.0000000007018109798926_REALK, &
     -0.0000032906611595218069_REALK,     -0.0000535674459489697926_REALK, &
     -0.0001207657087420467134_REALK,     -0.0001259275368625990430_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
     -0.0003336345951001555392_REALK,     -0.0000000031594298083185_REALK, &
     -0.0000148139787696964969_REALK,     -0.0002411512363522173613_REALK, &
     -0.0005436660168535897538_REALK,     -0.0005669036607445199737_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0002058193164808537197_REALK,      0.0000000019490535249264_REALK, &
      0.0000091387494867716424_REALK,      0.0001487662951128506522_REALK, &
      0.0003353877853976239233_REALK,      0.0003497230973002190147_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
     -0.0002658194788517256706_REALK,     -0.0000000025172389118213_REALK, &
     -0.0000118028650928694218_REALK,     -0.0001921344396325339369_REALK, &
     -0.0004331595685574253334_REALK,     -0.0004516738907516653365_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
     -0.0001433404143578456626_REALK,     -0.0000000013573951398024_REALK, &
     -0.0000063645733575655612_REALK,     -0.0001036065163783655199_REALK, &
     -0.0002335768330759555110_REALK,     -0.0002435604905052041408_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0004815110537008214639_REALK,      0.0000000045597800660946_REALK, &
      0.0000213799606865014462_REALK,      0.0003480364075624538646_REALK, &
      0.0007846344488284106179_REALK,      0.0008181716855531794615_REALK, &
     -0.0101008193515395826301_REALK,     -0.0000000956520403350703_REALK, &
     -0.0004484946274391233547_REALK,     -0.0073008767992508482106_REALK, &
     -0.0164595407804177459865_REALK,     -0.0171630626769622429695_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
     -0.0454721152304820158685_REALK,     -0.0000004306086911141595_REALK, &
     -0.0020190440665643046361_REALK,     -0.0328672654707447206568_REALK, &
     -0.0740979626463564577943_REALK,     -0.0772650947010408567506_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0280517662530376590824_REALK,      0.0000002656426754822167_REALK, &
      0.0012455491002071887759_REALK,      0.0202758293448337054143_REALK, &
      0.0457110630866043726739_REALK,      0.0476648681304250310120_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
     -0.0362293783389692342078_REALK,     -0.0000003430824606981418_REALK, &
     -0.0016086498505698341565_REALK,     -0.0261866110619971663442_REALK, &
     -0.0590366889522240673172_REALK,     -0.0615600645391534626039_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
     -0.0195363188787671254731_REALK,     -0.0000001850036810237381_REALK, &
     -0.0008674478527060342859_REALK,     -0.0141208601283436419732_REALK, &
     -0.0318349260681808354789_REALK,     -0.0331956303467889901615_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0000000000000000000000_REALK,      0.0000000000000000000000_REALK, &
      0.0656266659399085683368_REALK,      0.0000006214668611594881_REALK, &
      0.0029139425294547839293_REALK,      0.0474349838461207429230_REALK, &
      0.1069403131297596010185_REALK,      0.1115112092995736298162_REALK/
