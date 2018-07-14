#ifndef TEST1_H
#define TEST1_H

std::string inputCPPCommentTest("DOAUTOPHASE = 1;\n"
"OPTION, ECHO=FALSE;\n"
"OPTION, PSDUMPFREQ=200000;\n"
"OPTION, STATDUMPFREQ=50;\n"
"IF (DOAUTOPHASE == 1)\n"
" {OPTION, AUTOPHASE=4;}\n"
"ELSE\n"
" {OPTION, AUTOPHASE=0;}\n"
"OPTION, CZERO=TRUE;\n"
"\n"
"TITLE, STRING=\"SRF-Injector77pC, BK, 2014-10-16\";\n"
"\n"
"//////////////////////////////\n"
"// Begin Content\n"
"//////////////////////////////\n"
"\n"
"//////////////////////////////\n"
"// GAF_1\n"
"//////////////////////////////\n"
"MARK1L1GAF_1: MARKER, ELEMEDGE=0, L=0;\n"
"CATHODEVGAF_1: RFCAVITY, ELEMEDGE=0.0, L=0;\n"
"\n"
"IF (DOAUTOPHASE == 1) {\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  LAG=-5.9/180*PI, FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"} ELSE {\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  LAG=(44.68+90-5.9)/180*PI, FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"}\n");

std::string outputCPPCommentTest("DOAUTOPHASE = 1;\n"
"OPTION, ECHO=FALSE;\n"
"OPTION, PSDUMPFREQ=200000;\n"
"OPTION, STATDUMPFREQ=50;\n"
"IF (DOAUTOPHASE == 1)\n"
" {OPTION, AUTOPHASE=4;}\n"
"ELSE\n"
" {OPTION, AUTOPHASE=0;}\n"
"OPTION, CZERO=TRUE;\n"
"\n"
"TITLE, STRING=\"SRF-Injector77pC, BK, 2014-10-16\";\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"MARK1L1GAF_1: MARKER, ELEMEDGE=0, L=0;\n"
"CATHODEVGAF_1: RFCAVITY, ELEMEDGE=0.0, L=0;\n"
"\n"
"IF (DOAUTOPHASE == 1) {\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  LAG=-5.9/180*PI, FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"} ELSE {\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  LAG=(44.68+90-5.9)/180*PI, FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"}\n");

std::string inputCCommentTest("DOAUTOPHASE = 1;\n"
"OPTION, ECHO=FALSE;\n"
"OPTION, PSDUMPFREQ=200000;\n"
"OPTION, STATDUMPFREQ=50;\n"
"IF (DOAUTOPHASE == 1)\n"
" {OPTION, AUTOPHASE=4;}\n"
"ELSE\n"
" {OPTION, AUTOPHASE=0;}\n"
"OPTION, CZERO=TRUE;\n"
"\n"
"TITLE, STRING=\"SRF-Injector77pC, BK, 2014-10-16\";\n"
"\n"
"//////////////////////////////\n"
"// Begin Content\n"
"//////////////////////////////\n"
"\n"
"//////////////////////////////\n"
"// GAF_1\n"
"//////////////////////////////\n"
"MARK1L1GAF_1: MARKER, ELEMEDGE=0, L=0;\n"
"CATHODEVGAF_1: RFCAVITY, ELEMEDGE=0.0, L=0;\n"
"\n"
"IF (DOAUTOPHASE == 1) {\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  LAG=-5.9/180*PI, FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"/*} ELSE {*/\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  /*LAG=(44.68+90-5.9)/180*PI, */FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"}\n");

std::string outputCCommentTest("DOAUTOPHASE = 1;\n"
"OPTION, ECHO=FALSE;\n"
"OPTION, PSDUMPFREQ=200000;\n"
"OPTION, STATDUMPFREQ=50;\n"
"IF (DOAUTOPHASE == 1)\n"
" {OPTION, AUTOPHASE=4;}\n"
"ELSE\n"
" {OPTION, AUTOPHASE=0;}\n"
"OPTION, CZERO=TRUE;\n"
"\n"
"TITLE, STRING=\"SRF-Injector77pC, BK, 2014-10-16\";\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"MARK1L1GAF_1: MARKER, ELEMEDGE=0, L=0;\n"
"CATHODEVGAF_1: RFCAVITY, ELEMEDGE=0.0, L=0;\n"
"\n"
"IF (DOAUTOPHASE == 1) {\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  LAG=-5.9/180*PI, FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"}\n");

std::string inputTGaussTest("DOAUTOPHASE = 1;\n"
"OPTION, ECHO=FALSE;\n"
"OPTION, PSDUMPFREQ=200000;\n"
"OPTION, STATDUMPFREQ=50;\n"
"IF (DOAUTOPHASE == 1)\n"
" {OPTION, AUTOPHASE=4;}\n"
"ELSE\n"
" {OPTION, AUTOPHASE=0;}\n"
"OPTION, CZERO=TRUE;\n"
"\n"
"TITLE, STRING=\"SRF-Injector77pC, BK, 2014-10-16\";\n"
"\n"
"//////////////////////////////\n"
"// Begin Content\n"
"//////////////////////////////\n"
"\n"
"//////////////////////////////\n"
"// GAF_1\n"
"//////////////////////////////\n"
"MARK1L1GAF_1: MARKER, ELEMEDGE=0, L=0;\n"
"CATHODEVGAF_1: RFCAVITY, ELEMEDGE=0.0 + TGAUSS(-3.141e-3), L=0;\n"
"\n"
"IF (DOAUTOPHASE == 1) {\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  LAG=-5.9/180*PI, FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"/*} ELSE {*/\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  /*LAG=(44.68+90-5.9)/180*PI, */FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"}\n");

std::string outputTGaussTest1("DOAUTOPHASE = 1;\n"
"OPTION, ECHO=FALSE;\n"
"OPTION, PSDUMPFREQ=200000;\n"
"OPTION, STATDUMPFREQ=50;\n"
"IF (DOAUTOPHASE == 1)\n"
" {OPTION, AUTOPHASE=4;}\n"
"ELSE\n"
" {OPTION, AUTOPHASE=0;}\n"
"OPTION, CZERO=TRUE;\n"
"\n"
"TITLE, STRING=\"SRF-Injector77pC, BK, 2014-10-16\";\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"MARK1L1GAF_1: MARKER, ELEMEDGE=0, L=0;\n"
"CATHODEVGAF_1: RFCAVITY, ELEMEDGE=0.0 + ");

std::string outputTGaussTest2(", L=0;\n"
"\n"
"IF (DOAUTOPHASE == 1) {\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  LAG=-5.9/180*PI, FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"}\n");

std::string inputGaussTest("DOAUTOPHASE = 1;\n"
"OPTION, ECHO=FALSE;\n"
"OPTION, PSDUMPFREQ=200000;\n"
"OPTION, STATDUMPFREQ=50;\n"
"IF (DOAUTOPHASE == 1)\n"
" {OPTION, AUTOPHASE=4;}\n"
"ELSE\n"
" {OPTION, AUTOPHASE=0;}\n"
"OPTION, CZERO=TRUE;\n"
"\n"
"TITLE, STRING=\"SRF-Injector77pC, BK, 2014-10-16\";\n"
"\n"
"//////////////////////////////\n"
"// Begin Content\n"
"//////////////////////////////\n"
"\n"
"//////////////////////////////\n"
"// GAF_1\n"
"//////////////////////////////\n"
"MARK1L1GAF_1: MARKER, ELEMEDGE=0 + GAUSS(), L=0;\n"
"CATHODEVGAF_1: RFCAVITY, ELEMEDGE=0.0 + TGAUSS(-3.141e-3), L=0;\n"
"\n"
"IF (DOAUTOPHASE == 1) {\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  LAG=-5.9/180*PI, FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"/*} ELSE {*/\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  /*LAG=(44.68+90-5.9)/180*PI, */FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"}\n");

std::string outputGaussTest1("DOAUTOPHASE = 1;\n"
"OPTION, ECHO=FALSE;\n"
"OPTION, PSDUMPFREQ=200000;\n"
"OPTION, STATDUMPFREQ=50;\n"
"IF (DOAUTOPHASE == 1)\n"
" {OPTION, AUTOPHASE=4;}\n"
"ELSE\n"
" {OPTION, AUTOPHASE=0;}\n"
"OPTION, CZERO=TRUE;\n"
"\n"
"TITLE, STRING=\"SRF-Injector77pC, BK, 2014-10-16\";\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"MARK1L1GAF_1: MARKER, ELEMEDGE=0 + ");

std::string outputGaussTest2(", L=0;\n"
"CATHODEVGAF_1: RFCAVITY, ELEMEDGE=0.0 + ");

std::string inputCallTest("DOAUTOPHASE = 1;\n"
"OPTION, ECHO=FALSE;\n"
"OPTION, PSDUMPFREQ=200000;\n"
"OPTION, STATDUMPFREQ=50;\n"
"IF (DOAUTOPHASE == 1)\n"
" {OPTION, AUTOPHASE=4;}\n"
"ELSE\n"
" {OPTION, AUTOPHASE=0;}\n"
"OPTION, CZERO=TRUE;\n"
"\n"
"TITLE, STRING=\"SRF-Injector77pC, BK, 2014-10-16\";\n"
"\n"
"//////////////////////////////\n"
"// Begin Content\n"
"//////////////////////////////\n"
"\n"
"//////////////////////////////\n"
"// GAF_1\n"
"//////////////////////////////\n"
"MARK1L1GAF_1: MARKER, ELEMEDGE=0, L=0;\n"
"CATHODEVGAF_1: RFCAVITY, ELEMEDGE=0.0, L=0;\n"
"\n"
" call , 	file= 	\"path/to/some/OPAL_input-file.ver3.141.txt\" ;     \n"
"IF (DOAUTOPHASE == 1) {\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  LAG=-5.9/180*PI, FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"} ELSE {\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  LAG=(44.68+90-5.9)/180*PI, FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"}\n");

std::string outputCallTest1("DOAUTOPHASE = 1;\n"
"OPTION, ECHO=FALSE;\n"
"OPTION, PSDUMPFREQ=200000;\n"
"OPTION, STATDUMPFREQ=50;\n"
"IF (DOAUTOPHASE == 1)\n"
" {OPTION, AUTOPHASE=4;}\n"
"ELSE\n"
" {OPTION, AUTOPHASE=0;}\n"
"OPTION, CZERO=TRUE;\n"
"\n"
"TITLE, STRING=\"SRF-Injector77pC, BK, 2014-10-16\";\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"MARK1L1GAF_1: MARKER, ELEMEDGE=0, L=0;\n"
"CATHODEVGAF_1: RFCAVITY, ELEMEDGE=0.0, L=0;\n"
"\n"
" ");

std::string inputSubFileTest1("DOAUTOPHASE = 1;\n"
"OPTION, ECHO=FALSE;\n"
"OPTION, PSDUMPFREQ=200000;\n"
"OPTION, STATDUMPFREQ=50;\n"
"IF (DOAUTOPHASE == 1)\n"
" {OPTION, AUTOPHASE=4;}\n"
"ELSE\n"
" {OPTION, AUTOPHASE=0;}\n"
"OPTION, CZERO=TRUE;\n"
"\n"
"TITLE, STRING=\"SRF-Injector77pC, BK, 2014-10-16\";\n"
"\n"
"//////////////////////////////\n"
"// Begin Content\n"
"//////////////////////////////\n"
"\n"
"//////////////////////////////\n"
"// GAF_1\n"
"//////////////////////////////\n"
"MARK1L1GAF_1: MARKER, ELEMEDGE=0, L=0;\n"
"CATHODEVGAF_1: RFCAVITY, ELEMEDGE=0.0, L=0;\n"
"\n"
" call , 	file= 	\"SubFileTest.in\" ;     \n"
"IF (DOAUTOPHASE == 1) {\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  LAG=-5.9/180*PI, FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"} ELSE {\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384 + 1.234*TGauss(+2.1e-3), TYPE=\"STANDING\",\n"
"  LAG=(44.68+90-5.9)/180*PI, FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"}\n"
"TRACK, LINE=BerlinPro, BEAM=beam1, MAXSTEPS=180000, DT=1.0e-13, ZSTOP=15.0;\n"
"RUN, METHOD = \"STATISTICAL-ERRORS(PARALLEL-T, 4, 1000)\", BEAM = beam1, FIELDSOLVER = Fs1, DISTRIBUTION = Dist1;\n"
"ENDTRACK;\n");

std::string inputSubFileTest2("//////////////////////////////\n"
"// BAF_1\n"
"//////////////////////////////\n"
"MARK2L1BAF_1: MARKER, ELEMEDGE=2.10471, L=0;\n"
"IF (DOAUTOPHASE == 1) {\n"
"  CAVH1BAF_1: RFCAVITY, ELEMEDGE=2.892001 + 3.414e-3*GAUSS(), TYPE=\"STANDING\", LAG=-90./180*PI,\n"
"  FMAPFN=\"CBFp.dat\", VOLT=4, FREQ=1300;\n"
"  CAVH2BAF_1: RFCAVITY, ELEMEDGE=3.681949 + 1e-3*TGAUSS(.31e-3), TYPE=\"STANDING\", LAG=0/180*PI,\n"
"  FMAPFN=\"CBFp.dat\", VOLT=19, FREQ=1300;\n"
"  CAVH3BAF_1: RFCAVITY, ELEMEDGE=4.471991, TYPE=\"STANDING\", LAG=0/180*PI,\n"
"  FMAPFN=\"CBFp.dat\", VOLT=19, FREQ=1300;\n"
"\n"
"// fixed ASTRA autophase minus 180° seeks correct results\n"
"} ELSE {\n"
"  CAVH1BAF_1: RFCAVITY, ELEMEDGE=2.892001, TYPE=\"STANDING\", LAG=(234.63-180)/180*PI,\n"
"  FMAPFN=\"CBFp.dat\", VOLT=4, FREQ=1300;\n"
"  CAVH2BAF_1: RFCAVITY, ELEMEDGE=3.681949, TYPE=\"STANDING\", LAG=(60.62-90)/180*PI,\n"
"  FMAPFN=\"CBFp.dat\", VOLT=19, FREQ=1300;\n"
"  CAVH3BAF_1: RFCAVITY, ELEMEDGE=4.471991, TYPE=\"STANDING\", LAG=(260.13-90)/180*PI,\n"
"  FMAPFN=\"CBFp.dat\", VOLT=19, FREQ=1300;\n"
"}\n"
"\n"
"BAF_1: LINE = (MARK2L1BAF_1, CAVH1BAF_1, CAVH2BAF_1, CAVH3BAF_1);");

std::string outputSubFileTest1("     \n"
"IF (DOAUTOPHASE == 1) {\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  LAG=-5.9/180*PI, FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"} ELSE {\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384 + 1.234*");

std::string outputASTTest("DOAUTOPHASE = 1;\n"
"OPTION, ECHO=FALSE;\n"
"OPTION, PSDUMPFREQ=200000;\n"
"OPTION, STATDUMPFREQ=50;\n"
"IF (DOAUTOPHASE == 1)\n"
" {OPTION, AUTOPHASE=4;}\n"
"ELSE\n"
" {OPTION, AUTOPHASE=0;}\n"
"OPTION, CZERO=TRUE;\n"
"\n"
"TITLE, STRING=\"SRF-Injector77pC, BK, 2014-10-16\";\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"MARK1L1GAF_1: MARKER, ELEMEDGE=0, L=0;\n"
"CATHODEVGAF_1: RFCAVITY, ELEMEDGE=0.0, L=0;\n"
"\n"
" \n"
"\n"
"\n"
"MARK2L1BAF_1: MARKER, ELEMEDGE=2.10471, L=0;\n"
"IF (DOAUTOPHASE == 1) {\n"
"  CAVH1BAF_1: RFCAVITY, ELEMEDGE=2.892001 + 3.414e-3*0.0, TYPE=\"STANDING\", LAG=-90./180*PI,\n"
"  FMAPFN=\"CBFp.dat\", VOLT=4, FREQ=1300;\n"
"  CAVH2BAF_1: RFCAVITY, ELEMEDGE=3.681949 + 1e-3*0.0, TYPE=\"STANDING\", LAG=0/180*PI,\n"
"  FMAPFN=\"CBFp.dat\", VOLT=19, FREQ=1300;\n"
"  CAVH3BAF_1: RFCAVITY, ELEMEDGE=4.471991, TYPE=\"STANDING\", LAG=0/180*PI,\n"
"  FMAPFN=\"CBFp.dat\", VOLT=19, FREQ=1300;\n"
"\n"
"\n"
"} ELSE {\n"
"  CAVH1BAF_1: RFCAVITY, ELEMEDGE=2.892001, TYPE=\"STANDING\", LAG=(234.63-180)/180*PI,\n"
"  FMAPFN=\"CBFp.dat\", VOLT=4, FREQ=1300;\n"
"  CAVH2BAF_1: RFCAVITY, ELEMEDGE=3.681949, TYPE=\"STANDING\", LAG=(60.62-90)/180*PI,\n"
"  FMAPFN=\"CBFp.dat\", VOLT=19, FREQ=1300;\n"
"  CAVH3BAF_1: RFCAVITY, ELEMEDGE=4.471991, TYPE=\"STANDING\", LAG=(260.13-90)/180*PI,\n"
"  FMAPFN=\"CBFp.dat\", VOLT=19, FREQ=1300;\n"
"}\n"
"\n"
"BAF_1: LINE = (MARK2L1BAF_1, CAVH1BAF_1, CAVH2BAF_1, CAVH3BAF_1);\n"
"     \n"
"IF (DOAUTOPHASE == 1) {\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384, TYPE=\"STANDING\",\n"
"  LAG=-5.9/180*PI, FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"} ELSE {\n"
"  CAVH1GAF_1: RFCAVITY, ELEMEDGE=0.021384 + 1.234*0.0, TYPE=\"STANDING\",\n"
"  LAG=(44.68+90-5.9)/180*PI, FMAPFN=\"Ez_1_4cell_cathode1_5mm.txt\", VOLT=30, FREQ=1300;\n"
"}\n");

std::string originalTrack("TRACK, LINE=BerlinPro, BEAM=beam1, MAXSTEPS=180000, DT=1.0e-13, ZSTOP=15.0;\n"
"RUN, METHOD = \"STATISTICAL-ERRORS(PARALLEL-T, 4, 1000)\", BEAM = beam1, FIELDSOLVER = Fs1, DISTRIBUTION = Dist1;\n"
"ENDTRACK;\n");

std::string replacedTrack("TRACK, LINE=BerlinPro, BEAM=beam1, MAXSTEPS=180000, DT=1.0e-13, ZSTOP=15.0;\n"
"RUN, METHOD = \"PARALLEL-T\", BEAM = beam1, FIELDSOLVER = Fs1, DISTRIBUTION = Dist1;\n"
"ENDTRACK;\n");

#endif //TEST1_H