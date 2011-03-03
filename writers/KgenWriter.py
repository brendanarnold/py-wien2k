# tetcnt.f contains code for writign to .kgen
# See the WRITE(15,...) commands
# 1234 FORMAT(2i10,e20.12,2I10) NKP,NTT,V,MWRIT,NREC
#      6768    122648  0.171228185529E-05       101      1215
# 1235 FORMAT(6i10)(ITTFL(J),J=1,%*MWRIT)
#         4         1         1         2      1038         4
# 1235 FORMAT(6i10)ITTFL
#         0
## NKP         NUMBER OF IRREDUCIBLE K-POINTS                **
## MWRIT       INFORMATION FOR MWRIT TETRAHEDRA ARE WRITTEN  **
##             AT ONE TIME. 
# NTT - not declared
