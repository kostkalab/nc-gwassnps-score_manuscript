
## data:
#E999-sco.bw     intergrated score
#E083-sco.bw     Fetal Heart Score
#E105-sco.bw     Right Ventricle Score
#E065-sco.bw     Aorta Score
#E095-sco.bw     Left Ventricle Score 
#E104-sco.bw     Right Atrium Score

## It contains 127 tissue specific score (001-129 without 060 and 064)and 1 integreted score.
# H4 is heart specific score
mkdir data
cd data

## get the data
# Integreted score
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E999-sco.bw
# H1
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E001-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E002-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E003-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E004-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E005-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E006-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E007-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E008-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E009-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E010-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E011-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E012-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E013-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E014-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E015-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E016-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E018-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E019-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E020-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E021-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E022-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H1/E024-sco.bw
# H2
wget http://compgen.cshl.edu/fitCons2/hg19/H2/E053-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H2/E054-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H2/E067-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H2/E068-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H2/E069-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H2/E070-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H2/E071-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H2/E072-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H2/E073-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H2/E074-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H2/E081-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H2/E082-sco.bw
# H3
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E029-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E030-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E031-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E032-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E033-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E034-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E035-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E036-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E037-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E038-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E039-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E040-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E041-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E042-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E043-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E044-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E045-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E046-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E047-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E048-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E050-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E051-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H3/E062-sco.bw
# H4 Heart specific
wget http://compgen.cshl.edu/fitCons2/hg19/H4/E083-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H4/E105-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H4/E065-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H4/E095-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H4/E104-sco.bw
# H5
wget http://compgen.cshl.edu/fitCons2/hg19/H5/E075-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H5/E077-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H5/E079-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H5/E084-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H5/E085-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H5/E092-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H5/E094-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H5/E101-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H5/E102-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H5/E106-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H5/E109-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H5/E110-sco.bw
# H6
wget http://compgen.cshl.edu/fitCons2/hg19/H6/E076-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H6/E078-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H6/E089-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H6/E090-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H6/E100-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H6/E103-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H6/E107-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H6/E108-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H6/E111-sco.bw
# H7
wget http://compgen.cshl.edu/fitCons2/hg19/H7/E023-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H7/E025-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H7/E026-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H7/E027-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H7/E028-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H7/E049-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H7/E055-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H7/E056-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H7/E057-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H7/E058-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H7/E059-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H7/E061-sco.bw
# H8 
wget http://compgen.cshl.edu/fitCons2/hg19/H8/E017-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H8/E052-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H8/E063-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H8/E066-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H8/E080-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H8/E086-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H8/E087-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H8/E088-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H8/E091-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H8/E093-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H8/E096-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H8/E097-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H8/E098-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H8/E099-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H8/E112-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H8/E113-sco.bw
# H9
wget http://compgen.cshl.edu/fitCons2/hg19/H9/E114-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H9/E115-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H9/E116-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H9/E117-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H9/E118-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H9/E119-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H9/E120-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H9/E121-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H9/E122-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H9/E123-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H9/E124-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H9/E125-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H9/E126-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H9/E127-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H9/E128-sco.bw
wget http://compgen.cshl.edu/fitCons2/hg19/H9/E129-sco.bw



