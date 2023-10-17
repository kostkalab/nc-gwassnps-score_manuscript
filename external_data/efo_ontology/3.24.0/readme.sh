# to get the efo.owl file in 3.24.0 version, using the following command. 
wget http://www.ebi.ac.uk/efo/releases/v3.24.0/efo.owl

# to convert efo.owl into efo_3_24_0.obo file, use robot.  
# to install robot, please see: http://robot.obolibrary.org 
robot convert --input efo.owl --check false --output efo_3_24_0.obo
