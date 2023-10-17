# First install gdown if you don't have it installed https://github.com/wkentaro/gdown 
pip install gdown

gdown https://drive.google.com/u/0/uc?id=0B85dCRJ7-QypVm1mZlhYZnh3UHc&export=download
# download file from https://sites.google.com/site/emorydivan/download 

tar -xvf AllTSSscore.tar.gz
cd allTSSscore/

mkdir DIVANtoolkit
cd DIVANtoolkit
gdown --fuzzy https://drive.google.com/file/d/0B6ffiUfBhxaPVnpLOFVZSXVWYWs/view?usp=sharing&resourcekey=0-f9p685AuwpAJIJdDJcelCg
# download DIVANtoolkit.zip from https://drive.google.com/file/d/0B6ffiUfBhxaPVnpLOFVZSXVWYWs/view?resourcekey=0-f9p685AuwpAJIJdDJcelCg

unzip DIVANtoolkit.zip
rm DIVANtoolkit.zip

cd ..

# download scoreTSS file and ensembl, cosmic and 1KG files
gdown --fuzzy https://drive.google.com/file/d/0B6ffiUfBhxaPOEhjcUZ3YnBIU2M/view?usp=sharing&resourcekey=0-LYRytCqlRpSmKHqWAgzLbA
gdown https://drive.google.com/u/0/uc?id=0B6ffiUfBhxaPSW5taHc3TkJtZ0k&export=download
gdown --fuzzy https://drive.google.com/file/d/0B6ffiUfBhxaPOUlkczY3N1J4TGs/view?usp=sharing&resourcekey=0-GiBO-TC2wjOhHl86GcDlTw
gdown https://drive.google.com/u/0/uc?id=0B6ffiUfBhxaPeU5ZaXFUQVFhUVk&export=download

for i in *.tar.gz; do
  tar -xvf $i
  rm $i
done


##### important, please read ######
cd DIVANtoolkit
# for line 36 in scoreDIVAN.console.genome.R
# change it into the following to make it work
# load(file.path(distribution,paste(basename(disease),".score.dist.rda",sep="")))


