# create vectors of names to be assigned to each folder and to each batch file. A batch file and corresponding folder is needed 
#for adjustment parameters and running simulations within each chain. See Joseph et al., 2013. 
folder_name<-c("LREW1","LREW2","LREW3","LREW4","LREW5","LREW6","LREW7","LREW8")
bfname<-c("runLREW1.bat","runLREW2.bat","runLREW3.bat","runLREW4.bat",
"runLREW5.bat","runLREW6.bat","runLREW7.bat","runLREW8.bat")
#Define a template for the text of each batch file, and write one batch file for each folder to be craeted.
batch.text.template<-'@echo off
E:
cd "E:\\DATA-project\\docs\\my papers docs\\swat dream-parallel\\SWAT_DREAM\\!FOLDER_NAME!"
start /w SWAT_Edit.exe
start /w swat2009.exe
if %errorlevel%==0 exit 0
echo.'
for(i.batch in 1:length(folder_name)){
this.folder.name<-folder_name[i.batch]
this.bfname<-bfname[i.batch]
#Generate a uniquely named folder for each chain.
dir.create(this.folder.name)
file.copy(dir("TxtInOut",full.names=T),this.folder.name,recursive=TRUE)

#Replace this.folder.name into batch.text.template
batch.text<-gsub("!FOLDER_NAME!",this.folder.name,batch.text.template)
write(batch.text,this.bfname)
}
#
######### E N D #######
