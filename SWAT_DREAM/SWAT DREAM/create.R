# create vectors of names to be assigned to each folder and to each batch file. A batch file and corresponding folder is needed 
#for adjustment parameters and running simulations within each chain.
#library(coda)
#library(snow)
#library(dream)
folder_name<-c("LREW1")
bfname<-c("runLREW1.bat")
#Define a template for the text of each batch file, and write one batch file for each folder to be craeted.
batch.text.template<-'@echo off
E:
cd "E:\\Dr Pourreza\\SWAT DREAM\\!FOLDER_NAME!"
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