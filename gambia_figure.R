gambia_figure=function()
{
  #####Creates the Gambia paper Figure (see github readme for more info). Requires ggplot2, seqinr and ShortRead. Please install these packages. Apologies that the files haven't been renamed as the samples.  The order of files is seen in new_order which matches up to the labels.  LAIV53 and LAIV54 are the vaccine 2017 and 2018 respectively.
  new_order=c(15, 16, 17, 3, 4, 18, 21, 5, 6, 22, 1, 2, 7, 8, 9, 10, 11, 12, 13, 14, 19, 20)
  all_a=list.files(pattern="_ALL a")
  filenames=all_a[new_order]
  library(ggplot2)
  library(seqinr)
  library(ShortRead)
  
  fig_labels=c('B027J Day 2','B099K Day 2','B156H Day 2','B344K Day 2','D107J Day 2','D414H Day 2','A399H Day 7','A440E Day 7','A457G Day 7','A484G Day 7','B023C Day 7','B027J Day 7','B099K Day 7','B183G Day 7','B204D Day 7','B221F Day 7','B344K Day 7','D117A Day 7','D414H Day 7','D438H Day 7','Vaccine 2017','Vaccine 2018')
  print(paste("There are",length(filenames),"files",sep=" "))
  freq_results=vector("list", length(filenames))
  figure_data=array(dim=c(850*22,3))
  #figure_data[,1]=rep(29:878,22)
  #figure_data[,3]=rep(fig_labels,each=850)
  for(file_number in 1:length(filenames))
  { 
    print(paste("Reading in",strsplit(filenames[file_number],split="[.]")[[1]][1],"file number",file_number,sep=" "))
    dataf=readDNAStringSet(filenames[file_number])
    
    #Remove Sequences with N
    dataf=dataf[c(1,which(vcountPattern("N",dataf)==0)),]
    datam=as.matrix(dataf)
    subamp_locations=c(1,469,386,897) #Where the subamplicons are alternating start and finish
    s=subamp_locations*0-1
    for(a in 1:length(subamp_locations))
    {
      s[a]=which(datam[1,]!="-")[subamp_locations[a]] #A better way for finding the start and finish of the amplicon.
    }
    
    print("Found these subamplicon splits")
    print(s)
    if(-1%in%s)
    {
      print("Not found one of the starts/ends of the amplicons.  Will crash soon.  Sorry")
    }
    subamplist=array(dim=length(dataf))
    
    subamplist[which(substr(dataf[],s[2]+200,s[2]+210)!="-----------")]=2
    subamplist[which(substr(dataf[],s[1]+100,s[1]+110)!="-----------")]=1
    print(table(subamplist))
    
    frequency_array=array(dim=c(5,850))
    tempsub=which(subamplist[]==1)
    for(a in 29:447)   #Takes into account where the primer is hence the -28 and uses subamp 1 for the overlap.
    {
      base_location=which(datam[1,]!="-")[a]
      frequency_array[1,a-28]=length(which(datam[tempsub,base_location]%in%c("A","a")))
      frequency_array[2,a-28]=length(which(datam[tempsub,base_location]%in%c("C","c")))
      frequency_array[3,a-28]=length(which(datam[tempsub,base_location]%in%c("G","g")))
      frequency_array[4,a-28]=length(which(datam[tempsub,base_location]%in%c("T","t")))
      frequency_array[5,a-28]=length(which(datam[tempsub,base_location]=="-"))
    }
    frequency_array=frequency_array/table(subamplist)[[1]]
    
    tempsub=which(subamplist[]==2)
    for(a in 448:878) 
    {
      base_location=which(datam[1,]!="-")[a]
      frequency_array[1,a-28]=length(which(datam[tempsub,base_location]%in%c("A","a")))
      frequency_array[2,a-28]=length(which(datam[tempsub,base_location]%in%c("C","c")))
      frequency_array[3,a-28]=length(which(datam[tempsub,base_location]%in%c("G","g")))
      frequency_array[4,a-28]=length(which(datam[tempsub,base_location]%in%c("T","t")))
      frequency_array[5,a-28]=length(which(datam[tempsub,base_location]=="-"))
    }
    frequency_array[,420:850]=frequency_array[,420:850]/table(subamplist)[[2]]
    freq_results[[file_number]]=frequency_array
    figure_data[(file_number*850-849):(file_number*850),2]=t(colSums(frequency_array)-colMaxs(frequency_array))
    #plot(1:850,colSums(frequency_array)-colMaxs(frequency_array))
  }
   figure_df=data.frame(base=rep(9:858,22),samplen=rep(fig_labels,each=850),mutfreq=figure_data[,2]*100) #Makes the data frame has the figure labels provided and the numbers for 
  
   ggplot(figure_df)+geom_point(aes(base,mutfreq))+facet_wrap(figure_df$samplen)+ labs( x = "Nucleotide Position", y = "Mutation Percentage")
  
}