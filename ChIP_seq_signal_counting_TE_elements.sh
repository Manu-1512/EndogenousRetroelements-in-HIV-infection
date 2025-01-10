mkdir extendend_windows

## Let's extend the co-ordinates 5 KB upstream and 5 KB downstream (Total 10 KB)
# Will get a new file in "extendend_windows" directory with prefix "sorted_" name


for i in *bed; \
do awk '{\
         if ($6=="+") print $1,$2-5000,$2+5000,$4,$5,$6; \
         else print $1,$3-5000,$3+5000,$4,$5,$6}' OFS="\t" "$i" | \
         sort-bed - > extendend_windows/sorted_"$i"; done
   
   #Let's go to directory "extendend_windows"      
   cd extendend_windows   
   
   # We remove coloumn 5 (It has no info in my bed file and and putting coloumn 6 (strand information) in coloumn 4), resulting file would be prefixed "formatted_"

for i in *; do awk '{print $1,$2,$3,$6,$4}' OFS="\t" "$i" > formated_"$i"; done


# formatted_*bed file should look like this (it's for demo )

#chr12	4822413	4832413	+	LTR12
#chr11	23525631	23535631	-	LTR12
#chr13	49303395	49313395	-	LTR12
#chr3	4796707	4806707	+	LTR12
#chr2	132514018	132524018	+	LTR12
#chr11	95566068	95576068	-	LTR12
#chr4	11411104	11421104	+	LTR12
#chr2	91276799	91286799	-	LTR12




### Now we divide them into 100 bins of 100 base pairs (total was 10 KB)

for i in formated_*; \
           do \
              awk 'BEGIN { binNum=100; } \
                           { chr=$1; start=$2; stop=$3; \
                            binIdx=1; binSize=(stop-start)/binNum; \
                  for (i=start;i<stop;i+=binSize) \
                          { print chr"\t"i"\t"(i+binSize)"\t"binIdx"\t"$4"\t"$5; binIdx++; \
                           } \
                            }' "$i" \
           > bins_"$i"; \
             done
           
###### here you get the bins as prefix "bins_formatted_sorted_" names
###### Once you generate these files then you can use them for life time

#make a directory to throw binned files
mkdir bins

#Sort them just to be sure that file would be Okay
for i in bins_*; do sort -k1,1 -k2,2n "$i" > bins/"$i"; done


#### Now Let's count the ChIP-seq signal in every bins

#### We take raw signal file in "bedGraph" which is obtained by running "macs2" in -B mode
#### Since we are going to run "bedops" tool to count the signal in every bin so we have to put the value in 5th coloumn
#### Because bedops assume that counts are in 5th coloumn whereas the bedGraph file we get from macs2 has signal in 4th coulmn (first 3 coloumns are bed coordinate of reads)

for i in *.bedGraph; \
    do \
       awk '{print $1,$2,$3,"ID",$4}' OFS="\t" "$i" | \
          sort-bed - \
          > Signal_"$i"; \
    done 
    
#### Now we got the signal file which is ready for bedops
### Now we run bedops
### for example, we are counting signals of all ChIP-seq data which is in bedGraph format over the LTR7-up which is this file "/workdir/Manu/Data/for_Tommy/LTR_files/Tommy_list/bins_LTR7_hESC_inactive_2KB_ext.bed"
$FILE = /workdir/Manu/Data/for_Tommy/LTR_files/Tommy_list/bins_LTR7_hESC_active_2KB_ext.bed
for i in Signal*.bedGraph; \
         do bedmap \
            --sum \
               --prec 0 \
                  $FILE "$i" \
                  > /workdir/Manu/Data/for_Tommy/LTR_files/Tommy_list/bins_LTR7_active_2KB_ext_"$i"_count; \
        done
        
        
# We got the one file of each ChIP-seq data. File name ending with ChIP_seq_bedGraph_count.

### Now we go to R do the fun
R
### First import the bins which we counted the signal
bins_LTR7 <- read.delim("bins_LTR7_hESC_active_2KB_ext.bed",stringsAsFactors=F, header=F)
#### Second import the count which we got the output from bedops
FOXA1_active_count=read.delim("bins_LTR7_hESC_active_2KB_ext_Signal_GSM1505633_FoxA1_021113_h64.bw.bedGraph_count",stringsAsFactors=F, header=F) 

head(bins_LTR7)

## putting the counts in 7th coloumn of bineed file
bins_LTR7[,7]=FOXA1_active_count[,1]
Just checking if things are Okay
head(bins_LTR7)
head(bins_LTR7, 200)

## assigning the coloumn names to our file
colnames(bins_LTR7)=c("chr","start","stop","bin_num","strand","LTR5s","tag_count")
### make an empty dataframe , here we took all 262 LTR7-up but I realized later that it was 256 (6 were duplicated) but it doesn't matter here
FOXA1_active_clust=data.frame(1:262)
### make an empty list
jo=list();
### fill the embty list (jo) with 100 couloumns and 262 rows
for( i in 1:262){jo[[i]]=rep(i,100)}

# filled jo is transferred to ho because we don't want to change jo (so that if anything goes wrong then we can use it in future)
ho=list();
#### putting the list of LTR7, 262*100 times in "ho"
ho=rep("LTR7", 26200)
### Now we repeat jo (1 to 100) and combine with ho (LTR7) using _ between them and put them all to df named "go"
go=paste(unlist(ho),unlist(jo),sep="_")
## Now we have an unique name of every entry and we can just put this in 8th couloumn so that we have an unique ID for every bin
bins_LTR7[,8]=go; head(bins_LTR7)
# just checking if everything worked fine
length(unique(bins_LTR7[,8]))
# Assigning the unique IDs as rownames
rownames(FOXA1_active_clust)=(unique(bins_LTR7[,8]))
### the empty dataframe we created before , should have numeric strings so just putting zero into it
FOXA1_active_clust[,1:100]=0 
### 1 to 100 as coloumn names 
colnames(FOXA1_active_clust)=1:100
# just checking the dimensions if everything is Okay
dim(FOXA1_active_clust)
head(FOXA1_active_clust)
head(bins_LTR7)

####

### Now you have a dataframe with 100 coloumns which is every bin and 262 rows which is every locus
#### Now we put values into it, previously all jo, ho, go was for this action
### Here we also reverse the coloumns if it's on -ve strand while preserving the bin numbers
for(i in 1:262){
 cat("\r", i , " ")
if(head(bins_LTR7[(which(rownames(FOXA1_active_clust[i,]) == bins_LTR7[,8])),5],1) == "+"){FOXA1_active_clust[i,]=(bins_LTR7[(which(rownames(FOXA1_active_clust[i,]) == bins_LTR7[,8])),7])}
 else{FOXA1_active_clust[i,]=rev(bins_LTR7[(which(rownames(FOXA1_active_clust[i,]) == bins_LTR7[,8])),7])}
 }
 
 ### As the bins without signal would have "NA" values so we are putting 0 into that
FOXA1_active_clust[(is.na(FOXA1_active_clust))]=0

# preserving the obtained data as it is and transferring everything to new df with name FOXA1_active_LTR7_H9
 FOXA1_active_LTR7_H9=FOXA1_active_clust

# Now to make sure that everything is numeric or integer
for (i in 1:100){FOXA1_active_LTR7_H9[,i]=as.integer(FOXA1_active_clust[,i])}

# just checking
str(FOXA1_active_LTR7_H9)
summary(FOXA1_active_LTR7_H9[,1])
head(is.na(FOXA1_active_LTR7_H9))

# previous command can sometimes conver 0 to "NA" values so we are putting back 0 into that
FOXA1_active_LTR7_H9[(is.na(FOXA1_active_LTR7_H9))]=0
## Just checking
summary(FOXA1_active_LTR7_H9[,1])

## Our data is ready Now so we preserve the original one and put a df into "cho" name for further analysis

cho=(FOXA1_active_LTR7_H9)
cho[cho == "NAN"]=0
for(i in 1:100){cho[,i]=as.numeric(cho[,i])}
### We found that sum of entire genomewide signal was 6.95 million so we just calculate signal per million 
# We make final df LA_FOXP1 for plotting and statistics
cho <- cho/6.954196
LA_FOXP1=(cho)

plot(apply(LA_FOXP1, 2, mean), type="n", xlab="Distance from TSS (kb)", ylab="Average tag count", ylim=c(0,20))
lines(apply(LA_FOXP1, 2, mean), col="hotpink4",lwd=5)
dev.off()

### If plot is very spiky then we can smooth it and just save it
par(mfrow = c(1, 1), xaxt="n")
png("FOXA1_hESC_LTR7_Ave_profile.png",  width = 14.7, height = 12.8, units = "cm", res = 600, pointsize = 12)
plot(smooth.spline(apply(LA_FOXP1/6, 2, mean), spar=0.40), type="n", xlab="Distance from TSS (kb)", ylab="Average tag count", ylim=c(0,4))
lines(smooth.spline(apply(LA_FOXP1/6, 2, mean), spar=0.40), col="hotpink4",lwd=5)
axis(1,at=c(0,20,40,60,80,100),labels=c(-1,-0.6,-0.2,0.2,0.6,1))
legend(x = "topright", c("FOXP1 7-up","FOXP1 7-down"), bty = "n", cex=1, col=c("hotpink4", "midnightblue"),lwd=c(rep(5,2)),lty=c(rep("solid",3)))
dev.off()


################################ next segment is IRRELEVANT for the current analysis but nice to have it

#### Clustering the loci to see if different loci is having different bins where ChIP-seq signal is enriched
## here we use K-means clustering and asking to give top 5 clusters
# putting the clustering in df "cl"
cl=kmeans(cho,5)

## plotting them

plot(smooth.spline(apply(cluster1,2,mean)),type="n",ylim=c(0,3.5),xlab="Distance from LTR5 Left Boundary (kb)", ylab="Average tag count", xaxt="n")
axis(1,at=c(0,20,40,50,60,80,100),labels=c(-5,-3,-1,0,1,3,5))
lines(smooth.spline(apply(cluster1,2,mean)),col="darkred",lwd=2.5)
lines(smooth.spline(apply(cluster2,2,mean)),col="darkgrey",lwd=2.5)
lines(smooth.spline(apply(cluster4,2,mean)),col="orange",lwd=2.5)
lines(smooth.spline(apply(cluster5,2,mean)),col="darkblue",lwd=2.5)
lines(smooth.spline(apply(cluster3,2,mean)),col="green",lwd=2.5)
legend(x = "topleft", c(paste("Cluster 1","","(n=",nrow(cluster1),")",sep=""), paste("Cluster 2","","(n=",nrow(cluster2),")",sep=""), paste("Cluster 3","","(n=",nrow(cluster3),")",sep=""), paste("Cluster 4","","(n=",nrow(cluster4),")",sep=""), paste("Cluster 5","","(n=",nrow(cluster5),")",sep="")), bty = "n", cex=1,col=c("darkred","darkgrey","green","orange","darkblue"),lwd=c(rep(5,3)),lty=c(rep("solid",2)))

dev.off()

### Just renaming the clusters according their position
cluster1=cho[(match((names(which(cl$cluster==4))),rownames(cho))),]
cluster2=cho[(match((names(which(cl$cluster==5))),rownames(cho))),]
cluster3=cho[(match((names(which(cl$cluster==3))),rownames(cho))),]
cluster4=cho[(match((names(which(cl$cluster==2))),rownames(cho))),]
cluster5=cho[(match((names(which(cl$cluster==1))),rownames(cho))),]

########################################### keep plotting till it looks good