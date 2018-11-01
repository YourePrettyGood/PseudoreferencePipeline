#!/bin/awk -f
#Apply this script to the allSites VCF produced by GenotypeGVCFs
BEGIN{
#Set output field separator to tab to produce TSVs
   FS="\t";
   OFS=FS;
#Debug:
#   for (i=0;i<=ARGC;i=i+1) {
#      print i, ARGV[i];
#   }
#End Debug
   indelwindow=ARGV[ARGC-1];
   delete ARGV[ARGC-1];
#   print "indelwindow", indelwindow;
}
!/^#/{
#For non-header lines of a GATK VCF
#Check if ref is of length > 1 (implying a deletion)
   split($5, alts, ",");
   minlen=0;
   maxlen=0;
   for (elem in alts) {
      indellength=length(alts[elem])-length($4);
      if (indellength < minlen) {
         minlen=indellength;
      }
      if (indellength > maxlen) {
         maxlen=indellength;
      }
   }
#Not a SNP:
   if (maxlen > 0 || minlen < 0) {
#Handle the typical case of indelwindow > 0:
      if (indelwindow > 0) {
#Insertion:
         if (maxlen > -minlen) {
            start=$2-indelwindow+1;
            end=$2+indelwindow;
#Deletion (and uncaught edge case of maxlen==-minlen):
         } else {
            start=$2-indelwindow+1;
            end=$2-minlen+indelwindow;
         }
#Print the masking interval:
         if (start < 1) {
            BEDstart=0;
         } else {
            BEDstart=start-1;
         }
         print $1, BEDstart, end;
#Handle the special case of indelwindow == 0:
      } else {
#Deletion (and uncaught edge case of maxlen==-minlen):
         if (maxlen <= -minlen) {
            start=$2+1;
            end=$2-minlen;
#Print the masking interval:
            if (start < 1) {
               BEDstart=0;
            } else {
               BEDstart=start-1;
            }
            print $1, BEDstart, end;
         }
      }
   }
}
