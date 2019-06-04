#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
#If we don't pass in genomewide as a non-zero value, output per-scaffold:
   if (length(genomewide)==0) {
      genomewide=0;
   };
#Account for weights in column 4 if indicated, otherwise consider it omit:
   if (length(weighted)==0) {
      weighted=0;
   };
}
{
#Only include sites without the omission column if that column is present:
   if (weighted == 0 && (NF == 3 || $4 == 0)) {
      stat[$1]+=$3;
      count[$1]+=1;
   } else if (weighted == 1) {
#If weighted is selected, but we don't have a weight column, revert:
      if (weighted != 0 && NF > 3) {
         stat[$1]+=$3*$4;
         count[$1]+=$4;
      } else {
         stat[$1]+=$3;
         count[$1]+=1;
      };
   };
}
END{
   if (genomewide != 0) {
      for (scaf in stat) {
         stat_gw+=stat[scaf];
         count_gw+=count[scaf];
      };
      if (count_gw > 0) {
         print "Genome-wide", stat_gw/count_gw;
      } else {
         print "Genome-wide", "NA";
      };
   } else {
      for (scaf in stat) {
         if (count[scaf] > 0) {
            print scaf,stat[scaf]/count[scaf];
         } else {
            print scaf,"NA";
         };
      };
   };
}
