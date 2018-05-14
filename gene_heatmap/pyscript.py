
# coding: utf-8

# In[1]:


import subprocess


# In[ ]:


subprocess.call("sort -k1,1 -k2,2n reads.bed > readssort.bed",
                shell=True)
subprocess.call("bedtools genomecov -bg -i readssort.bed -g chromsizes.genome > Covreads.bedGraph",
                shell=True)
subprocess.call("bedGraphToBigWig Covreads.bedGraph http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes BigWig.bw",
                shell=True)
subprocess.call("computeMatrix reference-point -S BigWig.bw -R genome.bed -b 3000 -a 3000 --referencePoint TSS -o matrixTSS.gz", 
                shell=True)
subprocess.call("plotHeatmap -m matrixTSS.gz -out TSS.png",
                shell=True)
subprocess.call("computeMatrix reference-point -S BigWig.bw -R genome.bed -b 3000 -a 3000 --referencePoint TES -o matrixTES.gz",
                shell=True)
subprocess.call("plotHeatmap -m matrixTES.gz -out TES.png",
                shell=True)
subprocess.call("computeMatrix reference-point -S BigWig.bw -R genome.bed -b 3000 -a 3000 --referencePoint center -o matrixcenter.gz",
                shell=True)
subprocess.call("plotHeatmap -m matrixcenter.gz -out center.png", 
                shell=True)

