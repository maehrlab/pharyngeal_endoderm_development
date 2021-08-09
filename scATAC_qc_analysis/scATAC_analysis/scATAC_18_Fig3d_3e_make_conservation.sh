# Authors: Jack Huey

# Downlaod the EuarchontoGlire phastCon bigiw
wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/phastCons60way/mm10.60way.phastCons60wayEuarchontoGlire.bw

# Create bed file from TableS3
cat TableS3.csv | column -s, -t | tail +2 | awk 'OFS = "\t" { print $2,$3+1,$4,$8,$9,$15,$16,$18 }' | awk '$5 == "Distal" || $5 == "Promoter"' > TableS3.bed

# Get average conservation over each peak
bigwigaverageoverbed -a mm10.60way.phastCons60wayEuarchontoGlire.bw TableS3.bed TableS3.conservation_euar.bed

# Subset into conserved, distal peaks by cluster
mkdir cluster_conservation_euar_top
while read cluster; do
cat TableS3.conservation_euar.bed | awk -v cluster=$cluster '$8 == cluster' | awk '$13 > 0.5' | grep Distal | cut -f1-3 > cluster_conservation_euar_top/TableS3.conservation_euar.0_5.distal.$cluster.bed
done < <(cat TableS3.conservation_euar.bed | cut -f8 | sort | uniq)
