/*SNP Calling Pipeline
Does the following:
1. Align FASTQ reads to a reference Genome to create an alignment file
2. Processing the alignment file
3. Calling the variants 
*/
#!/bin/bash

VERBOSE=0
realign=0
gunzip=0
index=0
usage_info=0
while getopts "a:b:r:eo:f:zvih" option;
do
case $option in
a)reads1=$OPTARG;; #location of first reads file
b)reads2=$OPTARG;; #location of second reads file
r)ref=$OPTARG;; #location of reference genome
e)realign=1;; #user wants realignment
o)output=$OPTARG;; #output vcf file name that the user wants
f)millsFile=$OPTARG;; #location of Mills Indel File
z)gunzip=1;; #user wants gunzipped vcf file
v)VERBOSE=1;; #user wants verbose mode
i)index=1;; #user wants indexed file
h)usage_info=1;; #user wants usage information
esac
done
u=1
#if user wants usage information
if [ $usage_info -eq $u ];
then
echo "Welcome to the SNP calling pipeline. Here is the usage information." 
echo " The options and their meanings are mentioned below:
-a: Input reads file location - pai 1
-b: Input reads file location - pair 2
-r Reference genome file location
-e Perform read re-alignment
-o Output VCF file name in this format - outputfilename.vcf.gz
-f Mills file location
-z Output VCF file should be gunzipped (*.vcf.gz)
-v Verbose mode; print each instruction/command to tell the user 
what your script is doing right now
-i Index your output BAM file (using samtools index)
-h Print usage information"
exit
fi
#provided usage info and exit program

#checking whether input reads files exist
zero=0
test -e $reads1
temp=$(echo $?)
if [ $temp -ne $zero ]; then
echo "First input file doesn't exist. Please enter a valid file."
exit
fi
test -e $reads2
temp2=$(echo $?)
if [ $temp2 -ne $zero ]; then
echo "Second input file doesn't exist. Please enter a valid file."
exit
fi
#checking of input files done
#checking whether reference genome exists
test -e $ref
temp3=$(echo $?)
if [ $temp3 -ne $zero ]; then
echo "Reference genome file doesn't exist. Please enter a valid file."
exit
fi
#checking of reference genome done

#checking vcf output file
test -e $output
temp=$(echo $?)
if [ $temp -eq $zero ]; then
read -p "File exists. Do you want to overwrite file and continue? [Y/N]:" confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1 
echo "You have chosen to continue"
fi

#output vcf file check done

if [ $VERBOSE -eq $u ];
then
echo "Pipeline started -- Beginning mapping"
fi
#for mapping 
bwa index $ref #ref.fa is the reference file input by user
if [ $VERBOSE -eq $u ];
then
echo "Reference indexed"
fi
bwa mem -R '@RG\tID:test\tSM:test' $ref $reads1 $reads2 > lane.sam #reads mapped to ref
#lane.sam is the sam file which needs to be sorted further
if [ $VERBOSE -eq $u ];
then
echo "Reads mapped to reference"
fi
samtools fixmate -O bam lane.sam lane_fixmate.bam
samtools sort -O bam -o lane_sorted.bam -T lane_temp lane_fixmate.bam
samtools index lane_sorted.bam
if [ $VERBOSE -eq $u ];
then
echo "mapping is done"
fi
#mapping done

#for improvement
if [ $VERBOSE -eq $u ];
then
echo "Improvement started"
fi
java -Xmx2g -jar $HOME/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -I lane_sorted.bam -o lane.intervals --known $millsFile --log_to_file stambe6.log

java -Xmx4g -jar $HOME/bin/GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I lane_sorted.bam -targetIntervals lane.intervals -o lane_realigned.bam --log_to_file stambe6.log

#if user wants read realignment
rea=1
if [ $realign -eq $rea ]; then
java -Xmx2g -jar $HOME/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -I lane_realigned.bam -o lane1.intervals --known $millsFile --log_to_file stambe6.log
if [ $VERBOSE -eq $u ];
then
echo "RealignerTargetCreator done."
fi
java -Xmx4g -jar GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I lane_realigned.bam -targetIntervals lane1.intervals -o lane_realigned1.bam
mv lane_realigned1.bam lane_realigned.bam
fi
#read realignment done
if [ $VERBOSE -eq $u ];
then
echo "read realignment done"
fi
#if user wants indexed
if [ $index -eq 1 ]; then
samtools index lane_realigned.bam
fi
#indexing done
if [ $VERBOSE -eq $u ];
then
echo "Improvement done"
fi
#improvement done

#for variant calling
if [ $VERBOSE -eq $u ];
then
echo "Beginning variant calling"
fi
bcftools mpileup -Ou -f $ref lane_realigned.bam | bcftools call -vmO z -o output.vcf.gz #output.vcf.gz is the output file
zcat output.vcf.gz | gzip -c > $output #changed to file name that user wants
if [ $VERBOSE -eq $u ];
then
echo "Variant calling done"
fi
#variant calling done


#if user wants gunzipped vcf file
if [ $gunzip -eq 1 ]; then
gunzip $output
fi
#done gunzip
if [ $VERBOSE -eq $u ];
then
echo "Beginning conversion to bed file"
fi
#preparing file for conversion to bed file
zcat output.vcf.gz | sed '/##/d' | gzip -c > outputfile.gz #file without headers
gunzip outputfile.gz #now this file is outputfile


#for vcf to bed file
#for the test.vcf file, find the number of lines in its txt format
wc -l outputfile | cat > num_lines_temp
num_lines=$( awk '{print $1 }' num_lines_temp)
rm num_lines_temp
x=2

while [ $x -le $num_lines ] 
do
if [ $VERBOSE -eq $u ]; 
then
echo "Loop started. Loop number $x . Total loops to run $num_lines"
fi
head -n $x outputfile | tail +$x | awk '{print $1" "$2" "$3" "$4" "$5}' |cat > currentline #extract one line
alt_temp=$( awk '{print $5}' currentline )
alt=$(echo -n "$alt_temp" | wc -c)
ref_temp=$( awk '{print $4}' currentline )
reference=$(echo -n "$ref_temp" | wc -c)
pos=$(awk '{print $2}' currentline)
clmn3=$(($alt - $reference + $pos)) #value for column 3
clmn4=$(($alt -$reference)) #value of column 4
if [ $clmn4 -ne $zero ]; then
awk '{print $1" " $2" "'$clmn3'" "'$clmn4'}' currentline | sed 's/chr//g' |cat >> indels.txt
else
awk '{print $1" " $2" "'$clmn3'" "'$clmn4'}' currentline | sed 's/chr//g' |cat >> snps.txt
fi
x=$((x + 1))
done
if [ $VERBOSE -eq $u ];
then
echo "Converted to bed file"
fi
#converted from vcf to indels.txt and snps.txt
