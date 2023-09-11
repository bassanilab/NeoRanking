source configure.sh

echo "Mutation NCI CD8 SNV"
grep -cE "[[:space:]]NCI[[:space:]].*[[:space:]]CD8[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Mutation_data.txt 
echo "Mutation TESLA CD8 SNV"
grep -cE "[[:space:]]TESLA[[:space:]].*[[:space:]]CD8[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Mutation_data.txt 
echo "Mutation HiTIDE CD8 SNV"
grep -cE "[[:space:]]HiTIDE[[:space:]].*[[:space:]]CD8[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Mutation_data.txt 
echo "Mutation NCI negative SNV"
grep -cE "[[:space:]]NCI[[:space:]].*[[:space:]]negative[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Mutation_data.txt 
echo "Mutation TESLA negative SNV"
grep -cE "[[:space:]]TESLA[[:space:]].*[[:space:]]negative[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Mutation_data.txt 
echo "Mutation HiTIDE negative SNV"
grep -cE "[[:space:]]HiTIDE[[:space:]].*[[:space:]]negative[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Mutation_data.txt 
echo "Mutation NCI not screened SNV"
grep -cE "[[:space:]]NCI[[:space:]].*[[:space:]]not_tested[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Mutation_data.txt 
echo "Mutation TESLA not screened SNV"
grep -cE "[[:space:]]TESLA[[:space:]].*[[:space:]]not_tested[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Mutation_data.txt 
echo "Mutation HiTIDE not screened SNV"
grep -cE "[[:space:]]HiTIDE[[:space:]].*[[:space:]]not_tested[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Mutation_data.txt 

echo "Neopep NCI CD8 SNV"
grep -cE "[[:space:]]NCI[[:space:]].*[[:space:]]CD8[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Neopep_data.txt 
echo "Neopep TESLA CD8 SNV"
grep -cE "[[:space:]]TESLA[[:space:]].*[[:space:]]CD8[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Neopep_data.txt 
echo "Neopep HiTIDE CD8 SNV"
grep -cE "[[:space:]]HiTIDE[[:space:]].*[[:space:]]CD8[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Neopep_data.txt 
echo "Neopep NCI negative SNV"
grep -cE "[[:space:]]NCI[[:space:]].*[[:space:]]negative[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Neopep_data.txt 
echo "Neopep TESLA negative SNV"
grep -cE "[[:space:]]TESLA[[:space:]].*[[:space:]]negative[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Neopep_data.txt 
echo "Neopep HiTIDE negative SNV"
grep -cE "[[:space:]]HiTIDE[[:space:]].*[[:space:]]negative[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Neopep_data.txt 
echo "Neopep NCI not screened  SNV"
grep -cE "[[:space:]]NCI[[:space:]].*[[:space:]]not_tested[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Neopep_data.txt 
echo "Neopep TESLA not screened  SNV"
grep -cE "[[:space:]]TESLA[[:space:]].*[[:space:]]not_tested[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Neopep_data.txt 
echo "Neopep HiTIDE not screened SNV"
grep -cE "[[:space:]]HiTIDE[[:space:]].*[[:space:]]not_tested[[:space:]].*[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Neopep_data.txt 

echo "Mutation SNV"
grep -cE "[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Mutation_data.txt 
echo "Neopep SNV"
grep -cE "[[:space:]]SNV[[:space:]]" $NEORANKING_RESOURCE/data/Neopep_data.txt 

