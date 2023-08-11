setwd("/restricted/projectnb/pulmseq/Allegro/Bronch_microRNA/edrn243/Analysis/LGRC_paper/mir34c/")
frequency.table = read.table("mir34c_sequence_counts.txt", header=T, stringsAsFactor=F)
frequency.table=frequency.table[grepl("GGCAGTG", rownames(frequency.table)),]
frequencies=rowSums(frequency.table)
# frequencies=sort(frequencies, decreasing=T)
sequence.and.frequencies=cbind(rownames(frequency.table), frequencies)
colnames(sequence.and.frequencies)=c("Sequence", "Summing")
sequence.and.frequencies.ordered=sequence.and.frequencies[order(as.numeric(sequence.and.frequencies[,"Summing"]), decreasing = T),]
write.csv(sequence.and.frequencies.ordered, "MostRelevantSequencies.fin.csv")