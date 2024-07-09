library('data.table')
library('tidyverse')
juncs.long <- fread("/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/SplicingAnalysis/CombinedJuncTables/All.tsv.gz")

Annotations <- read_tsv('/project2/yangili1/qhauck/nmd_splicing_rules/testing_new_rules/testing_coordinate_system/Leaf2_junction_classifications.txt') %>%
  separate(Intron_coord, into=c("chrom", "start", "end"), sep="[:-]", convert=T, remove = F) %>%
  add_count(Intron_coord, name="NumberOfEntriesWithSameCoords")

Annotations.filtered <- Annotations %>%
  filter(NumberOfEntriesWithSameCoords==1)


juncs.long.summary <- juncs.long %>%
  dplyr::select(chrom, start, stop, strand, Dataset, Count) %>%
  group_by(Dataset, chrom, start, stop) %>%
  summarise(Sum=sum(Count)) %>%
  ungroup()

juncs.long.summary.joined <- juncs.long.summary %>%
  group_by(Dataset) %>%
  mutate(DatasetSum = sum(Sum)) %>%
  ungroup() %>%
  mutate(RPM = Sum/DatasetSum*1E6) %>%
  mutate(stop = stop+1) %>%
  inner_join(Annotations.filtered, by=c("chrom", "start", "stop"="end"))


to_plot <- juncs.long.summary.joined %>% 
  dplyr::select(-DatasetSum, -Sum) %>%
  pivot_wider(names_from="Dataset", values_from="RPM")


ggplot(to_plot, aes(x=Expression.Splicing, y=chRNA.Expression.Splicing)) +
  # geom_point(alpha=0.01, color='black') +
  # geom_density2d(color='red') +
  geom_hex(bins=70) +
  geom_abline(slope = 1, color='red') +
  geom_text(data = . %>%
              replace_na(list(Expression.Splicing=1E-5, chRNA.Expression.Splicing=1E-5)) %>%
              group_by(Coding) %>%
              summarise(med = round(median(chRNA.Expression.Splicing/Expression.Splicing, na.rm=F), 3),
                        n = n()) %>%
              ungroup(),
            aes(label=str_glue(" median FC={med}\n n={n}")),
            x=-Inf, y=Inf, vjust=1.1, hjust=0.) +
  scale_color_identity() +
  scale_fill_viridis_c(trans='log10', option = 'B') + 
  scale_x_continuous(trans="log10", breaks=c(1E-1, 10, 1000), labels=c("0.1", "10", "1000")) +
  scale_y_continuous(trans="log10", breaks=c(1E-1, 10, 1000), labels=c("0.1", "10", "1000")) +
  facet_wrap(~Coding) +
  theme_bw() +
  coord_fixed() +
  labs(x="Junction RPM, steady-state", y="Junction RPM, naRNA", title='Non-Coding vs Coding Splice Junctions')

