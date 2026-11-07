# Original Figure 3b script: only input/output paths adapted; no statistical or logical edits.
# Supply the original CSV files; these are not replaced by the 1003 DEG/ORA results.
library(ggplot2)
dir.create('results/Figure_3b', recursive = TRUE, showWarnings = FALSE)

preg_go = read.csv('data/legacy_Figure3b/preg sp vc preg bm.csv')
preg_go$logp = -log10(preg_go$pvalue)
ggplot(preg_go,aes(x=logp,y=reorder(Description,logp)))+geom_col()+
  theme_classic()+theme(text = element_text(size=14))+ylab('GO terms')+
  xlab('-log10 pvalue')+scale_x_continuous(expand =c(0,0))
ggsave('results/Figure_3b/preg sp vc preg bm.pdf',
       width = 10,height = 8)
preg_go = read.csv('data/legacy_Figure3b/preg sp vs non-preg sp.csv')
preg_go$logp = -log10(preg_go$pvalue)
ggplot(preg_go,aes(x=logp,y=reorder(Description,logp)))+geom_col()+
  theme_classic()+theme(text = element_text(size=14))+ylab('GO terms')+
  xlab('-log10 pvalue')+scale_x_continuous(expand =c(0,0))
ggsave('data/legacy_Figure3b/preg sp vs non-preg sp.csv.pdf',
       width = 10,height = 8)


### adj
preg_go$logadj = -log10(preg_go$p.adjust)
ggplot(preg_go,aes(x=logadj,y=reorder(Description,logadj)))+geom_col()+
  theme_classic()+theme(text = element_text(size=14))+ylab('GO terms')+
  xlab('-log10 adj pvalue')+scale_x_continuous(expand =c(0,0))
ggsave('results/Figure_3b/preg sp vc preg bm.pdf',
       width = 10,height = 8)
preg_go = read.csv('data/legacy_Figure3b/preg sp vs non-preg sp.csv')
preg_go$logadj = -log10(preg_go$p.adjust)
ggplot(preg_go,aes(x=logadj,y=reorder(Description,logadj)))+geom_col()+
  theme_classic()+theme(text = element_text(size=14))+ylab('GO terms')+
  xlab('-log10 adj pvalue')+scale_x_continuous(expand =c(0,0))
ggsave('data/legacy_Figure3b/preg sp vs non-preg sp.csv.pdf',
       width = 10,height = 8)
