library(ggplot2)
raw.sim.results<-read.csv('raw_canonical_sim.csv', header=TRUE)
partitioned.sim.results<-split(raw.sim.results, raw.sim.results$total_energy)

mean_energies<-as.vector(unlist(lapply(partitioned.sim.results, function(df) { mean(df$mean_energy)})))
std_energies<-as.vector(unlist(lapply(partitioned.sim.results, function(df) { mean(df$std_energy)})))
betas<-as.vector(unlist(lapply(partitioned.sim.results, function(df) { mean(df$beta)})))

stat.sim.results<-data.frame(mean_energy=mean_energies, std_energy=std_energies, beta=betas)

betafcn<-function(avg.e) { log((1+avg.e)/avg.e) }
stdefcn<-function(avg.e) { sqrt(avg.e*(1+avg.e)) }
avg.es<-seq(0, 0.01, 0.0001)
model<-data.frame(mean_energy=avg.es, beta=betafcn(avg.es), std_energy=stdefcn(avg.es))

g<-ggplot(stat.sim.results, aes(x=log(mean_energy), y=beta))+geom_point()
g<-g+geom_line(data=model, aes(x=log(mean_energy), y=beta))
g<-g+xlab(expression(textstyle(log)*(bar(epsilon))))+ylab(expression(beta))

g1<-ggplot(stat.sim.results, aes(x=mean_energy, y=std_energy))+geom_point()
g1<-g1+geom_line(data=model, aes(x=mean_energy, y=std_energy))
g1<-g1+xlab(expression(textstyle(log)*(bar(epsilon))))+ylab(expression(sqrt(bar((Delta*epsilon)^2))))

