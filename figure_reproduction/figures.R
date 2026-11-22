#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(svglite);library(ragg)})
root <- commandArgs(trailingOnly=TRUE)[1]
d <- file.path(root,"data"); out <- file.path(root,"figures")
theme_set(theme_classic(base_size=9,base_family="Helvetica")+
 theme(axis.line=element_line(linewidth=.35),axis.ticks=element_line(linewidth=.3),
 axis.text=element_text(colour="black"),plot.title=element_text(size=10,face="bold"),
 legend.title=element_blank(),legend.text=element_text(size=8),legend.key.width=unit(9,"mm"),
 plot.margin=margin(3,4,3,3,"mm"),strip.background=element_blank()))
cols<-c("#333333","#347D9E","#885080")
savef<-function(p,name,h=120){
 p<-p+plot_annotation(tag_levels="A") & theme(plot.tag=element_text(face="bold",size=11))
 svglite(file.path(out,paste0(name,".svg")),width=6.5,height=h/25.4);print(p);dev.off()
 pdf(file.path(out,paste0(name,".pdf")),width=6.5,height=h/25.4,family="Helvetica",useDingbats=FALSE);print(p);dev.off()
 agg_png(file.path(out,paste0(name,".png")),width=165,height=h,units="mm",res=300);print(p);dev.off()
 agg_tiff(file.path(out,paste0(name,".tiff")),width=165,height=h,units="mm",res=600,compression="lzw");print(p);dev.off()
}
f<-read.csv(file.path(d,"drop/fusion_verified.csv"));f$potential<-factor(f$potential_mV,levels=c(0,10,25),labels=c("0 mV","10 mV","25 mV"));neut<-subset(f,potential_mV==0)
a<-ggplot(neut,aes(radius_um,encounter_s))+geom_line(colour=cols[2],linewidth=.65)+scale_x_log10()+scale_y_log10()+labs(x="Droplet radius (µm)",y=expression(nu[enc]~(s^{-1})),title="Brownian encounters")
b<-ggplot(neut,aes(radius_um,capillary_kBT))+geom_line(colour=cols[2],linewidth=.65)+scale_x_log10()+labs(x="Droplet radius (µm)",y=expression(Delta*G[gamma]/k[B]*T),title="Capillary free energy")
c<-ggplot(f,aes(radius_um,barrier_kBT,colour=potential))+geom_line(linewidth=.65)+scale_x_log10()+scale_colour_manual(values=cols)+labs(x="Droplet radius (µm)",y=expression(U[el](0)/k[B]*T),title="Formal HHF contact energy")+theme(legend.position="top")
e<-ggplot(f,aes(radius_um,log10_contact_rate_s,colour=potential))+geom_line(linewidth=.65)+scale_x_log10()+scale_colour_manual(values=cols)+labs(x="Droplet radius (µm)",y=expression(log[10]*nu[contact]~(s^{-1})),title="Formal contact arrival rate")+theme(legend.position="top")
savef((a+b)/(c+e),"Fig_S1",125)
g<-read.csv(file.path(d,"drop/state_grid.csv"));g$state<-factor(g$state,levels=c("Dispersed","Coarsening droplets","Nanocondensate window","Gel/solid"),labels=c("Dispersed","Coarsening","Candidate window","Network arrest"))
a<-ggplot(g,aes(sticker_fraction,charged_fraction,fill=state))+geom_raster()+scale_fill_manual(values=c("#E0E3E5","#DFC18C","#86B5AA","#9B7C9F"),drop=FALSE)+coord_fixed()+labs(x="Sticker fraction",y="Charged fraction",title="Illustrative state map")+theme(legend.position="bottom",legend.text=element_text(size=7))+guides(fill=guide_legend(nrow=2))
bud<-data.frame(diameter_nm=exp(seq(log(50),log(2000),length.out=120)));bud$charge<-6*.03/(5*.005*bud$diameter_nm);bud$coverage<-6/(10*.005*bud$diameter_nm)
bb<-rbind(data.frame(diameter_nm=bud$diameter_nm,fraction=bud$charge,criterion="Charge budget"),data.frame(diameter_nm=bud$diameter_nm,fraction=bud$coverage,criterion="Complete coverage"))
b<-ggplot(bb,aes(diameter_nm,fraction,colour=criterion))+geom_hline(yintercept=1,linetype=2,linewidth=.3)+geom_line(linewidth=.7)+scale_x_log10()+scale_y_log10()+scale_colour_manual(values=cols[2:3])+labs(x="Droplet diameter (nm)",y="Required molecular fraction",title="Interfacial molecule budgets")+theme(legend.position="bottom",legend.text=element_text(size=7))
savef(a+b,"Fig_S2",98)
dia<-read.csv(file.path(d,"membrane/diameter.csv"));cv<-read.csv(file.path(d,"membrane/convergence.csv"));
stress<-rbind(data.frame(diameter=dia$diameter_nm,tension=dia$centre_tension_pN_per_nm,site="Centre"),data.frame(diameter=dia$diameter_nm,tension=dia$maximum_tension_pN_per_nm,site="Spatial maximum"))
a<-ggplot(stress,aes(diameter,tension,colour=site))+geom_line(linewidth=.65)+geom_hline(yintercept=8,linetype=2,linewidth=.3)+scale_colour_manual(values=cols[2:3])+labs(x="Diameter (nm)",y=expression(tau~(pN~nm^{-1})),title="Distinct local stress observables")+theme(legend.position="top")
b<-ggplot(dia,aes(diameter_nm,maximum_tension_position_fraction))+geom_line(colour=cols[3],linewidth=.65)+labs(x="Diameter (nm)",y=expression(r[max]/a[c]),title="Peak position within the contact")+ylim(0,1.15)
ref<-subset(cv,membrane_elements==3840,c(diameter_nm,centre_tension_pN_per_nm,maximum_tension_pN_per_nm));names(ref)[2:3]<-c("refc","refm");cv<-merge(cv,ref,by="diameter_nm");cv$centre_error<-100*abs(cv$centre_tension_pN_per_nm/cv$refc-1)
c<-ggplot(subset(cv,membrane_elements<3840),aes(membrane_elements,pmax(centre_error,1e-5),group=factor(diameter_nm)))+geom_line(colour=cols[2],alpha=.6,linewidth=.45)+geom_point(size=1.2,colour=cols[2])+scale_x_log10(breaks=c(240,480,960,1920))+scale_y_log10()+labs(x="Membrane elements",y="Centre error versus 3,840 (%)",title="Seven diameter refinement series")
res<-do.call(rbind,lapply(c("diameter","convergence","cohesive_grid","ensemble"),function(n)read.csv(file.path(d,paste0("membrane/",n,".csv")))))
e<-ggplot(res,aes(force_residual))+stat_ecdf(geom="step",colour=cols[2],linewidth=.6)+geom_vline(xintercept=1e-6,linetype=2,linewidth=.3)+scale_x_log10()+labs(x="Force residual (pN)",y="Cumulative fraction",title="All 3,838 retained solutions")
savef((a+b)/(c+e),"Fig_S3",125)
r<-read.csv(file.path(d,"membrane/lipid_replicates.csv"));fit<-lm(force_over_2pi_pN~radius_nm,r);boot<-read.csv(file.path(d,"membrane/bootstrap_verified.csv"))
a<-ggplot(r,aes(radius_nm,force_over_2pi_pN))+geom_point(colour=cols[2],alpha=.7,size=1.8)+geom_abline(intercept=coef(fit)[1],slope=coef(fit)[2],colour=cols[3],linewidth=.6)+labs(x="Restrained pore radius (nm)",y=expression(F[hold]/(2*pi)~(pN)),title="24 independent trajectory means")
b<-ggplot(boot,aes(line_tension_pN))+geom_histogram(bins=45,fill=cols[2],colour="white",linewidth=.1)+geom_vline(xintercept=quantile(boot$line_tension_pN,c(.025,.975)),linetype=2,linewidth=.35)+labs(x=expression(lambda~(pN)),y="Bootstrap resamples",title="Stratified uncertainty in line tension")
savef(a+b,"Fig_S4",87)
ph<-read.csv(file.path(d,"membrane/phase_boundaries.csv"));cr<-subset(ph,classification=="small-particle pore; large-particle wrap")
classes<-as.data.frame(table(ph$classification));names(classes)<-c("outcome","count");classes$label<-c("Pore at both limits","Small wrap / large pore","Small pore / large wrap","Wrap at both limits")[match(classes$outcome,c("pore at both diameter limits","small-particle wrap; large-particle pore","small-particle pore; large-particle wrap","wrap at both diameter limits"))]
a<-ggplot(classes,aes(reorder(label,count),count))+geom_col(fill=cols[2],width=.65)+geom_text(aes(label=count),hjust=-.2,size=3)+coord_flip(clip="off")+scale_y_continuous(limits=c(0,200),expand=expansion(mult=c(0,0)))+labs(x=NULL,y="Parameter combinations",title="All 336 conditional outcomes")+theme(axis.text.y=element_text(size=7))
b<-ggplot(cr,aes(first_50_percent_crossing_nm))+geom_histogram(binwidth=6,fill=cols[3],colour="white",linewidth=.2)+geom_vline(xintercept=median(cr$first_50_percent_crossing_nm),linetype=2,linewidth=.4)+labs(x="First 50% crossing diameter (nm)",y="Directional branches",title="Subset of 128 directional branches")
pc<-read.csv(file.path(d,"membrane/selected_probability_curves.csv"));pc$time<-factor(pc$lysis_reference_time_ms)
c<-ggplot(pc,aes(diameter_nm,two_stage_pore_probability_before_wrap,colour=time))+geom_line(linewidth=.65)+geom_hline(yintercept=.5,linetype=2,linewidth=.3)+scale_colour_manual(values=c(cols[2],cols[3],"#9A773C","#333333"))+labs(x="Diameter (nm)",y="Conditional pore probability",title="Sensitivity to defect reference time",colour="Reference time (ms)")+theme(legend.position="right")
savef((a+b)/c+plot_layout(heights=c(1,1)),"Fig_S5",130)
writeLines(capture.output(sessionInfo()),file.path(root,"R_session_info.txt"))
