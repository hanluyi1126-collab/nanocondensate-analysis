#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

d <- "data"
out <- file.path("outputs", "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

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
