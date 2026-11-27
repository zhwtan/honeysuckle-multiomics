#Author: liyan
#Mail: 383660766@qq.com
################################################################################
library(gtools)
library(circlize)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(scales)
library(argparser)
#######
#设置参数类型、默认值及说明
argv<-arg_parser('')
argv<-add_argument(argv,"--conf",help="config文件路径")
argv<-add_argument(argv,"--fai",help="染色体长度信息fai文件路径")
argv<-add_argument(argv,"--link",help="link文件路径")
argv<-add_argument(argv,"--link_col",default='#2969a6',help="link连接颜色")
argv<-add_argument(argv,"--width",type="double",default=8,help="输出pdf的宽度")
argv<-add_argument(argv,"--height",type="double",default=8,help="输出pdf的高度")
argv<-add_argument(argv,"--blank",type="logical",default=FALSE,help="是否显示大gap")
argv<-add_argument(argv,"--outputtype",default="pdf",help="输出文件类型:pdf/jpg/png")
argv<-parse_args(argv)
#将参数传入变量中
conf<-argv$conf
fai<-argv$fai
wid<-argv$width
hei<-argv$height
outputtype<-argv$outputtype
blank<-argv$blank
link<-argv$link
link.col<-argv$link_col
#######
FAI<-read.table(fai,stringsAsFactors=F,sep="\t")
names(FAI)[2]<-"len"
FAI<-FAI[mixedorder(FAI$V1),]
cols<-col_factor("#9EA0A1",FAI$V1)(FAI$V1)
bcols<-adjustcolor(cols,alpha.f=0.4)
maxlenchr<-max(nchar(FAI$V1))
gmar<-maxlenchr*0.15*5*0.5+1
chrcol<-data.frame(V1=mixedsort(unique(FAI$V1)),col="black",stringsAsFactors=F)
outter.track<-function(fai,height=0.05,mar=rep(gmar,4),cell.padding = c(0, 0, 0, 0),start.degree=270){
        par(mar=mar,xpd=NA)
        geno<-fai
        circos.clear()
        circos.par(cell.padding = cell.padding,start.degree = start.degree)
		if (!blank){
			ng<-data.frame(V1=geno$V1,len=0,stringsAsFactors=F)
			ng<-rbind(geno,ng)
			circos.initialize(factors=factor(ng$V1,levels=mixedsort(unique(ng$V1))),x=as.numeric(ng$len))
			circos.trackPlotRegion(factors=factor(geno$V1,levels=mixedsort(unique(geno$V1))),x=as.numeric(geno$len),ylim=c(0,1),track.height=height,bg.col=cols,bg.border=NA)
		}else{
			bfai<-data.frame(V1=c(geno$V1,"zzzzz",'zzzzz'),len=c(rep(0,length(geno$V1)),0,max(geno$len)*0.6),stringsAsFactors=F)
			bgeno<-rbind(geno,bfai)
			ncols<-c(cols,"white")
			circos.initialize(factors=factor(bgeno$V1,levels=mixedsort(unique(bgeno$V1))),x=as.numeric(bgeno$len))
			circos.trackPlotRegion(factors=factor(bgeno$V1,levels=mixedsort(unique(bgeno$V1))),x=as.numeric(bgeno$len),ylim=c(0,1),track.height=height,bg.col=ncols,bg.border=NA)
		}
        for (index in chrcol$V1){
			sect<-geno[geno$V1==index,]
			at<-pretty(range(c(0,sect$len)),n=3)
			lab<-paste(at/1000000,"M",sep="")
			circos.axis(major.at=at,labels=lab,labels.cex=0.3,major.tick.percentage=1,sector.index=index,track.index=1,lwd=0.4)
			nlen<-nchar(geno$V1)
			#tcex<-ifelse(any(nlen>=10),0.5,ifelse(any(nlen>=6),0.8,1.1))
			#circos.text(mean(sect$len),2,labels=index,adj=c(0,0.5),facing='clockwise',xpd=NA,col=cols[which(chrcol$V1==index)],font=2,sector.index=index,track.index=1,cex=tcex)
			circos.text(mean(c(0,sect$len)),2,labels=index,adj=c(0,0.5),facing='clockwise',xpd=NA,col=cols[which(chrcol$V1==index)],font=2,sector.index=index,track.index=1)
        }
}
############################
innter.track<-function(path,type,height=0.07,backgroud.color,border.color,border.lwd,main.color,point.cex,heatmap.palette,line.width){
			dat<-read.table(path,stringsAsFactors=F,sep="\t")
			dat$V2<-as.numeric(dat$V2)
			dat$V3<-as.numeric(dat$V3)
			if (type=="point"){
				if(!blank){
					circos.trackPlotRegion(ylim=range(dat$V4),track.height=height,bg.col=backgroud.color,bg.border=border.color,bg.lwd=border.lwd)
				}else{
					circos.trackPlotRegion(ylim=range(dat$V4),track.height=height,bg.col=backgroud.color,bg.border=c(rep(border.color,length(chrcol$V1)),'white'),bg.lwd=border.lwd)
				}
				for (index in chrcol$V1){
					pdat<-dat[dat$V1==index,]
					circos.points(x=(pdat$V2+pdat$V3)/2,y=pdat$V4,sector.index=index,pch=16,col=main.color,cex=point.cex)
				}
			}else if (type=="line"){
				if(!blank){
					circos.trackPlotRegion(ylim=range(dat$V4),track.height=height,bg.col=backgroud.color,bg.border=border.color,bg.lwd=border.lwd)
				}else{
					circos.trackPlotRegion(ylim=range(dat$V4),track.height=height,bg.col=backgroud.color,bg.border=c(rep(border.color,length(chrcol$V1)),'white'),bg.lwd=border.lwd)
				}
				for (index in chrcol$V1){
					pdat<-dat[dat$V1==index,]
					circos.lines(x=(pdat$V2+pdat$V3)/2,y=pdat$V4,sector.index=index,type="l",area=TRUE,col='white',border=main.color,lwd=line.width)
				}
			}else if (type=='hist'){
				if(!blank){
					circos.trackPlotRegion(ylim=range(dat$V4),track.height=height,bg.col=backgroud.color,bg.border=border.color,bg.lwd=border.lwd)
				}else{
					circos.trackPlotRegion(ylim=range(dat$V4),track.height=height,bg.col=backgroud.color,bg.border=c(rep(border.color,length(chrcol$V1)),'white'),bg.lwd=border.lwd)
				}
                for (index in chrcol$V1){
                    pdat<-dat[dat$V1==index,]
                    circos.lines(x=(pdat$V2+pdat$V3)/2,y=pdat$V4,sector.index=index,type="h",col=main.color,lwd=line.width)
				}
			}else if (type=="rect"){
				if(!blank){
					circos.trackPlotRegion(ylim=range(dat$V4),track.height=height,bg.col=backgroud.color,bg.border=border.color,bg.lwd=border.lwd)
				}else{
					circos.trackPlotRegion(ylim=range(dat$V4),track.height=height,bg.col=backgroud.color,bg.border=c(rep(border.color,length(chrcol$V1)),'white'),bg.lwd=border.lwd)
				}
				for (index in chrcol$V1){
					pdat<-dat[dat$V1==index,]
					circos.rect(xleft=pdat$V2,ybottom=rep(0,length(pdat$V1)),xright=pdat$V3,ytop=pdat$V4,sector.index=index,col=NA,border=main.color,lwd=line.width)
				}
			}else if (type=='heatmap'){
				rcols<-col_numeric(heatmap.palette,dat$V4)(dat$V4)
				if(!blank){
					circos.trackPlotRegion(ylim=c(0,1),track.height=height,bg.col=backgroud.color,bg.border=border.color,bg.lwd=border.lwd)
				}else{
					circos.trackPlotRegion(ylim=c(0,1),track.height=height,bg.col=backgroud.color,bg.border=c(rep(border.color,length(chrcol$V1)),'white'),bg.lwd=border.lwd)
				}
				for (index in chrcol$V1){
                    pdat<-dat[dat$V1==index,]
					icol<-rcols[dat$V1==index]
                    circos.rect(xleft=pdat$V2,ybottom=rep(0,length(pdat$V1)),xright=pdat$V3,ytop=rep(1,length(pdat$V1)),sector.index=index,col=icol,border=icol,lwd=line.width)
                }
			}else{
				print("Type should be assigned.")
			}
}
##########################
link.plot<-function(x){
	lkcol<-adjustcolor(link.col,alpha.f=0.4)
	x<-as.matrix(x)
    sector.index1<-x[1]
    sector.index2<-x[4]
    pos1<-as.numeric(c(x[2],x[3]))
    pos2<-as.numeric(c(x[5],x[6]))
	circos.link(sector.index1=sector.index1,point1=pos1,sector.index2=sector.index2,point2=pos2,col=lkcol,border=lkcol,lwd=0.0000000000001)
}


##########################
#如果要pdf格式，生成该格式图片
if(outputtype=="pdf"){
pdf("result.pdf",width=wid,height=hei)
outter.track(fai=FAI,start.degree=80)
CONF<-read.table(conf,stringsAsFactors=F)
for (i in seq(CONF$V1)){
	innter.track(path=CONF$V2[i],type=CONF$V3[i],height=CONF$V4[i],backgroud.color=CONF$V5[i],border.color=CONF$V6[i],border.lwd=CONF$V7[i],main.color=CONF$V8[i],point.cex=as.numeric(CONF$V9[i]),heatmap.palette=CONF$V10[i],line.width=CONF$V11[i])
}
#########################
if(!is.na(link)){
        lk<-read.table(link,stringsAsFactors=F,sep="\t")
        adply(lk,1,link.plot)
}
dev.off()

}
#如果要png格式，生成该格式图片
if(outputtype=="png"){
png("result.png",width=0.5*4000,height=0.5*4000,units="px",res=300)
outter.track(fai=FAI,start.degree=80)
CONF<-read.table(conf,stringsAsFactors=F)
head(CONF)
for (i in seq(CONF$V1)){
        innter.track(path=CONF$V2[i],type=CONF$V3[i],height=CONF$V4[i],backgroud.color=CONF$V5[i],border.color=CONF$V6[i],border.lwd=CONF$V7[i],main.color=CONF$V8[i],point.cex=as.numeric(CONF$V9[i]),heatmap.palette=CONF$V10[i],line.width=CONF$V11[i])
}
#########################
if(!is.na(link)){
        lk<-read.table(link,stringsAsFactors=F,sep="\t")
        adply(lk,1,link.plot)
}
dev.off()

}
#如果要jpg格式，生成该格式图片
if(outputtype=="jpg"){
jpeg("result.jpg",width=0.5*4000,height=0.5*4000,units="px",res=300)
outter.track(fai=FAI,start.degree=80)
CONF<-read.table(conf,stringsAsFactors=F)
for (i in seq(CONF$V1)){
        innter.track(path=CONF$V2[i],type=CONF$V3[i],height=CONF$V4[i],backgroud.color=CONF$V5[i],border.color=CONF$V6[i],border.lwd=CONF$V7[i],main.color=CONF$V8[i],point.cex=as.numeric(CONF$V9[i]),heatmap.palette=CONF$V10[i],line.width=CONF$V11[i])
}
#########################
if(!is.na(link)){
        lk<-read.table(link,stringsAsFactors=F,sep="\t")
        adply(lk,1,link.plot)
}
dev.off()

}
#########################
#if(!is.na(link)){
#	lk<-read.table(link,stringsAsFactors=F,sep="\t")
#	adply(lk,1,link.plot)
#}
#dev.off()
