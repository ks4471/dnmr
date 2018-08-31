

cmapr<-function(modg,bkg,de_thresh=0.1,n_genes=5,logfc_thresh=0){  ## combine the stuffs below to use with function rather than combined stuff as is atm
####===================================================================================
####  1. module info  ----------------------------------------------------------------
##  lmod$up   - character vector of genes up-regulated, if not directional supply all in either $up or $down
##  lmod$down - optional down-regulated module genes
####  2. module dataset background  ----------------------------------------------------------------
##  bkg           - character vector of all genes used to generate modules
####  3. cmap database  ----------------------------------------------------------------
##  degdb         - additional input - degdat$drug_name_etc[,c('logFC','FDR')]  ## column names require this format
##  logfc_thresh - use as absolute -> automatically *-1 for downreg
cat('\n\tNOTE: this function REQUIRES properly formatted database named "degdb", cmap version available from :
    \t\thttps://www.dropbox.com/s/18l1w2jbqld2mej/cmap.enrich.data.Rdata?dl=0\n\n')


cat('\t\tcommon background correction\n')
cmap_bkg=bgcommon(degdb)
bkg=intersect(cmap_bkg,bkg)

degdb=lapply(degdb,function(x){x[rownames(x)%in%bkg,]})


modg=lapply(modg,function(x){x[x%in%bkg]})

# Head(degdb)


mstat=list()
mlfcw=list()
sigen=list()
sumry=list()


cat('\t\tperform enrichment, genes considered significantly regulated by drug if:
      \t\t de_thresh=',de_thresh,'
      \t\t n_genes=',n_genes,'
      \t\t logfc_thresh=',logfc_thresh,'\n\n
')






k=1
calncol=0
# idru=names(degdb)[3]
for(idru in names(degdb)){
  holder=list()

  holder$bkg=(degdb[[idru]])
  holder$drug_sigup=(holder$bkg[holder$bkg$FDR<=de_thresh & holder$bkg$logFC>=logfc_thresh,])
  holder$drug_sigdown=(holder$bkg[holder$bkg$FDR<=de_thresh & holder$bkg$logFC<=(logfc_thresh*(-1)),])    ##

  if(nrow(holder$drug_sigup)>n_genes | nrow(holder$drug_sigdown)>n_genes){

    if('up'%in%names(modg)){
      sigen[[idru]]$up=holder$bkg[modg$up,]
      sigen[[idru]]$up_down=holder$bkg[intersect(rownames(holder$drug_sigdown),modg$up),]
      sigen[[idru]]$up_up=holder$bkg[intersect(rownames(holder$drug_sigup),modg$up),]

      # mstat[[paste(idru,'up_downreg',sep='.')]]=unlist(fet(sampl=modg$up,bkgrnd=rownames(holder$bkg),success=rownames(holder$drug_sigdown),counts=F,tail='greater'))
      # mstat[[paste(idru,'up_upreg',sep='.')]]  =unlist(fet(sampl=modg$up,bkgrnd=rownames(holder$bkg),success=rownames(holder$drug_sigup),counts=F,tail='greater'))

      mstat[[idru]]$up_downreg=unlist(fet(sampl=rownames(sigen[[idru]]$up),bkgrnd=rownames(holder$bkg),success=rownames(holder$drug_sigdown),counts=F,tail='greater'))
      mstat[[idru]]$up_upreg  =unlist(fet(sampl=rownames(sigen[[idru]]$up),bkgrnd=rownames(holder$bkg),success=rownames(holder$drug_sigup),counts=F,tail='greater'))

      mlfcw[[idru]]$up_downreg=unlist(fet(
                samp.success=round(sum(abs(sigen[[idru]]$up_down)$logFC))
                ,bkgrnd.success=round(sum(abs(holder$drug_sigdown[!rownames(holder$drug_sigdown)%in%rownames(sigen[[idru]]$up_down),]$logFC)))
                ,samp.fail=round(sum(abs(sigen[[idru]]$up[!sigen[[idru]]$up%in%sigen[[idru]]$up_down]$logFC)))
                ,bkgrnd.fail=round(sum(abs(holder$bkg[!rownames(holder$bkg)%in%rownames(holder$drug_sigdown),]$FDR)))
                ,counts=T,tail='greater'))
      mlfcw[[idru]]$up_upreg  =unlist(fet(
                samp.success=round(sum(abs(sigen[[idru]]$up_up)$logFC))
                ,bkgrnd.success=round(sum(abs(holder$drug_sigup[!rownames(holder$drug_sigup)%in%rownames(sigen[[idru]]$up_up),]$logFC)))
                ,samp.fail=round(sum(abs(sigen[[idru]]$up[!sigen[[idru]]$up%in%sigen[[idru]]$up_up]$logFC)))
                ,bkgrnd.fail=round(sum(abs(holder$bkg[!rownames(holder$bkg)%in%rownames(holder$drug_sigup),]$FDR)))
                ,counts=T,tail='greater'))


    }

    if('down'%in%names(modg)){
      sigen[[idru]]$down=holder$bkg[modg$down,]
      sigen[[idru]]$down_down=holder$bkg[intersect(rownames(holder$drug_sigdown),modg$down),]
      sigen[[idru]]$down_up=holder$bkg[intersect(rownames(holder$drug_sigup),modg$down),]

      # mstat[[paste(idru,'down_downreg',sep='.')]]=unlist(fet(sampl=modg$down,bkgrnd=rownames(holder$bkg),success=rownames(holder$drug_sigdown),counts=F,tail='greater'))
      # mstat[[paste(idru,'down_upreg',sep='.')]]  =unlist(fet(sampl=modg$down,bkgrnd=rownames(holder$bkg),success=rownames(holder$drug_sigup),counts=F,tail='greater'))

      mstat[[idru]]$down_downreg=unlist(fet(sampl=rownames(sigen[[idru]]$down),bkgrnd=rownames(holder$bkg),success=rownames(holder$drug_sigdown),counts=F,tail='greater'))
      mstat[[idru]]$down_upreg=unlist(fet(sampl=rownames(sigen[[idru]]$down),bkgrnd=rownames(holder$bkg),success=rownames(holder$drug_sigup),counts=F,tail='greater'))

      mlfcw[[idru]]$down_upreg=unlist(fet(
                samp.success=round(sum(abs(sigen[[idru]]$down_up)$logFC))
                ,bkgrnd.success=round(sum(abs(holder$drug_sigup[!rownames(holder$drug_sigup)%in%rownames(sigen[[idru]]$down_up),]$logFC)))
                ,samp.fail=round(sum(abs(sigen[[idru]]$down[!sigen[[idru]]$down%in%sigen[[idru]]$down_up]$logFC)))
                ,bkgrnd.fail=round(sum(abs(holder$bkg[!rownames(holder$bkg)%in%rownames(holder$drug_sigup),]$FDR)))
                ,counts=T,tail='greater'))

      mlfcw[[idru]]$down_downreg  =unlist(fet(
                samp.success=round(sum(abs(sigen[[idru]]$down_down)$logFC))
                ,bkgrnd.success=round(sum(abs(holder$drug_sigdown[!rownames(holder$drug_sigdown)%in%rownames(sigen[[idru]]$down_down),]$logFC)))
                ,samp.fail=round(sum(abs(sigen[[idru]]$down[!sigen[[idru]]$down%in%sigen[[idru]]$down_down]$logFC)))
                ,bkgrnd.fail=round(sum(abs(holder$bkg[!rownames(holder$bkg)%in%rownames(holder$drug_sigdown),]$FDR)))
                ,counts=T,tail='greater'))

    }

    if('down'%in%names(modg) & 'up'%in%names(modg)){
      # nsuccs=mstat[[paste(idru,'down_upreg',sep='.')]][1:2]+mstat[[paste(idru,'up_downreg',sep='.')]][1:2]
      nsuccs=mstat[[idru]]$down_upreg[1:2]+mstat[[idru]]$up_downreg[1:2]
      nfails=c(length(unique(unlist(modg)))-nsuccs[1],nrow(holder$bkg)-nsuccs[2])
      # mstat[[paste(idru,'reverse',sep='.')]]=unlist(fet(samp.success=nsuccs[1], bkgrnd.success=nsuccs[2], samp.fail=nfails[1], bkgrnd.fail=nfails[2], tail = "greater",counts=T))
      mstat[[idru]]$reverse=unlist(fet(samp.success=nsuccs[1], bkgrnd.success=nsuccs[2], samp.fail=nfails[1], bkgrnd.fail=nfails[2], tail = "greater",counts=T))


      nsuclf=mlfcw[[idru]]$down_upreg[1:2]+mlfcw[[idru]]$up_downreg[1:2]
      nfailf=c(length(unique(unlist(modg)))-nsuclf[1],nrow(holder$bkg)-nsuclf[2])

      mlfcw[[idru]]$reverse=unlist(fet(samp.success=nsuclf[1], bkgrnd.success=nsuclf[2], samp.fail=nfailf[1], bkgrnd.fail=nfailf[2], tail = "greater",counts=T))


      # nsuccs=mstat[[paste(idru,'down_downreg',sep='.')]][1:2]+mstat[[paste(idru,'up_upreg',sep='.')]][1:2]
      nsuccs=mstat[[idru]]$down_downreg[1:2]+mstat[[idru]]$up_upreg[1:2]
      nfails=c(length(unique(unlist(modg)))-nsuccs[1],nrow(holder$bkg)-nsuccs[2])
      # mstat[[paste(idru,'mimic',sep='.')]]  =unlist(fet(samp.success=nsuccs[1], bkgrnd.success=nsuccs[2], samp.fail=nfails[1], bkgrnd.fail=nfails[2], tail = "greater",counts=T))
      mstat[[idru]]$mimic=unlist(fet(samp.success=nsuccs[1], bkgrnd.success=nsuccs[2], samp.fail=nfails[1], bkgrnd.fail=nfails[2], tail = "greater",counts=T))

      nsuclf=mlfcw[[idru]]$down_downreg[1:2]+mlfcw[[idru]]$up_upreg[1:2]
      nfailf=c(length(unique(unlist(modg)))-nsuclf[1],nrow(holder$bkg)-nsuclf[2])

      mlfcw[[idru]]$mimic=unlist(fet(samp.success=nsuclf[1], bkgrnd.success=nsuclf[2], samp.fail=nfailf[1], bkgrnd.fail=nfailf[2], tail = "greater",counts=T))
    }
  }
  k=lcount(k,length(degdb))
}

  cat('\t\tcompile output\n')
  rm(holder)
  k=1
  mpval=''
  mlfpv=''
  msens=''
  mspec=''
  for(idru in names(mstat)){
    holder=as.data.frame(t(as.data.frame(mstat[[idru]])))
    holder$sensitivity=holder$samp.success/(holder$samp.success+holder$samp.fail)
    holder$specificity=holder$samp.success/(holder$samp.success+holder$bkgrnd.success)

    dummy=holder[,'FETp',drop=F]
      colnames(dummy)=idru
    mpval=rmerge(mpval,dummy,verbose=F)

    dummy=holder[,'sensitivity',drop=F]
      colnames(dummy)=idru
    msens=rmerge(msens,dummy,verbose=F)

    dummy=holder[,'specificity',drop=F]
      colnames(dummy)=idru
    mspec=rmerge(mspec,dummy,verbose=F)


    holder=as.data.frame(t(as.data.frame(mlfcw[[idru]])))
    dummy=holder[,'FETp',drop=F]
      colnames(dummy)=idru
    mlfpv=rmerge(mlfpv,dummy,verbose=F)

    k=lcount(k,length(mstat))
  }

  mpval=t(mpval[rownames(mpval)!='1',colnames(mpval)!='x'])
  msens=t(msens[rownames(msens)!='1',colnames(msens)!='x'])
  mspec=t(mspec[rownames(mspec)!='1',colnames(mspec)!='x'])
  mlfpv=t(mlfpv[rownames(mlfpv)!='1',colnames(mlfpv)!='x'])

  orcol=c('reverse','up_downreg','down_upreg','mimic','up_upreg','down_downreg')
  mpval=mpval[,orcol[orcol%in%colnames(mpval)]]
  msens=msens[,orcol[orcol%in%colnames(msens)]]
  mspec=mspec[,orcol[orcol%in%colnames(mspec)]]
  mlfpv=mlfpv[,orcol[orcol%in%colnames(mlfpv)]]

  cat('\n\n')

  return(list(pval=as.data.frame(mpval),pvlf=as.data.frame(mlfpv),sensitivity=as.data.frame(msens),specificity=as.data.frame(mspec),fullp=mstat,fullfp=mlfcw,sigen=sigen))

}



resid<-function(dat_mat,cov_mat){
   cat('\tcolnames(dat_mat) same order as rownames(cov_mat)\n')
   print(Table(colnames(dat_mat)==rownames(cov_mat)))
  residd=t(lm(t(dat_mat)~.,data=cov_mat)$residuals)
  return(invisible(residd))
}


