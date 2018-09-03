



fet<-function(sampl,bkgrnd,success,counts=F,samp.success,bkgrnd.success,samp.fail,bkgrnd.fail,tail='greater'){  # ,alternative='greater',...
# alternative ='greater'
# phyper(success_in_sample, success_in_bkgd, failure_in_bkgd, sample_size, lower.tail=TRUE)

#fisher.test(matrix(c(x, 13-x, 5-x, 34+x), 2, 2), alternative='less');
# Numerical parameters in order:
# (success-in-sample, success-in-left-part, failure-in-sample, failure-in-left-part).
  if(!counts){
      bkgrnd=bkgrnd[!(bkgrnd%in%sampl)]
      test_res=list(samp.success=sum(sampl%in%success),bkgrnd.success=sum(bkgrnd%in%success),samp.fail=sum(!(sampl%in%success)),bkgrnd.fail=sum(!(bkgrnd%in%success)))
      test_mat=matrix(unlist(test_res),nrow=2,dimnames=list(c('samp','bkgrnd'),c('success','fail')))
      test_out=fisher.test(test_mat,alternative=tail)
  
#   print(test_mat)

      test_res$n.genes=length(sampl)
      test_res$FETp=(test_out$p.value)
      test_res$fetOR=round(test_out$estimate) #,digits=3 
      test_res$lowerCI=round(test_out$conf.int[1])#,digits=3
      test_res$upperCI=round(test_out$conf.int[2])#,digits=3
#   test_res$samp.success=paste(sampl[sampl%in%success],collapse=' ')
      return(((test_res)))
  }

  if(counts){
    test_res=list(samp.success=samp.success,bkgrnd.success=bkgrnd.success,samp.fail=samp.fail,bkgrnd.fail=bkgrnd.fail)
    test_mat=matrix(unlist(test_res),nrow=2,dimnames=list(c('samp','bkgrnd'),c('success','fail')))
    test_out=fisher.test(test_mat,alternative=tail)
    
#      print(test_mat)

    test_res$n.genes=sum(samp.success,samp.fail)
    test_res$FETp=(test_out$p.value)
    test_res$fetOR=round(test_out$estimate,digits=3)
    test_res$lowerCI=round(test_out$conf.int[1],digits=3)
    test_res$upperCI=round(test_out$conf.int[2],digits=3)
#   test_res$samp.success=paste(sampl[sampl%in%success],collapse=' ')
    return(((test_res)))

  }
}




dnmr<-function(dat_lis,dtb='default',phen=''){ #bkg,
##  dtb - required list - DNM in 'conrols' ie healthy parents & offspring
##   dnmDB - pre-made dataset, part of "adds" package, downloaded from http://denovo-db.gs.washington.edu/denovo-db/  ## mapped to HUGO gene ids, genes with DNM in controls removed from DNM in
  if(dtb[1]=='default'){
    cat('\tLoad(~/Dropbox/PROJ/ednm/dtb/denovo-db.variants.v.1.5__frameshift_missense_stopgain.Rdata)\n')
    Load('~/Dropbox/PROJ/ednm/dtb/denovo-db.variants.v.1.5__frameshift_missense_stopgain.Rdata')
    dtb=dnmd
  }      ##  load pre-made dataset, part of "adds" package
  if(!('control'%in%names(dtb))){stop('"control" - list of DNM in healthy controls && offspring (named "control") is required')}


   cat('\tcheck DNM dtb and background compatibility\n')
  # bkg=overlap(unlist(dtb),bkg)

  # bkg=c(bkg$inter,bkg$inb)  ##  background is specified as overap && 'expressed' - ie have significant signal in the dataset used to derive the list
  ## the above definition is the same as the "bkg" input by definiton..

  dtb_contr=dtb$control
  dtb=dtb[names(dtb)!='control']

   cat('\n')
  # str(bkg)
  if(length(intersect(unlist(dtb),unlist(dat_lis)))==0){stop('check that dat_lis and bkg IDs are HUGO gene names OR match the provided dtb')}
  # dtb=lapply(dtb,function(x){x[x%in%bkg]})
  dnmen=list()
  fetp=list()
  cat('\n\nperform FET enrichment\n')
  for(idnm in names(dtb)){
    cat('\t',idnm)
    holder=list()
    for(idat in names(dat_lis)){
      holder[[idat]]=unlist(fet(
              # sampl=dat_lis[[idat]]
              # ,bkgrnd=bkg
              # ,success=dtb[[idnm]]
              # ,
              counts=T
              ,samp.success=sum(dtb[[idnm]]%in%dat_lis[[idat]])
              ,bkgrnd.success=sum(dtb_contr%in%dat_lis[[idat]])
              ,samp.fail=sum(!(dtb[[idnm]]%in%dat_lis[[idat]]))
              ,bkgrnd.fail=sum(!(dtb_contr%in%dat_lis[[idat]]))
              ,tail='greater'
              ))

    }
    dnmen[[idnm]]=as.data.frame(t(as.data.frame(holder)))
    fetp[[idnm]]=dnmen[[idnm]]$FETp

  }
    cat('\n\n')

  fetp=as.data.frame(fetp)
    rownames(fetp)=names(dat_lis)

  return(list(fetp=fetp,dnmen=dnmen))
}

