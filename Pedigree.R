library(MASS)
library(fda)

fourier.expansion<- function(x,n_of_basis,pos){
        frange <- c(pos[1], pos[length(pos)])
        rlt=list();
        rlt$fbasis<-create.fourier.basis(frange,nbasis=n_of_basis)
        rlt$phi = eval.basis(pos,rlt$fbasis);
        rlt$coef<-ginv(t(rlt$phi)%*%rlt$phi)%*%t(rlt$phi)%*%t(x)
        return(rlt)
}

fourier.expansion.smoothed<- function(x,n_of_basis,pos,lambda){
        frange <- c(pos[1], pos[length(pos)])
        rlt=list();
        rlt$fbasis<-create.fourier.basis(frange,nbasis=n_of_basis)

        rlt$phi = eval.basis(pos,rlt$fbasis) + lambda* eval.basis(pos,rlt$fbasis,2);
        rlt$coef<-ginv(t(rlt$phi)%*%rlt$phi)%*%t(rlt$phi)%*%t(x)
        
        return(rlt)
}

mp<-function(A,B){
  if (is.vector(A)) A<-matrix(A,nrow=1)
  if (is.vector(B)) B<-matrix(B,ncol=1)

  C<-array(NA,dim=c(nrow(A),ncol(B)))
  if (ncol(A)!=nrow(B)) stop("Can not do matrix product because of mismatched dimension in A and B.")
  for (i in 1:nrow(A)){
      for (j in 1:ncol(B)){
           temp<-NULL
           for (k in 1:ncol(A)){
                temp[k]<-A[i,k]*B[k,j]
            }
           C[i,j]<-sum(temp,na.rm=TRUE)
      }
   }
  return(C)
}

fpca.genotype <- function(x,cc,phi,pos=NULL,percentage=0.8,nbasis=37,lambda=NULL){
        nsnps <- dim(x)[2]
        ninds <- dim(x)[1]
        if ( ninds != length(cc) ){
                stop("Individual in data are not matched with those in affected status")
        }
        if ( is.null(pos) ){
                pos <- (0:( nsnps-1) )/(nsnps-1)
        }else {
                idx<-order(pos)
                x<-x[,idx]
                pos<-pos[idx]
                pos<- (pos-pos[1])/(pos[nsnps]-pos[1])
        }
        dataca <- x[cc==1,]
        dataco <- x[cc==2,]
        nA<-length(cc[cc==1]);
        nG<-length(cc[cc==2]);

        if( is.null(lambda)){
                expanded<-fourier.expansion(x,nbasis,pos)
        }else if( lambda > 0 ){
                expanded<-fourier.expansion.smoothed(x,nbasis,pos,lambda)
        }
        coef<-t(expanded$coef-rowMeans(expanded$coef))/sqrt(nA+nG)
        pca.rlt<-prcomp(coef)
        pca.rlt$scores<-coef%*%pca.rlt$rotation
        xi<-pca.rlt$scores[cc==1,]
        xa<-pca.rlt$scores[cc==2,]

        xi.bar <- colMeans(xi)
        xa.bar <- colMeans(xa)
        v0<-diag(var(pca.rlt$scores))
        vec<-(xi.bar-xa.bar)^2/(1/nA+1/nG)/v0;
        df<-1:length(vec)
        rlt<-list()
		dd=2-cc
		adjustment = as.numeric( nA*nG/(nA+nG)/( t(dd-nA/(nA+nG))%*%phi%*%(dd-nA/(nA+nG) ) ) )
        rlt$stat<-vec * adjustment
        rlt$pv.all<-1-pchisq(rlt$stat,df)
        rlt$prop<-cumsum(v0)/sum(v0)
        rlt$pv<-rlt$pv.all[rlt$prop>percentage][1]
		rlt$min.pv <- min(rlt$pv.all)
  
  
        return(rlt)
}


chi2.test<- function(x,cc,phi){
  cc<-as.vector(cc)
  x<-as.matrix(x)
  if( dim(x)[2] != 1) stop('Error in Chi-square: multiple dimension data')
  nA<-sum(cc==1)
  nG<-sum(cc==2)
  AA<-sum(x[cc==1,]>0)
  Aa<-nA-AA
  aA<-sum(x[cc==2,]>0)
  aa<-nG-aA
  rlt<-1
  dd=2-cc
  if( AA > 5 && Aa > 5 && aA > 5 && aa > 5){
    rlt<-chisq.test(matrix(c(AA,Aa,aA,aa),2,2))$statistic
    rlt = rlt *nA*nG/(nA+nG)/( t(dd-nA/(nA+nG))%*%phi%*%(dd-nA/(nA+nG) ) )
    rlt = pchisq(rlt,2,lower.tail=F)
  }else{
    rlt<-fisher.test( matrix(c(AA,Aa,aA,aa),2,2) )$p.value
    rlt<-qchisq(rlt,2,lower.tail=F)
     rlt = rlt *nA*nG/(nA+nG)/( t(dd-nA/(nA+nG))%*%phi%*%(dd-nA/(nA+nG) ) )
     rlt = pchisq(rlt,2,lower.tail=F)
  }
  return(rlt  )
}


chi2.min<- function(x,cc,phi){
  chi2.pv=c()
  for ( i in 1:dim(x)[2] ){
    chi2.pv[i]<- chi2.test(x[,i],cc,phi)
  }
  stat<-min( na.omit(chi2.pv) )
  rlt= stat
  return(rlt)
}

chi2.permutation<- function(x,cc,phi,times=5000){
  chi2.pv=c()
  for ( i in 1:dim(x)[2] ){
    chi2.pv[i]<- chi2.test(x[,i],cc,phi)
  }
  stat<-min( na.omit(chi2.pv) )
  if( is.na(stat) ){
    return(NA)
  }
  nv<-1;
  t<-1;
  while ( (t< times) && ( nv < 50 ) ){
    for (j in 1:100){
          pheno<-sample(cc)
          for ( i in 1:dim(x)[2] ){
                chi2.pv[i]<- chi2.test(x[,i],pheno,phi)
          }
          temp.stat<-min( na.omit(chi2.pv) )
          if( !is.na( temp.stat) ){
            if (stat >= temp.stat ) {
              nv<- nv + 1
            }
          }
        }
        t <- t+100;
  }
  pv<-nv/t
  rlt=pv
  return(rlt)
}

pedT2<-function(x,cc,phi){
    if( dim(x)[2]>1 ){
  ca<-as.matrix( x[cc==1,] )
  co<-as.matrix( x[cc==2,] )
  nA<-dim(ca)[1]
  nG<-dim(co)[1]
  s0<-1/(nA+nG-2)*((nA-1)*cov(ca)+(nG-1)*cov(co) );
  m1<-colMeans(ca)
  m2<-colMeans(co)
  rnk<-qr(s0)$rank
  t2<-(nA*nG)/(nA+nG)*t(m1-m2)%*%ginv(s0)%*%(m1-m2)
  dd=2-cc
  t2=t2*nA*nG/(nA+nG)/( t(dd-nA/(nA+nG))%*%phi%*%(dd-nA/(nA+nG) ) )
  pv<- pchisq(t2,rnk,lower.tail=F)
 }else{
  pv<-chi2.test(x,cc,phi)
 }
 pv
}

pedCMC <- function(x,cc,maf,phi,level){
 
 ic1<- (maf<level)
 ic2<- (maf>=level) 
 
 snps<- sum(ic1)
 if( snps > 0 ){
  m<-apply( as.matrix(x[,ic1]),1,max);
  reduced.x<-cbind(x[,ic2],m);
 }else{
  reduced.x<-x;
 }
 pv<-pedT2(reduced.x,cc,phi);
 pv  
}

estimatePhi <- function(x,maf){
 inds=dim(x)[1]
 x=x[,maf>0.05]
 maf=maf[maf>0.05]
 snps=dim(x)[2]
 rlt=matrix(0,inds,inds)
 
 for( i in 1:inds ){
  for( j in i:inds){
   if( j == i ){
    xi=as.vector( x[i,])
    temp=1+sum( ( xi^2-(1+2*maf)*xi +2*maf^2 )/maf/(1-maf)/2 ,na.rm=T)/sum ( !is.na(xi ) ) 
	rlt[i,i]=max(temp,1)
	rlt[i,i]=min(rlt[i,i],2)
   }else{
    xi=as.vector( x[i,] )
    xj=as.vector( x[j,] )
    temp = sum( (xi-2*maf)*(xj-2*maf)/maf/(1-maf)/2,na.rm=T )/(sum( !is.na(xi) & !is.na(xj) ) )
	rlt[i,j]= max(temp,0)
	rlt[i,j] = min(rlt[i,j],0.5)
    rlt[j,i]=rlt[i,j]
   }
  }
 }
 rlt[rlt == "NaN"] <- 0
 kinship = matrix(nrow=0,ncol=3)
 colnames(kinship)=c("id1","id2","kinship")
 for( i in 1:inds){
   for( j in i:inds){
     if( rlt[i,j] != 0 ){
        kinship =rbind(kinship, matrix(nrow=1,ncol=3,c(i,j,rlt[i,j])) )
     }
   }
 }
 
 write.csv(rlt,"kinship_matrix.csv",row.names=F)
 write.table(kinship,"kinship.txt",row.names=F)
 rlt
 
}

TestPedigree<-function(pedData,kinship=NULL,map,lambda=1E-5,basis=33,permu_times=5000){
  columns = dim(pedData)[2]
  inds = dim(pedData)[1]
  temp = as.matrix(pedData[6:columns])
  temp[temp==0]=NA
  d1=dim(temp)[1]
  d2=dim(temp)[2]
  ref = apply(temp,2,function(x) min(x,na.rm=T) )
  ref[ seq(2,d2,by=2) ] = ref[ seq(1,d2,by=2) ]
  for ( i in 1:d2){
    temp[ temp[,i] != ref[i], i] = 0
	temp[ temp[,i] == ref[i], i] = 1
  }
  temp <- apply(temp, 2 ,as.numeric)
  genodata=matrix(0,d1,d2/2)
  for( i in 1:(d2/2) ){
    genodata[,i] = temp[,i*2-1] + temp[,i*2]
  }
  
  maf = colMeans(genodata,na.rm=T)/2
  genodata[,maf>0.5] = 2-genodata[,maf>0.5]
  
  pheno = pedData[,5] 
  maf = colMeans(genodata,na.rm=T)/2
  dphi = matrix(0,d1,d1)
  if( is.null(kinship) ){
    dphi = estimatePhi(genodata,maf)
  }else{
	for( i in 1:length(pheno) ){
	   dphi[i,i] = 1
	}
	for( i in 1:dim(kinship)[1] ){
	 dphi[ kinship[i,1], kinship[i,2] ] = kinship[i,3]
	 dphi[ kinship[i,2], kinship[i,1] ] = kinship[i,3]
	}
  }
  gene.group = map[,4]
  rlt = pedigree.test(gene.group,pheno,genodata,dphi,permu_times,basis,lambda)
  return(rlt)
}



pedigree.test<- function( gene.group, pheno, genodata,phi=NULL,permu=5000,basis=33,dlambda=1E-6){
	genodata = as.matrix(genodata)
	genodata = genodata[!is.na(pheno),]
	pheno = pheno[!is.na(pheno) ]
	if( !is.matrix(genodata) ) stop('Data cannot be convert to matrix')
	genodata[is.na(genodata)] = 0
	wmaf = colMeans(genodata)/2
	genodata[,wmaf>0.5]=2-genodata[,wmaf>0.5]
	genodata = genodata[,wmaf>0]
	gene.group = gene.group[wmaf>0]
	wmaf = wmaf[wmaf>0]
	
	if( is.null(phi) ){
	  phi = estimatePhi(genodata,wmaf)
	}
	
	print("phi estimated")
	
	if( length(pheno) != dim(genodata)[1] ) stop('The dimension of phenotype is not match with sample size')
	rlt<-data.frame( )
	
	
	genes = levels(gene.group)
	gene_num = length(genes)
	for( i in 1:gene_num){
		idx = ( gene.group == genes[i] )
		x = genodata[,idx]
		snps = dim(x)[2]
		if ( snps >= 3 ){
			#print(length(pos))
			#print(dim(x))
			rlt[i,"GENE"]<- genes[i]
			maf=colMeans(x)/2
			rlt[i,"Snp_Number"]<-snps			#record how many snps is in the test
			pos=( 0:(snps-1) )/(snps -1 )
			
			rlt[i,"FPCA"]<-as.numeric( try( fpca.genotype(x,pheno,phi,pos,0.8,nbasis=basis,lambda=dlambda)$pv ) )			
			rlt[i,"Chi_permutation"]<- as.numeric( try( chi2.permutation(x,pheno,phi,permu) ) )
			rlt[i,"Chi_min"]<- as.numeric( try(chi2.min(x,pheno,phi)) )
			rlt[i,"T2"]<- as.numeric( try( pedT2(x,pheno,phi) ) )
			rlt[i,"CMC"]<- as.numeric( try( pedCMC(x,pheno,maf,phi,0.05)) )
		}
	}
	rlt
}


 









 



 